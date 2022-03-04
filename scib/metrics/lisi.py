import itertools
import logging
import multiprocessing as mp
import os
import pathlib
import subprocess
import tempfile

import anndata2ri
import numpy as np
import pandas as pd
import rpy2.rinterface_lib.callbacks
import rpy2.robjects as ro
import scanpy as sc
import scipy.sparse
from scipy.io import mmwrite

from ..utils import check_adata, check_batch

rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)  # Ignore R warning messages


# Main LISI

def lisi(
        adata,
        batch_key,
        label_key,
        k0=90,
        type_=None,
        scale=True,
        verbose=False
):
    """
    Compute lisi score (after integration)
    params:
        matrix: matrix from adata to calculate on
        covariate_key: variable to compute iLISI on
        cluster_key: variable to compute cLISI on
    return:
        pd.DataFrame with median cLISI and median iLISI scores (following the harmony paper)
    """

    check_adata(adata)
    check_batch(batch_key, adata.obs)
    check_batch(label_key, adata.obs)

    # if type_ != 'knn':
    #    if verbose:
    #        print("recompute kNN graph with {k0} nearest neighbors.")
    # recompute neighbours
    if (type_ == 'embed'):
        adata_tmp = sc.pp.neighbors(adata, n_neighbors=k0, use_rep='X_emb', copy=True)
    elif (type_ == 'full'):
        if 'X_pca' not in adata.obsm.keys():
            sc.pp.pca(adata, svd_solver='arpack')
        adata_tmp = sc.pp.neighbors(adata, n_neighbors=k0, copy=True)
    else:
        adata_tmp = adata.copy()
    # if knn - do not compute a new neighbourhood graph (it exists already)

    # lisi_score = lisi_knn(adata=adata, batch_key=batch_key, label_key=label_key, verbose=verbose)
    lisi_score = lisi_knn_py(adata=adata_tmp, batch_key=batch_key, label_key=label_key, verbose=verbose)

    # iLISI: nbatches good, 1 bad
    ilisi_score = np.nanmedian(lisi_score[batch_key])
    # cLISI: 1 good, nbatches bad
    clisi_score = np.nanmedian(lisi_score[label_key])

    if scale:
        # get number of batches
        nbatches = len(np.unique(adata.obs[batch_key]))
        ilisi_score, clisi_score = scale_lisi(ilisi_score, clisi_score, nbatches)

    return ilisi_score, clisi_score


def lisi_knn_py(
        adata,
        batch_key,
        label_key,
        perplexity=None,
        verbose=False
):
    """
    Compute LISI score on kNN graph provided in the adata object. By default, perplexity
    is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """

    if 'neighbors' not in adata.uns:
        raise AttributeError(f"key 'neighbors' not found. Please make sure that a " +
                             "kNN graph has been computed")
    elif verbose:
        print("using precomputed kNN graph")

    # get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")
    dist_mat = scipy.sparse.find(adata.obsp['distances'])
    # get number of nearest neighbours parameter
    if 'params' not in adata.uns['neighbors']:
        # estimate the number of nearest neighbors as the median
        # of the distance matrix
        _, e = np.unique(dist_mat[0], return_counts=True)
        n_nn = np.nanmedian(e)
        n_nn = n_nn.astype('int')
    else:
        n_nn = adata.uns['neighbors']['params']['n_neighbors'] - 1
    # initialise index and fill it with NaN values
    nn_index = np.empty(shape=(adata.obsp['distances'].shape[0],
                               n_nn))
    nn_index[:] = np.NaN
    nn_dists = np.empty(shape=(adata.obsp['distances'].shape[0],
                               n_nn))
    nn_dists[:] = np.NaN
    index_out = []
    for cell_id in np.arange(np.min(dist_mat[0]), np.max(dist_mat[0]) + 1):
        get_idx = dist_mat[0] == cell_id
        num_idx = get_idx.sum()
        # in case that get_idx contains more than n_nn neighbours, cut away the outlying ones
        fin_idx = np.min([num_idx, n_nn])
        nn_index[cell_id, :fin_idx] = dist_mat[1][get_idx][np.argsort(dist_mat[2][get_idx])][:fin_idx]
        nn_dists[cell_id, :fin_idx] = np.sort(dist_mat[2][get_idx])[:fin_idx]
        if num_idx < n_nn:
            index_out.append(cell_id)

    out_cells = len(index_out)

    if out_cells > 0:
        if verbose:
            print(f"{out_cells} had less than {n_nn} neighbors.")

    if perplexity is None:
        # use LISI default
        perplexity = np.floor(nn_index.shape[1] / 3)

    # run LISI in python
    if verbose:
        print("importing knn-graph")

    batch = adata.obs[batch_key].cat.codes.values
    n_batches = len(np.unique(adata.obs[batch_key]))
    label = adata.obs[label_key].cat.codes.values
    n_labels = len(np.unique(adata.obs[label_key]))

    if verbose:
        print("LISI score estimation")

    simpson_estimate_batch = compute_simpson_index(D=nn_dists,
                                                   knn_idx=nn_index,
                                                   batch_labels=batch,
                                                   n_batches=n_batches,
                                                   perplexity=perplexity,
                                                   )
    simpson_estimate_label = compute_simpson_index(D=nn_dists,
                                                   knn_idx=nn_index,
                                                   batch_labels=label,
                                                   n_batches=n_labels,
                                                   perplexity=perplexity
                                                   )
    simpson_est_batch = 1 / simpson_estimate_batch
    simpson_est_label = 1 / simpson_estimate_label
    # extract results
    d = {batch_key: simpson_est_batch, label_key: simpson_est_label}
    lisi_estimate = pd.DataFrame(data=d, index=np.arange(0, len(simpson_est_label)))

    return lisi_estimate


# Graph LISI (analoguous to lisi function)
def lisi_graph(
        adata,
        batch_key,
        label_key,
        **kwargs
):
    """
    Compute cLISI and iLISI scores on precomputed kNN graph

    :param adata: adata object to calculate on
    :param batch_key: batch column name in adata.obs
    :param label_key: label column name in adata.obs
    :param **kwargs: arguments to be passed to iLISI and cLISI functions
    :return:
        Median cLISI and iLISI scores
    """
    ilisi = ilisi_graph(adata, batch_key=batch_key, **kwargs)
    clisi = clisi_graph(adata, batch_key=batch_key, label_key=label_key, **kwargs)
    return ilisi, clisi


def ilisi_graph(
        adata,
        batch_key,
        k0=90,
        type_=None,
        subsample=None,
        scale=True,
        multiprocessing=None,
        nodes=None,
        verbose=False
):
    """
    Compute iLISI score adapted from Harmony paper (Korsunsky et al, Nat Meth, 2019)

    :param adata: adata object to calculate on
    :param batch_key: batch column name in adata.obs
    :param k0: number of nearest neighbors to compute lisi score
        Please note that the initial neighborhood size that is
        used to compute shortest paths is 15.
    :param type_: type of data integration, either knn, full or embed
    :param subsample: Percentage of observations (integer between 0 and 100)
        to which lisi scoring should be subsampled
    :param scale: scale output values between 0 and 1 (True/False)
    :param multiprocessing: parallel computation of LISI scores, if None, no parallisation
        via multiprocessing is performed
    :param nodes: number of nodes (i.e. CPUs to use for multiprocessing); ignored, if
        multiprocessing is set to None
    :return: Median of iLISI score
    """

    check_adata(adata)
    check_batch(batch_key, adata.obs)

    adata_tmp = recompute_knn(adata, type_)
    ilisi_score = lisi_graph_py(
        adata=adata_tmp,
        batch_key=batch_key,
        n_neighbors=k0,
        perplexity=None,
        subsample=subsample,
        multiprocessing=multiprocessing,
        nodes=nodes,
        verbose=verbose
    )

    # iLISI: nbatches good, 1 bad
    ilisi = np.nanmedian(ilisi_score)

    if scale:
        nbatches = adata.obs[batch_key].nunique()
        ilisi = (ilisi - 1) / (nbatches - 1)

    return ilisi


def clisi_graph(
        adata,
        batch_key,
        label_key,
        k0=90,
        type_=None,
        subsample=None,
        scale=True,
        multiprocessing=None,
        nodes=None,
        verbose=False
):
    """
    Compute cLISI score adapted from Harmony paper (Korsunsky et al, Nat Meth, 2019)

    :params adata: adata object to calculate on
    :param batch_key: batch column name in adata.obs
    :param label_key: label column name in adata.obs
    :param k0: number of nearest neighbors to compute lisi score
        Please note that the initial neighborhood size that is
        used to compute shortest paths is 15.
    :param type_: type of data integration, either knn, full or embed
    :param subsample: Percentage of observations (integer between 0 and 100)
        to which lisi scoring should be subsampled
    :param scale: scale output values between 0 and 1 (True/False)
    :param multiprocessing: parallel computation of LISI scores, if None, no parallisation
        via multiprocessing is performed
    :param nodes: number of nodes (i.e. CPUs to use for multiprocessing); ignored, if
        multiprocessing is set to None
    :return: Median of cLISI score
    """

    check_adata(adata)
    check_batch(batch_key, adata.obs)
    check_batch(label_key, adata.obs)

    adata_tmp = recompute_knn(adata, type_)

    scores = lisi_graph_py(
        adata=adata_tmp,
        batch_key=label_key,
        n_neighbors=k0,
        perplexity=None,
        subsample=subsample,
        multiprocessing=multiprocessing,
        nodes=nodes,
        verbose=verbose
    )

    # cLISI: 1 good, nlabs bad
    clisi = np.nanmedian(scores)

    if scale:
        nlabs = adata.obs[label_key].nunique()
        clisi = (nlabs - clisi) / (nlabs - 1)

    return clisi


def recompute_knn(adata, type_):
    """
    Recompute neighbours
    """
    if type_ == 'embed':
        return sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_emb', copy=True)
    elif type_ == 'full':
        if 'X_pca' not in adata.obsm.keys():
            sc.pp.pca(adata, svd_solver='arpack')
        return sc.pp.neighbors(adata, n_neighbors=15, copy=True)
    else:
        # if knn - do not compute a new neighbourhood graph (it exists already)
        return adata.copy()


def lisi_graph_py(
        adata,
        batch_key,
        n_neighbors=90,
        perplexity=None,
        subsample=None,
        multiprocessing=None,
        nodes=None,
        verbose=False
):
    """
    Function to prepare call of compute_simpson_index
    Compute LISI score on shortes path based on kNN graph provided in the adata object.
    By default, perplexity is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """

    if 'neighbors' not in adata.uns:
        raise AttributeError(f"key 'neighbors' not found. Please make sure that a " +
                             "kNN graph has been computed")
    elif verbose:
        print("using precomputed kNN graph")

    # get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")

    batch = adata.obs[batch_key].cat.codes.values
    n_batches = len(np.unique(adata.obs[batch_key]))

    if perplexity is None or perplexity >= n_neighbors:
        # use LISI default
        perplexity = np.floor(n_neighbors / 3)

    # setup subsampling
    subset = 100  # default, no subsampling
    if subsample is not None:
        subset = subsample  # do not use subsampling
        if isinstance(subsample, int) == False:  # need to set as integer
            subset = int(subsample)

    # run LISI in python
    if verbose:
        print("Compute knn on shortest paths")

        # set connectivities to 3e-308 if they are lower than 3e-308 (because cpp can't handle double values smaller than that).
    connectivities = adata.obsp['connectivities']  # csr matrix format
    large_enough = connectivities.data >= 3e-308
    if verbose:
        n_too_small = np.sum(large_enough == False)
        if n_too_small:
            print(f"{n_too_small} connectivities are smaller than 3e-308 and will be set to 3e-308")
            print(connectivities.data[large_enough == False])
    connectivities.data[large_enough == False] = 3e-308

    # define number of chunks
    n_chunks = 1

    if multiprocessing is not None:
        # set up multiprocessing
        if nodes is None:
            # take all but one CPU and 1 CPU, if there's only 1 CPU.
            n_cpu = mp.cpu_count()
            n_processes = np.max([n_cpu, np.ceil(n_cpu / 2)]).astype('int')
        else:
            n_processes = nodes
        # update numbr of chunks
        n_chunks = n_processes

    # temporary file
    tmpdir = tempfile.TemporaryDirectory(prefix="lisi_")
    dir_path = tmpdir.name + '/'
    mtx_file_path = dir_path + 'input.mtx'
    print(mtx_file_path, dir_path)
    mmwrite(
        mtx_file_path,
        connectivities,
        symmetry='general'
    )
    # call knn-graph computation in Cpp

    root = pathlib.Path(__file__).parent.parent  # get current root directory
    cpp_file_path = root / 'knn_graph/knn_graph.o'  # create POSIX path to file to execute compiled cpp-code
    # comment: POSIX path needs to be converted to string - done below with 'as_posix()'
    # create evenly split chunks if n_obs is divisible by n_chunks (doesn't really make sense on 2nd thought)
    n_splits = n_chunks - 1
    args_int = [cpp_file_path.as_posix(), mtx_file_path, dir_path, str(n_neighbors), str(n_splits), str(subset)]
    try:
        subprocess.run(args_int)
    except Exception as e:
        print(e)
        print("Couldn't compute LISI, returning NaN")
        return np.nan

    if verbose:
        print("LISI score estimation")

    # do the simpson call
    if multiprocessing is not None:

        if verbose:
            print(f"{n_processes} processes started.")
        pool = mp.Pool(processes=n_processes)
        count = np.arange(0, n_processes)

        # create argument list for each worker
        results = pool.starmap(
            compute_simpson_index_graph,
            zip(itertools.repeat(dir_path),
                itertools.repeat(batch),
                itertools.repeat(n_batches),
                itertools.repeat(n_neighbors),
                itertools.repeat(perplexity),
                count)
        )
        pool.close()
        pool.join()

        simpson_est_batch = 1 / np.concatenate(results)

    else:
        simpson_estimate_batch = compute_simpson_index_graph(
            input_path=dir_path,
            batch_labels=batch,
            n_batches=n_batches,
            perplexity=perplexity,
            n_neighbors=n_neighbors,
            chunk_no=None
        )
        simpson_est_batch = 1 / simpson_estimate_batch

    tmpdir.cleanup()

    # extract results
    d = {batch_key: simpson_est_batch}

    lisi_estimate = pd.DataFrame(data=d, index=np.arange(0, len(simpson_est_batch)))

    return lisi_estimate


# LISI core functions

def compute_simpson_index(
        D=None,
        knn_idx=None,
        batch_labels=None,
        n_batches=None,
        perplexity=15,
        tol=1e-5
):
    """
    Simpson index of batch labels subsetted for each group.
    params:
        D: distance matrix n_cells x n_nearest_neighbors
        knn_idx: index of n_nearest_neighbors of each cell
        batch_labels: a vector of length n_cells with batch info
        n_batches: number of unique batch labels
        perplexity: effective neighborhood size
        tol: a tolerance for testing effective neighborhood size
    returns:
        simpson: the simpson index for the neighborhood of each cell
    """
    n = D.shape[0]
    P = np.zeros(D.shape[1])
    simpson = np.zeros(n)
    logU = np.log(perplexity)

    # loop over all cells
    for i in np.arange(0, n, 1):
        beta = 1
        # negative infinity
        betamin = -np.inf
        # positive infinity
        betamax = np.inf
        # get active row of D
        D_act = D[i, :]
        H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0
        # first get neighbor probabilities
        while (np.logical_and(np.abs(Hdiff) > tol, tries < 50)):
            if (Hdiff > 0):
                betamin = beta
                if (betamax == np.inf):
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if (betamin == -np.inf):
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2

            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1

        if (H == 0):
            simpson[i] = -1
            continue

            # then compute Simpson's Index
        non_nan_knn = knn_idx[i][np.invert(np.isnan(knn_idx[i]))].astype('int')
        batch = batch_labels[non_nan_knn]
        # convertToOneHot omits all nan entries.
        # Therefore, we run into errors in np.matmul.
        if len(batch) == len(P):
            B = convertToOneHot(batch, n_batches)
            sumP = np.matmul(P, B)  # sum P per batch
            simpson[i] = np.dot(sumP, sumP)  # sum squares
        else:  # assign worst possible score
            simpson[i] = 1

    return simpson


def compute_simpson_index_graph(
        input_path=None,
        batch_labels=None,
        n_batches=None,
        n_neighbors=90,
        perplexity=30,
        chunk_no=0,
        tol=1e-5
):
    """
    Simpson index of batch labels subsetted for each group.
    params:
        input_path: file_path to pre-computed index and distance files
        batch_labels: a vector of length n_cells with batch info
        n_batches: number of unique batch labels
        n_neighbors: number of nearest neighbors
        perplexity: effective neighborhood size
        chunk_no: for parallelisation, chunk id to evaluate
        tol: a tolerance for testing effective neighborhood size
    returns:
        simpson: the simpson index for the neighborhood of each cell
    """

    # initialize
    P = np.zeros(n_neighbors)
    logU = np.log(perplexity)

    if chunk_no is None:
        chunk_no = 0
    # check if the target file is not empty
    if os.stat(input_path + '_indices_' + str(chunk_no) + '.txt').st_size == 0:
        print("File has no entries. Doing nothing.")
        lists = np.zeros(0)
        return lists

    # read distances and indices with nan value handling
    indices = pd.read_table(input_path + '_indices_' + str(chunk_no) + '.txt', index_col=0, header=None, sep=',')
    indices = indices.T

    distances = pd.read_table(input_path + '_distances_' + str(chunk_no) + '.txt', index_col=0, header=None, sep=',')
    distances = distances.T

    # get cell ids
    chunk_ids = indices.columns.values.astype('int')

    # define result vector
    simpson = np.zeros(len(chunk_ids))

    # loop over all cells in chunk
    for i in enumerate(chunk_ids):
        # get neighbors and distances
        # read line i from indices matrix
        get_col = indices[i[1]]

        if get_col.isnull().sum() > 0:
            # not enough neighbors
            print(i[1] + " has not enough neighbors.")
            simpson[i[0]] = 1  # np.nan #set nan for testing
            continue
        else:
            knn_idx = get_col.astype('int') - 1  # get 0-based indexing

        # read line i from distances matrix
        D_act = distances[i[1]].values.astype('float')

        # start lisi estimation
        beta = 1
        # negative infinity
        betamin = -np.inf
        # positive infinity
        betamax = np.inf

        H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0
        # first get neighbor probabilities
        while (np.logical_and(np.abs(Hdiff) > tol, tries < 50)):
            if (Hdiff > 0):
                betamin = beta
                if (betamax == np.inf):
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if (betamin == -np.inf):
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2

            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1

        if (H == 0):
            simpson[i[0]] = -1
            continue
            # then compute Simpson's Index
        batch = batch_labels[knn_idx]
        B = convertToOneHot(batch, n_batches)
        sumP = np.matmul(P, B)  # sum P per batch
        simpson[i[0]] = np.dot(sumP, sumP)  # sum squares

    return simpson


def Hbeta(D_row, beta):
    """
    Helper function for simpson index computation
    """
    P = np.exp(- D_row * beta)
    sumP = np.nansum(P)
    if (sumP == 0):
        H = 0
        P = np.zeros(len(D_row))
    else:
        H = np.log(sumP) + beta * np.nansum(D_row * P) / sumP
        P /= sumP
    return H, P


def convertToOneHot(vector, num_classes=None):
    """
    Converts an input 1-D vector of integers into an output
    2-D array of one-hot vectors, where an i'th input value
    of j will set a '1' in the i'th row, j'th column of the
    output array.

    Example:
        v = np.array((1, 0, 4))
        one_hot_v = convertToOneHot(v)
        print one_hot_v

        [[0 1 0 0 0]
         [1 0 0 0 0]
         [0 0 0 0 1]]
    """

    # assert isinstance(vector, np.ndarray)
    # assert len(vector) > 0

    if num_classes is None:
        num_classes = np.max(vector) + 1
    # else:
    #    assert num_classes > 0
    #    assert num_classes >= np.max(vector)

    result = np.zeros(shape=(len(vector), num_classes))
    result[np.arange(len(vector)), vector] = 1
    return result.astype(int)


# DEPRECATED
# This code scales clisi incorrectly!
def scale_lisi(ilisi_score, clisi_score, nbatches):
    # scale iLISI score to 0 bad 1 good
    ilisi_score = (ilisi_score - 1) / (nbatches - 1)
    # scale clisi score to 0 bad 1 good
    clisi_score = (nbatches - clisi_score) / (nbatches - 1) # Scaled incorrectly by n_batches
    return ilisi_score, clisi_score


def lisi_knn(
        adata,
        batch_key,
        label_key,
        perplexity=None,
        verbose=False
):
    """
    Deprecated
    Compute LISI score on kNN graph provided in the adata object. By default, perplexity
    is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """

    if 'neighbors' not in adata.uns:
        raise AttributeError(
            "key 'neighbors' not found. Please make sure that a "
            "kNN graph has been computed"
        )
    elif verbose:
        print("using precomputed kNN graph")

    # get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")
    dist_mat = scipy.sparse.find(adata.obsp['distances'])
    # get number of nearest neighbours parameter
    if 'params' not in adata.uns['neighbors']:
        # estimate the number of nearest neighbors as the median
        # of the distance matrix
        _, e = np.unique(dist_mat[0], return_counts=True)
        n_nn = np.nanmin(e)
        n_nn = n_nn.astype('int')
    else:
        n_nn = adata.uns['neighbors']['params']['n_neighbors'] - 1
    nn_index = np.empty(shape=(adata.obsp['distances'].shape[0],
                               n_nn))
    nn_dists = np.empty(shape=(adata.obsp['distances'].shape[0],
                               n_nn))
    index_out = []
    for cell_id in np.arange(np.min(dist_mat[0]), np.max(dist_mat[0]) + 1):
        get_idx = dist_mat[0] == cell_id
        num_idx = get_idx.sum()
        # in case that get_idx contains more than n_nn neighbours, cut away the outlying ones
        # potential enhancement: handle case where less than n_nn neighbours are reported
        if num_idx >= n_nn:
            nn_index[cell_id, :] = dist_mat[1][get_idx][np.argsort(dist_mat[2][get_idx])][:n_nn]
            nn_dists[cell_id, :] = np.sort(dist_mat[2][get_idx])[:n_nn]
        else:
            index_out.append(cell_id)

    out_cells = len(index_out)

    if out_cells > 0:
        # remove all indexes in nn_index and nn_dists, which are 0
        # COMMENT: Terrible idea and commented out
        # nn_dists = np.delete(nn_dists, index_out, 0)
        # nn_index = np.delete(nn_index, index_out, 0)
        if verbose:
            print(f"{out_cells} had less than {n_nn} neighbors and were omitted in LISI score.")

    if perplexity is None:
        # use LISI default
        perplexity = np.floor(nn_index.shape[1] / 3)

    # run LISI in R
    anndata2ri.activate()
    ro.r("library(lisi)")

    if verbose:
        print("importing knn-graph")
    ro.globalenv['nn_indx'] = nn_index.astype('int').T
    ro.globalenv['nn_dst'] = nn_dists.T
    ro.globalenv['perplexity'] = perplexity
    ro.globalenv['batch'] = adata.obs[batch_key].cat.codes.values
    ro.globalenv['n_batches'] = len(np.unique(adata.obs[batch_key]))
    ro.globalenv['label'] = adata.obs[label_key].cat.codes.values
    ro.globalenv['n_labels'] = len(np.unique(adata.obs[label_key]))

    if verbose:
        print("LISI score estimation")
    simpson_estimate_batch = ro.r(
        f"simpson.estimate_batch <- compute_simpson_index(nn_dst, nn_indx, batch, n_batches, perplexity)")  # batch_label_keys)")
    simpson_estimate_label = ro.r(
        f"simpson.estimate_label <- compute_simpson_index(nn_dst, nn_indx, label, n_labels, perplexity)")  # batch_label_keys)")
    simpson_est_batch = 1 / np.squeeze(ro.r("simpson.estimate_batch"))
    simpson_est_label = 1 / np.squeeze(ro.r("simpson.estimate_label"))

    anndata2ri.deactivate()

    # extract results
    d = {batch_key: simpson_est_batch, label_key: simpson_est_label}
    lisi_estimate = pd.DataFrame(data=d, index=np.arange(0, len(simpson_est_label)))

    return lisi_estimate


def lisi_matrix(
        adata,
        batch_key,
        label_key,
        matrix=None,
        verbose=False
):
    """
    Deprecated
    Computes the LISI scores for a given data matrix in adata.X. The scoring function of the
    LISI R package is called with default parameters. This function takes a data matrix and
    recomputes nearest neighbours.
    """

    if matrix is None:
        matrix = adata.X

    # lisi score runs only on dense matrices (knn search)
    if scipy.sparse.issparse(matrix):
        matrix = matrix.todense()

    # run LISI in R
    anndata2ri.activate()
    ro.r("library(lisi)")

    if verbose:
        print("importing expression matrix")
    ro.globalenv['data_mtrx'] = matrix

    if verbose:
        print(f"covariates: {batch_key} and {label_key}")
    metadata = adata.obs[[batch_key, label_key]]
    ro.globalenv['metadata'] = metadata
    batch_label_keys = ro.StrVector([batch_key, label_key])
    ro.globalenv['batch_label_keys'] = batch_label_keys

    if verbose:
        print("LISI score estimation")
    lisi_estimate = ro.r(f"lisi.estimate <- compute_lisi(data_mtrx, metadata, batch_label_keys)")  # batch_label_keys)")
    anndata2ri.deactivate()

    return lisi_estimate
