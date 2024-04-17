import itertools
import multiprocessing as mp
import os
import pathlib
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from numba import njit
from scipy.io import mmwrite

import scib

from ..utils import check_adata, check_batch


# Graph LISI (analogous to lisi function)
def lisi_graph(adata, batch_key, label_key, **kwargs):
    """cLISI and iLISI scores

    This is a reimplementation of the LISI (Local Inverse Simpson’s Index) metrics
    https://doi.org/10.1038/s41592-019-0619-0

    see :func:`~scib.metrics.clisi_graph` and :func:`~scib.metrics.ilisi_graph`

    :param adata: adata object to calculate on
    :param batch_key: batch column name in ``adata.obs``
    :param label_key: label column name in ``adata.obs``
    :params \\**kwargs: arguments to be passed to :func:`~scib.metrics.clisi_graph` and :func:`~scib.metrics.ilisi_graph`
    :return: Overall cLISI and iLISI scores
    """
    ilisi = ilisi_graph(adata, batch_key=batch_key, **kwargs)
    clisi = clisi_graph(adata, label_key=label_key, **kwargs)
    return ilisi, clisi


def ilisi_graph(
    adata,
    batch_key,
    type_,
    use_rep="X_emb",
    k0=90,
    subsample=None,
    scale=True,
    n_cores=1,
    verbose=False,
):
    """Integration LISI (iLISI) score

    Local Inverse Simpson’s Index metrics adapted from `Korsunsky et al. 2019`_ to run on all full
    feature, embedding and kNN integration outputs via shortest path-based distance computation on single-cell kNN
    graphs.
    By default, this function returns a value scaled between 0 and 1 instead of the original LISI range of 0 to the
    number of batches.

    .. _Korsunsky et al. 2019: https://doi.org/10.1038/s41592-019-0619-0

    :param adata: adata object to calculate on
    :param batch_key: batch column name in ``adata.obs``
    :param `type_`: type of data integration, one of 'knn', 'embed' or 'full'
    :param use_rep: embedding slot in ``.obsm``, only used for embedding input
    :param k0: number of nearest neighbors to compute lisi score
        Please note that the initial neighborhood size that is
        used to compute shortest paths is 15.
    :param subsample: Percentage of observations (integer between 0 and 100)
        to which lisi scoring should be subsampled
    :param scale: scale output values between 0 and 1 (True/False)
    :param n_cores: number of cores (i.e. CPUs or CPU cores to use for multiprocessing)
    :return: Median of iLISI scores per batch labels

    This function can be applied to all integration output types and recomputes the kNN graph for feature and embedding
    output with specific parameters.
    Thus, no preprocessing is required, but the correct output type must be specified in ``type_``.

    **Examples**

    .. code-block:: python

        # feature output or unintegrated object
        scib.me.ilisi_graph(adata, batch_key="batch", type="full")

        # embeding output
        scib.me.ilisi_graph(adata, batch_key="batch", type="embed", use_rep="X_emb")

        # knn output
        scib.me.ilisi_graph(adata, batch_key="batch", type="knn")

    """

    check_adata(adata)
    check_batch(batch_key, adata.obs)

    adata_tmp = recompute_knn(adata, type_, use_rep)
    ilisi_score = lisi_graph_py(
        adata=adata_tmp,
        obs_key=batch_key,
        n_neighbors=k0,
        perplexity=None,
        subsample=subsample,
        n_cores=n_cores,
        verbose=verbose,
    )

    # iLISI: nbatches good, 1 bad
    ilisi = np.nanmedian(ilisi_score)

    if scale:
        nbatches = adata.obs[batch_key].nunique()
        ilisi = (ilisi - 1) / (nbatches - 1)

    return ilisi


def clisi_graph(
    adata,
    label_key,
    type_,
    use_rep="X_emb",
    batch_key=None,
    k0=90,
    subsample=None,
    scale=True,
    n_cores=1,
    verbose=False,
):
    """Cell-type LISI (cLISI) score

    Local Inverse Simpson’s Index metrics adapted from `Korsunsky et al. 2019`_ to run on all full
    feature, embedding and kNN integration outputs via shortest path-based distance computation on single-cell kNN
    graphs.
    By default, this function returns a value scaled between 0 and 1 instead of the original LISI range of 0 to the
    number of labels.

    .. _Korsunsky et al. 2019: https://doi.org/10.1038/s41592-019-0619-0

    :param adata: adata object to calculate on
    :param label_key: label column name in ``adata.obs``
    :param `type_`: type of data integration, one of 'knn', 'embed' or 'full'
    :param use_rep: embedding slot in ``.obsm``, only used for embedding input
    :param batch_key: deprecated, not used
    :param k0: number of nearest neighbors to compute lisi score
        Please note that the initial neighborhood size that is
        used to compute shortest paths is 15.
    :param subsample: Percentage of observations (integer between 0 and 100)
        to which lisi scoring should be subsampled
    :param scale: scale output values between 0 and 1 (True/False)
    :param n_cores: number of cores (i.e. CPUs or CPU cores to use for multiprocessing)
    :return: Median of cLISI scores per cell type labels

    This function can be applied to all integration output types and recomputes the kNN graph for feature and embedding
    output with specific parameters.
    Thus, no preprocessing is required, but the correct output type must be specified in ``type_``.

    **Examples**

    .. code-block:: python

        # feature output or unintegrated object
        scib.me.clisi_graph(adata, label_key="celltype", type="full")

        # embeding output
        scib.me.clisi_graph(adata, label_key="celltype", type="embed", use_rep="X_emb")

        # knn output
        scib.me.clisi_graph(adata, label_key="celltype", type="knn")

    """
    if batch_key is not None:
        warnings.warn("'batch_key' is deprecated and will be ignore", stacklevel=1)

    check_adata(adata)
    check_batch(label_key, adata.obs)

    adata_tmp = recompute_knn(adata, type_, use_rep)

    scores = lisi_graph_py(
        adata=adata_tmp,
        obs_key=label_key,
        n_neighbors=k0,
        perplexity=None,
        subsample=subsample,
        n_cores=n_cores,
        verbose=verbose,
    )

    # cLISI: 1 good, nlabs bad
    clisi = np.nanmedian(scores)

    if scale:
        nlabs = adata.obs[label_key].nunique()
        clisi = (nlabs - clisi) / (nlabs - 1)

    return clisi


def recompute_knn(adata, type_, use_rep="X_emb"):
    """Recompute neighbours"""
    if type_ == "embed":
        return sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep, copy=True)
    elif type_ == "full":
        if "X_pca" not in adata.obsm.keys():
            sc.pp.pca(adata, svd_solver="arpack")
        return sc.pp.neighbors(adata, n_neighbors=15, copy=True)
    else:
        # if knn - do not compute a new neighbourhood graph (it exists already)
        return adata.copy()


def lisi_graph_py(
    adata,
    obs_key,
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
):
    """
    Function to prepare call of compute_simpson_index
    Compute LISI score on shortes path based on kNN graph provided in the adata object.
    By default, perplexity is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """
    # use no more than the available cores
    n_cores = max(1, min(n_cores, mp.cpu_count()))

    if "neighbors" not in adata.uns:
        raise AttributeError(
            "Key 'neighbors' not found. Please make sure that a kNN graph has been computed"
        )
    elif verbose:
        print("using precomputed kNN graph")

    # get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")

    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    batch_labels = adata.obs[obs_key].cat.codes.values
    n_batches = adata.obs[obs_key].nunique()

    if perplexity is None or perplexity >= n_neighbors:
        # use LISI default
        perplexity = np.floor(n_neighbors / 3)

    # setup subsampling
    subset = 100  # default, no subsampling
    if subsample is not None:
        subset = subsample  # do not use subsampling
        if isinstance(subsample, int) is False:  # need to set as integer
            subset = int(subsample)

    # run LISI in python
    if verbose:
        print("Compute knn on shortest paths")

    # set connectivities to 3e-308 if they are lower than 3e-308 (because cpp can't handle double values smaller than that).
    connectivities = adata.obsp["connectivities"]  # csr matrix format
    large_enough = connectivities.data >= 3e-308
    if verbose:
        n_too_small = np.sum(large_enough is False)
        if n_too_small:
            print(
                f"{n_too_small} connectivities are smaller than 3e-308 and will be set to 3e-308"
            )
            print(connectivities.data[large_enough is False])
    connectivities.data[large_enough is False] = 3e-308

    # temporary file
    with tempfile.TemporaryDirectory(prefix="lisi_") as tmpdir:
        prefix = f"{tmpdir}/graph_lisi"
        mtx_file_path = prefix + "_input.mtx"
        mmwrite(mtx_file_path, connectivities, symmetry="general")

        # call knn-graph computation in Cpp
        root = pathlib.Path(scib.__file__).parent  # get current root directory
        cpp_file_path = (
            root / "knn_graph/knn_graph.o"
        )  # create POSIX path to file to execute compiled cpp-code
        # comment: POSIX path needs to be converted to string - done below with 'as_posix()'
        # create evenly split chunks if n_obs is divisible by n_chunks (doesn't really make sense on 2nd thought)
        args_int = [
            cpp_file_path.as_posix(),
            mtx_file_path,
            prefix,
            str(n_neighbors),
            str(n_cores),  # number of splits
            str(subset),
        ]

        try:
            if verbose:
                print(f'call {" ".join(args_int)}')
            subprocess.run(args_int)
        except RuntimeError as ex:
            print(f"Error computing LISI kNN graph {ex}\nSetting value to np.nan")
            return np.nan

        if verbose:
            print("LISI score estimation")

        if n_cores > 1:
            if verbose:
                print(f"{n_cores} processes started.")
            pool = mp.Pool(processes=n_cores)
            chunk_no = np.arange(0, n_cores)

            # create argument list for each worker
            results = pool.starmap(
                compute_simpson_index_graph,
                zip(
                    itertools.repeat(prefix),
                    itertools.repeat(batch_labels),
                    itertools.repeat(n_batches),
                    itertools.repeat(n_neighbors),
                    itertools.repeat(perplexity),
                    chunk_no,
                ),
            )
            pool.close()
            pool.join()

            simpson_estimate_batch = np.concatenate(results)

        else:
            simpson_estimate_batch = compute_simpson_index_graph(
                file_prefix=prefix,
                batch_labels=batch_labels,
                n_batches=n_batches,
                perplexity=perplexity,
                n_neighbors=n_neighbors,
            )

    return 1 / simpson_estimate_batch


# LISI core functions


def compute_simpson_index_graph(
    file_prefix=None,
    batch_labels=None,
    n_batches=None,
    n_neighbors=90,
    perplexity=30,
    chunk_no=0,
    tol=1e-5,
):
    """
    Simpson index of batch labels subset by group.

    :param file_prefix: file_path to pre-computed index and distance files
    :param batch_labels: a vector of length n_cells with batch info
    :param n_batches: number of unique batch labels
    :param n_neighbors: number of nearest neighbors
    :param perplexity: effective neighborhood size
    :param chunk_no: for parallelization, chunk id to evaluate
    :param tol: a tolerance for testing effective neighborhood size
    :returns: the simpson index for the neighborhood of each cell
    """
    index_file = file_prefix + "_indices_" + str(chunk_no) + ".txt"
    distance_file = file_prefix + "_distances_" + str(chunk_no) + ".txt"

    # check if the target file is not empty
    if os.stat(index_file).st_size == 0:
        print("File has no entries. Doing nothing.")
        lists = np.zeros(0)
        return lists

    # read distances and indices with nan value handling
    header = ["index"] + list(range(1, n_neighbors + 1))
    indices = pd.read_table(index_file, index_col=0, header=None, sep=",", names=header)
    indices = indices.T

    distances = pd.read_table(
        distance_file, index_col=0, header=None, sep=",", names=header
    )
    distances = distances.T

    # get cell ids
    chunk_ids = indices.columns.values.astype("int")

    # initialize
    logU = np.log(perplexity)
    simpson = np.zeros(len(chunk_ids))

    # loop over all cells in chunk
    for i, chunk_id in enumerate(chunk_ids):
        # get neighbors and distances
        # read line i from indices matrix
        get_col = indices[chunk_id]

        if get_col.isnull().sum() > 0:
            # not enough neighbors
            print(f"Chunk {chunk_id} does not have enough neighbors. Skipping...")
            simpson[i] = 1  # np.nan #set nan for testing
            continue

        knn_idx = get_col.astype("int") - 1  # get 0-based indexing

        # read line i from distances matrix
        D_act = distances[chunk_id].values.astype("float")

        # start lisi estimation
        beta = 1
        betamin = -np.inf
        betamax = np.inf

        H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0

        # first get neighbor probabilities
        while np.logical_and(np.abs(Hdiff) > tol, tries < 50):
            if Hdiff > 0:
                betamin = beta
                if betamax == np.inf:
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if betamin == -np.inf:
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2

            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1

        if H == 0:
            simpson[i] = -1
            continue
            # then compute Simpson's Index
        batch = batch_labels[knn_idx]
        B = convert_to_one_hot(batch, n_batches)
        sumP = np.matmul(P, B)  # sum P per batch
        simpson[i] = np.dot(sumP, sumP)  # sum squares

    return simpson


@njit
def Hbeta(D_row, beta):
    """
    Helper function for simpson index computation
    """
    D_row[np.isnan(D_row)] = 0
    P = np.exp(-D_row * beta)
    sumP = np.sum(P)
    if sumP == 0:
        H = 0
        P = np.zeros(len(D_row))
    else:
        H = np.log(sumP) + beta * np.sum(D_row * P) / sumP
        P /= sumP
    return H, P


@njit
def convert_to_one_hot(vector, num_classes=None):
    """
    Converts an input 1-D vector of integers into an output 2-D array of one-hot vectors,
    where an i'th input value of j will set a '1' in the i'th row, j'th column of the
    output array.

    Example:

    .. code-block:: python

        v = np.array((1, 0, 4))
        one_hot_v = convertToOneHot(v)
        print(one_hot_v)

    .. code-block::

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

    # result = np.zeros(shape=(len(vector), num_classes))
    # result[np.arange(len(vector)), vector] = 1
    # return result.astype(int)
    return np.eye(num_classes)[vector]
