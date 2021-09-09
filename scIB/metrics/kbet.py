import numpy as np
import scipy.sparse
import pandas as pd
import logging
import rpy2.robjects as ro
import rpy2.rinterface_lib.callbacks
import anndata2ri
import scanpy as sc

from scIB.utils import checkAdata, checkBatch
from .utils import diffusion_conn, diffusion_nn

rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)  # Ignore R warning messages


def kBET_single(
        matrix,
        batch,
        type_=None,
        k0=10,
        knn=None,
        subsample=0.5,
        heuristic=True,
        verbose=False
):
    """
    params:
        matrix: expression matrix (at the moment: a PCA matrix, so do.pca is set to FALSE
        batch: series or list of batch assignemnts
        subsample: fraction to be subsampled. No subsampling if `subsample=None`
    returns:
        kBET p-value
    """

    anndata2ri.activate()
    ro.r("library(kBET)")

    if verbose:
        print("importing expression matrix")
    ro.globalenv['data_mtrx'] = matrix
    ro.globalenv['batch'] = batch
    # print(matrix.shape)
    # print(len(batch))

    if verbose:
        print("kBET estimation")
    # k0 = len(batch) if len(batch) < 50 else 'NULL'

    ro.globalenv['knn_graph'] = knn
    ro.globalenv['k0'] = k0
    batch_estimate = ro.r(
        f"batch.estimate <- kBET(data_mtrx, batch, knn=knn_graph, k0=k0, plot=FALSE, do.pca=FALSE, heuristic=FALSE, adapt=FALSE, verbose={str(verbose).upper()})")

    anndata2ri.deactivate()
    try:
        ro.r("batch.estimate$average.pval")[0]
    except rpy2.rinterface_lib.embedded.RRuntimeError:
        return np.nan
    else:
        return ro.r("batch.estimate$average.pval")[0]


def kBET(adata, batch_key, label_key, embed='X_pca', type_=None,
         hvg=False, subsample=0.5, heuristic=False, verbose=False):
    """
    Compare the effect before and after integration
    params:
        matrix: matrix from adata to calculate on
    return:
        pd.DataFrame with kBET p-values per cluster for batch
    """

    checkAdata(adata)
    checkBatch(batch_key, adata.obs)
    checkBatch(label_key, adata.obs)
    # compute connectivities for non-knn type data integrations
    # and increase neighborhoods for knn type data integrations
    if type_ != 'knn':
        adata_tmp = sc.pp.neighbors(adata, n_neighbors=50, use_rep=embed, copy=True)
    else:
        # check if pre-computed neighbours are stored in input file
        adata_tmp = adata.copy()
        if 'diffusion_connectivities' not in adata.uns['neighbors']:
            if verbose:
                print(f"Compute: Diffusion neighbours.")
            adata_tmp = diffusion_conn(adata, min_k=50, copy=True)
        adata_tmp.obsp['connectivities'] = adata_tmp.uns['neighbors']['diffusion_connectivities']

    if verbose:
        print(f"batch: {batch_key}")

    # set upper bound for k0
    size_max = 2 ** 31 - 1

    kBET_scores = {'cluster': [], 'kBET': []}
    for clus in adata_tmp.obs[label_key].unique():

        adata_sub = adata_tmp[adata_tmp.obs[label_key] == clus, :].copy()
        # check if neighborhood size too small or only one batch in subset
        if np.logical_or(adata_sub.n_obs < 10,
                         len(adata_sub.obs[batch_key].cat.categories) == 1):
            print(f"{clus} consists of a single batch or is too small. Skip.")
            score = np.nan
        else:
            quarter_mean = np.floor(np.mean(adata_sub.obs[batch_key].value_counts()) / 4).astype('int')
            k0 = np.min([70, np.max([10, quarter_mean])])
            # check k0 for reasonability
            if (k0 * adata_sub.n_obs) >= size_max:
                k0 = np.floor(size_max / adata_sub.n_obs).astype('int')

            matrix = np.zeros(shape=(adata_sub.n_obs, k0 + 1))

            if verbose:
                print(f"Use {k0} nearest neighbors.")
            n_comp, labs = scipy.sparse.csgraph.connected_components(adata_sub.obsp['connectivities'],
                                                                     connection='strong')
            if n_comp > 1:
                # check the number of components where kBET can be computed upon
                comp_size = pd.value_counts(labs)
                # check which components are small
                comp_size_thresh = 3 * k0
                idx_nonan = np.flatnonzero(np.in1d(labs,
                                                   comp_size[comp_size >= comp_size_thresh].index))
                # check if 75% of all cells can be used for kBET run
                if len(idx_nonan) / len(labs) >= 0.75:
                    # create another subset of components, assume they are not visited in a diffusion process
                    adata_sub_sub = adata_sub[idx_nonan, :].copy()
                    nn_index_tmp = np.empty(shape=(adata_sub.n_obs, k0))
                    nn_index_tmp[:] = np.nan
                    nn_index_tmp[idx_nonan] = diffusion_nn(adata_sub_sub, k=k0).astype('float')
                    # need to check neighbors (k0 or k0-1) as input?
                    score = kBET_single(
                        matrix=matrix,
                        batch=adata_sub.obs[batch_key],
                        knn=nn_index_tmp + 1,  # nn_index in python is 0-based and 1-based in R
                        subsample=subsample,
                        verbose=verbose,
                        heuristic=False,
                        k0=k0,
                        type_=type_
                    )
                else:
                    # if there are too many too small connected components, set kBET score to 1
                    # (i.e. 100% rejection)
                    score = 1

            else:  # a single component to compute kBET on
                # need to check neighbors (k0 or k0-1) as input?
                nn_index_tmp = diffusion_nn(adata_sub, k=k0).astype('float')
                score = kBET_single(
                    matrix=matrix,
                    batch=adata_sub.obs[batch_key],
                    knn=nn_index_tmp + 1,  # nn_index in python is 0-based and 1-based in R
                    subsample=subsample,
                    verbose=verbose,
                    heuristic=False,
                    k0=k0,
                    type_=type_
                )

        kBET_scores['cluster'].append(clus)
        kBET_scores['kBET'].append(score)

    kBET_scores = pd.DataFrame.from_dict(kBET_scores)
    kBET_scores = kBET_scores.reset_index(drop=True)

    return kBET_scores
