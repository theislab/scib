import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
from sklearn.linear_model import LinearRegression

from scIB.utils import checkAdata, checkBatch


def pcr_comparison(
        adata_pre,
        adata_post,
        covariate,
        embed=None,
        n_comps=50,
        scale=True,
        verbose=False
):
    """
    Compare the effect before and after integration
    Return either the difference of variance contribution before and after integration
    or a score between 0 and 1 (`scaled=True`) with 0 if the variance contribution hasn't
    changed. The larger the score, the more different the variance contributions are before
    and after integration.
    params:
        adata_pre: uncorrected adata
        adata_post: integrated adata
        embed   : if `embed=None`, use the full expression matrix (`adata.X`), otherwise
                  use the embedding provided in `adata_post.obsm[embed]`
        scale: if True, return scaled score
    return:
        difference of R2Var value of PCR
    """

    if embed == 'X_pca':
        embed = None

    pcr_before = pcr(adata_pre, covariate=covariate, recompute_pca=True,
                     n_comps=n_comps, verbose=verbose)
    pcr_after = pcr(adata_post, covariate=covariate, embed=embed, recompute_pca=True,
                    n_comps=n_comps, verbose=verbose)

    if scale:
        score = (pcr_before - pcr_after) / pcr_before
        if score < 0:
            print("Variance contribution increased after integration!")
            print("Setting PCR comparison score to 0.")
            score = 0
        return score
    else:
        return pcr_after - pcr_before


def pcr(
        adata,
        covariate,
        embed=None,
        n_comps=50,
        recompute_pca=True,
        verbose=False
):
    """
    PCR for Adata object
    Checks whether to
        + compute PCA on embedding or expression data (set `embed` to name of embedding matrix e.g. `embed='X_emb'`)
        + use existing PCA (only if PCA entry exists)
        + recompute PCA on expression matrix (default)
    params:
        adata: Anndata object
        embed   : if `embed=None`, use the full expression matrix (`adata.X`), otherwise
                  use the embedding provided in `adata_post.obsm[embed]`
        n_comps: number of PCs if PCA should be computed
        covariate: key for adata.obs column to regress against
    return:
        R2Var of PCR
    """

    checkAdata(adata)
    checkBatch(covariate, adata.obs)

    if verbose:
        print(f"covariate: {covariate}")
    covariate_values = adata.obs[covariate]

    # use embedding for PCA
    if (embed is not None) and (embed in adata.obsm):
        if verbose:
            print(f"compute PCR on embedding n_comps: {n_comps}")
        return pc_regression(adata.obsm[embed], covariate_values, n_comps=n_comps)

    # use existing PCA computation
    elif (recompute_pca == False) and ('X_pca' in adata.obsm) and ('pca' in adata.uns):
        if verbose:
            print("using existing PCA")
        return pc_regression(adata.obsm['X_pca'], covariate_values, pca_var=adata.uns['pca']['variance'])

    # recompute PCA
    else:
        if verbose:
            print(f"compute PCA n_comps: {n_comps}")
        return pc_regression(adata.X, covariate_values, n_comps=n_comps)


def pc_regression(
        data,
        covariate,
        pca_var=None,
        n_comps=50,
        svd_solver='arpack',
        verbose=False
):
    """
    params:
        data: expression or PCA matrix. Will be assumed to be PCA values, if pca_sd is given
        covariate: series or list of batch assignments
        n_comps: number of PCA components for computing PCA, only when pca_sd is not given. If no pca_sd is given and n_comps=None, comute PCA and don't reduce data
        pca_var: iterable of variances for `n_comps` components. If `pca_sd` is not `None`, it is assumed that the matrix contains PCA values, else PCA is computed
    PCA is only computed, if variance contribution is not given (pca_sd).
    """

    if isinstance(data, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix)):
        matrix = data
    else:
        raise TypeError(f'invalid type: {data.__class__} is not a numpy array or sparse matrix')

    # perform PCA if no variance contributions are given
    if pca_var is None:

        if n_comps is None or n_comps > min(matrix.shape):
            n_comps = min(matrix.shape)

        if n_comps == min(matrix.shape):
            svd_solver = 'full'

        if verbose:
            print("compute PCA")
        pca = sc.tl.pca(matrix, n_comps=n_comps, use_highly_variable=False,
                        return_info=True, svd_solver=svd_solver, copy=True)
        X_pca = pca[0].copy()
        pca_var = pca[3].copy()
        del pca
    else:
        X_pca = matrix
        n_comps = matrix.shape[1]

    ## PC Regression
    if verbose:
        print("fit regression on PCs")

    # handle categorical values
    if pd.api.types.is_numeric_dtype(covariate):
        covariate = np.array(covariate).reshape(-1, 1)
    else:
        if verbose:
            print("one-hot encode categorical values")
        covariate = pd.get_dummies(covariate)

    # fit linear model for n_comps PCs
    r2 = []
    for i in range(n_comps):
        pc = X_pca[:, [i]]
        lm = LinearRegression()
        lm.fit(covariate, pc)
        r2_score = np.maximum(0, lm.score(covariate, pc))
        r2.append(r2_score)

    Var = pca_var / sum(pca_var) * 100
    R2Var = sum(r2 * Var) / 100

    return R2Var
