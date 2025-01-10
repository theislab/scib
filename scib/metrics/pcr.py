import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from ..utils import check_adata, check_batch


def pcr_comparison(
    adata_pre,
    adata_post,
    covariate,
    embed=None,
    n_comps=50,
    linreg_method="numpy",
    recompute_pca=False,
    scale=True,
    verbose=False,
    n_threads=1,
):
    """Principal component regression score

    Compare the explained variance before and after integration using :func:`~scib.metrics.pc_regression`.
    Return either the difference of variance contribution before and after integration
    or a score between 0 and 1 (``scaled=True``) with 0 if the variance contribution hasn't
    changed. The larger the score, the more different the variance contributions are before
    and after integration.

    :param adata_pre: anndata object before integration
    :param adata_post: anndata object after integration
    :param covariate: Key for ``adata_post.obs`` column to regress against
    :param embed: Embedding to use for principal component analysis.
        If None, use the full expression matrix (``adata_post.X``), otherwise use the embedding
        provided in ``adata_post.obsm[embed]``.
    :param n_comps: Number of principal components to compute
    :param recompute_pca: Whether to recompute PCA with default settings
    :param linreg_method: Method for linear regression, either 'sklearn' or 'numpy'
    :param scale: If True, scale score between 0 and 1 (default)
    :param verbose:
    :return:
        Difference of variance contribution of PCR (scaled between 0 and 1 by default)

    The function can be computed on full corrected feature spaces and latent embeddings for both integrated and
    unintegrated ``anndata.Anndata`` objects.
    No preprocessing is needed, as the function will perform PCA directly on the feature or embedding space.

    **Examples**

    .. code-block:: python

        # full feature output
        scib.me.pcr_comparison(adata_unintegrated, adata, covariate="batch")

        # embedding output
        scib.me.pcr_comparison(adata_unintegrated, adata, covariate="batch", embed="X_emb")

    """

    if embed == "X_pca":
        embed = None

    pcr_before = pcr(
        adata_pre,
        covariate=covariate,
        recompute_pca=recompute_pca,
        n_comps=n_comps,
        linreg_method=linreg_method,
        n_threads=n_threads,
        verbose=verbose,
    )

    pcr_after = pcr(
        adata_post,
        covariate=covariate,
        embed=embed,
        recompute_pca=recompute_pca,
        n_comps=n_comps,
        linreg_method=linreg_method,
        n_threads=n_threads,
        verbose=verbose,
    )

    if scale:
        score = (pcr_before - pcr_after) / pcr_before
        if score < 0:
            print(
                "Variance contribution increased after integration!\n"
                "Setting PCR comparison score to 0."
            )
            score = 0
        return score
    else:
        return pcr_after - pcr_before


def pcr(
    adata,
    covariate,
    embed=None,
    n_comps=50,
    recompute_pca=False,
    linreg_method="sklearn",
    verbose=False,
    n_threads=1,
):
    """Principal component regression for anndata object

    Wraps :func:`~scib.metrics.pc_regression` while checking whether to:

        + compute PCA on embedding or expression data (set ``embed`` to name of embedding matrix e.g. ``embed='X_emb'``)
        + use existing PCA (only if PCA entry exists)
        + recompute PCA on expression matrix (default)

    :param adata: Anndata object
    :param covariate: Key for ``adata.obs`` column to regress against
    :param embed: Embedding to use for principal component analysis.
        If None, use the full expression matrix (``adata.X``), otherwise use the embedding
        provided in ``adata.obsm[embed]``.
    :param n_comps: Number of PCs, if PCA is recomputed. The PCA will be recomputed if neither PCA loadings nor the
        principle components can be found.
    :param recompute_pca: whether to recompute a PCA on the
    :return:
        Variance contribution of regression

    The function can be computed on full corrected feature spaces and latent embeddings.
    No preprocessing is needed, as the function can perform PCA if ``recompute_pca=True``.
    Alternatively, you can also provide precomputed PCA, if the principle components are saved under ``.obsm["X_pca"]``
    and the PC loadings are saved in ``.uns["pca"]["variance"]``.

    **Examples**

    .. code-block:: python

        # full feature output
        scib.me.pcr(adata, covariate="batch", recompute_pca=True)

        # embedding output
        scib.me.pcr(adata, covariate="batch", embed="X_emb")
    """

    check_adata(adata)
    check_batch(covariate, adata.obs)

    if verbose:
        print(f"covariate: {covariate}")
    covariate_values = adata.obs[covariate]

    # use embedding for PCA
    if embed is not None:
        assert embed in adata.obsm
        if verbose:
            print(f"Compute PCR on embedding n_comps: {n_comps}")
        return pc_regression(
            adata.obsm[embed],
            covariate_values,
            n_comps=n_comps,
            linreg_method=linreg_method,
            n_threads=n_threads,
        )

    # use existing PCA computation
    elif (recompute_pca is False) and ("X_pca" in adata.obsm) and ("pca" in adata.uns):
        if verbose:
            print("using existing PCA")
        return pc_regression(
            adata.obsm["X_pca"],
            covariate_values,
            pca_var=adata.uns["pca"]["variance"],
            linreg_method=linreg_method,
            n_threads=n_threads,
        )

    # recompute PCA
    else:
        if verbose:
            print(f"compute PCA n_comps: {n_comps}")
        return pc_regression(
            adata.X,
            covariate_values,
            n_comps=n_comps,
            linreg_method=linreg_method,
            n_threads=n_threads,
        )


def pc_regression(
    data,
    covariate,
    pca_var=None,
    n_comps=50,
    svd_solver="arpack",
    linreg_method="sklearn",
    verbose=False,
    n_threads=1,
):
    """Principal component regression

    Compute the overall variance contribution given a covariate according to the following formula:

    .. math::

        Var(C|B) = \\sum^G_{i=1} Var(C|PC_i) \\cdot R^2(PC_i|B)

    for :math:`G` principal components (:math:`PC_i`), where :math:`Var(C|PC_i)` is the variance of the data matrix
    :math:`C` explained by the i-th principal component, and :math:`R^2(PC_i|B)` is the :math:`R^2` of the i-th
    principal component regressed against a covariate :math:`B`.


    :param data: Expression or PC matrix. Assumed to be PC, if pca_sd is given.
    :param covariate: series or list of batch assignments
    :param n_comps: number of PCA components for computing PCA, only when pca_sd is not given.
        If no pca_sd is not defined and n_comps=None, compute PCA and don't reduce data
    :param pca_var: Iterable of variances for ``n_comps`` components.
        If ``pca_sd`` is not ``None``, it is assumed that the matrix contains PC,
        otherwise PCA is computed on ``data``.
    :param svd_solver:
    :param verbose:
    :return:
        Variance contribution of regression
    """
    if isinstance(data, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix)):
        matrix = data
    else:
        raise TypeError(
            f"invalid type: {data.__class__} is not a numpy array or sparse matrix"
        )

    if linreg_method == "sklearn":
        linreg_method = linreg_multiple_sklearn
    elif linreg_method == "numpy":
        linreg_method = linreg_multiple_np
    elif linreg_method == "sequential":
        linreg_method = linreg_sklearn
    else:
        raise ValueError(f"invalid linreg_method: {linreg_method}")

    # perform PCA if no variance contributions are given
    if pca_var is None:

        if n_comps is None or n_comps > min(matrix.shape):
            n_comps = min(matrix.shape)

        if n_comps == min(matrix.shape):
            svd_solver = "full"
            # convert to dense bc 'full' is not available for sparse matrices
            if sparse.issparse(matrix):
                matrix = matrix.toarray()

        if verbose:
            print("compute PCA...")
        X_pca, _, _, pca_var = sc.tl.pca(
            matrix,
            n_comps=n_comps,
            use_highly_variable=False,
            return_info=True,
            svd_solver=svd_solver,
            copy=True,
        )
    else:
        X_pca = matrix
        n_comps = matrix.shape[1]

    # PC Regression
    if verbose:
        print("fit regression on PCs")

    # handle categorical values
    if pd.api.types.is_numeric_dtype(covariate):
        covariate = np.array(covariate).reshape(-1, 1)
    else:
        if verbose:
            print("one-hot encode categorical values")
        covariate = pd.get_dummies(covariate).to_numpy()

    # fit linear model for n_comps PCs
    if verbose:
        print(f"Use {n_threads} threads for regression...")

    r2 = linreg_method(X_pca, covariate, n_jobs=n_threads)
    r2 = [np.maximum(0, x) for x in r2]

    Var = pca_var / sum(pca_var)  # * 100
    R2Var = sum(r2 * Var)  # / 100

    return R2Var


def linreg_sklearn(X_pca, covariate, n_jobs=None):
    from concurrent.futures import ThreadPoolExecutor

    from sklearn.linear_model import LinearRegression

    if n_jobs is None:
        n_jobs = 1

    def compute_r2(pc):
        model = LinearRegression()
        model.fit(covariate, pc)
        return model.score(covariate, pc)

    # Parallelizing with ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        r2 = list(executor.map(compute_r2, X_pca.T))
    return r2


def linreg_np(X, y):
    _, residuals, _, _ = np.linalg.lstsq(X, y, rcond=None)
    if residuals.size == 0:
        residuals = np.array([0])
    rss = residuals[0] if len(residuals) > 0 else 0
    tss = np.sum((y - y.mean()) ** 2)
    r2_score = 1 - (rss / tss)
    return np.maximum(0, r2_score)


def linreg_multiple_sklearn(X_pca, covariate, n_jobs=None):
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score

    X = covariate
    Y = X_pca

    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
        # X = X.toarray()

    model = LinearRegression(fit_intercept=True, n_jobs=n_jobs)
    model.fit(X, Y)
    Y_pred = model.predict(X)
    return [r2_score(Y[:, i], Y_pred[:, i]) for i in range(Y.shape[1])]


def linreg_multiple_np(X_pca, covariate, n_jobs=None):
    """
    :param n_jobs: Number of threads ignored
    """
    X = covariate
    Y = X_pca

    X = np.hstack([np.ones((X.shape[0], 1)), X])  # fit with intercept
    beta = np.linalg.pinv(X.T @ X) @ X.T @ Y
    Y_pred = X @ beta

    # Compute r2 scores
    ss_total = np.sum((Y - Y.mean(axis=0)) ** 2, axis=0)
    ss_residual = np.sum((Y - Y_pred) ** 2, axis=0)

    return 1 - ss_residual / ss_total
