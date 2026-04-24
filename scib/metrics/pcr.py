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
    Return either the raw difference in variance contribution before and after integration
    or a scaled score (``scale=True``).

    With ``scale=True``, the score is computed as
    ``(pcr_before - pcr_after) / pcr_before`` and clipped to 0 if it becomes negative
    (i.e. when variance contribution increases after integration).

    :param adata_pre: Anndata object before integration
    :param adata_post: Anndata object after integration
    :param covariate: Key for ``adata_post.obs`` column to regress against
    :param embed: Embedding to use for principal component analysis.
        If None, use the full expression matrix (``adata_post.X``), otherwise use the embedding
        provided in ``adata_post.obsm[embed]``.
    :param n_comps: Number of principal components to compute
    :param linreg_method: Method for linear regression passed to
        :func:`~scib.metrics.pc_regression` (``'sklearn'``, ``'numpy'``, or
        ``'sequential'``).
    :param recompute_pca: Whether to recompute PCA. If ``False`` (default), use existing PCA
        stored in ``adata.obsm["X_pca"]`` and ``adata.uns["pca"]["variance"]`` if available.
    :param scale: If ``True`` (default), scale score between 0 and 1, where 0 means no change
        in variance contribution after integration. If ``False``, return the raw difference
        of variance contributions (post minus pre).
    :param verbose: If ``True``, print progress information
    :param n_threads: Number of threads to pass to the selected regression backend
        (for example BLAS/OpenMP thread limits for ``'numpy'`` and estimator-level
        parallelism where supported).
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
    :param linreg_method: Regression backend. One of ``'sklearn'``, ``'numpy'``, or ``'sequential'``.

        - ``'sequential'`` calls :func:`~scib.metrics.pcr.linreg_sklearn`, the
            original implementation that fits one model per PC and is typically
            much slower.
        - ``'sklearn'`` calls :func:`~scib.metrics.pcr.linreg_multiple_sklearn`,
            a multi-output linear regression backend.
        - ``'numpy'`` calls :func:`~scib.metrics.pcr.linreg_multiple_np`, a
            vectorized numpy backend with a categorical one-way ANOVA shortcut.
    :param svd_solver:
    :param n_threads: Number of threads passed to the selected regression backend.
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

    linreg_methods = {
        "sklearn": linreg_multiple_sklearn,
        "numpy": linreg_multiple_np,
        "sequential": linreg_sklearn,
    }
    if linreg_method not in linreg_methods:
        raise ValueError(f"invalid linreg_method: {linreg_method}")
    linreg_method = linreg_methods[linreg_method]

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

    # fit linear model for n_comps PCs
    if verbose:
        print(f"Use {n_threads} threads for regression...")

    r2 = linreg_method(X_pca, covariate, n_jobs=n_threads)
    r2 = [np.maximum(0, x) for x in r2]

    Var = pca_var / sum(pca_var)  # * 100
    R2Var = sum(r2 * Var)  # / 100

    return R2Var


def _to_covariate_2d(covariate):
    """Return covariate as a 2D ndarray with shape (n_samples, n_covariates)."""
    if isinstance(covariate, pd.Series):
        return covariate.to_numpy(copy=False)[:, None]
    if isinstance(covariate, pd.DataFrame):
        return covariate.to_numpy(copy=False)

    covariate_arr = np.asarray(covariate)
    if covariate_arr.ndim == 1:
        covariate_arr = covariate_arr[:, None]
    return covariate_arr


def _encode_covariate(covariate, as_sparse=False):
    """Encode covariates for regression backends.

    :param as_sparse: Return a scipy sparse matrix for numeric and categorical
        covariates when ``True``. Otherwise return a dense numpy array.
    """
    covariate_arr = _to_covariate_2d(covariate)

    # Use a single dtype check for all input containers (Series, DataFrame, ndarray).
    is_numeric = (
        pd.DataFrame(covariate_arr).dtypes.map(pd.api.types.is_numeric_dtype).all()
    )

    if is_numeric:
        return sparse.csr_matrix(covariate_arr) if as_sparse else covariate_arr

    if as_sparse:
        from sklearn.preprocessing import OneHotEncoder

        try:
            # Drop one reference level since downstream linear models fit an intercept.
            # This avoids a rank-deficient design matrix (dummy-variable trap) and
            # tends to scale better as category count increases.
            encoder = OneHotEncoder(
                sparse_output=True, handle_unknown="ignore", drop="first"
            )
        except TypeError:
            encoder = OneHotEncoder(sparse=True, handle_unknown="ignore", drop="first")
        return encoder.fit_transform(covariate_arr)
    return pd.get_dummies(pd.DataFrame(covariate_arr), drop_first=True).to_numpy()


def linreg_sklearn(X_pca, covariate, n_jobs=None):
    """Sequential sklearn regression backend for PCR.

    Fits one ``LinearRegression`` model per principal component and returns one
    :math:`R^2` value per component. This is the original implementation used by
    ``linreg_method='sequential'`` in :func:`~scib.metrics.pc_regression`.
    """
    from concurrent.futures import ThreadPoolExecutor

    from sklearn.linear_model import LinearRegression

    if n_jobs is None:
        n_jobs = 1

    X = _encode_covariate(covariate, as_sparse=False)

    def compute_r2(pc):
        model = LinearRegression()
        model.fit(X, pc)
        return model.score(X, pc)

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
    """Multi-output sklearn regression backend for PCR.

    Encodes the covariate once and fits a single multi-output
    ``LinearRegression`` model against all PCs.
    """
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score

    X = _encode_covariate(covariate, as_sparse=True)
    Y = X_pca

    model = LinearRegression(fit_intercept=True, n_jobs=n_jobs)
    model.fit(X, Y)
    Y_pred = model.predict(X)
    return [r2_score(Y[:, i], Y_pred[:, i]) for i in range(Y.shape[1])]


def linreg_multiple_np(X_pca, covariate, n_jobs=None):
    """Compute per-PC :math:`R^2` with a dense numpy regression backend.

    Execution has two paths:
        1) If exactly one non-numeric covariate is provided, use a fast
           one-way ANOVA formulation based on between-group variance. This
           is automatically used for one non-numeric covariate.
        2) Otherwise, run dense least-squares regression using an intercept
           and pseudo-inverse.

    One-way ANOVA path:
        For a single categorical covariate, compute per-PC
        ``R2 = SS_between / SS_total`` where
        ``SS_between = sum_g n_g * (mu_g - mu)^2`` and
        ``SS_total = sum_i (y_i - mu)^2``.
        This is algebraically equivalent to linear regression with categorical
        indicators and an intercept, but avoids explicit design-matrix solves.

    Thread control:
        If ``n_jobs`` is provided and ``threadpoolctl`` is installed,
        BLAS/LAPACK thread pools are limited for the full function scope
        (including shortcut checks and dense regression). This keeps runtime
        behavior reproducible across environments and avoids BLAS oversubscription.

    Used by ``linreg_method='numpy'`` in :func:`~scib.metrics.pc_regression`.

    :param n_jobs: Preferred number of BLAS/LAPACK threads for numpy linalg calls.
        If ``None``, keep the runtime default. If ``threadpoolctl`` is unavailable,
        this hint is ignored.
    """

    def _fit_numpy_regression(X, Y):
        X = np.hstack([np.ones((X.shape[0], 1)), X])  # fit with intercept
        beta = np.linalg.pinv(X.T @ X) @ X.T @ Y
        Y_pred = X @ beta

        ss_total = np.sum((Y - Y.mean(axis=0)) ** 2, axis=0)
        ss_residual = np.sum((Y - Y_pred) ** 2, axis=0)

        return 1 - ss_residual / ss_total

    def _fit_one_way_anova(covariate, Y):
        """Compute per-target R2 from one-way ANOVA decomposition."""
        if isinstance(covariate, pd.Series) and isinstance(
            covariate.dtype, pd.CategoricalDtype
        ):
            # Reuse existing categorical encoding directly when available.
            codes = covariate.cat.codes.to_numpy(copy=False)
            n_categories = len(covariate.cat.categories)
        else:
            categories = pd.Categorical(covariate)
            n_categories = len(categories.categories)
            codes = categories.codes

        # Pandas categorical missing values are encoded as -1 in cat.codes.
        if np.any(codes < 0):
            raise ValueError(
                "linreg_method='numpy' does not support missing covariate values"
            )

        if n_categories <= 1:
            return np.zeros(Y.shape[1], dtype=Y.dtype)

        mean_global = Y.mean(axis=0, keepdims=True)
        counts = np.bincount(codes, minlength=n_categories).astype(Y.dtype)

        # Sparse one_hot (n_categories × n_cells) @ Y  →  per-group sums for all PCs at once.
        # Equivalent to a one-hot encoding but stored as sparse to avoid the n×k dense matrix.
        one_hot = sparse.csr_matrix(
            (
                np.ones(len(codes), dtype=Y.dtype),  # one non-zero per cell
                (
                    codes,  # rows=group
                    np.arange(len(codes), dtype=np.int32),  # cols=cell
                ),
            ),
            shape=(n_categories, Y.shape[0]),
        )
        sums = one_hot @ Y  # shape: (n_categories, n_pcs)

        # Divide in-place: sums becomes per-group means, avoiding a second allocation.
        sums /= counts[:, None]
        ss_between = np.sum(counts[:, None] * (sums - mean_global) ** 2, axis=0)
        # Computational formula avoids a full (n_cells × n_pcs) copy:
        #   Σ(y - μ)² = Σy² - n·μ²
        ss_total = np.einsum("ij,ij->j", Y, Y) - Y.shape[0] * mean_global.ravel() ** 2
        return np.divide(
            ss_between,
            ss_total,
            out=np.zeros_like(ss_between),
            where=ss_total > 0,
        )

    from contextlib import nullcontext

    threadpool_context = nullcontext()
    if n_jobs is not None:
        try:
            from threadpoolctl import threadpool_limits

            threadpool_context = threadpool_limits(limits=n_jobs)
        except ImportError:
            pass

    with threadpool_context:
        # a Series[category] with numeric-valued categories must still use ANOVA.
        if isinstance(covariate, pd.Series) and not pd.api.types.is_numeric_dtype(
            covariate
        ):
            return _fit_one_way_anova(covariate, X_pca)

        covariate = _to_covariate_2d(covariate)

        # Fallback for non-numeric ndarray/list inputs (e.g. object dtype strings)
        if covariate.shape[1] == 1 and not np.issubdtype(covariate.dtype, np.number):
            return _fit_one_way_anova(pd.Series(covariate[:, 0]), X_pca)

        X = _encode_covariate(covariate, as_sparse=False)
        return _fit_numpy_regression(X, X_pca)
