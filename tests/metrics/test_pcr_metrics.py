import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.sparse import csr_matrix

import scib
from scib.metrics.pcr import (
    _encode_covariate,
    linreg_multiple_np,
    linreg_multiple_sklearn,
)
from tests.common import LOGGER, add_embed, assert_near_exact


@pytest.mark.parametrize("sparse", [False, True])
def test_pc_regression(adata, sparse):
    if sparse:
        adata.X = csr_matrix(adata.X)
    score = scib.me.pc_regression(
        adata.X,
        covariate=adata.obs["batch"],
        n_comps=adata.n_vars,
    )
    LOGGER.info(score)
    assert_near_exact(score, 0, diff=1e-3)


@pytest.mark.parametrize("linreg_method", ["sequential", "numpy", "sklearn"])
def test_pc_regression_numeric_covariate(adata_pca, linreg_method):
    X_pca = adata_pca.obsm["X_pca"]
    pca_var = np.asarray(adata_pca.uns["pca"]["variance"])

    # Use the first PC as covariate: PCs are orthogonal, so R²_0 ≈ 1 and
    # R²_i ≈ 0 for i > 0, making the expected score equal to var[0] / sum(var).
    covariate = pd.Series(X_pca[:, 0])

    score = scib.me.pc_regression(
        X_pca,
        covariate=covariate,
        pca_var=pca_var,
        linreg_method=linreg_method,
    )
    score_legacy = scib.me.pc_regression(
        X_pca,
        covariate=covariate,
        pca_var=pca_var,
        linreg_method="sequential",
    )

    assert score_legacy > pca_var[0] / pca_var.sum() * 0.9
    assert_near_exact(score, score_legacy, diff=1e-3)


@pytest.mark.parametrize("linreg_method", ["sequential", "numpy", "sklearn"])
def test_pcr_correctness_all_implementations(adata_pca, linreg_method):
    score = scib.me.pcr(
        adata_pca,
        covariate="celltype",
        linreg_method=linreg_method,
        verbose=False,
    )
    LOGGER.info(score)
    assert_near_exact(score, 0.33431789694498437, diff=1e-3)


def test_pcr_comparison_batch(adata):
    score = scib.me.pcr_comparison(
        adata, adata, covariate="batch", n_comps=50, scale=True
    )
    LOGGER.info(f"no PCA precomputed: {score}")
    assert_near_exact(score, 0, diff=1e-6)


def test_pcr_comparison_batch_precomputed(adata_pca):
    score = scib.me.pcr_comparison(adata_pca, adata_pca, covariate="batch", scale=True)
    LOGGER.info(f"precomputed PCA: {score}")
    assert_near_exact(score, 0, diff=1e-6)


def test_pcr_comparison_batch_embedding(adata):
    # use different embedding
    score = scib.me.pcr_comparison(
        adata_pre=adata,
        adata_post=add_embed(adata, type_="full"),
        covariate="batch",
        embed="X_emb",
        n_comps=50,
        scale=True,
    )
    LOGGER.info(f"using X_emb: {score}")
    assert_near_exact(score, 0, diff=1e-6)


@pytest.mark.parametrize(
    "covariate,expected_shape",
    [
        pytest.param(pd.Series([0.0, 1.0, 2.0, 3.0], name="num"), (4, 1), id="numeric"),
        pytest.param(
            pd.Series(["a", "b", "a", "c"], name="cat").astype("category"),
            (4, 2),
            id="categorical",
        ),
    ],
)
def test_covariate_encoding_helpers_dense_and_sparse(covariate, expected_shape):
    dense_cov = _encode_covariate(covariate, as_sparse=False)
    sparse_cov = _encode_covariate(covariate, as_sparse=True)

    assert isinstance(dense_cov, np.ndarray)
    assert dense_cov.shape == expected_shape

    assert sparse.issparse(sparse_cov)
    assert sparse_cov.shape == expected_shape


def test_linreg_multiple_sklearn_uses_sparse_one_hot(monkeypatch):
    from sklearn.linear_model import LinearRegression as SklearnLinearRegression

    captured = {}

    class RecordingLinearRegression(SklearnLinearRegression):
        def fit(self, X, y, sample_weight=None):
            captured["is_sparse"] = sparse.issparse(X)
            captured["shape"] = X.shape
            return super().fit(X, y, sample_weight=sample_weight)

    # Monkeypatching is needed because linreg_multiple_sklearn only returns R2 scores;
    # intercepting fit() lets us assert the internal contract that sklearn receives
    # sparse one-hot encoded covariates instead of a dense matrix.
    monkeypatch.setattr(
        "sklearn.linear_model.LinearRegression", RecordingLinearRegression
    )

    rng = np.random.default_rng(0)
    X_pca = rng.standard_normal((120, 8))
    covariate = pd.Series([f"batch_{i % 40}" for i in range(120)], dtype="category")

    scores = linreg_multiple_sklearn(X_pca, covariate, n_jobs=1)

    assert len(scores) == X_pca.shape[1]
    assert captured["is_sparse"] is True
    assert captured["shape"] == (120, 39)


@pytest.mark.slow
@pytest.mark.parametrize(
    "impl,runner",
    [
        pytest.param("numpy", linreg_multiple_np, id="numpy"),
        # pytest.param("sklearn", linreg_multiple_sklearn, id="sklearn"),
    ],
)
@pytest.mark.benchmark(group="pcr_large_impl_compare")
def test_large_matrix_runtime_benchmark(benchmark, impl, runner):
    n_cells = 500_000
    n_pcs = 50
    n_categories = 3000
    rng = np.random.default_rng(0)

    # Build data once; do not include setup cost in measured runtime.
    X_pca = rng.standard_normal((n_cells, n_pcs))
    covariate = pd.Series(
        [f"batch_{i % n_categories}" for i in range(n_cells)],
        dtype="category",
    )

    benchmark.extra_info["implementation"] = impl
    benchmark.extra_info["n_cells"] = n_cells
    benchmark.extra_info["n_pcs"] = n_pcs
    benchmark.extra_info["n_categories"] = n_categories

    def run():
        return runner(X_pca, covariate, n_jobs=10)

    scores = benchmark.pedantic(
        run,
        rounds=100,
        iterations=1,
        warmup_rounds=2,
    )
    assert len(scores) == n_pcs


@pytest.mark.slow
@pytest.mark.parametrize(
    "linreg_method",
    [
        pytest.param("sequential", id="sequential"),
        pytest.param("numpy", id="numpy"),
        pytest.param("sklearn", id="sklearn"),
    ],
)
@pytest.mark.benchmark(group="pc_regression_numeric_impl_compare")
def test_numeric_pc_regression_runtime_benchmark(benchmark, linreg_method):
    n_cells = 100_000
    n_pcs = 50
    rng = np.random.default_rng(0)

    # Build data once; do not include setup cost in measured runtime.
    X_pca = rng.standard_normal((n_cells, n_pcs))
    covariate = pd.Series(rng.standard_normal(n_cells), name="numeric_covariate")
    pca_var = np.asarray(rng.uniform(0.5, 2.0, size=n_pcs), dtype=np.float32)

    benchmark.extra_info["linreg_method"] = linreg_method
    benchmark.extra_info["n_cells"] = n_cells
    benchmark.extra_info["n_pcs"] = n_pcs

    def run():
        return scib.me.pc_regression(
            X_pca,
            covariate=covariate,
            pca_var=pca_var,
            linreg_method=linreg_method,
            n_threads=10,
        )

    score = benchmark.pedantic(
        run,
        rounds=100,
        iterations=1,
        warmup_rounds=2,
    )
    assert np.isfinite(score)
