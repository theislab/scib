import matplotlib
import scanpy as sc

import scib

matplotlib.use("Agg")


def test_summarize_counts_creates_obs_columns():
    adata = sc.datasets.blobs()

    # should not raise and should add obs columns
    scib.pp.summarize_counts(adata)

    for col in ("n_counts", "log_counts", "n_genes"):
        assert col in adata.obs

    # counts should be positive
    assert (adata.obs["n_counts"] > 0).sum() > 0


def test_plot_functions_run_without_error():
    adata = sc.datasets.blobs()

    # prepare summary used by plotting
    scib.pp.summarize_counts(adata)

    # plotting functions should run without raising
    scib.pp.plot_qc(adata)
    scib.pp.plot_scatter(adata)
    scib.pp.plot_count_filter(adata)
