import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

from scib.preprocessing import hvg_batch
from scib.utils import check_sanity


def adata_scale_on_copy(adata):
    adata_scale = adata.copy()
    # PP each sample: scaling, pca, neighbours
    sc.pp.scale(adata_scale, max_value=10)
    sc.pp.pca(
        adata_scale,
        n_comps=15,
        use_highly_variable=True,
        svd_solver="arpack",
    )
    sc.pp.neighbors(adata_scale, n_pcs=15)
    return adata_scale


def morans_i_score(adata, hvgs):
    """
    Inplace function for computing Moran's I on non-constant HVGS
    and adding them to .var
    :param data: adata, result is added to var
    :param hvgs: Genes for which to compute Moran's I
    """
    # Make sure that genes used for computation are non-constant
    hvg_data = adata[:, hvgs].X
    hvg_data = hvg_data.A if issparse(hvg_data) else hvg_data
    std_hvgs = np.asarray(hvg_data.std(axis=0)).squeeze()
    hvgs_used = np.array(hvgs)[std_hvgs > 1e-8]

    if len(hvgs) > len(hvgs_used):
        print(f"Using not-constant subset of hvgs: {hvgs_used}/{hvgs}.")

    # Compute Moran's I
    morans_i = sc.metrics.morans_i(adata[:, hvgs_used])
    morans_i = pd.Series(morans_i, index=hvgs_used)

    # Check that values are in correct range
    if not morans_i.between(-1, 1).all():
        raise ValueError("Computed Moran's I values are not in [-1,1].")

    # We consider only values > 0 as biologically meaningful
    morans_i = np.maximum(morans_i, 0)
    adata.var["morans_i"] = morans_i


def morans_i(
    adata_pre,
    adata_post,
    batch_key,
    n_hvg=100,
    embed="X_pca",
    rescale=True,
    compare_pre=False,
    hvgs=None,
):
    """
    Compute Moran's I on HVGs (defined across batches) on integrated data as a measure of bio conservation. The metric can be computed either as

        1. mean difference between Morans's I on embedding of integrated data and max Moran's I of individual non-integrated batches or
        2. mean Morans's I on embedding of integrated data.

    :param adata_pre: Non-integrated data.
        Embedding is computed on scaled (m=0, s=1) X (expression)
        per batch (recomputing HVGs, scaling, pca, neighbours).
        For Moran's I unscaled X is used.
    :param adata_post: Integrate data. Should have precomputed embedding and expression in X.
    :param batch_key: Batch obs col in adata_pre.
    :param n_hvg: How many top HVGs to use for Moran's I computation.
    :param embed: Embedding to use for computation of neighbours in adata_post
    :param rescale: Whether to rescale result to [0,1] so that 1 is better bio conservation
    :param compare_pre: Whether to compute metric under scenario A instead of B.
    :param hvgs: Custom list of genes to use for Moran's I computation (instead of hvgs
        computed across batches).
    :return: The resulting metric.
    """
    # Get integrated connectivities for Moran's I
    check_sanity(adata_pre, batch_key, hvgs)
    adata_post = adata_post.copy()
    sc.pp.neighbors(adata_post, use_rep=embed)
    # Prepare pre data
    check_sanity(adata_post, batch_key, hvgs)
    adata_pre = adata_pre.copy()
    adata_pre.obs[batch_key] = adata_pre.obs[batch_key].astype("category")

    # Get HVGs across samples
    if hvgs is None:
        hvgs = hvg_batch(adata_pre, batch_key, n_hvg, flavor="cell_ranger")
    else:
        assert adata_pre.var_names.isin(hvgs).sum() == len(hvgs)

    # Moran's I on integrated data
    morans_i_score(adata_post, hvgs)
    score = adata_post.var["morans_i"].mean()

    if compare_pre:
        # Get connectivities per batch
        adatas_batches = {}
        for batch in adata_pre.obs[batch_key].unique():
            # Subset only to cells from 1 sample
            adata_batch = adata_pre[adata_pre.obs[batch_key] == batch, :].copy()

            # Compute HVG for each sample
            sc.pp.highly_variable_genes(
                adata_batch,
                flavor="cell_ranger",
                batch_key=batch_key,
                n_top_genes=n_hvg,
            )

            adata_scale = adata_scale_on_copy(adata_batch)

            # Add computed info back to unscaled data
            adata_batch.uns["neighbors"] = adata_scale.uns["neighbors"]
            adata_batch.uns["pca"] = adata_scale.uns["pca"]
            adata_batch.obsm["X_pca"] = adata_scale.obsm["X_pca"]
            adata_batch.obsp["connectivities"] = adata_scale.obsp["connectivities"]
            adata_batch.obsp["distances"] = adata_scale.obsp["distances"]

            # Save adata of sample
            adatas_batches[batch] = adata_batch
        del adata_batch
        del adata_scale

        # Compute differences in Moran's I of integrated vs. max in unintegrated per hvg
        batch_scores = []
        for batch, adata_batch in adatas_batches.items():
            morans_i_score(adata_batch, hvgs)
            morans_i_batch = adata_batch.var["morans_i"]
            morans_i_batch.name = batch
            batch_scores.append(morans_i_batch)

        batch_scores = pd.concat(batch_scores, axis=1)

        # Difference with integrated data, values between [-1,1]
        morans_i_difference = batch_scores.max(axis=1) - adata_post.var["morans_i"]
        # We only consider when the integration made the Moran's I score worse
        morans_i_difference = np.maximum(morans_i_difference, 0)

        score = 1 - morans_i_difference.mean()
        return score
    return score
