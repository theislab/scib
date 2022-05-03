import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

from scib.preprocessing import hvg_batch


def morans_i(
    adata_pre,
    adata_post,
    batch_key,
    n_hvg=1000,
    embed="X_pca",
    rescale=True,
    compare_pre=False,
    hvgs=None,
):
    """
    Compute Moran's I on HVGs (defined across batches) on integrated data
        as a measure of bio conservation. The metric can be computed either as
        A.) mean difference between Morans's I on embedding of integrated data and
        max Moran's I of individual non-integrated batches or
        B.) mean Morans's I on embedding of integrated data.
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
    # Prepare adatas
    # Get integrated connectivities for Moran's I
    adata_post = adata_post.copy()
    sc.pp.neighbors(adata_post, use_rep=embed)
    # Prepare pre data
    adata_pre = adata_pre.copy()
    adata_pre.obs[batch_key] = adata_pre.obs[batch_key].astype("category")

    # Get HVGs across samples
    if hvgs is None:
        hvgs = hvg_batch(
            adata_pre, batch_key=batch_key, target_genes=n_hvg, flavor="cell_ranger"
        )
    else:
        assert adata_pre.var_names.isin(hvgs).sum() == len(hvgs)

    def compute_mi(data, hvgs):
        """
        Helper function for computing Moran's I on HVGS that are non-constant
        and adding them to var
        :param data: adata, result is added to var
        :param hvgs: Genes for which to compute Moran's I
        """
        # Make sure that genes used for computation are non-constant
        hvg_data_0 = data[0, hvgs].X
        hvg_data_0 = hvg_data_0.A if issparse(hvg_data_0) else hvg_data_0
        hvg_data_0 = hvg_data_0.squeeze()
        mean_hvgs = np.asarray(data[:, hvgs].X.mean(axis=0)).squeeze()

        hvgs_used = np.array(hvgs)[mean_hvgs != hvg_data_0]

        if len(hvgs) > len(hvgs_used):
            print(
                "Using subset (%i) of hvgs (%i) that are not constant in the data"
                % (len(hvgs_used), len(hvgs))
            )
        # Compute Moran's I
        morans_i = sc.metrics.morans_i(data[:, hvgs_used])
        morans_i = pd.Series(morans_i, index=hvgs_used)
        # This should not happen, just in case
        if (morans_i > 1).any() or (morans_i < -1).any():
            raise ValueError(
                "Problems in computing Moran's I, value not between -1 and 1"
            )
        data.var["morans_i"] = morans_i

    if compare_pre:
        # Get per sample connectivities
        adatas_sample = dict()
        for sample in adata_pre.obs[batch_key].unique():
            # Subset only to cells from 1 sample
            adata_sample = adata_pre[adata_pre.obs[batch_key] == sample, :].copy()
            # Compute HVG for each sample
            sc.pp.highly_variable_genes(
                adata_sample,
                flavor="cell_ranger",
                batch_key="study_sample",
                n_top_genes=2000,
            )
            adata_sample_scl = adata_sample.copy()
            # PP each sample: scaling, pca, neighbours
            sc.pp.scale(adata_sample_scl, max_value=10)
            sc.pp.pca(
                adata_sample_scl,
                n_comps=15,
                use_highly_variable=True,
                svd_solver="arpack",
            )
            sc.pp.neighbors(adata_sample_scl, n_pcs=15)
            # Add computed info back to unscaled data
            adata_sample.uns["neighbors"] = adata_sample_scl.uns["neighbors"]
            adata_sample.obsp["connectivities"] = adata_sample_scl.obsp[
                "connectivities"
            ]
            adata_sample.obsp["distances"] = adata_sample_scl.obsp["distances"]
            adata_sample.uns["pca"] = adata_sample_scl.uns["pca"]
            adata_sample.obsm["X_pca"] = adata_sample_scl.obsm["X_pca"]
            # Save adata of sample
            adatas_sample[sample] = adata_sample
        del adata_sample
        del adata_sample_scl

        # Compute Moran's I across samples and on integrated data
        for data in list(adatas_sample.values()) + [adata_post]:
            compute_mi(data, hvgs)

        # Compute differences in Moran's I
        # Table of Moran's I-s across batches
        batch_mis = []
        for sample, data in adatas_sample.items():
            mi = data.var["morans_i"]
            mi.name = sample
            batch_mis.append(mi)
        batch_mis = pd.concat(batch_mis, axis=1)
        # Difference with integrated data
        mi_diffs = adata_post.var["morans_i"] - batch_mis.max(axis=1)
        avg_mi_diff = mi_diffs.mean()

        # Rescale so that it is between [0,1] where 1 is better
        # Moran's I will be in [-1,1] and thus difference can be in [-2,2]
        if rescale:
            res = (avg_mi_diff + 2) / 4
        else:
            res = avg_mi_diff

    else:
        # Moran's I on integrated data
        compute_mi(adata_post, hvgs)
        res = adata_post.var["morans_i"].mean()
        if rescale:
            # Moran's I will be in [-1,1]
            res = (res + 1) / 2

    return res
