import numpy as np
import scanpy as sc

try:
    from scanpy._utils import renamed_arg
except ImportError:
    from .._package_tools import renamed_arg

from ..utils import split_batches


def precompute_hvg_batch(adata, batch, features, n_hvg=500, save_hvg=False):
    """
    Compute HVGs per batch

    :param adata: anndata object
    :param batch: key in adata.obs
    :param features: features to subset to
    :param n_hvg: maximum number of HVGs to compute
    :save_hvg: whether to add hvg per batch information to adata object
    :return:
        dictionary of batch to HVG list
    """
    adata_list = split_batches(adata, batch, hvg=features)
    hvg_dir = {}
    for i in adata_list:
        sc.pp.filter_genes(i, min_cells=1)
        n_hvg_tmp = np.minimum(n_hvg, int(0.5 * i.n_vars))
        if n_hvg_tmp < n_hvg:
            print(i.obs[batch][0] + " has less than the specified number of genes")
            print("Number of genes: " + str(i.n_vars))
        hvg = sc.pp.highly_variable_genes(
            i, flavor="cell_ranger", n_top_genes=n_hvg_tmp, inplace=False
        )
        hvg_dir[i.obs[batch][0]] = i.var.index[hvg["highly_variable"]]

    if save_hvg:
        adata.uns["hvg_before"] = hvg_dir
    else:
        return hvg_dir


@renamed_arg("batch", "batch_key")
def hvg_overlap(adata_pre, adata_post, batch_key, n_hvg=500, verbose=False):
    """Highly variable gene overlap

    Metric that computes the average percentage of overlapping highly variable genes per batch pre post integration.

    :param adata_pre: Unintegrated anndata object
    :param adata_post: Integrated anndata object
    :param batch_key: Batch variable in ``adata_post.obs``
    :param n_hvg: Number of hvgs to compute per batch
    :return:
        Average percentage of overlapping highly variable genes

    The score can only be computed on feature spaces.
    No preprocessing is needed, as the function will perform highly variable gene selection.

    **Example**

    .. code-block:: python

        # full feature output
        scib.me.hvg_overlap(adata_unintegrated, adata, batch_key="batch")
    """
    hvg_post = adata_post.var_names

    adata_post_list = split_batches(adata_post, batch_key)
    overlap = []

    hvg_pre_list = precompute_hvg_batch(adata_pre, batch_key, hvg_post, n_hvg=n_hvg)

    for ad_post in adata_post_list:  # range(len(adata_pre_list)):
        # remove genes unexpressed (otherwise hvg might break)
        sc.pp.filter_genes(ad_post, min_cells=1)
        batch_var = ad_post.obs[batch_key][0]
        n_hvg_tmp = len(hvg_pre_list[batch_var])

        if verbose:
            print(n_hvg_tmp)

        tmp_pre = hvg_pre_list[batch_var]
        hvg_post = sc.pp.highly_variable_genes(
            ad_post, flavor="cell_ranger", n_top_genes=n_hvg_tmp, inplace=False
        )
        tmp_post = ad_post.var.index[hvg_post["highly_variable"]]
        n_hvg_real = np.minimum(len(tmp_pre), len(tmp_post))
        overlap.append((len(set(tmp_pre).intersection(set(tmp_post)))) / n_hvg_real)
    return np.mean(overlap)
