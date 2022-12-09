import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse.csgraph import connected_components

from ..utils import check_batch
from .utils import RootCellError


def trajectory_conservation(
    adata_pre, adata_post, label_key, pseudotime_key="dpt_pseudotime", batch_key=None
):
    """Trajectory conservation score

    Trajectory conservation is measured by  spearman’s rank correlation coefficient :math:`s`, between the pseudotime
    values before and after integration.
    The final score was scaled to a value between 0 and 1 using the equation

     .. math::

        trajectory \\, conservation = \\frac {s + 1} {2}

    :param adata_pre: unintegrated adata
    :param adata_post: integrated adata
    :param label_key: column in ``adata_pre.obs`` of the groups used to precompute the trajectory
    :param pseudotime_key: column in ``adata_pre.obs`` in which the pseudotime is saved in.
        Column can contain empty entries, the dataset will be subset to the cells with scores.
    :param batch_key: set to batch key if you want to compute the trajectory metric by batch. By default the batch
        information will be ignored (``batch_key=None``)

    This function requires pseudotime values in ``.obs`` of the unintegrated object (``adata_pre``) computed per batch
    and can be applied to all integration output types.
    The input trajectories should be curated manually as the quality of the metric depends on the quality of the metric
    depends on the quality of the annotation.
    The integrated object (``adata_post``) needs to have a kNN graph based on the integration output.
    See :ref:`preprocessing` for more information on preprocessing.

    **Examples**

    .. code-block:: python

        # feature output
        scib.pp.reduce_data(
            adata, n_top_genes=2000, batch_key="batch", pca=True, neighbors=True
        )
        scib.me.trajectory_conservation(adata_unintegrated, adata, label_key="cell_type")

        # embedding output
        sc.pp.neighbors(adata, use_rep="X_emb")
        scib.me.trajectory_conservation(adata_unintegrated, adata, label_key="celltype")

        # knn output
        scib.me.trajectory_conservation(adata_unintegrated, adata, label_key="celltype")

    """
    # subset to cells for which pseudotime has been computed
    cell_subset = adata_pre.obs.index[adata_pre.obs[pseudotime_key].notnull()]
    adata_pre_ti = adata_pre[cell_subset]
    adata_post_ti = adata_post[cell_subset]
    try:
        iroot, adata_post_ti2 = get_root(
            adata_pre_ti, adata_post_ti, label_key, pseudotime_key
        )
    except RootCellError:
        print("No root cell found, setting trajectory conservation metric to 0.")
        return 0  # failure to find root cell means no TI conservation

    adata_post_ti2.uns["iroot"] = iroot

    sc.tl.dpt(adata_post_ti2)  # stored in 'dpt_pseudotime'
    adata_post_ti2.obs.loc[
        adata_post_ti2.obs["dpt_pseudotime"] > 1, "dpt_pseudotime"
    ] = 0
    adata_post_ti.obs["dpt_pseudotime"] = 0
    adata_post_ti.obs["dpt_pseudotime"] = adata_post_ti2.obs["dpt_pseudotime"]
    adata_post_ti.obs["dpt_pseudotime"].fillna(0, inplace=True)

    if batch_key is None:
        pseudotime_before = adata_pre_ti.obs[pseudotime_key]
        pseudotime_after = adata_post_ti.obs["dpt_pseudotime"]
        correlation = pseudotime_before.corr(pseudotime_after, "spearman")
        return (correlation + 1) / 2  # scaled
    else:
        check_batch(batch_key, adata_pre.obs)
        check_batch(batch_key, adata_post.obs)

        # check if batches match
        if not np.array_equal(
            adata_post_ti.obs[batch_key], adata_pre_ti.obs[batch_key]
        ):
            raise ValueError(
                "Batch columns do not match\n"
                f"adata_post_ti.obs['batch']:\n {adata_post_ti.obs[batch_key]}\n"
                f"adata_pre_ti.obs['batch']:\n {adata_pre_ti.obs[batch_key]}\n"
            )

        corr = pd.Series()
        for i in adata_pre_ti.obs[batch_key].unique():
            pseudotime_before = adata_pre_ti.obs[adata_pre_ti.obs[batch_key] == i][
                pseudotime_key
            ]
            pseudotime_after = adata_post_ti.obs[adata_post_ti.obs[batch_key] == i][
                "dpt_pseudotime"
            ]
            corr[i] = pseudotime_before.corr(pseudotime_after, "spearman")

        return (corr.mean() + 1) / 2  # scaled


def get_root(adata_pre, adata_post, ct_key, pseudotime_key="dpt_pseudotime", dpt_dim=3):
    """Determine root cell for integrated adata based on unintegrated adata

    :param adata_pre: unintegrated adata
    :param adata_post: integrated adata
    :param label_key: column in ``adata_pre.obs`` of the groups used to precompute the trajectory
    :param pseudotime_key: column in ``adata_pre.obs`` in which the pseudotime is saved in.
        Column can contain empty entries, the dataset will be subset to the cells with scores.
    :param dpt_dim: number of diffmap dimensions used to determine root
    """
    n_components, adata_post.obs["neighborhood"] = connected_components(
        csgraph=adata_post.obsp["connectivities"], directed=False, return_labels=True
    )

    start_clust = adata_pre.obs.groupby([ct_key]).mean()[pseudotime_key].idxmin()
    min_dpt = adata_pre.obs[adata_pre.obs[ct_key] == start_clust].index
    which_max_neigh = (
        adata_post.obs["neighborhood"]
        == adata_post.obs["neighborhood"].value_counts().idxmax()
    )
    min_dpt = [
        value for value in min_dpt if value in adata_post.obs[which_max_neigh].index
    ]

    adata_post_ti = adata_post[which_max_neigh]

    min_dpt = [adata_post_ti.obs_names.get_loc(i) for i in min_dpt]

    # compute Diffmap for adata_post
    sc.tl.diffmap(adata_post_ti)

    # determine most extreme cell in adata_post Diffmap
    min_dpt_cell = np.zeros(len(min_dpt))
    for dim in np.arange(dpt_dim):

        diffmap_mean = adata_post_ti.obsm["X_diffmap"][:, dim].mean()
        diffmap_min_dpt = adata_post_ti.obsm["X_diffmap"][min_dpt, dim]

        # count opt cell
        if len(diffmap_min_dpt) == 0:
            raise RootCellError("No root cell in largest component")

        # choose optimum function
        if len(diffmap_min_dpt) > 0 and diffmap_min_dpt.mean() < diffmap_mean:
            opt = np.argmin
        else:
            opt = np.argmax

        min_dpt_cell[opt(diffmap_min_dpt)] += 1

    # root cell is cell with max vote
    return min_dpt[np.argmax(min_dpt_cell)], adata_post_ti
