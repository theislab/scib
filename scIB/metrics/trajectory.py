import numpy as np
from scipy.sparse.csgraph import connected_components
import scanpy as sc

from .utils import RootCellError


def get_root(
        adata_pre,
        adata_post,
        ct_key,
        pseudotime_key="dpt_pseudotime",
        dpt_dim=3
):
    """
    Determine root cell for integrated adata based on unintegrated adata

    :param adata_pre: unintegrated adata
    :param adata_post: integrated adata
    :param label_key: column in `adata_pre.obs` of the groups used to precompute the trajectory
    :param pseudotime_key: column in `adata_pre.obs` in which the pseudotime is saved in.
        Column can contain empty entries, the dataset will be subset to the cells with scores.
    :param dpt_dim: number of diffmap dimensions used to determine root
    """
    n_components, adata_post.obs['neighborhood'] = connected_components(
        csgraph=adata_post.obsp['connectivities'],
        directed=False,
        return_labels=True
    )

    start_clust = adata_pre.obs.groupby([ct_key]).mean()[pseudotime_key].idxmin()
    min_dpt = adata_pre.obs[adata_pre.obs[ct_key] == start_clust].index
    which_max_neigh = adata_post.obs['neighborhood'] == adata_post.obs['neighborhood'].value_counts().idxmax()
    min_dpt = [value for value in min_dpt if value in adata_post.obs[which_max_neigh].index]

    adata_post_sub = adata_post[which_max_neigh]

    min_dpt = [adata_post_sub.obs_names.get_loc(i) for i in min_dpt]

    # compute Diffmap for adata_post
    sc.tl.diffmap(adata_post_sub)

    # determine most extreme cell in adata_post Diffmap
    min_dpt_cell = np.zeros(len(min_dpt))
    for dim in np.arange(dpt_dim):

        diffmap_mean = adata_post_sub.obsm["X_diffmap"][:, dim].mean()
        diffmap_min_dpt = adata_post_sub.obsm["X_diffmap"][min_dpt, dim]

        # count opt cell
        if len(diffmap_min_dpt) == 0:
            raise RootCellError('No root cell in largest component')

        # choose optimum function
        if len(diffmap_min_dpt) > 0 and diffmap_min_dpt.mean() < diffmap_mean:
            opt = np.argmin
        else:
            opt = np.argmax

        min_dpt_cell[opt(diffmap_min_dpt)] += 1

    # root cell is cell with max vote
    return min_dpt[np.argmax(min_dpt_cell)], adata_post_sub


def trajectory_conservation(
        adata_pre,
        adata_post,
        label_key,
        pseudotime_key="dpt_pseudotime",
	batch_key=False
):
    """
    :param adata_pre: unintegrated adata
    :param adata_post: integrated adata
    :param label_key: column in `adata_pre.obs` of the groups used to precompute the trajectory
    :param pseudotime_key: column in `adata_pre.obs` in which the pseudotime is saved in.
        Column can contain empty entries, the dataset will be subset to the cells with scores.
    :param batch_key: set to batch key if if you want to compute the trajectory metric by batch
    """
    # subset to cells for which pseudotime has been computed
    cell_subset = adata_pre.obs.index[adata_pre.obs[pseudotime_key].notnull()]
    adata_pre_sub = adata_pre[cell_subset]
    adata_post_sub = adata_post[cell_subset]
    if batch_key==False:
        try:
            iroot, adata_post_sub2 = get_root(adata_pre_sub, adata_post_sub, label_key, pseudotime_key)
        except RootCellError:
            print('No root cell found, setting trajectory conservation metric to 0.')
            return 0  # failure to find root cell means no TI conservation

        adata_post_sub2.uns['iroot'] = iroot

        sc.tl.dpt(adata_post_sub2)  # stored in 'dpt_pseudotime'
        adata_post_sub2.obs.loc[adata_post_sub2.obs['dpt_pseudotime'] > 1, 'dpt_pseudotime'] = 0
        adata_post_sub.obs['dpt_pseudotime'] = 0
        adata_post_sub.obs['dpt_pseudotime'] = adata_post_sub2.obs['dpt_pseudotime']
        adata_post_sub.obs['dpt_pseudotime'].fillna(0, inplace=True)

        pseudotime_before = adata_post_sub.obs['dpt_pseudotime']
        pseudotime_after = adata_pre_sub.obs[pseudotime_key]
        correlation = pseudotime_before.corr(pseudotime_after, 'spearman')

        return (correlation + 1) / 2  # scaled
    else:
        corr = pd.Series()
        for i in batch_key:
            adata_pre_sub_tmp = adata_pre_sub[adata_pre_sub.obs[batch_key]==i]
            adata_post_sub_tmp = adata_post_sub[adata_post_sub.obs[batch_key]==i]
            try:
                iroot, adata_post_sub2_tmp = get_root(adata_pre_sub_tmp, adata_post_sub_tmp, label_key, pseudotime_key)
            except RootCellError:
                print('No root cell found, setting trajectory conservation metric to 0.')
                return 0  # failure to find root cell means no TI conservation
            
            adata_post_sub2_tmp.uns['iroot'] = iroot
            
            sc.tl.dpt(adata_post_sub2_tmp)  # stored in 'dpt_pseudotime'
            adata_post_sub2_tmp.obs.loc[adata_post_sub2_tmp.obs['dpt_pseudotime'] > 1, 'dpt_pseudotime'] = 0
            adata_post_sub_tmp.obs['dpt_pseudotime'] = 0
            adata_post_sub_tmp.obs['dpt_pseudotime'] = adata_post_sub2_tmp.obs['dpt_pseudotime']
            adata_post_sub_tmp.obs['dpt_pseudotime'].fillna(0, inplace=True)
            
            pseudotime_before = adata_post_sub_tmp.obs['dpt_pseudotime']
            pseudotime_after = adata_pre_sub_tmp.obs[pseudotime_key]
            corr[i] = pseudotime_before.corr(pseudotime_after, 'spearman')
        return (corr.mean() + 1) / 2 # scaled
