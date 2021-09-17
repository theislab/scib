import numpy as np
import scanpy as sc

from scIB.utils import splitBatches


def precompute_hvg_batch(adata, batch, features, n_hvg=500, save_hvg=False):
    adata_list = splitBatches(adata, batch, hvg=features)
    hvg_dir = {}
    for i in adata_list:
        sc.pp.filter_genes(i, min_cells=1)
        n_hvg_tmp = np.minimum(n_hvg, int(0.5 * i.n_vars))
        if n_hvg_tmp < n_hvg:
            print(i.obs[batch][0] + ' has less than the specified number of genes')
            print('Number of genes: ' + str(i.n_vars))
        hvg = sc.pp.highly_variable_genes(i, flavor='cell_ranger', n_top_genes=n_hvg_tmp, inplace=False)
        hvg_dir[i.obs[batch][0]] = i.var.index[hvg['highly_variable']]

    if save_hvg:
        adata.uns['hvg_before'] = hvg_dir
    else:
        return hvg_dir


def hvg_overlap(adata_pre, adata_post, batch, n_hvg=500, verbose=False):
    hvg_post = adata_post.var_names

    adata_post_list = splitBatches(adata_post, batch)
    overlap = []

    if ('hvg_before' in adata_pre.uns_keys()) and (set(hvg_post) == set(adata_pre.var_names)):
        print('Using precomputed hvgs per batch')
        hvg_pre_list = adata_pre.uns['hvg_before']
    else:
        hvg_pre_list = precompute_hvg_batch(adata_pre, batch, hvg_post)

        for i in range(len(adata_post_list)):  # range(len(adata_pre_list)):
            sc.pp.filter_genes(adata_post_list[i], min_cells=1)  # remove genes unexpressed (otherwise hvg might break)

            # ov = list(set(adata_pre_list[i].var_names).intersection(set(hvg_pre_list[i])))
            # adata_pre_list[i] = adata_pre_list[i][:,ov]
            # adata_post_list[i] = adata_post_list[i][:,ov]
            batch_var = adata_post_list[i].obs[batch][0]

            n_hvg_tmp = len(hvg_pre_list[batch_var])
            # adata_pre.uns['n_hvg'][hvg_post]#np.minimum(n_hvg, int(0.5*adata_post_list[i].n_vars))
            if verbose:
                print(n_hvg_tmp)
            # if n_hvg_tmp<n_hvg:
            #    print(adata_post_list[i].obs[batch][0]+' has less than the specified number of genes')
            #    print('Number of genes: '+str(adata_post_list[i].n_vars))
            # hvg_pre = sc.pp.highly_variable_genes(adata_pre_list[i], flavor='cell_ranger', n_top_genes=n_hvg_tmp, inplace=False)
            tmp_pre = hvg_pre_list[batch_var]  # adata_pre_list[i].var.index[hvg_pre['highly_variable']]
            hvg_post = sc.pp.highly_variable_genes(adata_post_list[i], flavor='cell_ranger', n_top_genes=n_hvg_tmp,
                                                   inplace=False)
            tmp_post = adata_post_list[i].var.index[hvg_post['highly_variable']]
            n_hvg_real = np.minimum(len(tmp_pre), len(tmp_post))
            overlap.append((len(set(tmp_pre).intersection(set(tmp_post)))) / n_hvg_real)
    return np.mean(overlap)
