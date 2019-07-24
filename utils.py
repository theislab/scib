import numpy as np

def summarize_counts(adata, count_matrix=None):
    if count_matrix is None:
        count_matrix = adata.X
    adata.obs['n_counts'] = count_matrix.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes'] = (count_matrix > 0).sum(1)

    mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
    mt_count = count_matrix[:, mt_gene_mask].sum(1)
    if mt_count.ndim > 1:
        mt_count = np.squeeze(np.asarray(mt_count))
    adata.obs['mt_frac'] = mt_count/adata.obs['n_counts']
