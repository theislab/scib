import numpy as np
import anndata

# first, some checker functions for data sanity
def checkAdata(adata):
    if type(adata) is not anndata.AnnData:
        raise TypeError('Input is not a valid AnnData object')

def checkBatch(batch, obs):
    if batch not in obs:
        raise ValueError('Selected batch column is not in obs')
    else:
        nBatch = obs[batch].nunique()
        print('Object contains '+str(nBatch)+' batches.')

def checkHVG(hvg, adata_var):
    if type(hvg) is not list:
        raise TypeError('HVG list is not a list')
    else:
        if not all(i in adata_var.index for i in hvg):
            raise ValueError('Not all HVGs are in the adata object')

def checkSanity(adata, batch, hvg):
    checkAdata(adata)
    checkBatch(batch, adata.obs)
    if hvg is not None:
        checkHVG(hvg, adata.var)


def splitBatches(adata, batch):
    split = []
    for i in adata.obs[batch].unique():
        split.append(adata[adata.obs[batch]==i])
    return split

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
