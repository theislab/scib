import numpy as np
import anndata
import scanpy as sc


# checker functions for data sanity
def checkAdata(adata):
    if type(adata) is not anndata.AnnData:
        raise TypeError('Input is not a valid AnnData object')

def checkBatch(batch, obs, verbose=False):
    if batch not in obs:
        raise ValueError(f'column {batch} is not in obs')
    elif verbose:
        print(f'Object contains {obs[batch].nunique()} batches.')

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


def splitBatches(adata, batch, hvg= None):
    split = []
    if hvg is not None:
        adata = adata[:, hvg]
    for i in adata.obs[batch].unique():
        split.append(adata[adata.obs[batch]==i])
    return split

# remove duplicated columns
def merge_adata(adata_list):
    adata = adata_list[0].concatenate(*adata_list[1:], index_unique=None)
    columns_to_keep = [name.split('-')[1] == '0' for name in adata.var.columns.values]
    clean_var = adata.var.loc[:, columns_to_keep]
    adata.var = clean_var.rename(columns={name : name.split('-')[0] for name in clean_var.columns.values})
    return adata

def todense(adata):
    import scipy
    if isinstance(adata.X, scipy.sparse.csr_matrix):
        adata.X = adata.X.todense()
