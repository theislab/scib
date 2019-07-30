import numpy as np
import anndata
import scanpy as sc


def import_rpy2():
    global callbacks
    import rpy2.rinterface_lib.callbacks
    global logging
    import logging
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Ignore R warning messages
    global ro
    import rpy2.robjects as ro
    global pandas2ri
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    global anndata2ri
    import anndata2ri
    anndata2ri.activate()

# checker functions for data sanity
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
    print(adata.n_vars)
    columns_to_keep = [name.split('-')[1] == '0' for name in adata.var.columns.values]
    clean_var = adata.var.loc[:, columns_to_keep]
    adata.var = clean_var.rename(columns={name : name.split('-')[0] for name in clean_var.columns.values})
    return adata

def summarize_counts(adata, count_matrix=None, mito=True):
    if count_matrix is None:
        count_matrix = adata.X
    adata.obs['n_counts'] = count_matrix.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes'] = (count_matrix > 0).sum(1)

    if mito:
        mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
        mt_count = count_matrix[:, mt_gene_mask].sum(1)
        if mt_count.ndim > 1:
            mt_count = np.squeeze(np.asarray(mt_count))
        adata.obs['mt_frac'] = mt_count/adata.obs['n_counts']

def subsetHVG(adata, batch, number):
    ## does not work yet, use hvg_intersect
    import scanpy as sc
    hvg = sc.pp.highly_variable_genes(adata, n_top_genes=number, batch_key=batch, flavor='cell_ranger', inplace=False)
    return hvg

def hvg_intersect(adata, batch, num=4000):
    split = splitBatches(adata, 'sample')
    variable = []
    for i in split:
        print(i)
        tmp = sc.pp.highly_variable_genes(i, flavor='cell_ranger', n_top_genes=num, inplace=False)
        variable.append(set(i.var[[j[0] for j in tmp]].index))
    return list(variable[0].intersection(*variable[1:]))

