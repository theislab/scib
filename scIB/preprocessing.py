import scanpy as sc
import rpy2.rinterface_lib.callbacks
import logging
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
import anndata2ri

# Ignore R warning messages
#Note: this can be commented out to get more verbose R output
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

## Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()

def normalize(adata, color_col='batch', min_mean = 0.1):
    ro.r('library("scran")')
    
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15, svd_solver='arpack')
    sc.pp.neighbors(adata_pp)
    sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)
    
    ro.globalenv['data_mat'] = adata.X.T
    ro.globalenv['input_groups'] = adata_pp.obs['groups']
    size_factors = ro.r(f'computeSumFactors(data_mat, clusters = input_groups, min.mean = {min_mean})')
    adata.obs['size_factors'] = size_factors