import scanpy as sc
from scIB.preprocessing import reduce_data
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def create_adata_dummy(pca=False, n_top_genes=None, neighbors=False):
    adata = sc.datasets.paul15()
    adata.obs['celltype'] = adata.obs['paul15_clusters']
    np.random.seed(42)
    adata.obs['batch'] = np.random.randint(1, 5, adata.n_obs)
    adata.obs['batch'] = adata.obs['batch'].astype(str)
    adata.obs['batch'] = adata.obs['batch'].astype("category")
    adata.layers['counts'] = adata.X
    reduce_data(adata, pca=pca, n_top_genes=n_top_genes,
                umap=False, neighbors=neighbors)
    return adata

def add_emb(adata, type_='pca'):
    if type_ == 'pca':
        if 'X_pca' in adata.obsm:
            mtx = adata.obsm['X_pca']
        else:
            mtx = sc.tl.pca(adata, copy=True).obsm['X_pca']
    elif type_ == 'full':
        mtx = adata.X
    adata.obsm['X_emb'] = mtx
