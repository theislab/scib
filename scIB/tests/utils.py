import scanpy as sc
from scIB.preprocessing import reduce_data
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def create_adata_dummy(pca=False, n_top_genes=None, neighbors=False):
    adata = sc.datasets.paul15()
    adata.obs['celltype'] = adata.obs['paul15_clusters']
    adata.obs['batch'] = np.random.randint(1,5,adata.n_obs)
    reduce_data(adata, pca=pca, n_top_genes=n_top_genes,
                umap=False, neighbors=neighbors)
    return adata

def add_emb(adata):
    if 'X_pca' in adata.obsm:
        pca = adata.obsm['X_pca']
    else:
        pca = sc.tl.pca(adata, copy=True).obsm['X_pca']
    adata.obsm['X_emb'] = pca
