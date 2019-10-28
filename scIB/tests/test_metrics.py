import scanpy as sc
import scIB
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def create_adata_dummy(pca=False, n_top_genes=None, neighbors=False):
    adata = sc.datasets.paul15()
    adata.obs['celltype'] = adata.obs['paul15_clusters']
    adata.obs['batch'] = np.random.randint(1,5,adata.n_obs)
    scIB.pp.reduce_data(adata, pca=pca, n_top_genes=n_top_genes,
                        umap=False, neighbors=neighbors)
    return adata

def add_emb(adata):
    if 'X_pca' in adata.obsm:
        pca = adata.obsm['X_pca']
    else:
        pca = sc.tl.pca(adata, copy=True).obsm['X_pca']
    adata.obsm['X_emb'] = pca

def cluster(adata, label_key, cluster_key, verbose=False):
    return scIB.cl.opt_louvain(adata, label_key=label_key, cluster_key=cluster_key,
                plot=False, inplace=True, force=True, verbose=verbose)

# Metrics
def silhouette():
    adata = create_adata_dummy(pca=True, n_top_genes=2000)
    score = scIB.me.silhouette(adata, group_key='celltype', embed='X_pca', scale=True)
    print(score)
    assert score >= 0
    assert score <= 1
    
def silhouette_batch():
    adata = create_adata_dummy(pca=True, n_top_genes=2000)
    _, sil = scIB.me.silhouette_batch(adata, batch_key='batch', group_key='celltype', 
                                      embed='X_pca', scale=True, verbose=False)
    score = sil['silhouette_score'].mean()
    print(score)
    assert score >= 0
    assert score <= 1
    
def nmi():
    adata = create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
    _, _, nmi_all = cluster(adata, cluster_key='cluster', label_key='celltype', verbose=True)
    for score in nmi_all['NMI']:
        assert score >= 0
        assert score <= 1
    
def ari():
    adata = create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
    cluster(adata, cluster_key='cluster', label_key='celltype')
    score = scIB.me.ari(adata, group1='cluster', group2='celltype')
    print(score)
    assert score >= 0
    assert score <= 1

def pcr_comparison():
    # no PCA precomputed
    adata = create_adata_dummy()
    adata_int = adata.copy()
    score = scIB.me.pcr_comparison(adata, adata_int, covariate='batch', scale=True, verbose=False)
    print(score)
    assert score >= 0
    assert score <= 1
    
    # precomputed PCA
    adata = create_adata_dummy(pca=True, n_top_genes=2000)
    adata_int = adata.copy()
    score = scIB.me.pcr_comparison(adata, adata_int, covariate='batch', scale=True, verbose=False)
    print(score)
    assert score >= 0
    assert score <= 1
    
    # use different embedding
    add_emb(adata_int)
    score = scIB.me.pcr_comparison(adata, adata_int, covariate='batch', embed='X_emb',
                                   scale=True, verbose=False)
    print(score)
    assert score >= 0
    assert score <= 1

def cell_cycle():
    adata = create_adata_dummy(pca=True, n_top_genes=2000)
    adata_int = adata.copy()
    
    score = scIB.me.cell_cycle(adata, adata_int, batch_key='batch', organism='mouse')
    print(score)
    assert score >= 0
    assert score <= 1

def hvg_overlap():
    adata = create_adata_dummy()
    adata_int = adata.copy()
    score = scIB.me.hvg_overlap(adata_int, adata, batch='batch', n_hvg=500)
    print(score)
    assert score >= 0
    assert score <= 1
    
    
def metrics_all_methods():
    adata = sc.create_adata_dummy()
    
    methods = {
        'scanorama': scIB.integration.runScanorama,
        'trvae': scIB.integration.runTrVae,
        'seurat': scIB.integration.runSeurat,
        'harmony': scIB.integration.runHarmony,
        'mnn': scIB.integration.runMNN,
        'bbknn': scIB.integration.runBBKNN,
        'conos': scIB.integration.runConos,
        'scvi': scIB.integration.runScvi
    }
    # for name, func in methods.items():
      
