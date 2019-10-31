import scanpy as sc
from scIB import clustering as cl
from scIB.tests import utils
from scIB import metrics as me
from scIB.integration import *
import numpy as np
import warnings
warnings.filterwarnings('ignore')


def cluster(adata, label_key, cluster_key, verbose=False):
    return cl.opt_louvain(adata, label_key=label_key, cluster_key=cluster_key,
                plot=False, inplace=True, force=True, verbose=verbose)


def silhouette():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000)
    score = me.silhouette(adata, group_key='celltype', embed='X_pca', scale=True)
    print(score)
    assert score >= 0
    assert score <= 1
    
def silhouette_batch():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000)
    _, sil = me.silhouette_batch(adata, batch_key='batch', group_key='celltype', 
                                      embed='X_pca', scale=True, verbose=False)
    score = sil['silhouette_score'].mean()
    print(score)
    assert score >= 0
    assert score <= 1
    
def nmi():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
    _, _, nmi_all = cluster(adata, cluster_key='cluster', label_key='celltype', verbose=True)
    for score in nmi_all['score']:
        assert score >= 0
        assert score <= 1
    
def ari():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
    cluster(adata, cluster_key='cluster', label_key='celltype')
    score = me.ari(adata, group1='cluster', group2='celltype')
    print(score)
    assert score >= 0
    assert score <= 1

def pcr_comparison():
    # no PCA precomputed
    adata = utils.create_adata_dummy()
    adata_int = adata.copy()
    score = me.pcr_comparison(adata, adata_int, covariate='batch', scale=True, verbose=False)
    print(score)
    assert score >= 0
    assert score <= 1
    
    # precomputed PCA
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000)
    adata_int = adata.copy()
    score = me.pcr_comparison(adata, adata_int, covariate='batch', scale=True, verbose=False)
    print(score)
    assert score >= 0
    assert score <= 1
    
    # use different embedding
    utils.add_emb(adata_int)
    score = me.pcr_comparison(adata, adata_int, covariate='batch', embed='X_emb',
                                   scale=True, verbose=False)
    print(score)
    assert score >= 0
    assert score <= 1

def cell_cycle():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000)
    adata_int = adata.copy()
    
    score = me.cell_cycle(adata, adata_int, batch_key='batch', organism='mouse')
    print(score)
    assert score >= 0
    assert score <= 1

def hvg_overlap():
    adata = utils.create_adata_dummy()
    adata_int = adata.copy()
    score = me.hvg_overlap(adata_int, adata, batch='batch', n_hvg=500)
    print(score)
    assert score >= 0
    assert score <= 1
    
    
def metrics_all_methods():
    adata = utils.create_adata_dummy()
    
    methods = {
        'scanorama': runScanorama,
        'trvae': runTrVae,
        'seurat': runSeurat,
        'harmony': runHarmony,
        'mnn': runMNN,
        'bbknn': runBBKNN,
        'conos': runConos,
        'scvi': runScvi
    }
    # for name, func in methods.items():
      
