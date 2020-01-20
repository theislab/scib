import os
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
from scIB import clustering as cl
from scIB.tests import utils
from scIB import metrics as me
from scIB.integration import *
import warnings
warnings.filterwarnings('ignore')


def cluster(adata, label_key, cluster_key, verbose=False):
    return cl.opt_louvain(adata, label_key=label_key, cluster_key=cluster_key,
                plot=False, inplace=True, force=True, verbose=verbose)


def silhouette():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000)
    score = me.silhouette(adata, group_key='celltype', embed='X_pca', scale=True)
    print(f"score: {score}")
    assert score >= 0
    assert score <= 1
    
def silhouette_batch():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000)
    _, sil = me.silhouette_batch(adata, batch_key='batch', group_key='celltype', 
                                      embed='X_pca', scale=True, verbose=False)
    score = sil['silhouette_score'].mean()
    print(f"score: {score}")
    assert score >= 0
    assert score <= 1
    
def nmi():
    
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
    
    # trivial score
    score = scIB.me.nmi(adata, 'celltype', 'celltype')
    assert score == 1
    
    # on cell type 
    _, _, nmi_all = cluster(adata, cluster_key='cluster', label_key='celltype', verbose=True)
    for score in nmi_all['score']:
        print(score)
        assert score >= 0
        assert score <= 1
    
def ari():
    
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
    
    # trivial score
    score = scIB.me.ari(adata, 'celltype', 'celltype')
    assert score == 1
    
    # on cell type 
    cluster(adata, cluster_key='cluster', label_key='celltype')
    score = me.ari(adata, group1='cluster', group2='celltype')
    print(f"score: {score}")
    assert score >= 0
    assert score <= 1

def pcr_comparison():
    
    verbose = True
    
    # no PCA precomputed
    adata = utils.create_adata_dummy()
    adata_int = adata.copy()
    score = me.pcr_comparison(adata, adata_int, covariate='batch', n_comps=50,
                              scale=True, verbose=verbose)
    print(f"no PCA precomputed: {score}")
    assert score < 1e-6
    
    # use different embedding
    adata = utils.create_adata_dummy()
    adata_int = adata.copy()
    utils.add_emb(adata_int, type_='full')
    score = me.pcr_comparison(adata, adata_int, covariate='batch', embed='X_emb', n_comps=50,
                                   scale=True, verbose=verbose)
    print(f"using embedding: {score}")
    assert score >= 0
    assert score <= 1
    assert score < 1e-6
    
    # precomputed PCA
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000)
    adata_int = adata.copy()
    score = me.pcr_comparison(adata, adata_int, covariate='batch', scale=True, verbose=verbose)
    print(f"precomputed PCA: {score}")
    assert score == 0 # same PCA values -> difference should be 0

def cell_cycle():
    adata = utils.create_adata_dummy()
    adata_int = adata.copy()
    
    # only final score implementation
    score = me.cell_cycle(adata, adata_int, batch_key='batch', organism='mouse',
                          agg_func=np.mean, verbose=True)
    print(f"score: {score}")    
    assert score == 1
    
    # get all intermediate scores
    scores_df = me.cell_cycle(adata, adata_int, batch_key='batch', organism='mouse',
                          agg_func=None, verbose=True)
    print(f"score: {scores_df}")
    assert isinstance(scores_df, pd.DataFrame)
    for i in scores_df['score']:
        assert i == 1

def hvg_overlap():
    adata = utils.create_adata_dummy()
    adata_int = adata.copy()
    score = me.hvg_overlap(adata_int, adata, batch='batch', n_hvg=500)
    print(f"score: {score}")
    assert score == 1
    
def isolated_labels():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
    
    # test 2 different implementations of score
    for impl in [True, False]:
        score = me.isolated_labels(adata, label_key='celltype', 
                                   batch_key='batch', cluster=impl,
                                   n=4, verbose=True)
        print(f"score: {score}")
        assert score <= 1
        assert score >= 0
    
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

def all_metrics(adata_u, adata_i, script, type_, method, dir_="./data", verbose=False):
    """
    params:
        adata_u: unintegrated anndata
        adata_i: integrated anndata
        script: path to metrics.py script
        dir_: directory to test output
        type_: one of 'full', 'embed', 'knn'
        method: name of method, for saving file
    """
    
    #script = os.path.join(os.path.dirname(scIB.__file__), "scripts", "metrics.py")
    
    unintegrated = os.path.join(dir_, "unintegrated.h5ad")
    integrated = os.path.join(dir_, f"{method}.h5ad")
    metrics_dir = os.path.join(dir_, "metrics_out")
    
    if not os.path.exists(metrics_dir):
        os.makedirs(metrics_dir)
    
    adata_u.write(unintegrated)
    adata_i.write(integrated)
    
    cmd = ["python", script, "-u", unintegrated, "-i",  integrated,
         "-o", metrics_dir, "-b", "batch", "-l", "celltype", "--type", type_, 
        "--organism", "mouse"]
    if verbose:
        cmd.append("-v")
        
    call = subprocess.Popen(cmd, 
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    
    stdout, stderr = call.communicate()
    print(stdout.decode())
    print(stderr)
    
    metrics_file = os.path.join(metrics_dir, f"{method}_{type_}_metrics.csv")
    metrics = pd.read_csv(metrics_file, index_col=0)
    
    for metric, value in metrics.iterrows():
        value = value[0]
        print(f'{metric}: {value}')
        if np.isnan(value):
            continue
        assert value >= 0
        assert value <= 1
