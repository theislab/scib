import os
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
from scIB.tests import utils
import warnings
warnings.filterwarnings('ignore')


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

def test_all_metrics():
    adata = utils.create_adata_dummy()
    adata_int = adata.copy()
    
    for ot in ["full", "embed", "knn"]:
        all_metrics(adata, adata_int, script="metrics.py", type_=ot, method="orig")
   

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
