from .test_pipeline import *
from scIB.integration import *

import os
import subprocess
import pandas as pd


def metrics_all_methods(adata_factory):
    adata = adata_factory()

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


def test_all_metrics(adata_factory, test_metrics):
    adata = adata_factory()
    adata_int = adata.copy()

    script = os.path.join(os.path.dirname(scIB.__file__), "scripts", "metrics.py")

    for ot in ["full", "embed", "knn"]:
        all_metrics(adata, adata_int, script=script, type_=ot, pipeline_dir=test_metrics, method="orig")


def all_metrics(adata_u, adata_i, script, type_, method, pipeline_dir, verbose=False):
    """
    params:
        adata_u: unintegrated anndata
        adata_i: integrated anndata
        script: path to metrics.py script
        pipeline_dir: directory to test output
        type_: one of 'full', 'embed', 'knn'
        method: name of method, for saving file
    """
    unintegrated = pipeline_dir["adata_raw"]
    integrated = os.path.join(pipeline_dir["output_dir"], "integrated", f"{method}.h5ad")
    metrics_dir = os.path.join(pipeline_dir["output_dir"], "metrics")

    if not os.path.exists(metrics_dir):
        os.makedirs(metrics_dir)

    adata_u.write(unintegrated)
    adata_i.write(integrated)

    cmd = ["python", script, "-u", unintegrated, "-i", integrated,
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
