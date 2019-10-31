import scanpy as sc
from scIB import clustering as cl
from scIB.tests import utils
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def test_cluster():
    adata = utils.create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
    _, _, score_all = cl.opt_louvain(adata,
                                  label_key='celltype',
                                  cluster_key='cluster',
                                  plot=False, inplace=True)
    assert isinstance(score_all, pd.DataFrame)
    _, _, _, clustering = cl.opt_louvain(adata,
                                      label_key='celltype',
                                      cluster_key='cluster',
                                      plot=False, inplace=False)
    assert isinstance(clustering, pd.Series)
    
    # see if this runs through
    cl.opt_louvain(adata, label_key='celltype', cluster_key='cluster',
                plot=True, inplace=True)
