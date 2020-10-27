from scIB import clustering as cl
import pandas as pd
import warnings
warnings.filterwarnings('ignore')


def test_cluster(adata_factory):
    adata = adata_factory(pca=True, n_top_genes=2000, neighbors=True)
    
    _, _, score_all, clustering = cl.opt_louvain(
        adata,
        label_key='celltype',
        cluster_key='cluster',
        plot=False,
        inplace=False
    )
    assert isinstance(score_all, pd.DataFrame)
    assert isinstance(clustering, pd.Series)
