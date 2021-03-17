from tests.common import *


def test_cluster(adata):
    scIB.pp.reduce_data(adata, pca=True, n_top_genes=200, neighbors=True, umap=False)
    _, _, score_all, clustering = scIB.cl.opt_louvain(
        adata,
        label_key='celltype',
        cluster_key='cluster',
        plot=False,
        inplace=False
    )
    assert isinstance(score_all, pd.DataFrame)
    assert isinstance(clustering, pd.Series)
