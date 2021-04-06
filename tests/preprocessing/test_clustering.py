from tests.common import *


def test_cluster(adata_neighbors):
    _, _, score_all, clustering = scIB.cl.opt_louvain(
        adata_neighbors,
        label_key='celltype',
        cluster_key='cluster',
        plot=False,
        inplace=False
    )
    assert isinstance(score_all, pd.DataFrame)
    assert isinstance(clustering, pd.Series)
