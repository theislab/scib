from tests.common import *


def test_graph_connectivity(adata_neighbors):
    score = scIB.me.graph_connectivity(adata_neighbors, label_key='celltype')
    LOGGER.info(f"score: {score}")
    assert score == 0.9670013350457753
