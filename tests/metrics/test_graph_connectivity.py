import scib
from tests.common import LOGGER, assert_near_exact


def test_graph_connectivity(adata_neighbors):
    score = scib.me.graph_connectivity(adata_neighbors, label_key="celltype")
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.9670013350457753, diff=1e-3)
