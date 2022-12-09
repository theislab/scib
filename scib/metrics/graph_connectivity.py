import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components


def graph_connectivity(adata, label_key):
    """Graph Connectivity

    Quantify the connectivity of the subgraph per cell type label.
    The final score is the average for all cell type labels :math:`C`, according to the equation:

    .. math::

        GC = \\frac {1} {|C|} \\sum_{c \\in C} \\frac {|{LCC(subgraph_c)}|} {|c|}

    where :math:`|LCC(subgraph_c)|` stands for all cells in the largest connected component and :math:`|c|` stands for all cells of
    cell type :math:`c`.

    :param adata: integrated adata with computed neighborhood graph
    :param label_key: name in adata.obs containing the cell identity labels

    This function can be applied to all integration output types.
    The integrated object (``adata``) needs to have a kNN graph based on the integration output.
    See :ref:`preprocessing` for more information on preprocessing.

    **Examples**

    .. code-block:: python

        # feature output
        scib.pp.reduce_data(
            adata, n_top_genes=2000, batch_key="batch", pca=True, neighbors=True
        )
        scib.me.graph_connectivity(adata, label_key="celltype")

        # embedding output
        sc.pp.neighbors(adata, use_rep="X_emb")
        scib.me.graph_connectivity(adata, label_key="celltype")

        # knn output
        scib.me.graph_connectivity(adata, label_key="celltype")

    """
    if "neighbors" not in adata.uns:
        raise KeyError(
            "Please compute the neighborhood graph before running this function!"
        )

    clust_res = []

    for label in adata.obs[label_key].cat.categories:
        adata_sub = adata[adata.obs[label_key].isin([label])]
        _, labels = connected_components(
            adata_sub.obsp["connectivities"], connection="strong"
        )
        tab = pd.value_counts(labels)
        clust_res.append(tab.max() / sum(tab))

    return np.mean(clust_res)
