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

    This function can be applied to all integration output types.
    See below for examples of preproceassing and function calls.

    :param adata: adata with computed neighborhood graph
    :param label_key: name in adata.obs containing the cell identity labels

    **Preprocessing: Feature output**

    Feature output requires processing of the count matrix in the following steps:

        1. Highly variable gene selection (skip, if working on feature space subset)
        2. PCA
        3. kNN graph

    .. code-block:: python

        scib.pp.reduce_data(adata, n_top_genes=2000, pca=True, neighbors=True)
        scib.me.graph_connectivity(adata, label_key="celltype")

    **Preprocessing Embedding output**

    The embedding should be stored in ``adata.obsm``, by default under key ``'X_emb'``.
    KNN graph computation must use

    .. code-block:: python

        scib.pp.reduce_data(adata, pca=False, neighbors=True, use_rep="X_emb")
        scib.me.graph_connectivity(adata, label_key="celltype")

    **Preprocessing: kNN graph output**

    No preprocessing required.
    The kNN graph is stored under ``adata.uns['neighbors']`` and will be used if ``embed`` is set to ``None``.

    .. code-block:: python

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
