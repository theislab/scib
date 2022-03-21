import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components


def graph_connectivity(adata, label_key):
    """"
    Quantify how connected the subgraph corresponding to each batch cluster is.
    Calculate per label: #cells_in_largest_connected_component/#all_cells
    Final score: Average over labels

    :param adata: adata with computed neighborhood graph
    :param label_key: name in adata.obs containing the cell identity labels
    """
    if 'neighbors' not in adata.uns:
        raise KeyError(
            'Please compute the neighborhood graph before running this function!'
        )

    clust_res = []

    for label in adata.obs[label_key].cat.categories:
        adata_sub = adata[adata.obs[label_key].isin([label])]
        _, labels = connected_components(
            adata_sub.obsp['connectivities'],
            connection='strong'
        )
        tab = pd.value_counts(labels)
        clust_res.append(tab.max() / sum(tab))

    return np.mean(clust_res)
