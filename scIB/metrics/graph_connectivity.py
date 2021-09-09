import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components


def graph_connectivity(adata_post, label_key):
    """"
    Metric that quantifies how connected the subgraph corresponding to each batch cluster is.
    """
    if 'neighbors' not in adata_post.uns:
        raise KeyError('Please compute the neighborhood graph before running this '
                       'function!')

    clust_res = []

    for ct in adata_post.obs[label_key].cat.categories:
        adata_post_sub = adata_post[adata_post.obs[label_key].isin([ct]),]
        _,labs = connected_components(
            adata_post_sub.obsp['connectivities'],
            connection='strong'
        )
        tab = pd.value_counts(labs)
        clust_res.append(tab[0]/sum(tab))

    return np.mean(clust_res)


