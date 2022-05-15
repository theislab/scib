import numpy as np
import pandas as pd
import scipy.special
from sklearn.metrics.cluster import adjusted_rand_score

from ..utils import check_adata, check_batch


def ari(adata, cluster_key, label_key, implementation=None):
    """Adjusted Rand Index

    The adjusted rand index is a chance-adjusted rand index, which evaluates the pair-wise accuracy of clustering vs.
    ground truth label assignments.
    The score ranges between

    This metric is useful for evaluating the clustering performance of cells with respect to annotated cell identity
    labels and can be applied to all integration output types.
    For the regular integration benchmark use-case, the metric is applied to the integrated data.
    A larger value indicates better conservation of data-driven cell identity discovery and after integration, based on
    previous annotation.
    The ``adata`` must contain cluster assignments that are based off the knn graph given or derived from the integration
    method output.
    Preprocessing the ``anndata`` object is required and differs, depending on the output type of the integration method.
    See below for more information on preprocessing.

    :param adata: anndata object with cluster assignments in ``adata.obs[cluster_key]``
    :param cluster_key: string of column in adata.obs containing cluster assignments
    :param label_key: string of column in adata.obs containing labels
    :param implementation: if set to 'sklearn', uses sklearn's implementation,
        otherwise native implementation is taken

    **Preprocessing: Feature output**

    Feature output requires processing of the count matrix in the following steps:

        1. Highly variable gene selection (skip, if working on feature space subset)
        2. PCA
        3. kNN graph
        4. Clustering (optimised resolution recommended)

    .. code-block:: python

        scib.pp.reduce_data(adata, n_top_genes=2000, pca=True, neighbors=True)
        scib.me.opt_louvain(adata, cluster_key="cluster", label_key="celltype")

    **Preprocessing Embedding output**

    The embedding should be stored in ``adata.obsm``, by default under key ``'X_embed'``

        1. kNN graph
        2. Clustering (optimised resolution recommended)

    .. code-block:: python

        scib.pp.reduce_data(
            adata, pca=False, use_rep="X_embed", neighbors=True, use_rep="X_emb"
        )
        scib.me.opt_louvain(adata, cluster_key="cluster", label_key="celltype")


    **Preprocessing: kNN graph output**

    KNN graph output only requires clustering on the (optimised resolution recommended)

    .. code-block:: python

        scib.me.opt_louvain(adata, cluster_key="cluster", label_key="celltype")

    **Call function**

    Once the ``adata`` has been preprocessed according to its integration output type, call the metric as follows:

    .. code-block:: python

        scib.me.ari(adata, cluster_key="cluster", label_key="celltype")
    """

    check_adata(adata)
    check_batch(cluster_key, adata.obs)
    check_batch(label_key, adata.obs)

    cluster_key = adata.obs[cluster_key].to_numpy()
    label_key = adata.obs[label_key].to_numpy()

    if len(cluster_key) != len(label_key):
        raise ValueError(
            f"different lengths in cluster_key ({len(cluster_key)}) and label_key ({len(label_key)})"
        )

    if implementation == "sklearn":
        return adjusted_rand_score(cluster_key, label_key)

    def binom_sum(x, k=2):
        return scipy.special.binom(x, k).sum()

    n = len(cluster_key)
    contingency = pd.crosstab(cluster_key, label_key)

    ai_sum = binom_sum(contingency.sum(axis=0))
    bi_sum = binom_sum(contingency.sum(axis=1))

    index = binom_sum(np.ravel(contingency))
    expected_index = ai_sum * bi_sum / binom_sum(n, 2)
    max_index = 0.5 * (ai_sum + bi_sum)

    return (index - expected_index) / (max_index - expected_index)
