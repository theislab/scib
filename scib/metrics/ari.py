import numpy as np
import pandas as pd
import scipy.special
from sklearn.metrics.cluster import adjusted_rand_score

try:
    from scanpy._utils import renamed_arg
except ImportError:
    from .._package_tools import renamed_arg

from ..utils import check_adata, check_batch


@renamed_arg("group1", "cluster_key")
@renamed_arg("group2", "label_key")
def ari(adata, cluster_key, label_key, implementation=None):
    """Adjusted Rand Index

    The adjusted rand index is a chance-adjusted rand index, which evaluates the pair-wise accuracy of clustering vs.
    ground truth label assignments.
    The score ranges between 0 and 1 with larger values indicating better conservation of the data-driven cell identity
    discovery after integration compared to annotated labels.

    :param adata: anndata object with cluster assignments in ``adata.obs[cluster_key]``
    :param cluster_key: string of column in adata.obs containing cluster assignments
    :param label_key: string of column in adata.obs containing labels
    :param implementation: if set to 'sklearn', uses sklearn's implementation,
        otherwise native implementation is taken

    This function can be applied to all integration output types.
    The ``adata`` must contain cluster assignments that are based off the knn graph given or derived from the integration
    method output.
    For this metric you need to include all steps that are needed for clustering.
    See :ref:`preprocessing` for more information on preprocessing.

    **Examples**

    .. code-block:: python

        # feature output
        scib.pp.reduce_data(
            adata, n_top_genes=2000, batch_key="batch", pca=True, neighbors=True
        )
        scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", label_key="celltype")
        scib.me.ari(adata, cluster_key="cluster", label_key="celltype")

        # embedding output
        sc.pp.neighbors(adata, use_rep="X_emb")
        scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", label_key="celltype")
        scib.me.ari(adata, cluster_key="cluster", label_key="celltype")

        # knn output
        scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", label_key="celltype")
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
