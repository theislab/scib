import numpy as np
from sklearn.metrics import f1_score

from scIB.clustering import opt_louvain
from .silhouette import silhouette


def isolated_labels(
        adata,
        label_key,
        batch_key,
        embed,
        cluster=True,
        n=None,
        all_=False,
        verbose=True
):
    """
    score how well labels of isolated labels are distiguished in the dataset by
        1. clustering-based approach F1 score
        2. average-width silhouette score on isolated-vs-rest label assignment
    params:
        cluster: if True, use clustering approach, otherwise use silhouette score approach
        embed: key in adata.obsm used for silhouette score if cluster=False, or
            as representation for clustering (if neighbors missing in adata)
        n: max number of batches per label for label to be considered as isolated.
            if n is integer, consider labels that are present for n batches as isolated
            if n=None, consider minimum number of batches that labels are present in
        all_: return scores for all isolated labels instead of aggregated mean
    return:
        by default, mean of scores for each isolated label
        retrieve dictionary of scores for each label if `all_` is specified
    """

    scores = {}
    isolated_labels = get_isolated_labels(
        adata,
        label_key,
        batch_key,
        n,
        verbose
    )
    for label in isolated_labels:
        score = score_isolated_label(
            adata,
            label_key,
            label,
            embed,
            cluster,
            verbose=verbose
        )
        scores[label] = score

    if all_:
        return scores
    return np.mean(list(scores.values()))


def get_isolated_labels(adata, label_key, batch_key, n, verbose):
    """
    get labels that are considered isolated by the number of batches
    """

    tmp = adata.obs[[label_key, batch_key]].drop_duplicates()
    batch_per_lab = tmp.groupby(label_key).agg({batch_key: "count"})

    # threshold for determining when label is considered isolated
    if n is None:
        n = batch_per_lab.min().tolist()[0]

    if verbose:
        print(f"isolated labels: no more than {n} batches per label")

    labels = batch_per_lab[batch_per_lab[batch_key] <= n].index.tolist()
    if len(labels) == 0 and verbose:
        print(f"no isolated labels with less than {n} batches")
    return labels


def score_isolated_label(
        adata,
        label_key,
        label,
        embed,
        cluster=True,
        iso_label_key='iso_label',
        verbose=False
):
    """
    compute label score for a single label
    params:
        adata: anndata object
        label_key: key in adata.obs of isolated label type (usually cell label)
        label: value of specific isolated label e.g. cell type/identity annotation
        embed: embedding to be passed to opt_louvain, if adata.uns['neighbors'] is missing
        cluster: if True, compute clustering-based F1 score, otherwise compute
            silhouette score on grouping of isolated label vs all other remaining labels
        iso_label_key: name of key to use for cluster assignment for F1 score or
            isolated-vs-rest assignment for silhouette score
    """
    adata_tmp = adata.copy()

    def max_f1(adata, label_key, cluster_key, label, argmax=False):
        """cluster optimizing over largest F1 score of isolated label"""
        obs = adata.obs
        max_cluster = None
        max_f1 = 0
        for cluster in obs[cluster_key].unique():
            y_pred = obs[cluster_key] == cluster
            y_true = obs[label_key] == label
            f1 = f1_score(y_pred, y_true)
            if f1 > max_f1:
                max_f1 = f1
                max_cluster = cluster
        if argmax:
            return max_cluster
        return max_f1

    if cluster:
        # F1-score on clustering
        opt_louvain(
            adata_tmp,
            label_key,
            cluster_key=iso_label_key,
            label=label,
            use_rep=embed,
            function=max_f1,
            verbose=False,
            inplace=True
        )
        score = max_f1(adata_tmp, label_key, iso_label_key, label, argmax=False)
    else:
        # AWS score between label
        adata_tmp.obs[iso_label_key] = adata_tmp.obs[label_key] == label
        score = silhouette(adata_tmp, iso_label_key, embed)

    del adata_tmp

    if verbose:
        print(f"{label}: {score}")

    return score
