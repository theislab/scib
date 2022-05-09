import pandas as pd
from sklearn.metrics import f1_score

from .clustering import opt_louvain
from .silhouette import silhouette


def isolated_labels(
    adata,
    label_key,
    batch_key,
    embed,
    cluster=True,
    iso_threshold=None,
    return_all=False,
    verbose=True,
):
    """Isolated label score

    Score how well labels of isolated labels are distiguished in the dataset by either

        1. clustering-based approach F1 score, or
        2. average-width silhouette score (ASW) on isolated label vs all other labels

    :param adata: anndata object
    :param label_key: column in adata.obs
    :param batch_key: column in adata.obs
    :param cluster: if True, use clustering approach, otherwise use silhouette score approach
    :param embed: key in adata.obsm used for silhouette score if cluster=False, or
        as representation for clustering (if neighbors missing in adata)
    :param iso_threshold: max number of batches per label for label to be considered as
        isolated, if iso_threshold is integer.
        If iso_threshold=None, consider minimum number of batches that labels are present in
    :param return_all: return scores for all isolated labels instead of aggregated mean
    :param verbose:
    :return:
        Mean of scores for each isolated label
        or dictionary of scores for each label if `return_all=True`
    """
    scores = {}
    isolated_labels = get_isolated_labels(
        adata, label_key, batch_key, iso_threshold, verbose
    )

    for label in isolated_labels:
        score = score_isolated_label(
            adata, label_key, label, embed, cluster, verbose=verbose
        )
        scores[label] = score
    scores = pd.Series(scores)

    if return_all:
        return scores

    return scores.mean()


def score_isolated_label(
    adata,
    label_key,
    isolated_label,
    embed,
    cluster=True,
    iso_label_key="iso_label",
    verbose=False,
):
    """
    Compute label score for a single label

    :param adata: anndata object
    :param label_key: key in adata.obs of isolated label type (usually cell label)
    :param isolated_label: value of specific isolated label e.g. cell type/identity annotation
    :param embed: embedding to be passed to opt_louvain, if adata.uns['neighbors'] is missing
    :param cluster: if True, compute clustering-based F1 score, otherwise compute
        silhouette score on grouping of isolated label vs all other remaining labels
    :param iso_label_key: name of key to use for cluster assignment for F1 score or
        isolated-vs-rest assignment for silhouette score
    :param verbose:
    :return:
        Isolated label score
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
            label=isolated_label,
            use_rep=embed,
            function=max_f1,
            verbose=False,
            inplace=True,
        )
        score = max_f1(
            adata_tmp, label_key, iso_label_key, isolated_label, argmax=False
        )
    else:
        # AWS score between label
        adata_tmp.obs[iso_label_key] = adata_tmp.obs[label_key] == isolated_label
        score = silhouette(adata_tmp, iso_label_key, embed)

    del adata_tmp

    if verbose:
        print(f"{isolated_label}: {score}")

    return score


def get_isolated_labels(adata, label_key, batch_key, iso_threshold, verbose):
    """
    Get labels that are isolated depending on the number of batches

    :param adata: anndata object
    :param label_key: column in adata.obs
    :param batch_key: column in adata.obs
    :param iso_threshold: Maximum number of batches per label for label to be considered as isolated.
    :param verbose:
    """

    tmp = adata.obs[[label_key, batch_key]].drop_duplicates()
    batch_per_lab = tmp.groupby(label_key).agg({batch_key: "count"})

    # threshold for determining when label is considered isolated
    if iso_threshold is None:
        iso_threshold = batch_per_lab.min().tolist()[0]

    if verbose:
        print(f"isolated labels: no more than {iso_threshold} batches per label")

    labels = batch_per_lab[batch_per_lab[batch_key] <= iso_threshold].index.tolist()
    if len(labels) == 0 and verbose:
        print(f"no isolated labels with less than {iso_threshold} batches")

    return labels
