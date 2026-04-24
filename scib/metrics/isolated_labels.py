import warnings

import pandas as pd
from sklearn.metrics import f1_score, silhouette_samples

from .clustering import cluster_optimal_resolution


def isolated_labels_f1(
    adata,
    label_key,
    batch_key,
    embed,
    cluster_key="iso_label",
    resolutions=None,
    iso_threshold=None,
    verbose=True,
    **kwargs,
):
    """Isolated label score F1

    Score how well isolated labels are distinguished from other labels by data-driven clustering.
    The F1 score is used to evaluate clustering with respect to the ground truth labels.

    :param adata: anndata object
    :param label_key: column in ``adata.obs``
    :param batch_key: column in ``adata.obs``
    :param embed: key in adata.obsm used for as representation for kNN graph computation.
        If ``embed=None``, use the existing kNN graph in ``adata.uns['neighbors']``.
    :param iso_threshold: max number of batches per label for label to be considered as
        isolated, if iso_threshold is integer.
        If ``iso_threshold=None``, consider minimum number of batches that labels are present in
    :param cluster_key: clustering key prefix to look or recompute for each resolution in resolutions.
        Is passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :param resolutions: list of resolutions to be passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :param verbose:
    :params \\**kwargs: additional arguments to be passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :return: Mean of F1 scores over all isolated labels

    This function performs clustering on a kNN graph and can be applied to all integration output types.
    For this metric the ``adata`` needs a kNN graph and can optionally make use of precomputed clustering (see example below).
    The precomputed clusters must be saved under ``adata.obs[cluster_key]`` as well as ``adata.obs[f"{cluster_key}_{resolution}"]`` for all resolutions.

    See :ref:`preprocessing` for more information on preprocessing.

    **Examples**

    .. code-block:: python

        # full feature output
        scib.pp.reduce_data(
            adata, n_top_genes=2000, batch_key="batch", pca=True, neighbors=True
        )
        scib.me.isolated_labels_f1(adata, label_key="celltype")

        # embedding output
        sc.pp.neighbors(adata, use_rep="X_emb")
        scib.me.isolated_labels_f1(adata, batch_key="batch", label_key="celltype")

        # knn output
        scib.me.isolated_labels_f1(adata, batch_key="batch", label_key="celltype")

        # use precomputed clustering
        scib.cl.cluster_optimal_resolution(adata, cluster_key="iso_label", label_key="celltype")
        scib.me.isolated_labels_f1(adata, batch_key="batch", label_key="celltype")

        # overwrite existing clustering
        scib.me.isolated_labels_f1(adata, batch_key="batch", label_key="celltype", force=True)

    """
    return isolated_labels(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed=embed,
        cluster=True,
        cluster_key=cluster_key,
        resolutions=resolutions,
        iso_threshold=iso_threshold,
        verbose=verbose,
        **kwargs,
    )


def isolated_labels_asw(
    adata,
    label_key,
    batch_key,
    embed,
    iso_threshold=None,
    scale=True,
    verbose=True,
):
    """Isolated label score ASW

    Score how well isolated labels are distinguished from all other labels using the average-width silhouette score
    (ASW) :func:`~scib.metrics.silhouette`.

    :param adata: anndata object
    :param label_key: column in ``adata.obs``
    :param batch_key: column in ``adata.obs``
    :param embed: key in ``adata.obsm`` used for `func:~scib.metrics.silhouette`
    :param iso_threshold: max number of batches per label for label to be considered as
        isolated, if iso_threshold is integer.
        If ``iso_threshold=None``, consider minimum number of batches that labels are present in
    :param scale: Whether to scale the score between 0 and 1. Only relevant for ASW scores.
    :param verbose:
    :params \\**kwargs: additional arguments to be passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :return: Mean of ASW over all isolated labels

    The function requires an embedding to be stored in ``adata.obsm`` and can only be applied to feature and embedding
    integration outputs.
    Please note, that the metric cannot be used to evaluate kNN graph outputs.
    See :ref:`preprocessing` for more information on preprocessing.

    **Examples**

    .. code-block:: python

        # full feature output
        scib.pp.reduce_data(
            adata, n_top_genes=2000, batch_key="batch", pca=True, neighbors=False
        )
        scib.me.isolated_labels_asw(adata, label_key="celltype", embed="X_pca")

        # embedding output
        scib.me.isolated_labels_asw(
            adata, batch_key="batch", label_key="celltype", embed="X_emb"
        )

    """
    return isolated_labels(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed=embed,
        cluster=False,
        iso_threshold=iso_threshold,
        scale=scale,
        verbose=verbose,
    )


def isolated_labels(
    adata,
    label_key,
    batch_key,
    embed,
    cluster=True,
    cluster_key="iso_label",
    resolutions=None,
    iso_threshold=None,
    scale=True,
    return_all=False,
    verbose=True,
    **kwargs,
):
    """Isolated label score

    Score how well labels of isolated labels are distinguished in the dataset by either

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
    :param cluster_key: name of key to be passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :param resolutions: list of resolutions to be passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :param scale: Whether to scale the score between 0 and 1. Only relevant for ASW scores.
    :param return_all: return scores for all isolated labels instead of aggregated mean
    :param verbose:
    :params \\**kwargs: additional arguments to be passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :return:
        Mean of scores for each isolated label
        or dictionary of scores for each label if `return_all=True`
    """

    # 1. determine isolated labels
    isolated_labels = get_isolated_labels(
        adata, label_key, batch_key, iso_threshold, verbose
    )
    if verbose:
        print(f"isolated labels: {isolated_labels}")

    # 2. compute isolated label score for each isolated label
    scores = {}
    if not cluster:
        adata.obs["silhouette_temp"] = silhouette_samples(
            adata.obsm[embed], adata.obs[label_key]
        )
    for label in isolated_labels:
        score = score_isolated_label(
            adata,
            label_key,
            label,
            embed,
            cluster_key=cluster_key,
            cluster=cluster,
            scale=scale,
            verbose=verbose,
            resolutions=resolutions,
            **kwargs,
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
    cluster_key,
    cluster=True,
    resolutions=None,
    scale=True,
    verbose=False,
    **kwargs,
):
    """
    Compute label score for a single label

    :param adata: anndata object
    :param label_key: key in adata.obs of isolated label type (usually cell label)
    :param isolated_label: value of specific isolated label e.g. cell type/identity annotation
    :param embed: embedding to be passed to opt_louvain, if adata.uns['neighbors'] is missing
    :param cluster: if True, compute clustering-based F1 score, otherwise compute
        silhouette score on grouping of isolated label vs all other remaining labels
    :param cluster_key: name of key to use for cluster assignment for F1 score or
        isolated-vs-rest assignment for silhouette score
    :param resolutions: list of resolutions to be passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :param scale: Whether to scale the score between 0 and 1. Only relevant for ASW scores.
    :param verbose:
    :params \\**kwargs: additional arguments to be passed to :func:`~scib.metrics.cluster_optimal_resolution`
    :return:
        Isolated label score
    """
    adata = adata.copy()

    if cluster:
        # F1-score on clustering
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

        cluster_optimal_resolution(
            adata,
            label_key,
            cluster_key=cluster_key,
            use_rep=embed,
            metric=max_f1,
            metric_kwargs={"label": isolated_label},
            resolutions=resolutions,
            force=False,
            verbose=verbose,
        )
        score = max_f1(adata, label_key, cluster_key, isolated_label, argmax=False)
    else:
        # AWS score between isolated label vs rest
        if "silhouette_temp" not in adata.obs:
            adata.obs["silhouette_temp"] = silhouette_samples(
                adata.obsm[embed], adata.obs[label_key]
            )
        # aggregate silhouette scores for isolated label only
        score = adata.obs[adata.obs[label_key] == isolated_label].silhouette_temp.mean()
        if scale:
            score = (score + 1) / 2

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

    if iso_threshold == adata.obs[batch_key].nunique():
        warnings.warn(
            "iso_threshold is equal to number of batches in data, no isolated labels will be found",
            stacklevel=2,
        )
        return []

    if verbose:
        print(f"isolated labels: no more than {iso_threshold} batches per label")

    labels = batch_per_lab[batch_per_lab[batch_key] <= iso_threshold].index.tolist()
    if len(labels) == 0 and verbose:
        print(f"no isolated labels with less than {iso_threshold} batches")

    return labels
