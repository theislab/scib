import warnings

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from deprecated import deprecated

from .nmi import nmi


def get_resolutions(n=20, min=0, max=2):
    """
    Get equally spaced resolutions for optimised clustering

    :param n: number of resolutions
    :param min: minimum resolution
    :param max: maximum resolution

    .. code-block:: python

        scib.cl.get_resolutions(n=10)

    Output:

    .. code-block::

        [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

    """
    res_range = max - min
    return [res_range * (x + 1) / n for x in range(n)]


def cluster_optimal_resolution(
    adata,
    label_key,
    cluster_key,
    cluster_function=None,
    metric=None,
    resolutions=None,
    use_rep=None,
    force=False,
    verbose=True,
    return_all=False,
    metric_kwargs=None,
    **kwargs,
):
    """Optimised clustering

    Leiden, louvain or any custom clustering algorithm with resolution optimised against a metric

    :param adata: anndata object
    :param label_key: name of column in adata.obs containing biological labels to be
        optimised against
    :param cluster_key: name and prefix of columns to be added to adata.obs during clustering.
        Each resolution will be saved under "{cluster_key}_{resolution}", while the optimal clustering will be under ``cluster_key``.
        If ``force=True`` and one of the keys already exists, it will be overwritten.
    :param cluster_function: a clustering function that takes an anndata.Anndata object. Default: Leiden clustering
    :param metric: function that computes the cost to be optimised over. Must take as
        arguments ``(adata, label_key, cluster_key, **metric_kwargs)`` and returns a number for maximising
        Default is :func:`~scib.metrics.nmi()`
    :param resolutions: list of resolutions to be optimised over. If ``resolutions=None``,
        by default 10 equally spaced resolutions ranging between 0 and 2 will be used (see :func:`~scib.metrics.get_resolutions`)
    :param use_rep: key of embedding to use only if ``adata.uns['neighbors']`` is not
        defined, otherwise will be ignored
    :param force: whether to overwrite the cluster assignments in the ``.obs[cluster_key]``
    :param verbose: whether to print out intermediate results
    :param return_all: whether to results for all resolutions
    :param metric_kwargs: arguments to be passed to metric
    :param kwargs: arguments to pass to clustering function
    :returns:
        Only if ``return_all=True``, return tuple of ``(res_max, score_max, score_all)``
        ``res_max``: resolution of maximum score;
        ``score_max``: maximum score;
        ``score_all``: ``pd.DataFrame`` containing all scores at resolutions. Can be used to plot the score profile.

    If you specify an embedding that was not used for the kNN graph (i.e. ``adata.uns["neighbors"]["params"]["use_rep"]`` is not the same as ``use_rep``),
    the neighbors will be recomputed in-place.
    """

    def call_cluster_function(adata, res, resolution_key, cluster_function, **kwargs):
        if resolution_key in adata.obs.columns:
            warnings.warn(
                f"Overwriting existing key {resolution_key} in adata.obs", stacklevel=2
            )

        # check or recompute neighbours
        knn_rep = adata.uns.get("neighbors", {}).get("params", {}).get("use_rep")
        if use_rep is not None and use_rep != knn_rep:
            print(f"Recompute neighbors on rep {use_rep} instead of {knn_rep}")
            sc.pp.neighbors(adata, use_rep=use_rep)

        # call clustering function
        print(f"Cluster for {resolution_key} with {cluster_function.__name__}")
        cluster_function(adata, resolution=res, key_added=resolution_key, **kwargs)

    if cluster_function is None:
        cluster_function = sc.tl.leiden

    if cluster_key is None:
        cluster_key = cluster_function.__name__

    if metric is None:
        metric = nmi

    if metric_kwargs is None:
        metric_kwargs = {}

    if resolutions is None:
        resolutions = get_resolutions(n=10, max=2)

    score_max = 0
    res_max = resolutions[0]
    clustering = None
    score_all = []

    for res in resolutions:
        resolution_key = f"{cluster_key}_{res}"

        # check if clustering exists
        if resolution_key not in adata.obs.columns or force:
            call_cluster_function(
                adata, res, resolution_key, cluster_function, **kwargs
            )

        # score cluster resolution
        score = metric(adata, label_key, resolution_key, **metric_kwargs)
        score_all.append(score)

        if verbose:
            print(f"resolution: {res}, {metric.__name__}: {score}", flush=True)

        # optimise score
        if score_max < score:
            score_max = score
            res_max = res
            clustering = adata.obs[resolution_key]

    if verbose:
        print(f"optimised clustering against {label_key}")
        print(f"optimal cluster resolution: {res_max}")
        print(f"optimal score: {score_max}")

    score_all = pd.DataFrame(
        zip(resolutions, score_all), columns=["resolution", "score"]
    )

    # save optimal clustering in adata.obs
    if cluster_key in adata.obs.columns:
        warnings.warn(
            f"Overwriting existing key {cluster_key} in adata.obs", stacklevel=2
        )
    adata.obs[cluster_key] = clustering

    if return_all:
        return res_max, score_max, score_all
    return res_max, score_max


@deprecated
def opt_louvain(
    adata,
    label_key,
    cluster_key,
    function=None,
    resolutions=None,
    use_rep=None,
    inplace=True,
    plot=False,
    force=True,
    verbose=True,
    **kwargs,
):
    """Optimised Louvain clustering

    DEPRECATED: Use :func:`~scib.metrics.cluster_optimal_resolution` instead

    Louvain clustering with resolution optimised against a metric

    :param adata: anndata object
    :param label_key: name of column in adata.obs containing biological labels to be
        optimised against
    :param cluster_key: name of column to be added to adata.obs during clustering.
        Will be overwritten if exists and ``force=True``
    :param function: function that computes the cost to be optimised over. Must take as
        arguments ``(adata, group1, group2, **kwargs)`` and returns a number for maximising
    :param resolutions: list of resolutions to be optimised over. If ``resolutions=None``,
        default resolutions of 20 values ranging between 0.1 and 2 will be used
    :param use_rep: key of embedding to use only if ``adata.uns['neighbors']`` is not
        defined, otherwise will be ignored
    :returns:
        Tuple of ``(res_max, score_max, score_all)`` or
        ``(res_max, score_max, score_all, clustering)`` if ``inplace=False``.
        ``res_max``: resolution of maximum score;
        ``score_max``: maximum score;
        ``score_all``: ``pd.DataFrame`` containing all scores at resolutions. Can be used to plot the score profile.
        ``clustering``: only if ``inplace=False``, return cluster assignment as ``pd.Series``
    """

    if verbose:
        print("Clustering...")

    if function is None:
        function = nmi

    if cluster_key in adata.obs.columns:
        if force:
            print(
                f"Warning: cluster key {cluster_key} already exists "
                "in adata.obs and will be overwritten"
            )
        else:
            raise ValueError(
                f"cluster key {cluster_key} already exists in "
                + "adata, please remove the key or choose a different name."
                + "If you want to force overwriting the key, specify `force=True`"
            )

    if resolutions is None:
        n = 20
        resolutions = [2 * x / n for x in range(1, n + 1)]

    score_max = 0
    res_max = resolutions[0]
    clustering = None
    score_all = []

    try:
        adata.uns["neighbors"]
    except KeyError:
        if verbose:
            print("computing neighbours for opt_cluster")
        sc.pp.neighbors(adata, use_rep=use_rep)

    for res in resolutions:
        sc.tl.louvain(adata, resolution=res, key_added=cluster_key)
        score = function(adata, label_key, cluster_key, **kwargs)
        if verbose:
            print(f"resolution: {res}, {function.__name__}: {score}")
        score_all.append(score)
        if score_max < score:
            score_max = score
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]

    if verbose:
        print(f"optimised clustering against {label_key}")
        print(f"optimal cluster resolution: {res_max}")
        print(f"optimal score: {score_max}")

    score_all = pd.DataFrame(
        zip(resolutions, score_all), columns=("resolution", "score")
    )
    if plot:
        # score vs. resolution profile
        sns.lineplot(data=score_all, x="resolution", y="score").set_title(
            "Optimal cluster resolution profile"
        )
        plt.show()

    if inplace:
        adata.obs[cluster_key] = clustering
        return res_max, score_max, score_all
    else:
        return res_max, score_max, score_all, clustering
