import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from deprecated import deprecated

from .nmi import nmi


def get_resolutions(n=20, min=0.1, max=2):
    min = np.max([1, int(min * 10)])
    max = np.max([min, max * 10])
    frac = n / 10
    return [frac * x / n for x in range(min, max + 1)]


def cluster_optimal_resolution(
    adata,
    label_key,
    cluster_key,
    cluster_function=None,
    metric=None,
    resolutions=None,
    use_rep=None,
    force=True,
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
    :param cluster_key: name of column to be added to adata.obs during clustering.
        Will be overwritten if exists and ``force=True``
    :param cluster_function: a clustering function that takes an anndata.Anndata object. Default: Leiden clustering
    :param metric: function that computes the cost to be optimised over. Must take as
        arguments ``(adata, group1, group2, **kwargs)`` and returns a number for maximising
        Default is :func:`~scib.metrics.nmi()`
    :param resolutions: list of resolutions to be optimised over. If ``resolutions=None``,
        default resolutions of 10 values ranging between 0.1 and 2 will be used
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
    """
    if cluster_key in adata.obs.columns:
        if force:
            print(
                f"WARNING: cluster key {cluster_key} already exists in adata.obs and will be overwritten because "
                "force=True "
            )
        else:
            raise ValueError(
                f"cluster key {cluster_key} already exists in adata, please remove the key or choose a different "
                "name. If you want to force overwriting the key, specify `force=True` "
            )

    if cluster_function is None:
        cluster_function = sc.tl.leiden

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

    if use_rep is None:
        try:
            adata.uns["neighbors"]
        except KeyError:
            raise RuntimeError(
                "Neighbours must be computed when setting use_rep to None"
            )
    else:
        print(f"Compute neighbors on rep {use_rep}")
        sc.pp.neighbors(adata, use_rep=use_rep)

    for res in resolutions:
        cluster_function(adata, resolution=res, key_added=cluster_key, **kwargs)
        score = metric(adata, label_key, cluster_key, **metric_kwargs)
        if verbose:
            print(f"resolution: {res}, {metric.__name__}: {score}")
        score_all.append(score)

        # optimise score
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
        zip(resolutions, score_all), columns=["resolution", "score"]
    )

    adata.obs[cluster_key] = clustering

    if return_all:
        return res_max, score_max, score_all


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
