import numpy as np
import pandas as pd

from ..preprocessing import score_cell_cycle
from ..utils import check_adata
from .pcr import pc_regression


def cell_cycle(
    adata_pre,
    adata_post,
    batch_key,
    embed=None,
    agg_func=np.mean,
    organism="mouse",
    n_comps=50,
    verbose=False,
    recompute_cc=True,
    precompute_pcr_key=None,
):
    """Cell cycle conservation score

    Compare the variance contribution of S-phase and G2/M-phase cell cycle scores before and
    after integration. Cell cycle scores are computed per batch on the unintegrated data set,
    eliminating the batch effect confounded by the ``batch_key`` variable.

    .. math::

        CC \\, conservation = 1 - \\frac { |Var_{after} - Var_{before}| } {Var_{before}}

    Variance contribution is obtained through principal component regression using :func:`~scib.metrics.pc_regression`.
    The score can be computed on full corrected feature spaces and latent embeddings.

    :param adata_pre: adata before integration
    :param adata_post: adata after integration
    :param embed: Name of embedding in adata_post.obsm.
        If ``embed=None``, use the full expression matrix (``adata.X``), otherwise use the
        embedding provided in ``adata_post.obsm[embed]``
    :param agg_func: any function that takes a list of numbers and aggregates them into
        a single value. If ``agg_func=None``, all results will be returned
    :param organism: 'mouse' or 'human' for choosing cell cycle genes
    :param recompute_cc: If True, force recompute cell cycle score, otherwise use
        precomputed scores if available as 'S_score' and 'G2M_score' in adata.obs
    :param precompute_pcr_key: Key in adata_pre for precomputed PCR values for cell
        cycle scores. Ignores cell cycle scores in adata_pre if present.

    :return:
        A score between 1 and 0. The larger the score, the stronger the cell cycle
        variance is conserved.
    """
    check_adata(adata_pre)
    check_adata(adata_post)

    if embed == "X_pca":
        embed = None

    batches = adata_pre.obs[batch_key].unique()
    scores_final = []
    scores_before = []
    scores_after = []

    recompute_cc = (
        recompute_cc
        or "S_score" not in adata_pre.obs_keys()
        or "G2M_score" not in adata_pre.obs_keys()
    )
    recompute_pcr = (
        precompute_pcr_key is None or precompute_pcr_key not in adata_pre.uns_keys()
    )

    for batch in batches:
        before, after = get_pcr_before_after(
            adata_pre,
            adata_post,
            batch_key=batch_key,
            batch=batch,
            embed=embed,
            organism=organism,
            recompute_cc=recompute_cc,
            recompute_pcr=recompute_pcr,
            pcr_key=precompute_pcr_key,
            n_comps=n_comps,
            verbose=verbose,
        )

        # scale result
        score = 1 - abs(after - before) / before

        if score < 0:
            # Here variance contribution becomes more than twice as large as before
            if verbose:
                print(
                    "Variance contrib more than twice as large after integration.\n"
                    "Setting cell cycle score to 0."
                )
            score = 0

        if verbose:
            print(
                f"batch: {batch}\t before: {before}\t after: {after}\t score: {score}"
            )

        scores_before.append(before)
        scores_after.append(after)
        scores_final.append(score)

    if agg_func is None:
        return pd.DataFrame(
            [batches, scores_before, scores_after, scores_final],
            columns=["batch", "before", "after", "score"],
        )
    else:
        return agg_func(scores_final)


def get_pcr_before_after(
    adata_pre,
    adata_post,
    batch_key,
    batch,
    embed,
    organism,
    recompute_cc,
    recompute_pcr,
    pcr_key,
    n_comps,
    verbose,
):
    """
    Principle component regression value on cell cycle scores for one batch
    before and after integration

    :param adata_pre: adata before integration
    :param adata_post: adata after integration
    :param embed: Name of embedding in adata_post.obsm.
        If ``embed=None``, use the full expression matrix (``adata.X``), otherwise use the
        embedding provided in ``adata_post.obsm[embed]``
    :param organism: 'mouse' or 'human' for choosing cell cycle genes
    :param recompute_cc: If True, force recompute cell cycle score, otherwise use
        precomputed scores if available as 'S_score' and 'G2M_score' in adata.obs
    :param recompute_pcr: Whether to recompute principal component regression.
    :param pcr_key: Key in adata_pre for precomputed PCR values for cell cycle scores.
    :param n_comps: Number of components for PCA
    :param verbose:

    :return:
        A score between 1 and 0. The larger the score, the stronger the cell cycle
        variance is conserved.
    """
    # subset adata objects
    raw_sub = adata_pre[adata_pre.obs[batch_key] == batch]
    int_sub = adata_post[adata_post.obs[batch_key] == batch].copy()
    int_sub = int_sub.obsm[embed] if embed is not None else int_sub.X

    # sanity check: subsets have same number of rows?
    if raw_sub.shape[0] != int_sub.shape[0]:
        raise ValueError(
            f'batch "{batch}" of batch_key "{batch_key}" has unequal number of '
            f"entries before and after integration.\n"
            f"before: {raw_sub.shape[0]} after: {int_sub.shape[0]}"
        )

    # compute cc scores if necessary
    if recompute_cc:
        if verbose:
            print("Score cell cycle...")
        score_cell_cycle(raw_sub, organism=organism)

    # regression variable
    covariate = raw_sub.obs[["S_score", "G2M_score"]]

    # PCR on adata before integration
    if recompute_pcr:
        before = pc_regression(
            raw_sub.X, covariate, pca_var=None, n_comps=n_comps, verbose=verbose
        )
    else:
        before = pd.Series(raw_sub.uns[pcr_key])

    # PCR on adata after integration
    after = pc_regression(
        int_sub, covariate, pca_var=None, n_comps=n_comps, verbose=verbose
    )

    return before, after
