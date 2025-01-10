import numpy as np
import pandas as pd
from tqdm import tqdm

from ..preprocessing import score_cell_cycle
from ..utils import check_adata
from .pcr import pc_regression


def cell_cycle(
    adata_pre,
    adata_post,
    batch_key,
    embed=None,
    agg_func=np.nanmean,
    organism="mouse",
    n_comps=50,
    recompute_cc=True,
    precompute_pcr_key=None,
    verbose=False,
    linreg_method="numpy",
    n_threads=1,
):
    """Cell cycle conservation score

    Compare the variance contribution of S-phase and G2/M-phase cell cycle scores before and
    after integration. Cell cycle scores are computed per batch on the unintegrated data set,
    eliminating the batch effect confounded by the ``batch_key`` variable.

    .. math::

        CC \\, conservation = 1 - \\frac { |Var_{after} - Var_{before}| } {Var_{before}}

    Variance contribution is obtained through principal component regression using :func:`~scib.metrics.pc_regression`.

    :param adata_pre: adata before integration
    :param adata_post: adata after integration
    :param batch_key: Batch key in ``adata_post.obs``
    :param embed: Name of embedding in ``adata_post.obsm``.
        If ``embed=None``, use the full expression matrix (``adata_post.X``), otherwise use the
        embedding provided in ``adata_post.obsm[embed]``
    :param agg_func: any function that takes a list of numbers and aggregates them into
        a single value. If ``agg_func=None``, all results will be returned
    :param organism: 'mouse' or 'human' for choosing cell cycle genes
    :param n_comps: number of principle components
    :param recompute_cc: If True, force recompute cell cycle score, otherwise use
        precomputed scores if available as 'S_score' and 'G2M_score' in ``adata_post.obs``
    :param precompute_pcr_key: Key in adata_pre for precomputed PCR values for cell
        cycle scores. Ignores cell cycle scores in adata_pre if present.
    :param n_threads: Number of threads for linear regressions per principle component

    :return:
        A score between 1 and 0. The larger the score, the stronger the cell cycle
        variance is conserved.

    The function can be computed on full corrected feature spaces and latent embeddings for both integrated and
    unintegrated ``anndata.Anndata`` objects.
    No preprocessing is needed, as the function will perform PCA directly on the feature or embedding space.

    **Examples**

    .. code-block:: python

        # full feature output
        scib.me.cell_cycle(adata_unintegrated, adata, batch_key="batch")

        # embedding output
        scib.me.cell_cycle(adata_unintegrated, adata, batch_key="batch", embed="X_emb")

    """
    check_adata(adata_pre)
    check_adata(adata_post)

    if embed == "X_pca":
        embed = None

    recompute_cc = (
        recompute_cc
        or "S_score" not in adata_pre.obs_keys()
        or "G2M_score" not in adata_pre.obs_keys()
    )
    recompute_pcr = (
        precompute_pcr_key is None or precompute_pcr_key not in adata_pre.uns_keys()
    )

    batches = adata_pre.obs[batch_key].unique()
    scores_before = []
    scores_after = []
    scores_final = []

    for batch in tqdm(batches):
        before, after = get_pcr_before_after(
            adata_pre,
            adata_post,
            batch_key=batch_key,
            batch=batch,
            embed=embed,
            organism=organism,
            pcr_key=precompute_pcr_key,
            recompute_cc=recompute_cc,
            recompute_pcr=recompute_pcr,
            n_comps=n_comps,
            verbose=verbose,
            n_threads=n_threads,
            linreg_method=linreg_method,
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
            {
                "batch": pd.Series(batches, dtype=str),
                "before": pd.Series(scores_before, dtype=float),
                "after": pd.Series(scores_after, dtype=float),
                "score": pd.Series(scores_final, dtype=float),
            }
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
    pcr_key,
    recompute_cc=False,
    recompute_pcr=False,
    n_comps=50,
    verbose=True,
    n_threads=1,
    linreg_method="numpy",
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
            raw_sub.X,
            covariate,
            pca_var=None,
            n_comps=n_comps,
            verbose=verbose,
            n_threads=n_threads,
            linreg_method=linreg_method,
        )
    else:
        before = pd.Series(raw_sub.uns[pcr_key])

    # PCR on adata after integration
    after = pc_regression(
        int_sub,
        covariate,
        pca_var=None,
        n_comps=n_comps,
        verbose=verbose,
        n_threads=n_threads,
        linreg_method=linreg_method,
    )

    return before, after
