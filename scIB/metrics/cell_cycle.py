import numpy as np
import pandas as pd

from .pcr import pc_regression
from scIB.utils import checkAdata
from scIB.preprocessing import score_cell_cycle


def precompute_cc_score(
        adata,
        batch_key,
        organism='mouse',
        n_comps=50,
        verbose=False
):
    batches = adata.obs[batch_key].cat.categories
    scores_before = {}
    s_score = []
    g2m_score = []

    for batch in batches:
        raw_sub = adata[adata.obs[batch_key] == batch].copy()
        # score cell cycle if not already done
        if (np.in1d(['S_score', 'G2M_score'], adata.obs_keys()).sum() < 2):
            score_cell_cycle(raw_sub, organism=organism)
            s_score.append(raw_sub.obs['S_score'])
            g2m_score.append(raw_sub.obs['G2M_score'])

        covariate = raw_sub.obs[['S_score', 'G2M_score']]

        before = pc_regression(raw_sub.X, covariate, pca_var=None, n_comps=n_comps, verbose=verbose)
        scores_before.update({batch: before})

    if (np.in1d(['S_score', 'G2M_score'], adata.obs_keys()).sum() < 2):
        adata.obs['S_score'] = pd.concat(s_score)
        adata.obs['G2M_score'] = pd.concat(g2m_score)
    adata.uns['scores_before'] = scores_before
    return


def cell_cycle(
        adata_pre,
        adata_post,
        batch_key,
        embed=None,
        agg_func=np.mean,
        organism='mouse',
        n_comps=50,
        verbose=False,
        recompute_cc=True,
        precompute_pcr_key=None
):
    """
    Compare the variance contribution of S-phase and G2/M-phase cell cycle scores before and
    after integration. Cell cycle scores are computed per batch on the unintegrated data set,
    eliminatimg the batch effect confounded by the `batch_key` variable. This function
    returns a score between 1 and 0. The larger the score, the stronger the cell cycle
    variance is conserved.
    This score can be calculated on full corrected feature spaces and latent embeddings as
    variance contributions of a fixed score can be obtained via PC regression here.
    params:
        adata_pre, adata_post: adatas before and after integration
        embed   : if `embed=None`, use the full expression matrix (`adata.X`), otherwise
                  use the embedding provided in `adata_post.obsm[embed]`
        agg_func: any function that takes a list of numbers and aggregates them into a single number.
                  If `agg_func=None`, all results will be returned
        organism: 'mouse' or 'human' for choosing cell cycle genes
        recompute_cc: whether to force recompute cell cycle score (True) or use precomputed scores (False)
        recompute_pcr_key: key in adata_pre for precomputed PCR values for cell cycle scores (does not require
            cell cycle scores to be present in adata_pre)
    """
    checkAdata(adata_pre)
    checkAdata(adata_post)

    if embed == 'X_pca':
        embed = None

    batches = adata_pre.obs[batch_key].unique()
    scores_final = []
    scores_before = []
    scores_after = []

    recompute_cc = recompute_cc \
                   or 'S_score' not in adata_pre.obs_keys() \
                   or 'G2M_score' not in adata_pre.obs_keys()
    recompute_pcr = precompute_pcr_key is None \
                    or precompute_pcr_key not in adata_pre.uns_keys()

    for batch in batches:
        # subset adata objects
        raw_sub = adata_pre[adata_pre.obs[batch_key] == batch]
        int_sub = adata_post[adata_post.obs[batch_key] == batch].copy()
        int_sub = int_sub.obsm[embed] if embed is not None else int_sub.X

        # sanity check: subsets have same number of rows?
        if raw_sub.shape[0] != int_sub.shape[0]:
            message = f'batch "{batch}" of batch_key "{batch_key}" '
            message += 'has unequal number of entries before and after integration.'
            message += f'before: {raw_sub.shape[0]} after: {int_sub.shape[0]}'
            raise ValueError(message)

        # compute cc scores if necessary
        if recompute_cc:
            if verbose:
                print("score cell cycle")
            score_cell_cycle(raw_sub, organism=organism)

        # regression variable
        covariate = raw_sub.obs[['S_score', 'G2M_score']]

        # PCR on adata before integration
        before = pc_regression(raw_sub.X, covariate, pca_var=None, n_comps=n_comps, verbose=verbose) \
            if recompute_pcr else pd.Series(raw_sub.uns[precompute_pcr_key])
        scores_before.append(before)

        # PCR on adata after integration
        after = pc_regression(int_sub, covariate, pca_var=None, n_comps=n_comps, verbose=verbose)
        scores_after.append(after)

        # scaled result
        score = 1 - abs(after - before) / before
        if score < 0:
            # Here variance contribution becomes more than twice as large as before
            if verbose:
                print("Variance contrib more than twice as large after integration.")
                print("Setting score to 0.")
            score = 0
        scores_final.append(score)

        if verbose:
            print(f"batch: {batch}\t before: {before}\t after: {after}\t score: {score}")

    if agg_func is None:
        return pd.DataFrame([batches, scores_before, scores_after, scores_final],
                            columns=['batch', 'before', 'after', 'score'])
    else:
        return agg_func(scores_final)
