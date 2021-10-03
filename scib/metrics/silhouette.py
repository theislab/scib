import pandas as pd
from sklearn.metrics.cluster import silhouette_samples, silhouette_score


def silhouette(
        adata,
        group_key,
        embed,
        metric='euclidean',
        scale=True
):
    """
    Wrapper for sklearn silhouette function values range from [-1, 1] with
        1 being an ideal fit
        0 indicating overlapping clusters and
        -1 indicating misclassified cells
    By default, the score is scaled between 0 and 1. This is controlled `scale=True`

    :param group_key: key in adata.obs of cell labels
    :param embed: embedding key in adata.obsm, default: 'X_pca'
    :param scale: default True, scale between 0 (worst) and 1 (best)
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    asw = silhouette_score(
        X=adata.obsm[embed],
        labels=adata.obs[group_key],
        metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw


def silhouette_batch(
        adata,
        batch_key,
        group_key,
        embed,
        metric='euclidean',
        return_all=False,
        scale=True,
        verbose=True
):
    """
    Absolute silhouette score of batch labels subsetted for each group.

    :param batch_key: batches to be compared against
    :param group_key: group labels to be subsetted by e.g. cell type
    :param embed: name of column in adata.obsm
    :param metric: see sklearn silhouette score
    :param scale: if True, scale between 0 and 1
    :param return_all: if True, return all silhouette scores and label means
        default False: return average width silhouette (ASW)
    :param verbose:
    :return:
        average width silhouette ASW
        mean silhouette per group in pd.DataFrame
        Absolute silhouette scores per group label
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')

    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])

    for group in adata.obs[group_key].unique():
        adata_group = adata[adata.obs[group_key] == group]
        n_batches = adata_group.obs[batch_key].nunique()

        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue

        sil_per_group = silhouette_samples(
            adata_group.obsm[embed],
            adata_group.obs[batch_key],
            metric=metric
        )

        # take only absolute value
        sil_per_group = [abs(i) for i in sil_per_group]

        if scale:
            # scale s.t. highest number is optimal
            sil_per_group = [1 - i for i in sil_per_group]

        sil_all = sil_all.append(
            pd.DataFrame({
                'group': [group] * len(sil_per_group),
                'silhouette_score': sil_per_group
            })
        )

    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby('group').mean()
    asw = sil_means['silhouette_score'].mean()

    if verbose:
        print(f'mean silhouette per cell: {sil_means}')

    if return_all:
        return asw, sil_means, sil_all

    return asw
