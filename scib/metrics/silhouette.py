import numpy as np
import pandas as pd
from sklearn.metrics.cluster import silhouette_samples, silhouette_score

try:
    from scanpy._utils import renamed_arg
except ImportError:
    from .._package_tools import renamed_arg


@renamed_arg("group_key", "label_key")
def silhouette(adata, label_key, embed, metric="euclidean", scale=True):
    """Average silhouette width (ASW)

    Wrapper for sklearn silhouette function values range from [-1, 1] with

        * 1 indicates distinct, compact clusters
        * 0 indicates overlapping clusters
        * -1 indicates core-periphery (non-cluster) structure

    By default, the score is scaled between 0 and 1 (``scale=True``).

    :param label_key: key in adata.obs of cell labels
    :param embed: embedding key in adata.obsm, default: 'X_pca'
    :param metric: type of distance metric to use for the silhouette scores
    :param scale: default True, scale between 0 (worst) and 1 (best)

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
        scib.me.silhouette(adata, label_key="celltype", embed="X_pca")

        # embedding output
        scib.me.silhouette(adata, label_key="celltype", embed="X_emb")
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f"{embed} not in obsm")
    asw = silhouette_score(
        X=adata.obsm[embed], labels=adata.obs[label_key], metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw


@renamed_arg("group_key", "label_key")
def silhouette_batch(
    adata,
    batch_key,
    label_key,
    embed,
    metric="euclidean",
    return_all=False,
    scale=True,
    verbose=True,
):
    """Batch ASW

    Modified average silhouette width (ASW) of batch

    This metric measures the silhouette of a given batch.
    It assumes that a silhouette width close to 0 represents perfect overlap of the batches, thus the absolute value of
    the silhouette width is used to measure how well batches are mixed.
    For all cells :math:`i` of a cell type :math:`C_j`, the batch ASW of that cell type is:

    .. math::

        batch \\, ASW_j = \\frac{1}{|C_j|} \\sum_{i \\in C_j} |silhouette(i)|

    The final score is the average of the absolute silhouette widths computed per cell type :math:`M`.

    .. math::

        batch \\, ASW = \\frac{1}{|M|} \\sum_{i \\in M} batch \\, ASW_j

    For a scaled metric (which is the default), the absolute ASW per group is subtracted from 1 before averaging, so that
    0 indicates suboptimal label representation and 1 indicates optimal label representation.

    .. math::

        batch \\, ASW_j = \\frac{1}{|C_j|} \\sum_{i \\in C_j} 1 - |silhouette(i)|

    :param batch_key: batch labels to be compared against
    :param label_key: group labels to be subset by e.g. cell type
    :param embed: name of column in adata.obsm
    :param metric: see sklearn silhouette score
    :param scale: if True, scale between 0 and 1
    :param return_all: if True, return all silhouette scores and label means
        default False: return average width silhouette (ASW)
    :param verbose: print silhouette score per group
    :return:
        Batch ASW  (always)
        Mean silhouette per group in pd.DataFrame (additionally, if return_all=True)
        Absolute silhouette scores per group label (additionally, if return_all=True)

    The function requires an embedding to be stored in ``adata.obsm`` and can only be applied to feature and embedding
    integration outputs.
    Please note, that the metric cannot be used to evaluate kNN graph outputs.
    See :ref:`preprocessing` for more information on preprocessing.

    **Examples**

    .. code-block:: python

        # feature output
        scib.pp.reduce_data(
            adata, n_top_genes=2000, batch_key="batch", pca=True, neighbors=False
        )
        scib.me.silhouette_batch(adata, batch_key="batch", label_key="celltype", embed="X_pca")

        # embedding output
        scib.me.silhouette_batch(adata, batch_key="batch", label_key="celltype", embed="X_emb")

    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f"{embed} not in obsm")

    sil_per_label = []
    for group in adata.obs[label_key].unique():
        adata_group = adata[adata.obs[label_key] == group]
        n_batches = adata_group.obs[batch_key].nunique()

        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue

        sil = silhouette_samples(
            adata_group.obsm[embed], adata_group.obs[batch_key], metric=metric
        )

        # take only absolute value
        sil = [abs(i) for i in sil]

        if scale:
            # scale s.t. highest number is optimal
            sil = [1 - i for i in sil]

        sil_per_label.extend([(group, score) for score in sil])

    sil_df = pd.DataFrame.from_records(
        sil_per_label, columns=["group", "silhouette_score"]
    )

    if len(sil_per_label) == 0:
        sil_means = np.nan
        asw = np.nan
    else:
        sil_means = sil_df.groupby("group").mean()
        asw = sil_means["silhouette_score"].mean()

    if verbose:
        print(f"mean silhouette per group: {sil_means}")

    if return_all:
        return asw, sil_means, sil_df

    return asw
