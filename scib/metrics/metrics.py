import numpy as np
import pandas as pd

from ..utils import check_adata, check_batch
from .ari import ari
from .cell_cycle import cell_cycle
from .clustering import cluster_optimal_resolution
from .graph_connectivity import graph_connectivity
from .highly_variable_genes import hvg_overlap
from .isolated_labels import isolated_labels
from .kbet import kBET
from .lisi import clisi_graph, ilisi_graph
from .nmi import nmi
from .pcr import pcr_comparison
from .silhouette import silhouette, silhouette_batch
from .trajectory import trajectory_conservation


def metrics_fast(adata, adata_int, batch_key, label_key, **kwargs):
    """Only metrics with minimal preprocessing and runtime


    :Biological conservation:
        + HVG overlap :func:`~scib.metrics.hvg_overlap`
        + Cell type ASW :func:`~scib.metrics.silhouette`
        + Isolated label ASW :func:`~scib.metrics.isolated_labels`

    :Batch correction:
        + Graph connectivity :func:`~scib.metrics.graph_connectivity`
        + Batch ASW :func:`~scib.metrics.silhouette_batch`
        + Principal component regression :func:`~scib.metrics.pcr_comparison`

    :param adata: unintegrated, preprocessed anndata object
    :param adata_int: integrated anndata object
    :param batch_key: name of batch column in adata.obs and adata_int.obs
    :param label_key: name of biological label (cell type) column in adata.obs and adata_int.obs
    :param kwargs:
        Parameters to pass on to :func:`~scib.metrics.metrics` function:

            + ``embed``
            + ``si_metric``
            + ``n_isolated``
    """
    return metrics(
        adata,
        adata_int,
        batch_key,
        label_key,
        isolated_labels_asw_=True,
        silhouette_=True,
        hvg_score_=True,
        graph_conn_=True,
        pcr_=True,
        **kwargs,
    )


def metrics_slim(adata, adata_int, batch_key, label_key, **kwargs):
    """All metrics apart from kBET and LISI scores

    :Biological conservation:
        + HVG overlap :func:`~scib.metrics.hvg_overlap`
        + Cell type ASW :func:`~scib.metrics.silhouette`
        + Isolated label ASW :func:`~scib.metrics.isolated_labels`
        + Isolated label F1 :func:`~scib.metrics.isolated_labels`
        + NMI cluster/label :func:`~scib.metrics.nmi`
        + ARI cluster/label :func:`~scib.metrics.ari`
        + Cell cycle conservation :func:`~scib.metrics.cell_cycle`
        + Trajectory conservation :func:`~scib.metrics.trajectory_conservation`

    :Batch correction:
        + Graph connectivity :func:`~scib.metrics.graph_connectivity`
        + Batch ASW :func:`~scib.metrics.silhouette_batch`
        + Principal component regression :func:`~scib.metrics.pcr_comparison`

    :param adata: unintegrated, preprocessed anndata object
    :param adata_int: integrated anndata object
    :param batch_key: name of batch column in adata.obs and adata_int.obs
    :param label_key: name of biological label (cell type) column in adata.obs and adata_int.obs
    :param kwargs:
        Parameters to pass on to :func:`~scib.metrics.metrics` function:

            + ``embed``
            + ``cluster_key``
            + ``cluster_nmi``
            + ``nmi_method``
            + ``nmi_dir``
            + ``si_metric``
            + ``organism``
            + ``n_isolated``
    """
    return metrics(
        adata,
        adata_int,
        batch_key,
        label_key,
        isolated_labels_asw_=True,
        silhouette_=True,
        hvg_score_=True,
        graph_conn_=True,
        pcr_=True,
        isolated_labels_f1_=True,
        trajectory_=True,
        nmi_=True,
        ari_=True,
        cell_cycle_=True,
        **kwargs,
    )


def metrics_all(adata, adata_int, batch_key, label_key, **kwargs):
    """All metrics

    :Biological conservation:
        + HVG overlap :func:`~scib.metrics.hvg_overlap`
        + Cell type ASW :func:`~scib.metrics.silhouette`
        + Isolated label ASW :func:`~scib.metrics.isolated_labels`
        + Isolated label F1 :func:`~scib.metrics.isolated_labels`
        + NMI cluster/label :func:`~scib.metrics.nmi`
        + ARI cluster/label :func:`~scib.metrics.ari`
        + Cell cycle conservation :func:`~scib.metrics.cell_cycle`
        + cLISI (cell type Local Inverse Simpson's Index) :func:`~scib.metrics.clisi_graph`
        + Trajectory conservation :func:`~scib.metrics.trajectory_conservation`

    :Batch correction:
        + Graph connectivity :func:`~scib.metrics.graph_connectivity`
        + Batch ASW :func:`~scib.metrics.silhouette_batch`
        + Principal component regression :func:`~scib.metrics.pcr_comparison`
        + kBET (k-nearest neighbour batch effect test) :func:`~scib.metrics.kBET`
        + iLISI (integration Local Inverse Simpson's Index) :func:`~scib.metrics.ilisi_graph`

    :param adata: unintegrated, preprocessed anndata object
    :param adata_int: integrated anndata object
    :param batch_key: name of batch column in adata.obs and adata_int.obs
    :param label_key: name of biological label (cell type) column in adata.obs and adata_int.obs
    :param kwargs:
        Parameters to pass on to :func:`~scib.metrics.metrics` function:

            + ``embed``
            + ``cluster_key``
            + ``cluster_nmi``
            + ``nmi_method``
            + ``nmi_dir``
            + ``si_metric``
            + ``organism``
            + ``n_isolated``
            + ``subsample``
            + ``type_``
    """
    return metrics(
        adata,
        adata_int,
        batch_key,
        label_key,
        isolated_labels_asw_=True,
        silhouette_=True,
        hvg_score_=True,
        graph_conn_=True,
        pcr_=True,
        isolated_labels_f1_=True,
        trajectory_=True,
        nmi_=True,
        ari_=True,
        cell_cycle_=True,
        kBET_=True,
        ilisi_=True,
        clisi_=True,
        **kwargs,
    )


def metrics(
    adata,
    adata_int,
    batch_key,
    label_key,
    embed="X_pca",
    cluster_key="cluster",
    cluster_nmi=None,
    ari_=False,
    nmi_=False,
    nmi_method="arithmetic",
    nmi_dir=None,
    silhouette_=False,
    si_metric="euclidean",
    pcr_=False,
    cell_cycle_=False,
    organism="mouse",
    hvg_score_=False,
    isolated_labels_=False,  # backwards compatibility
    isolated_labels_f1_=False,
    isolated_labels_asw_=False,
    n_isolated=None,
    graph_conn_=False,
    trajectory_=False,
    kBET_=False,
    lisi_graph_=False,
    ilisi_=False,
    clisi_=False,
    subsample=0.5,
    n_cores=1,
    type_=None,
    verbose=False,
):
    """Master metrics function

    Wrapper for all metrics used in the study.
    Compute of all metrics given unintegrated and integrated anndata object

    :param adata:
        unintegrated, preprocessed anndata object
    :param adata_int:
        integrated anndata object
    :param batch_key:
        name of batch column in adata.obs and adata_int.obs
    :param label_key:
        name of biological label (cell type) column in adata.obs and adata_int.obs
    :param embed:
        embedding representation of adata_int

        Used for:

            + silhouette scores (label ASW, batch ASW),
            + PC regression,
            + cell cycle conservation,
            + isolated label scores, and
            + kBET
    :param cluster_key:
        name of column to store cluster assignments. Will be overwritten if it exists
    :param cluster_nmi:
        Where to save cluster resolutions and NMI for optimal clustering
        If None, these results will not be saved
    :param `ari_`:
        whether to compute ARI using :func:`~scib.metrics.ari`
    :param `nmi_`:
        whether to compute NMI using :func:`~scib.metrics.nmi`
    :param nmi_method:
        which implementation of NMI to use
    :param nmi_dir:
        directory of NMI code for some implementations of NMI
    :param `silhouette_`:
        whether to compute the average silhouette width scores for labels and batch
        using :func:`~scib.metrics.silhouette` and :func:`~scib.metrics.silhouette_batch`
    :param si_metric:
        which distance metric to use for silhouette scores
    :param `pcr_`:
        whether to compute principal component regression using :func:`~scib.metrics.pc_comparison`
    :param `cell_cycle_`:
        whether to compute cell cycle score conservation using :func:`~scib.metrics.cell_cycle`
    :param organism:
        organism of the datasets, used for computing cell cycle scores on gene names
    :param `hvg_score_`:
        whether to compute highly variable gene conservation using :func:`~scib.metrics.hvg_overlap`
    :param `isolated_labels_`:
        whether to compute both isolated label scores using :func:`~scib.metrics.isolated_labels`
    :param `isolated_labels_f1_`:
        whether to compute isolated label score based on F1 score of clusters vs labels using
        :func:`~scib.metrics.isolated_labels`
    :param `isolated_labels_asw_`:
        whether to compute isolated label score based on ASW (average silhouette width) using
        :func:`~scib.metrics.isolated_labels`
    :param `n_isolated`:
        maximum number of batches per label for label to be considered as isolated
    :param `graph_conn_`:
        whether to compute graph connectivity score using :func:`~scib.metrics.graph_connectivity`
    :param `trajectory_`:
        whether to compute trajectory score using :func:`~scib.metrics.trajectory_conservation`
    :param `kBET_`:
        whether to compute kBET score using :func:`~scib.metrics.kBET`
    :param `lisi_graph_`:
        whether to compute both cLISI and iLISI using :func:`~scib.metrics.lisi_graph`
    :param `clisi_`:
        whether to compute cLISI using :func:`~scib.metrics.clisi_graph`
    :param `ilisi_`:
        whether to compute iLISI using :func:`~scib.metrics.ilisi_graph`
    :param subsample:
        subsample fraction for LISI scores
    :param n_cores: number of cores to be used for LISI functions
    :param `type_`:
        one of 'full', 'embed' or 'knn' (used for kBET and LISI scores)
    """

    check_adata(adata)
    check_batch(batch_key, adata.obs)
    check_batch(label_key, adata.obs)

    check_adata(adata_int)
    check_batch(batch_key, adata_int.obs)
    check_batch(label_key, adata_int.obs)

    # clustering
    if nmi_ or ari_:
        res_max, nmi_max, nmi_all = cluster_optimal_resolution(
            adata_int,
            label_key=label_key,
            cluster_key=cluster_key,
            metric=nmi,
            use_rep=embed,
            force=True,
            verbose=verbose,
            return_all=True,
        )
        if cluster_nmi is not None:
            nmi_all.to_csv(cluster_nmi, header=False)
            print(f"saved clustering NMI values to {cluster_nmi}")

    if nmi_:
        print("NMI...")
        nmi_score = nmi(
            adata_int,
            cluster_key=cluster_key,
            label_key=label_key,
            implementation=nmi_method,
            nmi_dir=nmi_dir,
        )
    else:
        nmi_score = np.nan

    if ari_:
        print("ARI...")
        ari_score = ari(adata_int, cluster_key=cluster_key, label_key=label_key)
    else:
        ari_score = np.nan

    if silhouette_:
        print("Silhouette score...")
        # global silhouette coefficient
        asw_label = silhouette(
            adata_int, label_key=label_key, embed=embed, metric=si_metric
        )
        # silhouette coefficient per batch
        asw_batch = silhouette_batch(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            embed=embed,
            metric=si_metric,
            return_all=False,
            verbose=False,
        )
    else:
        asw_label = np.nan
        asw_batch = np.nan

    if pcr_:
        print("PC regression...")
        pcr_score = pcr_comparison(
            adata, adata_int, embed=embed, covariate=batch_key, verbose=verbose
        )
    else:
        pcr_score = np.nan

    if cell_cycle_:
        print("cell cycle effect...")
        cc_score = cell_cycle(
            adata,
            adata_int,
            batch_key=batch_key,
            embed=embed,
            agg_func=np.mean,
            organism=organism,
        )
    else:
        cc_score = np.nan

    if isolated_labels_f1_ or isolated_labels_:
        print("Isolated labels F1...")
        il_score_f1 = isolated_labels(
            adata_int,
            label_key=label_key,
            batch_key=batch_key,
            embed=embed,
            cluster=True,
            iso_threshold=n_isolated,
            verbose=False,
        )
    else:
        il_score_f1 = np.nan

    if isolated_labels_asw_ or isolated_labels_:
        print("Isolated labels ASW...")
        il_score_asw = (
            isolated_labels(
                adata_int,
                label_key=label_key,
                batch_key=batch_key,
                embed=embed,
                cluster=False,
                iso_threshold=n_isolated,
                verbose=False,
            )
            if silhouette_
            else np.nan
        )
    else:
        il_score_asw = np.nan

    if graph_conn_:
        print("Graph connectivity...")
        graph_conn_score = graph_connectivity(adata_int, label_key=label_key)
    else:
        graph_conn_score = np.nan

    if kBET_:
        print("kBET...")
        kbet_score = kBET(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            type_=type_,
            embed=embed,
            scaled=True,
            verbose=verbose,
        )
    else:
        kbet_score = np.nan

    if lisi_graph_:
        clisi_ = True
        ilisi_ = True

    if clisi_:
        print("cLISI score...")
        clisi = clisi_graph(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            type_=type_,
            subsample=subsample * 100,
            scale=True,
            n_cores=n_cores,
            verbose=verbose,
        )
    else:
        clisi = np.nan

    if ilisi_:
        print("iLISI score...")
        ilisi = ilisi_graph(
            adata_int,
            batch_key=batch_key,
            type_=type_,
            subsample=subsample * 100,
            scale=True,
            n_cores=n_cores,
            verbose=verbose,
        )
    else:
        ilisi = np.nan

    if hvg_score_:
        hvg_score = hvg_overlap(adata, adata_int, batch_key)
    else:
        hvg_score = np.nan

    if trajectory_:
        print("Trajectory conservation score...")
        trajectory_score = trajectory_conservation(
            adata,
            adata_int,
            label_key=label_key,
            # batch_key=batch_key
        )
    else:
        trajectory_score = np.nan

    results = {
        "NMI_cluster/label": nmi_score,
        "ARI_cluster/label": ari_score,
        "ASW_label": asw_label,
        "ASW_label/batch": asw_batch,
        "PCR_batch": pcr_score,
        "cell_cycle_conservation": cc_score,
        "isolated_label_F1": il_score_f1,
        "isolated_label_silhouette": il_score_asw,
        "graph_conn": graph_conn_score,
        "kBET": kbet_score,
        "iLISI": ilisi,
        "cLISI": clisi,
        "hvg_overlap": hvg_score,
        "trajectory": trajectory_score,
    }

    return pd.DataFrame.from_dict(results, orient="index")
