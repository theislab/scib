import numpy as np
import pandas as pd

from ..utils import check_adata, check_batch
from .ari import ari
from .cell_cycle import cell_cycle
from .clustering import opt_louvain
from .graph_connectivity import graph_connectivity
from .highly_variable_genes import hvg_overlap
from .isolated_labels import isolated_labels
from .kbet import kBET
from .lisi import clisi_graph, ilisi_graph
from .nmi import nmi
from .pcr import pcr_comparison
from .silhouette import silhouette, silhouette_batch
from .trajectory import trajectory_conservation


def metrics_fast(
        adata,
        adata_int,
        batch_key,
        label_key,
        **kwargs
):
    """
    Only fast metrics:

    Biological conservation
        HVG overlap
        Cell type ASW
        Isolated label ASW

    Batch conservation
        Graph connectivity
        Batch ASW
        PC regression
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
        **kwargs
    )


def metrics_slim(
        adata,
        adata_int,
        batch_key,
        label_key,
        **kwargs
):
    """
    All metrics apart from kBET and LISI scores:

    Biological conservation
        HVG overlap
        Cell type ASW
        Isolated label ASW
        Isolated label F1
        NMI cluster/label
        ARI cluster/label
        Cell cycle conservation

    Batch conservation
        Graph connectivity
        Batch ASW
        PC regression
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
        **kwargs
    )


def metrics_all(
        adata,
        adata_int,
        batch_key,
        label_key,
        **kwargs
):
    """
    All metrics

    Biological conservation
        HVG overlap
        Cell type ASW
        Isolated label ASW
        Isolated label F1
        NMI cluster/label
        ARI cluster/label
        Cell cycle conservation
        cLISI

    Batch conservation
        Graph connectivity
        Batch ASW
        PC regression
        kBET
        iLISI
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
        **kwargs
    )


def metrics(
        adata,
        adata_int,
        batch_key,
        label_key,
        hvg_score_=False,
        cluster_key='cluster',
        cluster_nmi=None,
        ari_=False,
        nmi_=False,
        nmi_method='arithmetic',
        nmi_dir=None,
        silhouette_=False,
        embed='X_pca',
        si_metric='euclidean',
        pcr_=False,
        cell_cycle_=False,
        organism='mouse',
        isolated_labels_=False,  # backwards compatibility
        isolated_labels_f1_=False,
        isolated_labels_asw_=False,
        n_isolated=None,
        graph_conn_=False,
        kBET_=False,
        subsample=0.5,
        lisi_graph_=False,
        ilisi_=False,
        clisi_=False,
        trajectory_=False,
        type_=None,
        verbose=False,
):
    """
    Master metrics function: Wrapper for all metrics used in the study
    Compute of all metrics given unintegrate and integrated anndata object
    """

    check_adata(adata)
    check_batch(batch_key, adata.obs)
    check_batch(label_key, adata.obs)

    check_adata(adata_int)
    check_batch(batch_key, adata_int.obs)
    check_batch(label_key, adata_int.obs)

    # clustering
    if nmi_ or ari_:
        res_max, nmi_max, nmi_all = opt_louvain(
            adata_int,
            label_key=label_key,
            cluster_key=cluster_key,
            use_rep=embed,
            function=nmi,
            plot=False,
            verbose=verbose,
            inplace=True,
            force=True
        )
        if cluster_nmi is not None:
            nmi_all.to_csv(cluster_nmi, header=False)
            print(f'saved clustering NMI values to {cluster_nmi}')

    results = {}

    if nmi_:
        print('NMI...')
        nmi_score = nmi(
            adata_int,
            group1=cluster_key,
            group2=label_key,
            method=nmi_method,
            nmi_dir=nmi_dir
        )
    else:
        nmi_score = np.nan

    if ari_:
        print('ARI...')
        ari_score = ari(
            adata_int,
            group1=cluster_key,
            group2=label_key
        )
    else:
        ari_score = np.nan

    if silhouette_:
        print('Silhouette score...')
        # global silhouette coefficient
        asw_label = silhouette(
            adata_int,
            group_key=label_key,
            embed=embed,
            metric=si_metric
        )
        # silhouette coefficient per batch
        asw_batch = silhouette_batch(
            adata_int,
            batch_key=batch_key,
            group_key=label_key,
            embed=embed,
            metric=si_metric,
            return_all=False,
            verbose=False
        )
    else:
        asw_label = np.nan
        asw_batch = np.nan

    if pcr_:
        print('PC regression...')
        pcr_score = pcr_comparison(
            adata,
            adata_int,
            embed=embed,
            covariate=batch_key,
            verbose=verbose
        )
    else:
        pcr_score = np.nan

    if cell_cycle_:
        print('cell cycle effect...')
        cc_score = cell_cycle(
            adata,
            adata_int,
            batch_key=batch_key,
            embed=embed,
            agg_func=np.mean,
            organism=organism
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
            verbose=False
        )
    else:
        il_score_f1 = np.nan

    if isolated_labels_asw_ or isolated_labels_:
        print("Isolated labels ASW...")
        il_score_asw = isolated_labels(
            adata_int,
            label_key=label_key,
            batch_key=batch_key,
            embed=embed,
            cluster=False,
            iso_threshold=n_isolated,
            verbose=False
        ) if silhouette_ else np.nan
    else:
        il_score_asw = np.nan

    if graph_conn_:
        print('Graph connectivity...')
        graph_conn_score = graph_connectivity(
            adata_int,
            label_key=label_key
        )
    else:
        graph_conn_score = np.nan

    if kBET_:
        print('kBET...')
        kbet_score = kBET(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            type_=type_,
            embed=embed,
            scaled=True,
            verbose=verbose
        )
    else:
        kbet_score = np.nan

    if lisi_graph_:
        clisi_ = True
        ilisi_ = True

    if clisi_:
        print('cLISI score...')
        clisi = clisi_graph(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            type_=type_,
            subsample=subsample * 100,
            scale=True,
            multiprocessing=True,
            verbose=verbose
        )
    else:
        clisi = np.nan

    if ilisi_:
        print('iLISI score...')
        ilisi = ilisi_graph(
            adata_int,
            batch_key=batch_key,
            type_=type_,
            subsample=subsample * 100,
            scale=True,
            multiprocessing=True,
            verbose=verbose
        )
    else:
        ilisi = np.nan

    if hvg_score_:
        hvg_score = hvg_overlap(adata, adata_int, batch_key)
    else:
        hvg_score = np.nan

    if trajectory_:
        print('Trajectory conservation score...')
        trajectory_score = trajectory_conservation(
            adata,
            adata_int,
            label_key=label_key,
            # batch_key=batch_key
        )
    else:
        trajectory_score = np.nan

    results = {
        'NMI_cluster/label': nmi_score,
        'ARI_cluster/label': ari_score,
        'ASW_label': asw_label,
        'ASW_label/batch': asw_batch,
        'PCR_batch': pcr_score,
        'cell_cycle_conservation': cc_score,
        'isolated_label_F1': il_score_f1,
        'isolated_label_silhouette': il_score_asw,
        'graph_conn': graph_conn_score,
        'kBET': kbet_score,
        'iLISI': ilisi,
        'cLISI': clisi,
        'hvg_overlap': hvg_score,
        'trajectory': trajectory_score
    }

    return pd.DataFrame.from_dict(results, orient='index')


# Deprecated

def measureTM(*args, **kwargs):
    """
    Deprecated
    params:
        *args: function to be tested for time and memory
        **kwargs: list of function paramters
    returns:
        tuple : (memory (MB), time (s), list of *args function outputs)
    """
    import cProfile
    from pstats import Stats

    import memory_profiler

    prof = cProfile.Profile()
    out = memory_profiler.memory_usage((prof.runcall, args, kwargs), retval=True)
    mem = np.max(out[0]) - out[0][0]
    print(f'memory usage:{round(mem, 0)} MB')
    print(f'runtime: {round(Stats(prof).total_tt, 0)} s')
    return mem, Stats(prof).total_tt, out[1:]
