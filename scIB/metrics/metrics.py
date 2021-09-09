import pandas as pd
from scIB.utils import *
from scIB.clustering import opt_louvain
import cProfile
from pstats import Stats
import memory_profiler

from .ari import ari
from .cell_cycle import cell_cycle
from .graph_connectivity import graph_connectivity
from .highly_variable_genes import hvg_overlap
from .isolated_labels import isolated_labels
from .kbet import kBET
from .lisi import lisi_graph, scale_lisi
from .nmi import nmi
from .pcr import pcr_comparison
from .silhouette import silhouette, silhouette_batch
from .trajectory import trajectory_conservation
from .utils import NeighborsError, RootCellError


def measureTM(*args, **kwargs):
    """
    params:
        *args: function to be tested for time and memory
        **kwargs: list of function paramters
    returns:
        tuple : (memory (MB), time (s), list of *args function outputs)
    """
    prof = cProfile.Profile()
    out = memory_profiler.memory_usage((prof.runcall, args, kwargs), retval=True)
    mem = np.max(out[0]) - out[0][0]
    print(f'memory usage:{round(mem, 0)} MB')
    print(f'runtime: {round(Stats(prof).total_tt, 0)} s')
    return mem, Stats(prof).total_tt, out[1:]


def metrics(
        adata,
        adata_int,
        batch_key,
        label_key,
        hvg_score_=False,
        cluster_key='cluster',
        cluster_opt_out=None,
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
        isolated_labels_f1_=False,
        isolated_labels_asw_=False,
        n_isolated=None,
        graph_conn_=False,
        kBET_=False,
        kBET_sub=0.5,
        lisi_graph_=False,
        lisi_raw=False,
        trajectory_=False,
        type_=None,
        verbose=False,
):
    """
    Master metrics function: Wrapper for all metrics used in the study
    Compute of all metrics given unintegrate and integrated anndata object
    """

    checkAdata(adata)
    checkBatch(batch_key, adata.obs)
    checkBatch(label_key, adata.obs)

    checkAdata(adata_int)
    checkBatch(batch_key, adata_int.obs)
    checkBatch(label_key, adata_int.obs)

    # clustering
    if nmi_ or ari_:
        print('Clustering...')
        res_max, nmi_max, nmi_all = opt_louvain(
            adata_int,
            label_key=label_key,
            cluster_key=cluster_key,
            function=nmi,
            plot=False,
            verbose=verbose,
            inplace=True,
            force=True
        )
        if cluster_opt_out is not None:
            nmi_all.to_csv(cluster_opt_out, header=False)
            print(f'saved clustering NMI values to {cluster_opt_out}')

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
        sil_global = silhouette(
            adata_int,
            group_key=label_key,
            embed=embed,
            metric=si_metric
        )
        # silhouette coefficient per batch
        _, sil_clus = silhouette_batch(
            adata_int,
            batch_key=batch_key,
            group_key=label_key,
            embed=embed,
            metric=si_metric,
            verbose=False
        )
        sil_clus = sil_clus['silhouette_score'].mean()
    else:
        sil_global = np.nan
        sil_clus = np.nan

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

    if isolated_labels_f1_:
        print("Isolated labels F1...")
        il_score_f1 = isolated_labels(
            adata_int,
            label_key=label_key,
            batch_key=batch_key,
            embed=embed,
            cluster=True,
            n=n_isolated,
            verbose=False
        )
    else:
        il_score_f1 = np.nan

    if isolated_labels_asw_:
        print("Isolated labels ASW...")
        il_score_asw = isolated_labels(
            adata_int,
            label_key=label_key,
            batch_key=batch_key,
            embed=embed,
            cluster=False,
            n=n_isolated,
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
        try:
            kbet_results = kBET(
                adata_int,
                batch_key=batch_key,
                label_key=label_key,
                type_=type_,
                embed=embed,
                subsample=kBET_sub,
                heuristic=True,
                verbose=verbose
            )
            kbet_score = 1 - np.nanmean(kbet_results['kBET'])
        except NeighborsError:
            print('Not enough neighbours')
            kbet_score = 0
    else:
        kbet_score = np.nan

    # if lisi_:
    #    print('LISI score...')
    #    ilisi_score, clisi_score = lisi(adata_int, batch_key=batch_key, label_key=label_key,
    #                                    type_ = type_, verbose=verbose)
    # else:
    #    ilisi_score = np.nan
    #    clisi_score = np.nan
    # results['iLISI'] = ilisi_score
    # results['cLISI'] = clisi_score

    if lisi_graph_:
        print('LISI graph score...')
        ilisi_raw, clisi_raw = lisi_graph(
            adata_int,
            batch_key=batch_key,
            label_key=label_key,
            type_=type_,
            subsample=kBET_sub * 100,
            scale=False,
            multiprocessing=True,
            verbose=verbose
        )
        nbatches = len(np.unique(adata_int.obs[batch_key]))
        ilisi_scaled, clisi_scaled = scale_lisi(
            ilisi_raw,
            clisi_raw,
            nbatches
        )
        if lisi_raw:
            results['iLISI_raw'] = ilisi_raw
            results['cLISI_raw'] = clisi_raw
    else:
        ilisi_scaled = np.nan
        clisi_scaled = np.nan

    if hvg_score_:
        hvg_score = hvg_overlap(adata, adata_int, batch_key)
    else:
        hvg_score = np.nan

    if trajectory_:
        print('Trajectory conservation score...')
        try:
            trajectory_score = trajectory_conservation(adata, adata_int, label_key=label_key)
        except RootCellError:
            print('No cell of root cluster in largest connected component')
            trajectory_score = 0
    else:
        trajectory_score = np.nan

    results = {
        'NMI_cluster/label': nmi_score,
        'ARI_cluster/label': ari_score,
        'ASW_label': sil_global,
        'ASW_label/batch': sil_clus,
        'PCR_batch': pcr_score,
        'cell_cycle_conservation': cc_score,
        'isolated_label_F1': il_score_f1,
        'isolated_label_silhouette': il_score_asw,
        'graph_conn': graph_conn_score,
        'kBET': kbet_score,
        'iLISI': ilisi_scaled,
        'cLISI': clisi_scaled,
        'hvg_overlap': hvg_score,
        'trajectory': trajectory_score
    }

    return pd.DataFrame.from_dict(results, orient='index')
