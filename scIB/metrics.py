import numpy as np
from scipy import sparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata
from scIB.utils import *
from scIB.preprocessing import hvg_intersect

import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Ignore R warning messages
import rpy2.robjects as ro
import anndata2ri


### Silhouette score
def silhouette_score(adata, batch_key='study', group_key='cell_type', metric='euclidean', 
                     embed='X_pca', means=False, verbose=True):
    """
    Silhouette score of batch labels subsetted for each group.
    params:
        adata: 
        batch_key: batches to be compared against
        group_key: group labels to be subsetted by e.g. cell type
        metric: see sklearn silhouette score
        embed: name of column in adata.obsm
    returns:
        all scores: absolute silhouette scores per group label
        group means: if `mean=True`
    """
    import sklearn.metrics as scm

    
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    
    # ony look at group labels that are present for all batches
    n_batch = adata.obs[batch_key].nunique()
    groups = adata.obs.groupby(group_key)[batch_key].nunique()
    groups = groups[groups == n_batch]
    
    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])
    for group in groups.index:
        adata_group = adata[adata.obs[group_key] == group]
        sil_per_group = scm.silhouette_samples(adata_group.obsm[embed],
                                               adata_group.obs[batch_key],
                                               metric=metric)
        sil_per_group = [abs(i) for i in sil_per_group] # take only absolute value
        d = pd.DataFrame({'group' : [group]*len(sil_per_group), 'silhouette_score' : sil_per_group})
        sil_all = sil_all.append(d)
    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby('group').mean()
    
    if verbose:
        print(f'mean silhouette per cell: {sil_means}')
    
    if means:
        return sil_all, sil_means
    
    return sil_all

def plot_silhouette_score(adata_dict, batch_key='study', group_key='cell_type', metric='euclidean', 
                     embed='X_pca', palette='Dark2', per_group=False, verbose=True):
    """
    params:
        adata_dict: dictionary of adata objects, each labeled by e.g. integration method name
    """
    
    with sns.color_palette(palette):
        for label, adata in adata_dict.items():
            checkAdata(adata)
            sil_scores = silhouette_score(adata, 
                                          batch_key=batch_key,
                                          group_key=group_key,
                                          metric=metric,
                                          embed=embed,
                                          means=False,
                                          verbose=verbose)
            sns.distplot(sil_scores['silhouette_score'], label=label, hist=False)
        plt.title('Silhouette scores per cell for all groups')
        plt.show()
        
        if per_group:
            print()
            for data_set, adata in adata_dict.items():
                sil_scores = silhouette_score(adata,
                                              batch_key=batch_key,
                                              group_key=group_key,
                                              metric=metric,
                                              embed=embed,
                                              means=False,
                                              verbose=verbose)
                # plot for all groups
                for group in sil_scores['group'].unique():
                    group_scores = sil_scores[sil_scores['group'] == group]
                    sns.distplot(group_scores['silhouette_score'], label=group, hist=False)
                plt.title(f'Silhouette scores per cell for {data_set}')
                plt.show()

### Naive cluster overlap
def cluster_overlap(adata, group1='louvain', group2='louvain_post'):
    
    checkAdata(adata)
    checkBatch(group1, adata.obs)
    checkBatch(group2, adata.obs)
    
    cluster_ov = {}
    louv_post_sizes = adata.obs.groupby(group1).size()
    for i in adata.obs[group2].cat.categories:
        a = adata.obs[adata.obs[group2] == i]
        overlap = a.groupby(group1).size()
        cluster_ov[i] = (overlap / louv_post_sizes).sum() / len(overlap[overlap > 0])
    return cluster_ov

def plot_cluster_overlap(adata_dict, group1, group2, df=False):
    """
    params:
        adata_dict: dictionary of adata objects, each labeled by e.g. integration method name
        group1: column containing cluster assignments
        group2: column containing cluster assignments
    return:
        clust_df: dataframe with plotted data points
    """
    series = []
    dict_keys = list(adata_dict.keys())
    for i in dict_keys:
        c_ov = cluster_overlap(adata_dict[i], group1=group1, group2=group2)
        series.append(pd.Series(c_ov))
    clust_df = pd.DataFrame(series).transpose()
    clust_df.columns = dict_keys
    
    with sns.set_palette('Dark2'):
        sns.boxplot(data=clust_df)
        sns.swarmplot(data=clust_df, color=".25")
    
    if df:
        return clust_df
    return None


### NMI normalised mutual information
def nmi(adata, group1, group2, method="max", nmi_dir=None):
    """
    Normalized mutual information NMI based on 2 different cluster assignments `group1` and `group2`
    params:
        adata: Anndata object
        group1: column name of `adata.obs`
        group2: column name of `adata.obs`
        method: NMI implementation
            'max': scikit method with `average_method='max'`
            'min': scikit method with `average_method='min'`
            'geometric': scikit method with `average_method='geometric'`
            'arithmetic': scikit method with `average_method='arithmetic'`
            'Lancichinetti': implementation by A. Lancichinetti 2009 et al.
            'ONMI': implementation by Aaron F. McDaid et al. (https://github.com/aaronmcdaid/Overlapping-NMI) Hurley 2011
        nmi_dir: directory of compiled C code if 'Lancichinetti' or 'ONMI' are specified as `method`. Compilation should be done as specified in the corresponding README.
    return:
        normalized mutual information (NMI)
    """
    
    if method in ['max', 'min', 'geometric', 'arithmetic']:
        nmi_value = nmi_scikit(adata, group1, group2, average_method=method)
    elif method == "Lancichinetti":
        nmi_value = nmi_Lanc(adata, group1, group2, nmi_dir=nmi_dir)
    elif method == "ONMI":
        nmi_value = onmi(adata, group1, group2, nmi_dir=nmi_dir)
    else:
        raise ValueError(f"Method {method} not valid")
    
    return nmi_value


def nmi_scikit(adata, group1, group2, average_method='max'):
    """
    implementation from scikit-learn
    params:
        average_method: different ways of averaging for normalization
            'max': scikit method with `average_method='max'`
            'min': scikit method with `average_method='min'`
            'geometric': scikit method with `average_method='geometric'`
            'arithmetic': scikit method with `average_method='arithmetic'`
    """
    checkAdata(adata)
    checkBatch(group1, adata.obs)
    checkBatch(group2, adata.obs)
    
    from sklearn.metrics import normalized_mutual_info_score
    
    group1_list = adata.obs[group1].tolist()
    group2_list = adata.obs[group2].tolist()
    
    return normalized_mutual_info_score(group1_list, group2_list, average_method=average_method)
    

def onmi(adata, group1, group2, nmi_dir=None, verbose=True):
    """
    Based on implementation https://github.com/aaronmcdaid/Overlapping-NMI
    publication: Aaron F. McDaid, Derek Greene, Neil Hurley 2011
    params:
        nmi_dir: directory of compiled C code
    """
    
    checkAdata(adata)
    checkBatch(group1, adata.obs)
    checkBatch(group2, adata.obs)
    
    if nmi_dir is None:
        raise FileNotFoundError("Please provide the directory of the compiled C code from https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz")
    
    import subprocess
    import os
    
    group1_file = write_tmp_labels(adata, group1, to_int=False)
    group2_file = write_tmp_labels(adata, group2, to_int=False)
    
    nmi_call = subprocess.Popen(
        [nmi_dir+"onmi", group1_file, group2_file], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT)
    
    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)
    
    nmi_out = stdout.decode()
    if verbose:
        print(nmi_out)
    
    nmi_split = [x.strip().split('\t') for x in nmi_out.split('\n')]
    nmi_max = float(nmi_split[0][1])
    
    # remove temporary files
    os.remove(group1_file)
    os.remove(group2_file)
    
    return nmi_max


def nmi_Lanc(adata, group1, group2, nmi_dir="external/mutual3/", verbose=True):
    """
    paper by A. Lancichinetti 2009
    https://sites.google.com/site/andrealancichinetti/mutual
    recommended by Malte
    """
    
    checkAdata(adata)
    
    if nmi_dir is None:
        raise FileNotFoundError("Please provide the directory of the compiled C code from https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz")
    
    import subprocess
    import os
    
    group1_file = write_tmp_labels(adata, group1, to_int=False)
    group2_file = write_tmp_labels(adata, group2, to_int=False)
    
    nmi_call = subprocess.Popen(
        [nmi_dir+"mutual", group1_file, group2_file], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT)
    
    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)
    nmi_out = stdout.decode().strip()
    
    return float(nmi_out.split('\t')[1])

def write_tmp_labels(adata, group_name, to_int=False, delim='\n'):
    """
    write the values of a specific obs column into a temporary file in text format
    params:
        adata: anndata object
        group_name: name of column to be saved
        to_int: rename the unique column entries by integers in range(1,len(adata.obs[group_name])+1)
    """
    import tempfile
    
    checkAdata(adata)
    checkBatch(group_name, adata.obs)
    
    if to_int:
        label_map = {}
        i = 1
        for label in adata.obs[group_name].unique():
            label_map[label] = i
            i += 1
        labels = delim.join([str(label_map[name]) for name in adata.obs[group_name]])
    else:
        labels = delim.join([str(name) for name in adata.obs[group_name]])
        
    clusters = {label:[] for label in adata.obs[group_name].unique()}
    for i, label in enumerate(adata.obs[group_name]):
        clusters[label].append(str(i))
    
    output = '\n'.join([' '.join(c) for c in clusters.values()])
    output = str.encode(output)
    
    # write to file
    with tempfile.NamedTemporaryFile(delete=False) as f:
        f.write(output)
        filename = f.name
    
    return filename

### ARI adjusted rand index
def ari(adata, group1, group2):
    """
    params:
        adata: anndata object
        group1: "true" cluster assignments
        group2: "predicted" cluster assignments
    """
    from sklearn.metrics.cluster import adjusted_rand_score
    
    checkAdata(adata)
    checkBatch(group1, adata.obs)
    checkBatch(group2, adata.obs)
    
    group1_list = adata.obs[group1].tolist()
    group2_list = adata.obs[group2].tolist()
    
    return adjusted_rand_score(group1_list, group2_list)


### Cell cycle effect
def cell_cycle(adata, hvg=False, s_phase_key='S_score', g2m_phase_key='G2M_score'):
    """
    params:
        adata:
        s_phase_key: key of column containing S-phase score
        g2m_phase_key: key of column containing G2M-phase score
    
    """
    s_phase = pcr(adata, adata.X, hvg=hvg, covariate=s_phase_key)
    g2m_phase = pcr(adata, adata.X, hvg=hvg, covariate=g2m_phase_key)
    
    return s_phase, g2m_phase
    
    #returns:
    #    sum of variance difference of S-phase score and G2M-phase score
    #s_phase = pcr_comparison(adata, raw, corrected, hvg=hvg, covariate=s_phase_key)
    #g2m_phase = pcr_comparison(adata, raw, corrected, hvg=hvg, covariate=g2m_phase_key)
        
    #return s_phase + g2m_phase

### Highly Variable Genes conservation
def hvg_overlap(adata_post, adata_pre, batch, n_hvg=500):
    hvg_pre= set(hvg_intersect(adata_pre, batch=batch, max_genes=n_hvg))
    hvg_post= set(hvg_intersect(adata_post, batch=batch, max_genes=n_hvg))
    jaccard = len(hvg_pre.intersection(hvg_post))/len(hvg_pre.union(hvg_post))
    return jaccard

### PC Regression
def get_hvg_indices(adata):
    if "highly_variable" not in adata.var.columns:
        print("No highly variable genes computed, continuing with full matrix")
        return np.array(range(adata.n_vars))
    return np.where((adata.var["highly_variable"] == True))[0]
        
def pcr_comparison(adata, raw, corrected, hvg=False, covariate='sample', verbose=True):
    """
    Compare the effect before and after integration
    params:
        raw: count matrix before integration
        corrected: count matrix after correction
    return:
        difference of R2Var value of PCR
    """
    
    pcr_before = pcr(adata, matrix=raw, hvg=hvg, covariate=covariate, verbose=verbose)
    pcr_after = pcr(adata, matrix=corrected, hvg=hvg, covariate=covariate, verbose=verbose)
    
    return pcr_before - pcr_after

def pcr(adata, matrix=None, hvg=False, covariate='sample', verbose=True):
    """
    PCR for Adata object
    params:
        adata: Anndata object
        matrix: if None, take adata.X else specify count matrix
        covariate: key for adata.obs column
    return:
        R2Var of PCR
    """
    
    checkAdata(adata)
    checkBatch(covariate, adata.obs)
    
    if hvg:
        hvg_idx = get_hvg_indices(adata)
        if verbose:
            print(f"subsetting to {len(hvg_idx)} highly variable genes")
        raw = raw[:, hvg_idx]
        corrected = corrected[:, hvg_idx]
    
    if verbose:
        print(f"covariate: {covariate}")
    batch = adata.obs[covariate]
    
    return pc_regression(raw, batch, verbose)

def pc_regression(data, batch, pca_sd=None, n_comps=None, svd_solver='arpack', verbose=True):
    """
    params:
        data: Anndata or count matrix
        batch: series or list of batch assignemnts
        n_comps: number of PCA components, only when PCA is not yet computed
        pca_sd: iterable of variances for `n_comps` components. If `pca_sd` is not `None`, it is assumed that the matrix contains PCA values, else PCA is computed
    """
    
    if isinstance(data, anndata.AnnData):
        matrix = adata.X
        pca_sd = None
    elif isinstance(data, (np.matrix, np.ndarray, sparse.csr_matrix)):
        matrix = data
    else:
        raise TypeError(f'invalid type {data.__class__} for data')
    
    # perform PCA if necessary
    if pca_sd is None:
        if verbose:
            print("PCA")
            
        if n_comps is None or n_comps > min(matrix.shape):
            n_comps = min(matrix.shape)

        if n_comps == min(matrix.shape):
            svd_solver = 'full'
    
        pca = sc.tl.pca(matrix, n_comps=n_comps, use_highly_variable=False,
                        return_info=True, svd_solver=svd_solver, copy=True)
        X_pca = pca[0].copy()
        pca_sd = pca[3].copy()
        del pca
    else:
        X_pca = matrix
        n_comps = matrix.shape[1]
    
    ## PC Regression
    if verbose:
        print("PC regression")    
                
    batch = pd.get_dummies(batch) if 'category' == str(batch.dtype) else np.array(batch)
    
    # fit linear model for n_comps PCs
    from sklearn.linear_model import LinearRegression
    r2 = []
    for i in range(n_comps):
        lm = LinearRegression()
        lm.fit(X_pca[:, [i]], batch)
        r2.append(lm.score(X_pca[:, [i]], batch))
    
    Var = pca_sd**2 / sum(pca_sd**2) * 100
    R2Var = sum(r2*Var)/100
    
    return R2Var

### kBET
def kBET_single(matrix, batch, subsample=0.5, heuristic=True, verbose=False):
    """
    params:
        matrix: count matrix
        batch: series or list of batch assignemnts
        subsample: fraction to be subsampled. No subsampling if `subsample=None`
    returns:
        kBET p-value
    """
    if isinstance(subsample, float):
        matrix, indices = sc.pp.subsample(matrix, fraction=subsample, copy=True)
        batch = batch[indices]
    
    anndata2ri.activate()
    ro.r("library(kBET)")
    
    if verbose:
        print("importing count matrix")
    ro.globalenv['data_mtrx'] = matrix
    ro.globalenv['batch'] = batch
    
    if verbose:
        print("kBET estimation")
    k0 = 1 if len(batch) < 50 else 'NULL'
    batch_estimate = ro.r(f"batch.estimate <- kBET(data_mtrx, batch, plot=FALSE, k0={k0}, heuristic={str(heuristic).upper()}, verbose={str(verbose).upper()})")
    
    anndata2ri.deactivate()
    return ro.r("batch.estimate$average.pval")[0]

def kBET(adata, matrix, covariate_key='sample', cluster_key='louvain',
                    hvg=False, subsample=0.5, heuristic=False, verbose=False):
    """
    Compare the effect before and after integration
    params:
        matrix: matrix from adata to calculate on
    return:
        pd.DataFrame with kBET p-values per cluster for batch
    """
    
    checkAdata(adata)
    checkBatch(covariate_key, adata.obs)
    checkBatch(cluster_key, adata.obs)
    
    if hvg:
        hvg_idx = get_hvg_indices(adata)
        if verbose:
            print(f"subsetting to {len(hvg_idx)} highly variable genes")
        matrix = matrix[:, hvg_idx]
    
    if verbose:
        print(f"covariate: {covariate_key}")
    batch = adata.obs[covariate_key]
    
    kBET_scores = {'cluster': [], 'kBET': []}
    for clus in adata.obs[cluster_key].unique():
        if verbose:
            print(f'cluster {clus}')
        idx = np.where((adata.obs[cluster_key] == clus))[0]
        score = kBET_single(
            matrix[idx, :],
            batch[idx],
            subsample=subsample,
            verbose=verbose,
            heuristic=heuristic
        )
        kBET_scores['cluster'].append(clus)
        kBET_scores['kBET'].append(score)
    
    kBET_scores = pd.DataFrame.from_dict(kBET_scores)
    kBET_scores = kBET_scores.reset_index(drop=True)
    return kBET_scores


def kBET_comparison(adata, raw, corrected, covariate_key='sample', cluster_key='louvain', hvg=False, subsample=0.5, heuristic=False, verbose=False):
    """
    Compare the effect before and after integration
    params:
        raw: count matrix before integration
        corrected: count matrix after correction
    return:
        pd.DataFrame with difference of kBET p-values
    """
    
    checkAdata(adata)
    checkBatch(covariate_key, adata.obs)
    checkBatch(cluster_key, adata.obs)
    
    kBET_before = kBET(adata, raw, 
                       covariate_key=covariate_key,
                       cluster_key= cluster_key,
                       subsample=subsample,
                       heuristic=heuristic,
                       verbose=verbose)
    kBET_after = kBET(adata, corrected,
                      covariate_key=covariate_key,
                      cluster_key= cluster_key,
                      subsample=subsample,
                      heuristic=heuristic,
                      verbose=verbose)
    kBET_scores = kBET_before.merge(kBET_after, on='cluster', suffixes=('_before','_after'))
    kBET_scores['difference'] = kBET_scores['kBET_before'] - kBET_scores['kBET_after']
    return kBET_scores

### Time and Memory
def measureTM(*args, **kwargs):
    """
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
    mem = np.max(out[0])- out[0][0]
    print(f'memory usage:{round(mem,0) } MB')
    print(f'runtime: {round(Stats(prof).total_tt,0)} s')
    return mem, Stats(prof).total_tt, out[1:]


### All Metrics
def metrics(adata_dict):
    """
    summary of all metrics for all tools in a DataFrame
    params:
        adata_dict: list of adata results from different integration methods
            ["seurat", "scanorama", "mnn", "scGen", "Harmony", "CONOS"]
        
    """
    
    metrics = {}
    for tool, adata in adata_dict.items():
        metrics[tool] = metrics_per_tool(adata)
    
    return pd.DataFrame(metrics) #TODO: name columns

def metrics_per_tool(adata, matrix=None,
                     silhouette=True, si_batch='tissue', si_group='cell_type', si_embed='X_pca', 
                     nmi=True, ari=True, group1='cell_type', group2='louvain_post', nmi_method='max', nmi_dir=None, 
                     cell_cycle=True, s_phase_key='S_score', g2m_phase_key='G2M_score',
                     hvg = True,
                     kBET=True, kBET_batch='sample', kBET_cluster='louvain', kBET_sub=0.5,
                    ):
    """
    summary of all metrics for one Anndata object
    params:
        adata:
        silhouette: compute silhouette score on batch `si_batch`, `si_group` using the embedding `si_embed` (check `silhouette_score` function for details)
        nmi: compute normalized mutual information NMI
    """
    
    if not matrix:
        matrix = adata.X
    metrics = {}
    if silhouette:
         metrics['silhouette'] =  silhouette_score(adata, batch=si_batch, group=si_group,
                                         metric='euclidean', embed=si_embed, verbose=False)
    if nmi:
        metrics['NMI'] = nmi(adata, group1, group2, method=nmi_method, nmi_dir=nmi_dir)
    if ari:
        metrics['ARI'] = ari(adata, group1, group2, method=nmi_method)
    if cell_cycle:
        metrics['S-phase'], metrics['G2M-phase'] = cell_cycle(
            adata, hvg=hvg, s_phase_key='S_score', g2m_phase_key='G2M_score')
    if kBET:
        kbet_scores = kBET(adata, matrix, covariate_key=kBET_batch, cluster_key=kBET_cluster,
                           hvg=hvg, subsample=kBET_sub, heuristic=True, verbose=False)
        metrics['kBET'] = kbet_scores['kBET'].mean()
        
    metrics['HVG'] = hvg
    return metrics

