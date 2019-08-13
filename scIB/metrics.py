import scanpy as sc
import pandas as pd
import seaborn as sns
from scIB.utils import *
import cProfile
from pstats import Stats
from memory_profiler import profile
import memory_profiler
import scIB
import timeit
import numpy as np

sns.set_context('talk')
sns.set_palette('Dark2')

### Silhouette score
def silhouette_score(adata, batch='method', group='cell_ontology_class', metric='euclidean', embed='X_pca', verbose=True):
    """
    Silhouette score subsetted for each cluster (group) between batches.
    This results in 1 score per group label
    params:
        adata: 
        batch: batches to be compared against
        group: group labels to be subsetted by e.g. cell type
        metric: 
        embed: name of column in adata.obsm
    returns:
        per_group: scores per group label
        sil_means: mean silhouette score of group means
    """
    import sklearn.metrics as scm
    
    checkAdata(adata)
    checkBatch(batch, adata.obs)
    checkBatch(group, adata.obs)
    
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    
    # ony look at group labels that are present for all batches
    n_batch = adata.obs[batch].nunique()
    labels = adata.obs.groupby(group)[batch].nunique()
    labels = labels[labels == n_batch]
    
    sil_means = []
    per_group = []
    for j in labels.index:
        tmp_type = adata[adata.obs[group] == j]
        sil = scm.silhouette_score(tmp_type.obsm[embed], tmp_type.obs[batch], metric=metric)
        sil_means.append(sil)
        per_group.extend(scm.silhouette_samples(tmp_type.obsm[embed], tmp_type.obs[batch], metric=metric))
    per_group = [abs(i) for i in per_group] # take only absolute value
    
    if verbose:
        print(f'mean silhouette over label means: {np.mean(sil_means)}')
        print(f'mean silhouette per cell: {np.mean(per_group)}')
    
    return per_group, sil_means

def plot_silhouette_score(adata_dict, verbose=True):
    """
    params:
        adata_dict: dictionary of adata objects, each labeled by e.g. integration method name
    """
    for label, adata in adata_dict.items():
        checkAdata(adata)
        per_group, sil_means = silhouette_score(adata, verbose=verbose)
        sns.distplot(per_group, label=label, hist=False)

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
    sns.boxplot(data=clust_df)
    sns.swarmplot(data=clust_df, color=".25")
    
    if df:
        return clust_df
    return None


### NMI normalised mutual information
'''
def nmi(adata, labelColumn, res=0.5):
    import sklearn.metrics as scm
    sc.tl.louvain(adata, resolution=res, key_added='louvain_post')
    labels_pre = adata.obs[labelColumn]
    labels_post = adata.obs['louvain_post']
    
    return scm.normalized_mutual_info_score(labels_pre, labels_post)
'''

def nmi(adata, group1, group2, average_method='max'):
    """
    """
    from sklearn.metrics import normalized_mutual_info_score
    
    checkAdata(adata)
    checkBatch(group1, adata.obs)
    checkBatch(group2, adata.obs)
    
    group1_list = adata.obs[group1].tolist()
    group2_list = adata.obs[group2].tolist()
    
    return normalized_mutual_info_score(group1_list, group2_list, average_method=average_method)
    

def onmi(adata, group1, group2, onmi_dir="../../Overlapping-NMI/", verbose=False):
    """
    compute normalized mutual information based on 2 different cluster assignments
    runs the compiled onmi C code
    """
    import subprocess
    import os
    
    group1_file = write_tmp_labels(adata, group1, to_int=False)
    group2_file = write_tmp_labels(adata, group2, to_int=False)
    
    nmi_call = subprocess.Popen(
        [onmi_dir+"onmi", group1_file, group2_file], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT)
    
    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)
    
    nmi_out = stdout.decode()
    if verbose:
        print(nmi_out)
    
    nmi_split = [x.strip().split('\t') for x in nmi_out.split('\n')]
    nmi_max = nmi_split[0][1]
    
    # remove temporary files
    os.remove(group1_file)
    os.remove(group2_file)
    
    return nmi_max


def nmi_Lanc(adata, group1, group2, nmi_dir="../../mutual3/"):
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
    
    return nmi_out.split('\t')[1]

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
    

def measureTM(*args, **kwargs):
    prof = cProfile.Profile()
    out = memory_profiler.memory_usage((prof.runcall, args, kwargs), retval=True) 
    mem = np.max(out[0])- out[0][0]
    print(f'memory usage:{round(mem,0) } MB')
    print(f'runtime: {round(Stats(prof).total_tt,0)} s')
    return mem, Stats(prof).total_tt, out[1:]

def pcr_hvg(pre, post, n_hvg, batch):
    import rpy2.rinterface_lib.callbacks
    import logging
    rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Ignore R warning messages
    import rpy2.robjects as ro
    import anndata2ri
    anndata2ri.activate()
    ro.numpy2ri.activate()
    cons = []
    for x in pre.obs['batch'].unique():
        tmp = pre[pre.obs['batch']==x]
        tmp_post = post[post.obs['batch']==x]
        
        ro.globalenv['tmp'] = tmp.X
        ro.globalenv['tmp_post'] = tmp_post.X

        ro.globalenv['pca.data'] = ro.r("pca.data <- prcomp(tmp, center=TRUE)")
        ro.globalenv['pca_post.data'] = ro.r("pca_post.data <- prcomp(tmp_post, center=TRUE)")
        
        hvg = pd.DataFrame.from_records(sc.pp.highly_variable_genes(tmp, n_top_genes=n_hvg, inplace=False)).highly_variable
        hvg_1 = tmp[:,hvg]
        hvg_2 = tmp_post[:,hvg]
        pcr_pre_all = []
        pcr_post_all = []
        summ = 0
        for i in range(hvg_1.shape[1]):
            ro.globalenv['y_1']= hvg_1[:,i].X
            ro.globalenv['y_2']= hvg_2[:,i].X
            pcr_pre = ro.r("batch.pca <- pcRegression(pca.data, y_1, n_top=50)")
            pcr_post = ro.r("batch.pca <- pcRegression(pca_post.data, y_2, n_top=50)")
            pcr_preV = dict(zip(pcr_pre.names, list(pcr_pre)))
            pcr_postV = dict(zip(pcr_post.names, list(pcr_post)))
            pcr_pre_all.append(pcr_preV['pcRegscale'])
            pcr_post_all.append(pcr_postV['pcRegscale'])
            print(i)
        for i in range(len(pcr_pre_all)):
            diff = (pcr_pre_all[i]-pcr_post_all[i])**2
            summ = summ + diff
        cons.append(summ)
    anndata2ri.deactivate()
    ro.numpy2ri.deactivate()
    return np.mean(cons)