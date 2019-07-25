import scanpy as sc
import pandas as pd
import seaborn as sns
import sklearn.metrics as scm
from scIB.utils import *

sns.set_context('talk')
sns.set_palette('Dark2')

def silhouette_score(adata, batch='method', group='cell_ontology_class', metric='euclidean', embed='X_pca'):
    """
    params:
        adata: 
        batch: 
        group: e.g. cell type
        metric: 
        embed: name of column in adata.obsm
    returns:
        per_group: scores per group label
        sil_means: mean silhouette score of group means
    """
    checkAdata(adata)
    if embed not in adata.obsm:
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
    print(f'mean silhouette over label means: {np.mean(sil_means)}')
    per_group = [abs(i) for i in per_group] # take only absolute value
    print(f'mean silhouette per cell: {np.mean(per_group)}')
    
    return per_group, sil_means

def plot_silhouette_score(adata_dict):
    """
    params:
        adata_dict: dictionary of adata objects, each labeled by e.g. integration method name
    """
    for label, adata in adata_dict.items():
        checkAdata(adata)
        per_group, sil_means = silhouette_score(adata)
        sns.distplot(per_group, label=label, hist=False)

def cluster_overlap(adata, group1='louvain', group2='louvain_post'):
    checkAdata(adata)
    cluster_ov = {}
    louv_post_sizes = adata.obs.groupby(group1).size()
    for i in adata.obs[group2].cat.categories:
        a = adata.obs[adata.obs[group2] == i]
        overlap = a.groupby(group1).size()
        cluster_ov[i] = (overlap / louv_post_sizes).sum() / len(overlap[overlap > 0])
    return cluster_ov

def plot_cluster_overlap(adata_dict, group1='louvain', group2='louvain_post'):
    """
    params:
        adata_dict: dictionary of adata objects, each labeled by e.g. integration method name
        group1: 
        group2: 
    return:
        clust_df: dataframe with plotted data points
    """
    series = []
    method_names = list(adata_dict.keys())
    for i in method_names:
        c_ov = cluster_overlap(adata[i], group1=group1, group2=group2)
        series.append(pd.Series(c_ov))
    clust_df = pd.DataFrame(series).transpose()

    clust_df.columns = method_names
    sns.boxplot(data=clust_df)
    sns.swarmplot(data=clust_df, color=".25")
    
    return clust_df
