import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from scIB import utils
from scIB import metrics


def opt_louvain(adata, label_key, cluster_key, function=None, resolutions=None,
                inplace=True, plot=False, force=True, verbose=True, **kwargs):
    """
    params:
        label_key: name of column in adata.obs containing biological labels to be
            optimised against
        cluster_key: name of column to be added to adata.obs during clustering. 
            Will be overwritten if exists and `force=True`
        function: function that computes the cost to be optimised over. Must take as
            arguments (adata, group1, group2, **kwargs) and returns a number for maximising
        resolutions: list if resolutions to be optimised over. If `resolutions=None`,
            default resolutions of 20 values ranging between 0.1 and 2 will be used
    returns:
        res_max: resolution of maximum score
        score_max: maximum score
        score_all: `pd.DataFrame` containing all scores at resolutions. Can be used to plot the score profile.
        clustering: only if `inplace=False`, return cluster assingment as `pd.Series`
        plot: if `plot=True` plot the score profile over resolution
    """
    
    if function is None:
        function = metrics.nmi
    
    if cluster_key in adata.obs.columns:
        if force:
            if verbose:
                print(f"Warning: cluster key {cluster_key} already exists " +
                      "in adata.obs and will be overwritten")
        else:
            raise ValueError(f"cluster key {cluster_key} already exists in " +
                             "adata, please remove the key or choose a different name." +
                             "If you want to force overwriting the key, specify `force=True`")
    
    if resolutions is None:
        n = 20
        resolutions = [2*x/n for x in range(1,n+1)]
    
    score_max = 0
    res_max = resolutions[0]
    clustering = None
    score_all = []
    
    #maren's edit - recompute neighbors if not existing
    try:
        adata.uns['neighbors']
    except KeyError:
        print('computing neigbours for opt_cluster')
        sc.pp.neighbors(adata)

    for res in resolutions:
        sc.tl.louvain(adata, resolution=res, key_added=cluster_key)
        score = function(adata, label_key, cluster_key, **kwargs)
        score_all.append(score)
        if score_max < score:
            scor_max = score
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]
    
    if verbose:
        print(f'optimised clustering against {label_key}')
        print(f'optimal cluster resolution: {res_max}')
        print(f'optimal score: {score_max}')
    
    score_all = pd.DataFrame(zip(resolutions, score_all), columns=('resolution', 'score'))
    if plot:
        # score vs. resolution profile
        sns.lineplot(data= score_all, x='resolution', y='score').set_title('Optimal cluster resolution profile')
        plt.show()
    
    if inplace:
        adata.obs[cluster_key] = clustering
        return res_max, score_max, score_all
    else:
        return res_max, score_max, score_all, clustering
