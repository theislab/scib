import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from scIB import utils
from scIB import metrics


def opt_louvain(adata, label='cell_type', cluster_key='louvain', resolutions=None, nmi_method='max', nmi_dir=None, inplace=True, plot=True, force=False, verbose=True):
    """
    returns:
        res_max: resolution of maximum NMI
        nmi_max: maximum NMI score
        nmi_all: `pd.DataFrame` containing all NMI at resolutions. Can be used to plot the NMI profile.
        clustering: only if `inplace=False`, return cluster assingment as `pd.Series`
        plot: if `plot=True` plot the NMI profile over resolution
    """
    
    if cluster_key in adata.obs.columns:
        if force:
            if verbose:
                print(f"Warning: cluster key {cluster_key} already exists in adata.obs and will be overwritten")
        else:
            raise ValueError(f"cluster key {cluster_key} already exists in adata, please remove the key or choose a different name. If you want to force overwriting the key, specify `force=True`")
    
    if not resolutions:
        n = 20
        resolutions = [2*x/n for x in range(1,n+1)]
    
    nmi_max = 0
    res_max = resolutions[0]
    clustering = None
    nmi_all = []
    
    for res in resolutions:
        sc.tl.louvain(adata, resolution=res, key_added=cluster_key)
        nmi = metrics.nmi(adata, group1=label, group2=cluster_key, method=nmi_method, nmi_dir=nmi_dir)
        nmi_all.append(nmi)
        if nmi_max < nmi:
            nmi_max = nmi
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]
    
    nmi_all = pd.DataFrame(zip(resolutions, nmi_all), columns=('resolution', 'NMI'))
    if plot:
        # NMI vs. resolution profile
        sns.lineplot(data= nmi_all, x='resolution', y='NMI').set_title('NMI profile')
        plt.show()
    
    if inplace:
        adata.obs[cluster_key] = clustering
        return res_max, nmi_max, nmi_all
    else:
        return res_max, nmi_max, nmi_all, clustering
