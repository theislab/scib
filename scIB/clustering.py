import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from scIB import utils
from scIB import metrics


def opt_louvain(adata, label='cell_type', resolution=None, nmi_method='max', nmi_dir=None, inplace=True, plot=False):
    """
    returns:
        res_max: resolution of maximum NMI
        nmi_max: maximum NMI score
        nmi_all: `pd.DataFrame` containing all NMI at resolutions. Can be used to plot the NMI profile.
        clustering: only if `inplace=False`, return cluster assingment as `pd.Series`
        plot: if `plot=True` plot the NMI profile over resolution
    """
    
    if not resolution:
        n = 20
        resolution = [x/n for x in range(1,n+1)]
    
    nmi_max = 0
    res_max = resolution[0]
    clustering = None
    nmi_all = []
    
    for res in resolution:
        sc.tl.louvain(adata, resolution=res, key_added='louvain')
        nmi = metrics.nmi(adata, group1=label, group2='louvain', method=nmi_method, nmi_dir=nmi_dir)
        nmi_all.append(nmi)
        if nmi_max < nmi:
            nmi_max = nmi
            res_max = res
            clustering = adata.obs['louvain']
        del adata.obs['louvain']
    
    nmi_all = pd.DataFrame(zip(resolution, nmi_all), columns=('resolution', 'NMI'))
    if plot:
        # NMI vs. resolution profile
        sns.lineplot(data= nmi_all, x='resolution', y='NMI').set_title('NMI profile')
        plt.show()
    
    if inplace:
        adata.obs['louvain'] = clustering
        return res_max, nmi_max, nmi_all
    else:
        return res_max, nmi_max, nmi_all, clustering
