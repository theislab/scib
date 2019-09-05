#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import warnings
warnings.filterwarnings('ignore')

if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Compute all metrics')

    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-b', '--batch', required=True, help='Batch variable')
    parser.add_argument('-g', '--group', required=True, help='Group variable (cell type or cluster)')
    parser.add_argument('-v', '--hvgs', help='Number of highly variable genes', default=None)
    args = parser.parse_args()
    
    batch = parser.batch
    group = parser.group
    hvg = paser.hvg is not None
    
    adata = sc.read(file = args.input_file, cache=True)
    
    scIB.preprocessing.reduce_data(adata_int, pca=True, umap=True, hvg=hvg, n_top_genes=parser.hvg)
    scIB.me.metrics(adata, matrix=adata.X, batch_key=batch, group_key=group, cluster_key=None,
                     silhouette_=True,  si_embed='X_pca', si_metric='euclidean',
                     nmi_=True, ari_=True, nmi_method='max', nmi_dir=None, 
                     pcr_=True, kBET_=False, cell_cycle_=True, hvg=hvg, verbose=False
                    )

