#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import warnings
warnings.filterwarnings('ignore')

if __name__=='__main__':
    """
    read adata object, compute all metrics and output csv.
    """
    
    import argparse

    parser = argparse.ArgumentParser(description='Compute all metrics')

    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-m', '--metrics', required=True, help='output location of metrics csv')
    parser.add_argument('-n', '--nmi', required=True, help='output location of nmi values of optimised clustering for later plotting')
    parser.add_argument('-b', '--batch_key', required=True, help='Key of batch')
    parser.add_argument('-l', '--label_key', required=True, help='Key of annotated labels e.g. "cell_type"')
    parser.add_argument('-c', '--cluster_key', required=True, help='Name of key to be created for cluster assignment')
    parser.add_argument('-e', '--type', help='Type of result: full, embed, knn')
    parser.add_argument('-s', '--s_phase', default=None, help='S-phase marker genes')
    parser.add_argument('-g', '--g2m_phase', default=None, help='G2-/M-phase marker genes')
    parser.add_argument('-v', '--hvgs', help='Number of highly variable genes', default=None, type=int)
    args = parser.parse_args()
    
    result_types = [
        "full", # reconstructed expression data
        "embed", # embedded/latent space
        "knn" # only corrected neighbourhood graph as output
    ]
    if args.type not in result_types:
        raise ValueError(f'{args.type} is not a valid result type flag')
    
    batch_key = args.batch_key
    label_key = args.label_key
    cluster_key = args.cluster_key
    cc = (args.s_phase is not None) and (args.g2m_phase is not None)
    hvg = args.hvgs is not None
    
    ###
    
    print("reading file")
    adata = sc.read(args.input_file, cache=True)
    print(adata)
    
    if (args.type == "embed"):
    
    print("reducing data")
    scIB.preprocessing.reduce_data(adata,
                                   neighbors=True, pca=True, umap=False,
                                   hvg=hvg, n_top_genes=args.hvgs)
    
    print("clustering")
    res_max, nmi_max, nmi_all = scIB.cl.opt_louvain(adata, 
                        label_key=label_key, cluster_key=cluster_key, 
                        plot=False, force=True, inplace=True)
    # save data for NMI profile plot
    nmi_all.to_csv(args.nmi)
    
    if cc:
        print("scoring cell cycle genes")
        s_genes = open(args.s_phase).read().split('\n')
        s_genes = [gene for gene in s_genes if gene in adata.var_names]
        g2m_genes = open(args.g2m_phase).read().split('\n')
        g2m_genes = [gene for gene in g2m_genes if gene in adata.var_names]
        if len(s_genes)+len(g2m_genes) == 0:
            print('no cell cycle genes in adata, skipping cell cycle effect')
            cc = False
        else:
            sc.tl.score_genes_cell_cycle(adata, s_genes, g2m_genes)
    
    print("computing metrics")
    embed = (args.type == 'embed')
    si_embed = 'X_pca'
    
    results = scIB.me.metrics(adata, matrix=adata.X, hvg=hvg
                    batch_key=batch_key, group_key=label_key, cluster_key=cluster_key,
                    silhouette_=True,  si_embed=si_embed, si_metric='euclidean',
                    nmi_=True, ari_=True, nmi_method='max', nmi_dir=None, 
                    pcr_=True, kBET_=False, cell_cycle_=cc, verbose=False
                    )
    # save metrics' results
    results.to_csv(args.metrics)
    
    print("done")

