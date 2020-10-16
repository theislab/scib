#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import argparse
import os
import warnings
warnings.filterwarnings('ignore')

RESULT_TYPES = [
    "full",  # reconstructed expression data
    "embed", # embedded/latent space
    "knn"    # corrected neighbourhood graph
]

if __name__=='__main__':
    """
    Save embeddings for all scenarios, methods and settings
    """

    parser = argparse.ArgumentParser(description='Save embeddings')

    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--outfile', required=True,
                        help='Output coordinates CSV file')
    parser.add_argument('-m', '--method', required=True, help='Name of method')
    parser.add_argument('-r', '--result', required=True, choices=RESULT_TYPES,
                        help='Type of result: full, embed, knn')
    parser.add_argument('-b', '--batch_key', required=True, help='Key of batch')
    parser.add_argument('-l', '--label_key', required=True,
                        help='Key of annotated labels e.g. "cell_type"')

    args = parser.parse_args()

    outfile = args.outfile
    method = args.method
    result = args.result
    setup = f'{method}_{result}'
    batch_key = args.batch_key
    label_key = args.label_key
    adata = sc.read_h5ad(args.input)

    print('Preparing dataset...')
    if result == 'embed':
        scIB.pp.reduce_data(adata, n_top_genes=None, neighbors=True,
                            use_rep='X_emb', pca=False, umap=False)
    elif result == 'full':
        sc.pp.filter_genes(adata, min_cells=1)
        scIB.pp.reduce_data(adata, n_top_genes=2000, neighbors=True,
                            use_rep='X_pca', pca=True, umap=False)

    # Calculate embedding
    if args.method.startswith('conos'):
        print('Calculating graph embedding...')
        sc.tl.draw_graph(adata, key_added_ext='graph')
        basis = 'draw_graph_graph'
        label = 'Graph'
    else:
        print('Calculating UMAP...')
        sc.tl.umap(adata)
        basis = 'umap'
        label = 'UMAP'

    # Save embedding plots
    outdir = os.path.dirname(outfile)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    print(f'Saving embedding plot for labels "{label_key}"...')
    fig = sc.pl.embedding(adata, basis=basis, color=label_key, return_fig=True)
    fig.set_size_inches(10, 8)
    fig.savefig(os.path.join(outdir, f'{setup}_labels.png'),
                bbox_inches='tight')
    print(f'Saving embedding plot for batches "{batch_key}"...')
    fig = sc.pl.embedding(adata, basis=basis, color=batch_key, return_fig=True)
    fig.set_size_inches(10, 8)
    fig.savefig(os.path.join(outdir, f'{setup}_batch.png'),
                bbox_inches='tight')

    # Save embedding coordinates
    print('Saving embedding coordinates...')
    adata.obs[label + '1'] = adata.obsm['X_' + basis][:, 0]
    adata.obs[label + '2'] = adata.obsm['X_' + basis][:, 1]
    coords = adata.obs[[label_key, batch_key, label + '1', label + '2' ]]

    coords.to_csv(os.path.join(outfile), index_label='CellID')
