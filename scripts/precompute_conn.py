#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scipy.io as scio
from scIB.metrics import diffusion_conn 
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# types of integration output
RESULT_TYPES = [
    "full", # reconstructed expression data
    "embed", # embedded/latent space
    "knn" # only corrected neighbourhood graph as output
]

if __name__=='__main__':
    """
    read adata object, precompute diffusion connectivities for knn data integration methods and output connectivity matrix.
    """
    
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Precompute diffusion connectivities for knn data integration methods.')
    
    
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True, help='output directory')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-t', '--type', required=True, choices=RESULT_TYPES, help='Type of result: full, embed, knn\n full: scanorama, seurat, MNN\n embed: scanorama, Harmony\n knn: BBKNN')
    
    args = parser.parse_args()
    
    verbose = args.verbose
    type_ = args.type
        
    # set prefix for output and results column name
    base = os.path.basename(args.input).split('.h5ad')[0]
    
    
    if verbose:
        print('Options')
        print(f'    type:\t{type_}')
    
      
    
    ###
 
    print("reading adata input file")
    if os.stat(args.input).st_size>0:
        adata = sc.read(args.input, cache=True)
        print(adata)
        if (type_ == 'knn'):
            neighbors = adata.uns['neighbors']['connectivities']
            del adata
            diff_neighbors = diffusion_conn(neighbors, min_k=50, copy=False, max_iterations=20)
            scio.mmwrite(target = os.path.join(args.output, f'{base}_diffconn.mtx'), a = diff_neighbors)
            print("done")
        else:
            print('Wrong type chosen, doing nothing.')
    else:
        print("No file found. Doing nothing.")
    

    
