#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import argparse
import os
import sys
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

    args = parser.parse_args()

    outfile = args.outfile

    outdir = os.path.dirname(outfile)

    if os.stat(args.input).st_size == 0:
        print(f'{args.input} is empty, empty output files will be created')

        # Create empty output files
        for file in [labels_png, batches_png, outfile]:
            print(f'Creating empty {file}')
            open(file, "w").close()

        # End script
        sys.exit()

    adata = sc.read_h5ad(args.input, backed='r')

    print('Saving embedding coordinates...')
    adata.obs.select_dtypes(include='category').to_csv(os.path.join(outfile), index_label='CellID')
