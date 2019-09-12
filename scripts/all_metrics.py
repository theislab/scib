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
    import os
    import glob

    parser = argparse.ArgumentParser(description='Collect all metrics')

    parser.add_argument('-i', '--input', required=True, help='input directory')
    parser.add_argument('-i', '--input', required=True, help='file')
    args = parser.parse_args()
    
    ##
    
    adata_dict = {}
    files = glob.glob(f'{args.input}*metrics*')
    for file in files:
        adata = sc.read(file)
        nice_name = '_'.join(file.split('_')[2:3])
        adata_dict[nice_name] = adata
    
    print(adata_dict)
    
    # call all_metrics function
    results = scIB.me.all_metrics(adata_dict)
    results.to_csv(args.output)
    