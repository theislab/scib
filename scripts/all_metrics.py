#!/usr/bin/env python
# coding: utf-8

import pandas as pd
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
    parser.add_argument('-o', '--output', required=True, help='output file')
    args = parser.parse_args()
    
    ##
    
    res_dict = {}
    files = glob.glob(f'{args.input}*metrics*')
    print(files)
    for file in files:
        res = pd.read_csv(file)
        nice_name = '_'.join(file.split('_')[2:4])
        print(nice_name)
        res_dict[nice_name] = res
    
    print(res_dict)
    
    # call all_metrics function
    results = scIB.me.metrics_all(res_dict)
    results.to_csv(args.output)
    
