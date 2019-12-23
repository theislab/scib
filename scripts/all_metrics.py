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
    from functools import reduce

    parser = argparse.ArgumentParser(description='Collect all metrics')

    parser.add_argument('-i', '--input', required=True, help='input directory')
    parser.add_argument('-o', '--output', required=True, help='output file')
    args = parser.parse_args()
    
    ##
    
    res_list = []
    results = pd.DataFrame()
    files = glob.glob(os.path.join(args.input, '*metrics*'))
    for file_ in files:
        res = pd.read_csv(file_, index_col=0)
        res_list.append(res)
    
    results = reduce(lambda  left,right: pd.merge(left, right, left_index=True, right_index=True), res_list)
    results.to_csv(args.output)
    
