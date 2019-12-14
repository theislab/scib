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

    parser.add_argument('-i', '--input', nargs='+', required=True, help='input directory')
    parser.add_argument('-o', '--output', required=True, help='output file')
    parser.add_argument('-r', '--root', required=True, 
                        help='root directory for inferring column names from path')
    args = parser.parse_args()
    
    """
    TODO:
        1. remove root from file names
        2. read metric output files and add clean file name as column name
        3. merge columns together
    """
    
    
    res_list = []
    results = pd.DataFrame()
    for file in args.input:
        clean_name = file.replace(args.root, "")
        clean_name = clean_name.replace(".csv", "")
        res = pd.read_csv(file, index_col=0)
        res = res.rename(columns={res.columns[0]: clean_name})
        res_list.append(res)
    
    results = reduce(lambda  left,right: pd.merge(left, right, left_index=True, right_index=True), res_list)
    results = results.T
    print(results.head())
    results.to_csv(args.output)
    
