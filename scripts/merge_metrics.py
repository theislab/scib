#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import scIB
import warnings
warnings.filterwarnings('ignore')
import argparse
from functools import reduce

if __name__=='__main__':
    """
    Merge metrics output for all scenarios, methods and settings
    """
    
    parser = argparse.ArgumentParser(description='Collect all metrics')

    parser.add_argument('-i', '--input', nargs='+', required=True, help='input directory')
    parser.add_argument('-o', '--output', required=True, help='output file')
    parser.add_argument('-r', '--root', required=True, 
                        help='root directory for inferring column names from path')
    args = parser.parse_args()
    
    
    res_list = []
    for file in args.input:
        clean_name = file.replace(args.root, "").replace(".csv", "")
        res = pd.read_csv(file, index_col=0)
        res.rename(columns={res.columns[0]: clean_name}, inplace=True)
        res_list.append(res)
    
    results = reduce(lambda  left,right: pd.merge(left, right, left_index=True, right_index=True), res_list)
    results = results.T
    results.to_csv(args.output)
    
