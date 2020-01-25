#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import scIB
import warnings
warnings.filterwarnings('ignore')
import argparse

if __name__=='__main__':
    """
    Merge results for cell cycle variance for all scenarios, methods and settings
    """

    parser = argparse.ArgumentParser(description='Collect all metrics')

    parser.add_argument('-i', '--input', nargs='+', required=True, help='input directory')
    parser.add_argument('-o', '--output', required=True, help='output file')
    parser.add_argument('-r', '--root', required=True, 
                        help='root directory for inferring column names from path')
    args = parser.parse_args()
    
    res_list = []
    for file in args.input:
        res = pd.read_csv(file)
        clean_name = file.replace(args.root, "").replace(".csv", "")
        res['setting'] = clean_name
        res_list.append(res)
    
    results = pd.concat(res_list)
    print(results)
    results.to_csv(args.output, index=False)

