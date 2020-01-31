#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import scIB
import warnings
warnings.filterwarnings('ignore')
import argparse
from functools import reduce

if __name__=='__main__':
    """
    Scale LISI scores in the output of merge_metrics.py
    """
    
    parser = argparse.ArgumentParser(description='Scale LISI scores')

    parser.add_argument('-i', '--input', required=True, help='Merged metrics file in csv file format')
    parser.add_argument('-o', '--output', required=True, help='Merged and LISI-scaled metrics file in csv file format')
    
    args = parser.parse_args()
    
    #get input file
    file = args.input
    res = pd.read_csv(file, index_col=0)
    #get scenario info
    scenario= pd.Series([idx.split('/')[0] for idx in res.index], dtype='category')
    
    #select columns with LISI score, PCR and cc variance (halfopen support)
    columns_oi = ['cLISI', 'iLISI']#, 'PCR_batch', 'cell_cycle_conservation']
    #select only those columns, which are present in the input file
    columns = res.columns[np.in1d(res.columns, columns_oi)]
    
    #scale by scenario
    for scen in scenario.cat.categories:
        scen_idx = res.index[scenario == scen]
        res_tmp = res.loc[scen_idx]
        #scale selected metrics (assuming scale=True in metrics call)
        for column in columns:
            #scale to 1 good and 0 bad
            max_cl = res_tmp[column].max() #take the max of observed score
            min_cl = res_tmp[column].min() #take the min of observed score
            intval_len = max_cl-min_cl
            #divide by 1 if all scores are the same, otherwise divide by max-min
            divisor = [1 if intval_len < 1e-10 else intval_len]
            #scale score
            res.loc[scen_idx,column] = (res.loc[scen_idx,column]-min_cl)/divisor
    
    results = res #do not transpose the file
    #write output file
    results.to_csv(args.output)
