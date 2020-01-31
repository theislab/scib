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
    
    #scale by scenario
    for scen in scenario.cat.categories:
        scen_idx = res.index[scenario == scen]
        res_tmp = res.loc[scen_idx]
        #select columns with LISI score and scale (assuming scale=True in metrics call)
        #original cLISI in [0,1] with 1 good and 0 bad
        #scale to 1 good and 0 bad
        max_cl = res_tmp['cLISI'].max() #take the max of observed score
        min_cl = res_tmp['cLISI'].min() #take the min of observed score
        intval_clen = max_cl-min_cl
        #divide by 1 if all scores are the same, otherwise divide by max-min
        divisor = [1 if intval_clen < 1e-10 else intval_clen]
        #scale cLISI
        res.loc[scen_idx,'cLISI'] = (res.loc[scen_idx,'cLISI']-min_cl)/divisor
        #original iLISI in [0,max] with 0 bad and max good
        max_il = res_tmp['iLISI'].max()
        min_il = res_tmp['iLISI'].min()
        intval_ilen = max_il-min_il
        #divide by 1 if all scores are the same, otherwise divide by max-min
        divisori = [1 if intval_ilen < 1e-10 else intval_ilen]
        #scale iLISI
        res.loc[scen_idx,'iLISI']=(res.loc[scen_idx,'iLISI']-min_il)/divisori
    
    results = res #do not transpose the file
    #write output file
    results.to_csv(args.output)
