#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 15:18:30 2021

@author: polaris
"""

import numpy as np
import pandas as pd
import sys
sys.path.append('../../lib/')
from regimeshifts import ews


def significance_rob(ts,max_wL=100,max_bW=80,n=1000, res=5, trend='positive'):
    """
    Estimates the p-values associated to each combination of window length and
    bandwidth. 
    Parameters
    ------------
    ts: Pandas time-series
    max_wL: integer
            Maximum window length.
    max_bW: integer
            Maximum bandwidth.
    n: integer
       Number of surrogate series in the ensemble
    res: integer
        resolution of the significant values over the parameter space
    trend: string ('positive'|'negative')
        If a positive or negative trend is expected, the function estimates the
        p-value accordingly
            
    """
    bW_vs = np.arange(5,max_bW,res)    
    wL_vs = np.arange(20,max_wL,res)
    pval_ar1 = np.full([len(bW_vs),len(wL_vs)],np.nan)
    pval_var = np.full([len(bW_vs),len(wL_vs)],np.nan)
    for i,b in enumerate(bW_vs):
        for j,w in enumerate(wL_vs): 
            print('AR(1)')
            sig_test_ar1 = ts.significance(indicator='ar1',n=n, detrend=True,bW=b, wL=w,test=trend)
            print('Variance')
            sig_test_var = ts.significance(indicator='var',n=n, detrend=True,bW=b, wL=w,test=trend)
            pval_ar1[i][j] = sig_test_ar1.pvalue
            pval_var[i][j] = sig_test_var.pvalue
    return bW_vs, wL_vs, pval_ar1, pval_var

## The records are publicly available at https://www.ncei.noaa.gov/access/paleo-search/
rec_ids = {'GI':'NIS_PB_RCSars','d18O':'NIS_DR_d18O','d13C':'NIS_DR_d13C'}
data_path = 'data/bivalve_rec/'
md = pd.read_csv(data_path+'chronologies_metadata.csv',index_col=0) ## Records metadata 
series = []
for chron,rec_id in rec_ids.items():    
    df = pd.read_csv(data_path+md.loc[rec_id].file, comment='#', sep='\t', index_col=0) ## Reads the data-file for each record
    ts = df[md.loc[rec_id]['column_name']]
    ts.rename(chron,inplace=True)
    ts.index = ts.index.astype(int)
    ts.index.name = 'year_ce'
    series.append(ts.dropna())
records = pd.concat(series,axis=1,join='outer') ## Merges the series in a dataframe
records = ews.Ews(records[records.index>999]) ## Filters out the data prior to the year 1000 CE


## Episodes
episodes = {'1':[1110,1260,'positive'],
            #'int':[1190,1330,'negative'],
            '2':[1260,1380,'positive']}

model_size = 2000

#just_gi = {'GI':datasets['GI']}
for ep,years in episodes.items():
    print(ep)
    df_ep = records[(records.index>=years[0])&(records.index<=years[1])] # Selects the correspondent interval in each record
    for rec in rec_ids.keys():
        print(rec)
        ts = ews.Ews(df_ep[rec])
        bW_vs, wL_vs, par1, pvar = significance_rob(ts,max_wL=100,max_bW=80,n=model_size,trend=years[2])
        np.save('data/signif-tests/bW_vs.npy', bW_vs)
        np.save('data/signif-tests/wL_vs.npy', wL_vs)
        np.save(f'data/signif-tests/par1_{rec}_{ep}.npy', par1)
        np.save(f'data/signif-tests/pvar_{rec}_{ep}.npy', pvar)
        del bW_vs, wL_vs, par1, pvar





