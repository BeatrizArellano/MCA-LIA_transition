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
            kcb_ar1 = ts.bootstrap(detrend=True, method='ar1',bW=b, wL=w,n=n)
            kcb_var = ts.bootstrap(detrend=True, method='var',bW=b, wL=w,n=n)
            ar1_ts = ts.ar1(detrend = True, bW=b,wL=w)
            var_ts= ts.var(detrend = True, bW=b,wL=w)
            if trend=='positive':
                pval_ar1[i][j] = len(kcb_ar1[kcb_ar1>=ar1_ts.kendall])/len(kcb_ar1)
                pval_var[i][j] = len(kcb_var[kcb_var>=var_ts.kendall])/len(kcb_var)
            else:
                pval_ar1[i][j] = len(kcb_ar1[kcb_ar1<=ar1_ts.kendall])/len(kcb_ar1)
                pval_var[i][j] = len(kcb_var[kcb_var<=var_ts.kendall])/len(kcb_var)
    return bW_vs, wL_vs, pval_ar1, pval_var

## The records are publicly available at https://www.ncei.noaa.gov/access/paleo-search/
rec_ids = {'GI':'NIS_PB_RCSars','d18O':'NIS_DR_d18O','d13C':'NIS_DR_d13C'}
datasets={}
data_path = 'data/bivalve_rec/'
md = pd.read_csv(data_path+'chronologies_metadata.csv',index_col=0) ## Records metadata 
for chron,rec_id in rec_ids.items():    
    df = pd.read_csv(data_path+md.loc[rec_id].file, comment='#', sep='\t', index_col=0) ## Reads the data-file for each record
    ts = df[md.loc[rec_id]['column_name']]
    ts.rename(chron,inplace=True)
    ts.index = ts.index.astype(int)
    ts.index.name = 'year_ce'
    datasets[chron] = ews.Ews(ts.dropna())



## Episodes
episodes = {'1':[1110,1260,'positive'],
            'int':[1190,1330,'negative'],
            '2':[1260,1380,'positive']}

model_size = 2000

#just_gi = {'GI':datasets['GI']}
for ep,years in episodes.items():
    for rec,ds in datasets.items():
#    for rec,ds in just_gi.items():
        print(ep,rec)
        ts = ews.Ews(ds[(ds.index>=years[0])&(ds.index<=years[1])])
        bW_vs, wL_vs, par1, pvar = significance_rob(ts,max_wL=100,max_bW=80,n=model_size,trend=years[2])
        np.save('data/signif-tests/bW_vs.npy', bW_vs)
        np.save('data/signif-tests/wL_vs.npy', wL_vs)
        np.save(f'data/signif-tests/par1_{rec}_{ep}.npy', par1)
        np.save(f'data/signif-tests/pvar_{rec}_{ep}.npy', pvar)
        del bW_vs, wL_vs, par1, pvar





