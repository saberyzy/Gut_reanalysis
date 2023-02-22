#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""
This script formatting the data from the cellranger aggr and filter cells based on custom cutoffs

author: Zhiyuan Yao

Update Date: Aug-10-2022
"""
import scanpy as sc 
import anndata as ann
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
import logging
import os
import glob
import matplotlib
import math


# In[2]:


def read_data(path):
    """read mtx data and return the anndata object
    
    Parameters:
    ----------
    path: str
        The path where cellranger aggr filtered matrix located

    Returns:
        adata: scanpy.adata
        scanpy adata object
    """     
    adata = sc.read_10x_mtx(path)
    return(adata)

def add_s_info(adata, csv_path):
    """reformat the cell barcode info, add the sample info to the metadata
    
    Parameters:
    ----------
    adata: scanpy.adata
        scanpy adata object
    csv_path: str
        path of the aggr csv
        
    Returns:
        adata: scanpy.adata
        scanpy adata object
    """     
    sample_info = pd.read_csv(csv_path, sep = ',')
    adata.obs['cell_barcode'] = adata.obs.index
    #add sample id 
    adata.obs['sample_id'] = ''
    for ind in list(sample_info.index):
        n = str(ind + 1)
        s_name = sample_info.loc[ind, 'sample_id']
        adata.obs.loc[adata.obs['cell_barcode'].str.endswith(n), 'sample_id'] = s_name
    #reformat cell barcodes
    adata.obs.cell_barcode = adata.obs.cell_barcode.str.split('-', expand = True)[0].tolist()
    adata.obs.cell_barcode = adata.obs.cell_barcode+'-1-'+ adata.obs.sample_id
    adata.obs.set_index('cell_barcode', drop = True, inplace = True)
    
    return(adata)

def qc_metric(adata):
    """calculate the qc metric and store the results in the obs and var
    
    Parameters:
    ----------
    adata: scanpy.adata
        scanpy adata object
        
    Returns:
        adata: scanpy.adata
        scanpy adata object
    """      
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=[50], log1p=False, inplace=True)
    return(adata)


def filter_data(adata, min_counts = None, max_counts = None, mt_frac = None, min_genes = None, max_genes = None):
    """filter the data based on the custom filters, mark the filtered cells in the obs
    
    Parameters:
    ----------
    adata: scanpy.adata
        scanpy adata object
    min_counts: int
        any integers > 0 
    max_counts: int
        any integers > 0
    mt_frac: num
        ranger from 0 to 100
    min_genes: int
        any integers > 0 
    max_genes: int
        any integers > 0 
        
    Returns:
        adata: scanpy.adata
        scanpy adata object, the unfiltered object with a column in the obs marking filtered cells
        temp: scanpy.adata
        scanpy adata object, the anndata object only contains filtered cells
    """ 
    temp = adata.copy()
    print('Total number of cells: {:d}'.format(temp.n_obs))
    if min_counts is not None:
        sc.pp.filter_cells(temp, min_counts = min_counts)
        print('Number of cells after min count filter: {:d}'.format(temp.n_obs))
    if max_counts is not None:
        sc.pp.filter_cells(temp, max_counts = max_counts)
        print('Number of cells after max count filter: {:d}'.format(temp.n_obs))
    if min_genes is not None:
        sc.pp.filter_cells(temp, min_genes = min_genes)
        print('Number of cells after min genes filter: {:d}'.format(temp.n_obs))
    if max_genes is not None:
        sc.pp.filter_cells(temp, max_genes = max_genes)
        print('Number of cells after max genes filter: {:d}'.format(temp.n_obs))
    if mt_frac is not None:
        temp = temp[temp.obs['pct_counts_mt'] < mt_frac]
        print('Number of cells after MT filter: {:d}'.format(temp.n_obs))
        
    adata.obs['cell_quality'] = 'low'
    adata.obs.loc[temp.obs.index, 'cell_quality'] = 'high'
    
    return(adata, temp)    

def write_adata(adata, savepath):
    """write adata to the h5ad file 
    
    Parameters:
    ----------
    adata: scanpy.adata
        scanpy adata object
    savepath: str
        the path where the h5ad will be stored
        
    Returns:
        none
    """
    adata.write(savepath)
    


# In[3]:


if __name__ == '__main__':
    #load merged data
    #mtx path
    mtxpath = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/aggr/Gut_merge/outs/count/filtered_feature_bc_matrix'
    #aggregation csv path
    csv_path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/aggr/Gut_merge/outs/aggregation.csv'
    #h5ad save path
    savepath = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/aggr/processed_data/pp_2.h5ad'
    
    adata = read_data(mtxpath)
    adata = add_s_info(adata, csv_path)
    
    
    adata = qc_metric(adata)
#this is for the paper's cutoff
#     adata, adata_2 = filter_data(adata, 
#                              min_counts = None, 
#                              max_counts = None, 
#                              mt_frac = 10, 
#                              min_genes = 600, 
#                              max_genes = 10000)
    
#     write_adata(adata, savepath)
#this is for our groups's cutoff
    adata, adata_2 = filter_data(adata, 
                             min_counts = 400, 
                             max_counts = None, 
                             mt_frac = 25, 
                             min_genes = 200, 
                             max_genes = None)
    
    write_adata(adata, savepath)




# In[ ]:





# In[ ]:




