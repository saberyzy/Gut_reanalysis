#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""
This script rename and reorganize the original fastq files downloaded from the Genome Sequence Archive 
https://ngdc.cncb.ac.cn, accession number HRA001730. Make these files ready for cellranger count alignment

author: Zhiyuan Yao

Update Date: Aug-10-2022
"""
import numpy as np
import pandas as pd
from copy import copy
import os
import glob
import shutil


# In[7]:


def filename_reorg(raw_info, raw_info_2):
    """Rename fastq files to meet cellranger's standard
    
    Parameters:
    ----------
    raw_info: dataframe
        The pandas dataframe that stores the fastq original file name, Accession number, corresponding experimental
        accession number and Run title etc. 
    raw_info_2: dataframe
        The pandas dataframe that stores experimental accession number and its corresponding Sample name and Run title
        raw_info and raw_info_2 are extracted from the excel table file downloaded from the Genome Sequence Archive 
        https://ngdc.cncb.ac.cn, accession number HRA001730

    Returns:
    -------
    raw_info: dataframe
        The pandas dataframe that stores old names and new names for each fastq files
    """
    #get the lane numbers
    lane_number = raw_info['Run title'].str[-1]
    
    #add one column of sample names
    for ind in raw_info.index:
        exp_id = raw_info.loc[ind, 'Experiment accession']
        raw_info.loc[ind, 'sample_name'] = raw_info_2.loc[exp_id, 'BioSample name']
        
    
    raw_info['new_name_R1'] = raw_info['sample_name']+'_S1_L00'+lane_number+'_R1_001.fastq.gz'
    raw_info['new_name_R2'] = raw_info['sample_name']+'_S1_L00'+lane_number+'_R2_001.fastq.gz'
    
    return(raw_info)
    
def make_input_folder(raw_info,path):
    """Generate input folders to store renamed fastq files
    
    Parameters:
    ----------
    raw_info: dataframe
        The pandas dataframe that stores old names and new names for each fastq files, generated from the function 
        filename_reorg  
    path: str
        The root path where the input folders to be created

    Returns:
        None
    """    
    print("making input folders...")
    sample_list = list(raw_info['sample_name'].unique())
    for sample in sample_list:
        folder_dir = os.path.join(path, sample)
        if not os.path.isdir(folder_dir):
            os.mkdir(folder_dir)
        else:
            raise Exception("The folder already exists")
    print("{} input folders has been created".format(len(sample_list)))
        
def mv_input_file(raw_info, src_path, dst_path):
    """Rename and move the fastq files to the corresponding input folders
    
    Parameters:
    ----------
    raw_info: dataframe
        The pandas dataframe that stores old names and new names for each fastq files, generated from the function 
        filename_reorg  
    src_path: str
        The path where the original fastq files located
    dst_path: str
        The root path where the input folders to be created
        
    Returns:
        None
    """ 
    print("moving files from {} to {}".format(src_path, dst_path))
    for ind in raw_info.index:
        filename_1 = raw_info.loc[ind, 'File name 1']
        new_filename_1 = raw_info.loc[ind, 'new_name_R1']
        filename_2 = raw_info.loc[ind, 'File name 2']
        new_filename_2 = raw_info.loc[ind, 'new_name_R2']
        
        s_id = raw_info.loc[ind, 'sample_name']
        
        file_src_path_1 = os.path.join(src_path, filename_1)
        file_dst_path_1 = os.path.join(dst_path, s_id, new_filename_1)
        
        file_src_path_2 = os.path.join(src_path, filename_2)
        file_dst_path_2 = os.path.join(dst_path, s_id, new_filename_2)
        
        shutil.copyfile(file_src_path_1, file_dst_path_1)
        print("file {} has been sucessfully renamed to {} and moved".format(filename_1, new_filename_1))
        shutil.copyfile(file_src_path_2, file_dst_path_2)
        print("file {} has been sucessfully renamed to {} and moved".format(filename_2, new_filename_2))
        
def batch(input_path, batch_file_path, n_sub):
    """split samples into batches and record the information to a tsv file for each batch
    
    Parameters:
    ----------
    input_path: str
        The root path where the input folders to be created 
    batch_file_path: str
        The path where the tsv files will be stored
    n_sub: int
        The number of samples in each batch
        
    Returns:
        None
    """ 
    path_list = glob.glob(os.path.join(input_path,'*'))
    if len(path_list) == 0:
        raise IOError('No input subfolders')
    sample_list = [x.split('/')[-1] for x in path_list]
    
    df_temp = pd.DataFrame(sample_list, columns = ['name'])
    n_group = len(df_temp.index) // n_sub + 1
    
    for g_id in np.arange(n_group):
        f_name = 's_batch_'+str(g_id)+'.tsv'
        if g_id != (n_group-1):
            temp = test.iloc[g_id*n_sub : (g_id+1)*n_sub]
        else:
            temp = test.iloc[g_id*n_sub : ]
        
        temp.to_csv(os.path.join(batch_file_path, f_name), sep = '\t')
    


# In[7]:


if __name__ == "__main__":
    
    #load infomation metrics
    raw_info_path_1 = '/domino/edv/id-td-virology/Zhiyuan/public/Gut/processed_data/raw_info.csv'
    raw_info_path_2 = '/domino/edv/id-td-virology/Zhiyuan/public/Gut/processed_data/raw_info_2.csv'
    raw_info = pd.read_csv(raw_info_path_1, sep = ',', index_col = 0)
    raw_info_2 = pd.read_csv(raw_info_path_2, sep = ',', index_col = 0)
        
    #add new filenames that fit the cellranger formats
    info_test = filename_reorg(raw_info, raw_info_2)
    
    #make folders for mapping-ready files
    dst_path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/input/'
    make_input_folder(info_test,dst_path)
    
    #path where the source files located
    src_path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/rawdata/raw_fastq/'
    #move files to the folders 
    mv_input_file(info_test, src_path, dst_path)
    
    #split samples into batches
    input_path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/input/'
    batch_file_path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/batch_file/'
    batch(input_path, batch_file_path, n_sub = 8)
    
    
    
    
    
    
    
    


# In[ ]:




