#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""
This script do the preparation for the cellranger aggr, which integrate all samples into one big matrix.
The preparation includes aggr csv file generation which is required to run cellranger aggr, 
correction of the sample name due to the midleading information from the original metadata, collecting all 
the web_summary files from each sample etc.

After running this script, submit aggr.sh file to domino jobs

author:Zhiyuan Yao
Date:Aug-10-2022
"""

import numpy as np
import pandas as pd
from copy import copy
import os
import glob
import shutil
import subprocess


# In[2]:


def correct_metadata(temp1,temp2):
    """generate the dataframe stored original sample name and the corrected sample name, this is special for the Gut paper
    the Gut paper metadata is midleading resulting in the original fastq unmatched the sample names, this function
    together with c_meta_gen function aims to address the problem.
    
    Parameters:
    ----------
    temp1: dataframe
        The table extracted from the original metadata, contains fastq accession, 
        Run title and Experiment accession. 
    temp2: dataframe
        The table extracted from the original metadata, contains Experiment accession, 
        Experiment title and BioSample name.

    Returns:
        df: dataframe
        The table contains the Experimental accession, the original corresponding sample name and 
        the corrected sample name
    """     
    df = pd.DataFrame(index = temp2.index, columns = ['origin_sample_name', 'corrected_sample_name'])
    for ind in df.index:
        temp_run_title = temp1.loc[temp1['Experiment accession'] == ind, 'Run title'][0]


        correct_name = temp2.loc[temp2['Experiment title'] == temp_run_title, 'BioSample name'][0]
        wrong_name = temp2.loc[ind, 'BioSample name']
        df.loc[ind, 'origin_sample_name'] = wrong_name
        df.loc[ind, 'corrected_sample_name'] = correct_name
        
    return(df)

def c_meta_gen(raw_info, raw_info_2, savepath):
    """generate and save the dataframe stored original sample name and the corrected sample name.
    Detailed description please see function correct_metadata
    
    Parameters:
    ----------
    raw_info: dataframe
        The pandas dataframe that stores the fastq original file name, Accession number, corresponding experimental
        accession number and Run title etc.         
    raw_info_2: dataframe
        The pandas dataframe that stores experimental accession number and its corresponding Sample name and Run title
        raw_info and raw_info_2 are extracted from the excel table file downloaded from the Genome Sequence Archive 
        https://ngdc.cncb.ac.cn, accession number HRA001730
    savepath: str
        path to save the corrected table
        
    Returns:
        none
    """    
    temp_1 = raw_info.loc[:, ['Run title', 'Experiment accession']].copy()
    temp_1['Run title'] = temp_1['Run title'].str.split(':', expand = True)[0]
    temp_1  = temp_1.drop_duplicates(subset = 'Run title')
    
    temp_2 = raw_info_2.loc[:,['Experiment title', 'BioSample name']].copy()
    
    correct_meta = correct_metadata(temp1,temp2)
    correct_meta.to_csv(savepath, sep = '\t')

def get_output_sample(count_output_path):
    """Generate a table containing the sample_id and the path to its coressponding molecule_info.h5 file.
    
    Parameters:
    ----------
    count_output_path: str
        The root path where the output folders generated during the mapping process located

    Returns:
        df: dataframe
        The uncorrected csv contains the sample name and the path to the molecule_info.h5 file
    """     
    path_list = glob.glob(os.path.join(count_output_path,'*'))
    if len(path_list) == 0:
            raise IOError('No count output subfolders')
    sample_list = [x.split('/')[-1] for x in path_list]
    dic = {}
    dic['sample_id'] = sample_list
    dic['molecule_h5'] = path_list
    df = pd.DataFrame(dic)
    df['molecule_h5'] = df['molecule_h5'] + '/outs/molecule_info.h5'
    df.set_index('sample_id', inplace = True)
    return(df)
    
def correct_csv(df, c_meta, write_path):
    """Correct the table generated from function get_output_sample with the corrected metadata from function 
    c_meta_gen
    
    Parameters:
    ----------
    df: dataframe
        The output of the function get_output_sample
    c_meta: dataframe
        The corrected metatable generated from the funtion c_meta_gen
    write_path:str
        The path where the corrected csv will be stored
    Returns:
        none
    """     
    df_temp = df.copy()
    df_temp['sample_id'] = ''
    for ind in df_temp.index:
        df_temp.loc[ind, 'sample_id'] = c_meta.loc[c_meta.origin_sample_name == ind,
                                                 'corrected_sample_name'][0]
    df_temp = df_temp.set_index('sample_id')
    
    #write
    df_temp.to_csv(write_path)

def move_summary_files(src_path, dst_path, c_meta):
    """collect all the web_summary files together with the corrected sample names
    
    Parameters:
    ----------
    src_path: str
        root path where the output folders from mapping located
    dst_path: str
        path where all web_summary files will be moved to
    c_meta: dataframe
        The corrected metatable generated from the funtion c_meta_gen

    Returns:
        none
    """     
    for sample in list(c_meta.origin_sample_name):
        file_src_path = os.path.join(src_path,sample,'outs/web_summary.html')
        file_new_name = c_meta.loc[c_meta.origin_sample_name == sample, 
                                   'corrected_sample_name'][0]
        
        file_dst_path = os.path.join(dst_path,'{}.html'.format(file_new_name))
        
        shutil.copyfile(file_src_path, file_dst_path)
        print("web_summary for original name {} has been corrected to {} and copied to the new folder".format(sample, file_new_name))

# def cellranger_aggr(output_path, csv_path, name):
#     os.chdir(output_path)
#     cellranger_command = "cellranger aggr \
# --id={} \
# --csv={} \
# --normalize=none \
# --nosecondary".format(name, csv_path)
    
#     print("start aggregation...")
#     print(cellranger_command)
    
#     a = subprocess.run(cellranger_command, capture_output=True, check = True, text = True, shell = True)        
#     if a.returncode == 0:
#         print('aggregation completed'.format(sample))
            
#     else:
#         print("aggregation failed".format(sample))
    


# In[3]:


if __name__ == "__main__":
    # path for cellranger count outputs
    count_output_path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/output/'
    # path to write the aggr csv file
    write_path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/aggr/aggr_csv.csv'

    #generate the corrected metatable
    raw_info_path_1 = '/domino/edv/id-td-virology/Zhiyuan/public/Gut/processed_data/raw_info.csv'
    raw_info_path_2 = '/domino/edv/id-td-virology/Zhiyuan/public/Gut/processed_data/raw_info_2.csv'
    c_meta_path = '/domino/edv/id-td-virology/Zhiyuan/public/Gut/processed_data/corrected_sample_name.csv'
    raw_info = pd.read_csv(raw_info_path_1, sep = ',', index_col = 0)
    raw_info_2 = pd.read_csv(raw_info_path_2, sep = ',', index_col = 0)
    c_meta_gen(raw_info, raw_info_2, c_meta_path)

    # path to load corrected metatable
    c_meta_path = '/domino/edv/id-td-virology/Zhiyuan/public/Gut/processed_data/corrected_sample_name.csv'
    c_meta = pd.read_csv(c_meta_path, sep = '\t', index_col = 0)

    # generate the csv that necesseary for cellranger aggr
    csv = get_output_sample(count_output_path)
    correct_csv(csv, c_meta, write_path)

    #rename and move summary files together
    dst_path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/aggr/web_summary/'
    move_summary_files(count_output_path, dst_path, c_meta)

