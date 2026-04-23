#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 5 - Count table ############################
###############################################################################

#sbatch --export=All -t 24:00:00 --mem=100G  -e run_log/1-ct.err -o run_log/1-ct.out -J 1-ct --wrap="python make_count-table_allMuscles.py"

#combine all bed files into one file with chr start end n_lpnpI_sites counts_s1 counts_s2 .... counts_sx


import matplotlib.pyplot as plt
import numpy as np
import gzip
import shutil
import sys
import os #os is used to instruct pytho to run a bash command. NOTE!! it can not run commands that windows don't recongnize such as ls, cat, gzip, etc.
import subprocess
import pandas as pd

#build functions for opening and saving files
def accessData(list_files, path, list_labels):
    """opens each file and stores in dictionary with key a label"""
    #create empty dictionary to store the dataframes
    dic_dataframes = {}
    #iterate to get the file names
    #make a counter to access the correct label
    count = 0
    for file in list_files:
        #open each file
        data_frame = pd.read_csv(path+file, sep='\t', index_col=[0,1,2,3,4],  names=[list_labels[count]])
        #add the dataframe as an element of a dictionary with key the filename
        dic_dataframes[list_labels[count]] = data_frame
        count = count + 1
    return dic_dataframes

def saveCSV(dic_data, path, label):
    """Save dataframes inside dictionary to csv files"""
    for key in dic_data:
        dic_data[key].to_csv(path+key+'_'+label+'.bed', sep='\t', header=True)
    return print('Files saved')

###############################################################################
###############################################################################

#assign the arguments to variables
id_files=['TM6-13', 'TM6-14', 'TM6-15', 'TM6-16', 'TM6-27', 'TM6-28', 'TM6-29', 'TM6-30', 'TM6-38', 'TM6-39', 'TM6-40', 'TM6-41']
files_cpg=[]
labels_cpg=[]
path = 'processed/counttables/genomewide/'

#access each file
for filename in id_files:
    #combine all files into one counttable
    files_cpg.append(filename + '_cpg_coverage_filtered_sorted.bed')
    labels_cpg.append(filename)

#use the formulas to combine the datasets
#access data
data_cpg = accessData(files_cpg, path, labels_cpg)

#Construct the count CpG and DCM table that cotains each bin (chr, start, end, number of lpnp1CpG sites) and as columns each section
data_all={}
count=0
for label in data_cpg:
    if count == 0:
        data_all['Muscle_1']=data_cpg[label]
        count=1
    else:
        data_all['Muscle_1']= pd.concat([data_all['Muscle_1'], data_cpg[label]], axis=1).dropna()


saveCSV(data_all, path, label='cpg_coverage_table')

