import numpy as np
import pandas as pd
import sys


def accessData(list_files, path, list_labels):
    """opens each file and stores in dictionary with key a label"""
    #create empty dictionary to store the dataframes
    dic_dataframes = {}
    #iterate to get the file names
    #make a counter to access the correct label
    count = 0
    for file in list_files:
        #open each file
        data_frame = pd.read_csv(path+file, sep='\t', index_col=[0,1,2,3])
        #add the dataframe as an element of a dictionary with key the filename
        dic_dataframes[list_labels[count]] = data_frame
        count = count + 1
    return dic_dataframes

def saveCSV(dic_data, path, label):
    """Save dataframes inside dictionary to csv files"""
    for key in dic_data:
        dic_data[key].to_csv(path+key+label+'.bed', sep='\t', header=True)
    return print('Files saved')


#############################################
#############################################

#access de data
#Create a variable list with the file names and one with the labels to use as dictionary keys
path = 'processed/counttables/binned/DMR/split/'
id=sys.argv[1]
split_id=sys.argv[2:]
new_label = [id + '_DMR-ttest_pval_MvsF_v2',
             id + '_DMR-ttest_p0.05_MvsF_v2']
labels={}
labels[new_label[0]]  = []
labels[new_label[1]]  = []
file_names={}
file_names[new_label[0]]  = []
file_names[new_label[1]]  = []

for split in split_id:
    labels[new_label[0]].append(split + '_DMR-ttest_pval_MvsF_v2')
    labels[new_label[1]].append(split + '_DMR-ttest_p0.05_MvsF_v2')

    file_names[new_label[0]].append(split + '_DMR-ttest_pval_MvsF_v2.bed')
    file_names[new_label[1]].append(split + '_DMR-ttest_p0.05_MvsF_v2.bed')


outpath = 'processed/counttables/binned/DMR/'

for key in new_label:
    #use accessData() function to obtain a dictonary with each dataset with labels as key
    split_data = accessData(file_names[key], path, labels[key])
    #re-join dataframes into one
    rejoin_data = pd.concat([v for k,v in split_data.items()])

    #save re-join data to dic
    dic_data={}
    dic_data[key]=rejoin_data

    #save the re-join dic
    saveCSV(dic_data, outpath, label='')

