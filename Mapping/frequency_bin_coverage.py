import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics
import sys


#sbatch --export=All -t 48:00:00 --mem=100G -e run_log/hist_bin.err -o run_log/hist_bin.out -J hist_bin --wrap="python frequency_bin_coverage.py TM6_1000bp200ovl_filtered TM4_1000bp200ovl M3_1000bp200ovl"

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
        data_frame = pd.read_csv(path+file, sep='\t', index_col=[0,1,2,3],  header=None)
        #add the dataframe as an element of a dictionary with key the filename
        dic_dataframes[list_labels[count]] = data_frame
        count = count + 1

    return dic_dataframes


def saveCSV(dic_data, path, label):
    """Save dataframes inside dictionary to csv files"""
    for key in dic_data:
        dic_data[key].to_csv(path+key+label+'.bed', sep='\t', header=False)

    return print('Files saved')




    #############################################################
    #############################################################
    #############################################################
    #############################################################

path = 'processed/counttables/binned/'
id_files = sys.argv[1:]


#access de data
#Create a variable list with the file names and one with the labels to use as dictionary keys
file_names = []
labels = []
for id in id_files:
    file_names.append(id + '_cpg_coverage_table_genes.bed')
    labels.append(id)

#use accessData() function to obtain a dictonary with each dataset with labels as key
data = accessData(file_names, path, labels)



#plot for each muscle a histograme of coverage frequencies in bins in each section
for key in data:
    count=0
    #assing last 5 columns to index
    for column in data[key].columns[-5:]:
        data[key]=data[key].set_index(column, append=True)
    for section in data[key]:
        count=count+1
        coverage=data[key][section]
        #plot coverage of each bin
        fig, axs = plt.subplots(1, 1, figsize=(15,5))
        fig.suptitle('1kb-bin coverage frequency ' + key + ' section ' + str(count), fontsize=16)
        axs.hist(coverage, bins=np.arange(0,500,10),log=True)
        axs.grid(which='major', axis='y', linestyle='--')
        #axs.set_title('Chromosome ' + chromosome)
        axs.set_xlabel("Coverage(read counts)")
        axs.set_ylabel("Frequency")

        #save figure to pdf
        plt.savefig('processed/counttables/binned/' + key + '_section_' + str(count) + '_cpg_coverage_genes_histogram.pdf')
        plt.close()

#access de data
#Create a variable list with the file names and one with the labels to use as dictionary keys
file_names = []
labels = []
for id in id_files:
    file_names.append(id + '_cpg_coverage_table_TSS3kb.bed')
    labels.append(id)

#use accessData() function to obtain a dictonary with each dataset with labels as key
data = accessData(file_names, path, labels)

#plot for each muscle a histograme of coverage frequencies in bins in each section
for key in data:
    count=0
    #assing last 5 columns to index
    for column in data[key].columns[-5:]:
        data[key]=data[key].set_index(column, append=True)
    for section in data[key]:
        count=count+1
        coverage=data[key][section]
        #plot coverage of each bin
        fig, axs = plt.subplots(1, 1, figsize=(15,5))
        fig.suptitle('1kb-bin coverage frequency ' + key + ' section ' + str(count), fontsize=16)
        axs.hist(coverage, bins=np.arange(0,500,10),log=True)
        axs.grid(which='major', axis='y', linestyle='--')
        #axs.set_title('Chromosome ' + chromosome)
        axs.set_xlabel("Coverage(read counts)")
        axs.set_ylabel("Frequency")

        #save figure to pdf
        plt.savefig('processed/counttables/binned/' + key + '_section_' + str(count) + '_cpg_coverage_TSS3kb_histogram.pdf')
        plt.close()

