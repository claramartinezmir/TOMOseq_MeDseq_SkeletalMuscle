import numpy as np
import pandas as pd
import sys
from scipy import stats
from scipy.stats import chi2_contingency
from scipy.stats import ttest_ind
from itertools import combinations
from statsmodels.sandbox.stats.multicomp import multipletests
import multiprocessing

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


def findDMRsTtest(inp):
    data_meth_g1, data_meth_g2, index = inp
    g1 = pd.DataFrame(np.array([i for i in data_meth_g1]))
    g2 = pd.DataFrame(np.array([i for i in data_meth_g2]))
    #t-test
    #if all are zeros give p=1:
    if float(np.mean(g1))==0 and float(np.mean(g2))==0:
        p_val=1
        log2FC=0
    else:
        statistic, p = ttest_ind(g1, g2, equal_var=False)
        p_val=float(p)
        log2FC=float(np.log2((np.mean(g1)+0.1)/(np.mean(g2)+0.1)))
    return index, p_val, log2FC

#############################################
#############################################

#access de data
#Create a variable list with the file names and one with the labels to use as dictionary keys
id=sys.argv[1]
#id='M1M2M3_1000bp0ovl_cpg_coverage_table_split_10'
path = 'processed/counttables/binned/split/'
file_names = [id + '.bed']
labels = [id]
outpath = 'processed/counttables/binned/DMR/split/'
#use accessData() function to obtain a dictonary with each dataset with labels as key
data_normal = accessData(file_names, path, labels)
# Convert the any columns from string to float
data_normal[id] = data_normal[id].astype(float)

list_M_sections=['M3_s17', 'M3_s18','M3_s19', 'M3_s20', 'M3_s1', 'M3_s2', 'M3_s3','M3_s4','M3_s33','M3_s34', 'M3_s35', 'M3_s36']
list_F_sections=['M1_s10', 'M1_s11', 'M1_s12', 'M1_s13', 'M2_s8', 'M2_s9', 'M2_s10', 'M2_s11', 'M2_s22', 'M2_s23', 'M2_s24', 'M2_s25', 'M1_s35','M1_s36', 'M1_s38','M2_s36','M2_s37', 'M2_s38', 'M2_s39']


#we divide the table by groups
data_normal_M=data_normal[id][list_M_sections]
data_normal_F=data_normal[id][list_F_sections]

#create df to store pvalues
pvalues=pd.DataFrame(columns=['Log2FC_MvsF', 'p_val_MvsF', 'adj_p_val_MvsF'], index= data_normal_M.index)
Log2FC_MvsF={}
p_val_MvsF={}
pool = multiprocessing.Pool(8) #threads
for label, p, fc in pool.imap_unordered(findDMRsTtest, [(data_normal_M.loc[index], data_normal_F.loc[index], index) for index in data_normal_M.index]):
    #assing p_vales to dataframe by index
    Log2FC_MvsF[label] = fc
    p_val_MvsF[label] = p

#add the values to the dataframe
df_Log2FC_MvsF=pd.DataFrame.from_dict(Log2FC_MvsF, orient='index')
pvalues['Log2FC_MvsF']=df_Log2FC_MvsF
df_p_val_MvsF=pd.DataFrame.from_dict(p_val_MvsF, orient='index')
pvalues['p_val_MvsF']=df_p_val_MvsF

#Correct the p-values for multiple testing
pvalues['adj_p_val_MvsF']=multipletests(pvalues['p_val_MvsF'], method='fdr_bh')[1]

#search for bins with pvalues significant
list_index=[]
for index in pvalues.index:
    if pvalues.loc[index]['adj_p_val_MvsF']>=0.05:
        list_index.append(index)

data_DMR=data_normal[id].drop(list_index)

#save corrected p values
dic_pvalues_corr={}
dic_pvalues_corr[id]=pvalues
saveCSV(dic_pvalues_corr, outpath, label = '_DMR-ttest_pval_MvsF_v2')

#save DMRs
data_dmr={}
data_dmr[id]=data_DMR
saveCSV(data_dmr, outpath, label = '_DMR-ttest_p0.05_MvsF_v2')
