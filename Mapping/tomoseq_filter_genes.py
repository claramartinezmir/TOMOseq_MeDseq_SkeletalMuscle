#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 5 - TOMOseq Processing #####################
###############################################################################

#sbatch --export=All -t 48:00:00 --mem=100G -e run_log/filter_genes.err -o run_log/filter_genes.out -J filter_genes --wrap="python tomoseq_filter_genes.py TM6_1000bp200ovl TM4_1000bp200ovl M3_1000bp200ovl"


import numpy as np
import gzip
import sys
import os #os is used to instruct pytho to run a bash command. NOTE!! it can not run commands that windows don't recongnize such as ls, cat, gzip, etc.
import subprocess
import pandas as pd

#assign the arguments to variables
path = 'processed/counttables/binned/'
thr=int(sys.argv[1])
id_files = sys.argv[2:]

#access de data
#Create a variable list with the file names and one with the labels to use as dictionary keys
file_names = []
labels = []

for key in id_files:
    seq_count = 0
    iter_count = 0
    print('Filtering non methylated data ' + key)

    #create new file
    fout = open(path + key + '_cpg_coverage_table_genes_TOMOseqfilter'+ str(thr)+'.bed','w')
    #open file
    with open(path + key + '_cpg_coverage_table_genes.bed') as f:
        for idx, line in enumerate(f):
            line = line.rsplit()
            index1 = line[0]
            index2 = line[1]
            index3 = line[2]
            index4 = line[3]
            #only bins with more than 4 counts in at least 3 sections from the same group
            if sum([float(i)>thr for i in line[4:8]])>=3:
                fout.write('\n'.join(['\t'.join([index1, index2, index3, index4, line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20]]), '']))
            elif sum([float(i)>thr for i in line[8:12]])>=3:
                fout.write('\n'.join(['\t'.join([index1, index2, index3, index4, line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20]]), '']))
            elif sum([float(i)>thr for i in line[12:16]])>=3:
                fout.write('\n'.join(['\t'.join([index1, index2, index3, index4, line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18], line[19], line[20]]), '']))

