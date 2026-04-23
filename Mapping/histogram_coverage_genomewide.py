#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 3.3 - Histogram coverage #####################
###############################################################################

import matplotlib.pyplot as plt
import numpy as np
import sys
import os #os is used to instruct pytho to run a bash command. NOTE!! it can not run commands that windows don't recongnize such as ls, cat, gzip, etc.


#assign the arguments to variables
filename=sys.argv[1]

'''
#dcm
coverage = []

#open file
with open('processed/counttables/genomewide/' + filename + '_dcm_coverage_filtered.bed') as f: #if open more than file just add ", gzip.open(fq2) as f2:
    for line in f: #for index and file_line in file (each line will have a number (index) and the content of the line)
        line = line.rsplit()
        coverage.append(int(line[5]))

#plot coverage on chromosome position
fig, axs = plt.subplots(1, 1, figsize=(15,5))
fig.suptitle('Histogram filtered DCM coverage section ' + filename, fontsize=16)
axs.hist(coverage, bins=np.arange(0,150,1),log=True)
axs.grid(which='major', axis='y', linestyle='--')
#axs.set_title('Chromosome ' + chromosome)
axs.set_xlabel("Coverage")
axs.set_ylabel("Frequency")

#save figure to pdf
plt.savefig('processed/counttables/genomewide/' + filename + '_dcm_coverage_filtered_histogram.pdf')
plt.close()
'''

#cpg
coverage = []

#open file
with open('processed/counttables/genomewide/' + filename + '_cpg_coverage.bed') as f: #if open more than file just add ", gzip.open(fq2) as f2:
    for line in f: #for index and file_line in file (each line will have a number (index) and the content of the line)
        line = line.rsplit()
        coverage.append(int(line[5]))

#plot coverage on chromosome position
fig, axs = plt.subplots(1, 1, figsize=(15,5))
fig.suptitle('Filtered CpG coverage section ' + filename, fontsize=16)
axs.hist(coverage, bins=np.arange(0,200,1),log=True)
axs.grid(which='major', axis='y', linestyle='--')
#axs.set_title('Chromosome ' + chromosome)
axs.set_xlabel("Coverage")
axs.set_ylabel("Frequency")

#save figure to pdf
plt.savefig('processed/counttables/genomewide/figures/' + filename + '_cpg_coverage_histogram.pdf')
plt.close()
