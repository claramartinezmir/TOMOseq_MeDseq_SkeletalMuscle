#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 3.4 - Power log plots  #####################
###############################################################################

import matplotlib.pyplot as plt
import numpy as np
import gzip
import shutil
import sys
import os #os is used to instruct pytho to run a bash command. NOTE!! it can not run commands that windows don't recongnize such as ls, cat, gzip, etc.
import subprocess

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
axs.hist(coverage, bins=10**(np.arange(0,4,0.1)),log=True)
axs.set_xscale('log')
axs.grid(which='major', axis='y', linestyle='--')
#axs.set_title('Chromosome ' + chromosome)
axs.set_xlabel("Coverage")
axs.set_ylabel("Frequency")

#save figure to pdf
plt.savefig('processed/counttables/genomewide/' + filename + '_dcm_coverage_filtered_powerlog.pdf')
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
axs.hist(coverage, bins=10**(np.arange(0,4,0.1)),log=True)
axs.set_xscale('log')
axs.grid(which='major', axis='y', linestyle='--')
#axs.set_title('Chromosome ' + chromosome)
axs.set_xlabel("Coverage")
axs.set_ylabel("Frequency")

#save figure to pdf
plt.savefig('processed/counttables/genomewide/figures/' + filename + '_cpg_coverage_powerlog.pdf')
plt.close()
