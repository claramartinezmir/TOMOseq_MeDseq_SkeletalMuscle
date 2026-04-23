#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 3.5 - Power log plots  #####################
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

#cpg
coverage = []

threshold=dict({'M3-4': 200, 'M3-5': 500, 'M3-6': 200, 'M3-7': 300, 'M3-21': 200, 'M3-22': 200, 'M3-23': 200, 'M3-24': 200, 'M3-37': 200, 'M3-38': 200, 'M3-39': 200, 'M3-40': 200,
        'TM4-21': 600, 'TM4-22': 300, 'TM4-23': 400, 'TM4-24': 500, 'TM4-35': 400, 'TM4-36': 300, 'TM4-37': 700, 'TM4-38': 300, 'TM4-49': 600, 'TM4-50': 250, 'TM4-51': 250, 'TM4-52': 500,
        'TM6-13': 700, 'TM6-14': 600, 'TM6-15': 650, 'TM6-16': 700, 'TM6-27': 500, 'TM6-28': 1000, 'TM6-29': 250, 'TM6-30': 1000, 'TM6-38': 800, 'TM6-39': 800, 'TM6-40': 500, 'TM6-41': 700
        })

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
axs.axvline(x=threshold[filename], ls='-.', c='r')
axs.set_xlabel("Coverage")
axs.set_ylabel("Frequency")

#save figure to pdf
plt.savefig('processed/counttables/genomewide/figures/' + filename + '_cpg_coverage_powerlog_thr.pdf')
plt.close()

