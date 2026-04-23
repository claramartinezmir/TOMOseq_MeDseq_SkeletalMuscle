#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 4 - Genomewide coverage and binning ########
###############################################################################

#bed files contains chr start end n_lpnpI_sites counts

import matplotlib.pyplot as plt
import numpy as np
import gzip
import shutil
import sys
import os #os is used to instruct pytho to run a bash command. NOTE!! it can not run commands that windows don't recongnize such as ls, cat, gzip, etc.
import subprocess
import pandas as pd

#assign the arguments to variables
bins = sys.argv[1]
filename = sys.argv[2]

path = 'processed/counttables/binned/'



#binned coverage file with bedtools map
#we calculate the total coverage for each bin
os.system('bedtools map -a /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/Lpnp1CpG/mm10_GRCm38/binned/mm10_Lpnp1CpG_index_mC_' + str(bins) + '_sorted.bed -b processed/counttables/genomewide/' + filename + '_cpg_coverage_filtered_sorted.bed -o sum -c 6 > ' + path + filename + '_' + bins + '_cpg_coverage.bed')

#remove uninformative bins from fixed number of bases binning --> NOTE! not all bins have LpnPI-CpG/DCM sites. But also from 10 sites the bin might fall in a blacklisted region and there is no coverage
#these are annotated as '.', we need to eliminate those bins
#create new file
fout = open(path + filename + '_' + bins + '_cpg_coverage_filtered.bed','w')
#open file
with open(path + filename + '_' + bins + '_cpg_coverage.bed') as f:
    for line in f:
        line = line.rsplit()
        if line[4] != '.':
            fout.write('\n'.join(['\t'.join([line[0], line[1], line[2], line[3], line[4]]), '']))

