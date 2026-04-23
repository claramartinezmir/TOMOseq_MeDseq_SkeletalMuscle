#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 3.1 - Genomewide coverage ##################
###############################################################################

import matplotlib.pyplot as plt
import numpy as np
import gzip
import shutil
import sys
import os #os is used to instruct pytho to run a bash command. NOTE!! it can not run commands that windows don't recongnize such as ls, cat, gzip, etc.
import subprocess
import pandas as pd

#assign the arguments to variables
filename=sys.argv[1]

#1-coverage file with bedtools
os.system('bedtools coverage -a /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/Lpnp1CpG/mm10_GRCm38/MM10_Lpnp1CpG_index_mC.bed -b processed/mapped_filtered/' + filename + '_cpg_mapped_filtered.sort.bam -counts > processed/counttables/genomewide/' + filename + '_cpg_coverage.bed')
#os.system('bedtools coverage -a /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/dcm/MM10_dcm_index_mC.bed -b processed/mapped_filtered/' + filename + '_dcm_mapped_filtered.sort.bam -counts > processed/counttables/genomewide/' + filename + '_dcm_coverage.bed')


#2-Filter blacklisted regions
os.system('bedtools intersect -v -a processed/counttables/genomewide/' + filename + '_cpg_coverage.bed -b /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/blacklist/mm10-blacklist.v2.sorted.bed > processed/counttables/genomewide/' + filename + '_cpg_coverage_filtered.bed')
#os.system('bedtools intersect -v -a processed/counttables/genomewide/' + filename + '_dcm_coverage.bed -b /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/blacklist/mm10-blacklist.v2.sorted.bed > processed/counttables/genomewide/' + filename + '_dcm_coverage_filtered.bed')
