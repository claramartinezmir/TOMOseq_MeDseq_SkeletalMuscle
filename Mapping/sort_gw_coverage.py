#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 3.2 - Sort genomewide coverage  ############
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

#first we sort the coverage files (if I do this step on the methylation fraction script then the sort starts before the unsorted file is finished, then it's smaller)
os.system('sort -k 1,1 -k2,2n processed/counttables/genomewide/' + filename + '_cpg_coverage_filtered.bed > processed/counttables/genomewide/' + filename + '_cpg_coverage_filtered_sorted.bed')
os.system('rm processed/counttables/genomewide/' + filename + '_cpg_coverage_filtered.bed')
