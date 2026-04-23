#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 5 - TSS3kb and gene body annotation ########
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
bins = sys.argv[1]
lib = sys.argv[2]


path = 'processed/counttables/binned/DMR/'

#annotation of bins with bedtools and references
os.system('bedtools intersect -a ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_p0.05_MvsF_v2.bed -b /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/genes/ensemble-102/Mus_musculus.GRCm38.102_genes.bed -wb -wa > ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_p0.05_MvsF_genes_v2.bed')

os.system('bedtools intersect -a ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_p0.05_MvsF_v2.bed -b /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/genes/ensemble-102/Mus_musculus.GRCm38.102_TSS3kb.bed -wb -wa > ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_p0.05_MvsF_TSS3kb_v2.bed')

os.system('bedtools intersect -a ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_p0.05_MvsF_v2.bed -b /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/genes/ensemble-102/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff -wb -wa > ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_p0.05_MvsF_regulatory_v2.bed')

os.system('bedtools intersect -a ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_pval_MvsF_v2.bed -b /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/genes/ensemble-102/Mus_musculus.GRCm38.102_genes.bed -wb -wa > ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_pval_MvsF_genes_v2.bed')

os.system('bedtools intersect -a ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_pval_MvsF_v2.bed -b /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/genes/ensemble-102/Mus_musculus.GRCm38.102_TSS3kb.bed -wb -wa > ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_pval_MvsF_TSS3kb_v2.bed')

os.system('bedtools intersect -a ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_pval_MvsF_v2.bed -b /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/genes/ensemble-102/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff -wb -wa > ' + path + lib + '_' + bins + '_cpg_coverage_table_filtered_normal_DMR-ttest_pval_MvsF_regulatory_v2.bed')

