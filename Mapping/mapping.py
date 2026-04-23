#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 2 - Mapping ################################
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

### 1- Mapping using bwa
#step 1
#os.system('bwa aln /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/genome/ensemble-102/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz processed/trimmed/' + filename + '_trimmed_cpg.fq.gz > processed/mapped/' + filename + '_cpg.sai')
#os.system('~/bwa/bwa aln referencegenome/mm10.fa.gz processed/trimmed/' + filename + '_trimmed_dcm.fq.gz > processed/mapped/' + filename + '_dcm.sai')
#step 2
#os.system('bwa samse /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/genome/ensemble-102/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz processed/mapped/' + filename + '_cpg.sai processed/trimmed/' + filename + '_trimmed_cpg.fq.gz > processed/mapped/' + filename + '_cpg.sam')
#os.system('~/bwa/bwa samse referencegenome/mm10.fa.gz processed/mapped/' + filename + '_dcm.sai processed/trimmed/' + filename + '_trimmed_dcm.fq.gz > processed/mapped/' + filename + '_dcm.sam')
#convert sam to bam
#os.system('samtools view -S -b processed/mapped/' + filename + '_cpg.sam > processed/mapped/' + filename + '_cpg.bam')
#os.system('samtools view -S -b processed/mapped/' + filename + '_dcm.sam > processed/mapped/' + filename + '_dcm.bam')
#report
# Since the BAM file contains records for both mapped and unmapped reads, just counting records doesn't provide
# information about the mapping rate of our alignment. The samtools flagstat tool provides a simple analysis of
# mapping rate based on the the SAM flag fields.
# Here's how to run samtools flagstat and both see the output in the terminal and save it in a file – the samtools
# flagstat standard output is piped to tee, which both writes it to the specified file and sends it to its standard output:
#os.system('samtools flagstat processed/mapped/' + filename + '_cpg.bam | tee processed/mapped/' + filename + '_cpg.flagstat.txt')
#os.system('samtools flagstat processed/mapped/' + filename + '_dcm.bam | tee processed/mapped/' + filename + '_dcm.flagstat.txt')

### 2- Filtering with SAMTools
# The most common samtools view filtering options are:
    #-q N – only report alignment records with mapping quality of at least N (>= N). Threshold at 60 for removing double mapping and bad quality

# 1- Filter uniquely  mapped reads. NOTE!! we need the header in the output (-h)
#os.system('samtools view -h -q 30 processed/mapped/' + filename + '_cpg.bam > processed/mapped_filtered/' + filename + '_cpg_mapped_filtered.bam')
#os.system('samtools view -h -q 30 processed/mapped/' + filename + '_dcm.bam > processed/mapped_filtered/' + filename + '_dcm_mapped_filtered.bam')

# 2- sort files --> samtools sort re-orders entries in the SAM file either by locus (contig name + coordinate position) or by read name.
#os.system('samtools sort -O bam -T processed/mapped/' + filename + '_cpg_mapped_filtered.tmp processed/mapped_filtered/' + filename + '_cpg_mapped_filtered.bam  > processed/mapped_filtered/' + filename + '_cpg_mapped_filtered.sort.bam')
#os.system('samtools sort -O bam -T processed/mapped/' + filename + '_dcm_mapped_filtered.tmp processed/mapped_filtered/' + filename + '_dcm_mapped_filtered.bam  > processed/mapped_filtered/' + filename + '_dcm_mapped_filtered.sort.bam')

# 3- index files
#os.system('samtools index processed/mapped_filtered/' + filename + '_cpg_mapped_filtered.sort.bam')
#os.system('samtools index processed/mapped_filtered/' + filename + '_dcm_mapped_filtered.sort.bam')

#Report
# More information about the alignment is provided by the samtools idxstats report, which shows how many reads aligned to
# each contig in your reference. Note that samtools idxstats must be run on a sorted, indexed BAM file.
    #The output has four tab-delimited columns:
      #contig name
      #contig length
      #number of mapped reads
      #number of unmapped reads
#os.system('samtools idxstats processed/mapped_filtered/' + filename + '_cpg_mapped_filtered.sort.bam | tee processed/mapped_filtered/' + filename + '_cpg.idxstats.txt')
#os.system('samtools idxstats processed/mapped_filtered/' + filename + '_dcm_mapped_filtered.sort.bam | tee processed/mapped_filtered/' + filename + '_dcm.idxstats.txt')

### 3- Deeptools --> Create heatmap of coverage at TSS and TTS with deepTools
#make the bigwig -->

#NOTE!!! --binSize 100 for next time!!!!!

#os.system('bamCoverage -b processed/mapped_filtered/' + filename + '_dcm_mapped_filtered.sort.bam -o processed/counttables/deeptools/' + filename + '_dcm_coverage_100bp.bw -bs 100 -bl /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/blacklist/mm10-blacklist.v2.sorted.bed')
#os.system('bamCoverage -b processed/mapped_filtered/' + filename + '_cpg_mapped_filtered.sort.bam -o processed/counttables/deeptools/' + filename + '_cpg_coverage_100bp.bw -bs 100')

#os.system('computeMatrix scale-regions -S processed/counttables/deeptools/' + filename + '_dcm_coverage_100bp.bw -R /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/genes/ensembl/Mus_musculus.GRCm38.102.gtf.gz -o processed/counttables/deeptools/' + filename + '_dcm.matrix_100bp.gz  -b 3000 -a 3000 --regionBodyLength 5000 --skipZeros')
os.system('computeMatrix scale-regions -S processed/counttables/deeptools/' + filename + '_cpg_coverage_100bp.bw -R /exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/references/referencebed/genes/ensembl/Mus_musculus.GRCm38.102.gtf.gz -o processed/counttables/deeptools/' + filename + '_cpg.matrix_100bp.gz  -b 3000 -a 3000 --regionBodyLength 5000 --skipZeros')

#plot the heatmap
#os.system('plotHeatmap -m processed/counttables/deeptools/' + filename + '_dcm.matrix_100bp.gz -out processed/counttables/deeptools/' + filename + '_dcm_heatmap_100bp.pdf')
os.system('plotHeatmap -m processed/counttables/deeptools/' + filename + '_cpg.matrix_100bp.gz -out processed/counttables/deeptools/' + filename + '_cpg_heatmap_100bp.pdf')

#IGV browser with the bw files
#remove bam and sai files after mapping
#os.system('rm processed/mapped/*.sa*')
#remove unsorted bam files
#os.system('rm processed/mapped_filtered/*_mapped_filtered.bam')

