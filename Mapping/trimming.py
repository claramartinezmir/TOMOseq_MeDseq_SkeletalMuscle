#!/usr/bin/env python

###############################################################################
####### Workflow MeDseq data: step 1 - Trimming and LpnPI Filtering ###########
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

### 1 - FASTQ Quality Assurance tool fastQC
os.system('fastqc rawqdata/' + filename + '.fastq.gz -o processed/fastqc')
### 2- TrimGalore: adapter trimming, low quality bases
#remove adapter (universal illumina: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA), bad quality reads and reads smaller than 20 after trimming
os.system('trim_galore rawdata/' + filename + '.fastq.gz -q 20 --phred33 --gzip -o processed/trimmed')
### 3- LpnPI filtering
#count the total number of reads
# In the code below we pipe the output from wc -l (number of lines in the FASTQ file) to awk,
# which executes its body (the statements between the curly braces ( {  } ) for each line of input.
# Here the input is just one line, with one field – the line count. The awk body just divides the 1st input field ($1) by 4 and writes the result to standard output.
output = subprocess.check_output('zcat rawdata/' + filename + ".fastq.gz | wc -l | awk '{print $1 / 4}'", shell=True)
total = int(output)
#We start by opening the fastq file that we are going to create
#assing the variables to use
fout_cpg = open('processed/trimmed/' + filename + '_trimmed_cpg.fq','w')
fout_dcm = open('processed/trimmed/' + filename + '_trimmed_dcm.fq','w')
seq_count = 0
iter_count = 0
sizes = []
sizes_filtered =[]
cpg_count = 0
dcm_count = 0
cpg_positions = []
dcm_positions = []

print('Processing ' + filename + '_trimmed.fq.gz')
#use gzip.open(file.fq.gz) or just open(file.fq)
with gzip.open('processed/trimmed/' + filename + '_trimmed.fq.gz') as f: #if open more than file just add ", gzip.open(fq2) as f2:
    for idx, fl in enumerate(f): #for index and file_line in file (each line will have a number (index) and the content of the line)
        #we have to fix a bit the text to make sense to python
        fl = fl.rstrip().rsplit()[0] #the line has a start of line (b) and new line (\n)
        fl = str(fl, 'utf-8') #transform to string
        #now with mod we ask to divide the index by 4, and it gives the residue (only round numbers, not decimal division)
        l = np.mod(idx,4)
        if l == 0: #if fl = 0, l = 0 --> first line is the read name
            #first we reset the cpg or dcm read to False for the new read
            cpg = False
            dcm = False
            name = fl
            if '@' not in name:
                print (name)
                sys.exit('ERROR: corrupted fastq file %s' % ('processed/trimmed/' + filename + '_trimmed.fq.gz'))
        if l == 1: #fl = 1, l = 1 --> read sequence
            seq = fl
            seq_count = seq_count + 1
            sizes.append(len(seq))
            # reads were filtered based on LpnPI restriction site occurrence between 14 and 19 bp from either 5′ or 3′ end of the read for cpg and between 12 to 20 for dcm
            if "CCG" in seq[13:19]: ## if C*CG --> CpG methylation
                cpg = True
                cpg_count = cpg_count + 1
                sizes_filtered.append(len(seq))
                cpg_positions.append(seq.find('CCG', 13, 19)+1+1) #plus 1 because python positions start at 0, then it will make more sense for the plot, +1 for methylated C being second
            elif "CCG" in seq[-19:-13]: ## if C*CG --> CpG methylation
                cpg = True
                cpg_count = cpg_count + 1
                sizes_filtered.append(len(seq))
                cpg_positions.append(seq.find('CCG', -19, -13)+1+1)
            elif "CGG" in seq[13:19]: ## if *CGG --> CpG methylation
                cpg = True
                cpg_count = cpg_count + 1
                sizes_filtered.append(len(seq))
                cpg_positions.append(seq.find('CGG', 13, 19)+1)
            elif "CGG" in seq[-19:-13]: ## if *CGG --> CpG methylation
                cpg = True
                cpg_count = cpg_count + 1
                sizes_filtered.append(len(seq))
                cpg_positions.append(seq.find('CGG', -19, -13)+1)
            elif "GCGC" in seq[13:19]: ## if G*CGC --> CpG methylation
                cpg = True
                cpg_count = cpg_count + 1
                sizes_filtered.append(len(seq))
                cpg_positions.append(seq.find('GCGC', 13, 19)+1+1)
            elif "GCGC" in seq[-19:-13]: ## if G*CGC --> CpG methylation
                cpg = True
                cpg_count = cpg_count + 1
                sizes_filtered.append(len(seq))
                cpg_positions.append(seq.find('GCGC', -19, -13)+1+1)
            elif "CCAGG" in seq[11:20]: ## if C*C(A/T)GG --> DCM methylation
                dcm = True
                dcm_count = dcm_count + 1
                sizes_filtered.append(len(seq))
                dcm_positions.append(seq.find('CCAGG', 11, 20)+1+1)
            elif "CCAGG" in seq[-20:-11]: ## if C*C(A/T)GG --> DCM methylation
                dcm = True
                dcm_count = dcm_count + 1
                sizes_filtered.append(len(seq))
                dcm_positions.append(seq.find('CCAGG', -20, -11)+1+1)
            elif "CCTGG" in seq[11:20]: ## if C*C(A/T)GG --> DCM methylation
                dcm = True
                dcm_count = dcm_count + 1
                sizes_filtered.append(len(seq))
                dcm_positions.append(seq.find('CCTGG',11, 20)+1+1)
            elif "CCTGG" in seq[-20:-11]: ## if C*C(A/T)GG --> DCM methylation
                dcm = True
                dcm_count = dcm_count + 1
                sizes_filtered.append(len(seq))
                dcm_positions.append(seq.find('CCTGG', -20, -11)+1+1)
            if seq_count == 1000000:
                iter_count = iter_count + 1
                seq_proc = iter_count * seq_count
                print("Progress report: %i" % (seq_proc))
                seq_count = 0
        if l == 2: #fl = 2, l = 2 --> optional sequence description starting with + sign, it's relevant for paired end sequencing
            des = fl[0] #we only save the + sign, not the description
        if l == 3: #fl = 3, l = 3 --> sequence quality
            qua = fl
            if len(qua) != len(seq): #check point that lenght sequence and quality are equal
                sys.exit('ERROR: phred and read length not mathch!')
            #now we can write the file, we add new lines for each line using join,
            if cpg:
                fout_cpg.write('\n'.join([':'.join([name]), seq, des, qua, '']))
            if dcm:
                fout_dcm.write('\n'.join([':'.join([name]), seq, des, qua, '']))


print("Found %i CpG sequences, %i DCM sequences" % (cpg_count, dcm_count))

print("Percentage: %i CpG sequences, %i DCM sequences" % (cpg_count/total*100, dcm_count/total*100))

#gzip subfiles
os.system('gzip processed/trimmed/' + filename + '_trimmed_cpg.fq')
os.system('gzip processed/trimmed/' + filename + '_trimmed_dcm.fq')

# plot the histogram of read length distribution
fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
fig.suptitle('section ' + filename, fontsize=16)

# We can set the number of bins with the `bins` kwarg, for now I don't use bins because it only goes from ~20 to ~50
axs.hist([sizes, sizes_filtered], bins=max(sizes)-min(sizes), color=('blue','green'), edgecolor='white', linewidth=1, align='left', alpha=0.5)
axs.set_title("%i  total reads\n%i filtered reads" % (len(sizes), len(sizes_filtered)))
axs.set_xlabel("Sequence length (bp)")
axs.set_ylabel("Count (log10)")
labels_2 = [str(i) for i in np.arange(25,50,5)]
labels_3 = ['>49']
labels = labels_2 + labels_3
axs.set_xticks(ticks=np.arange(25,55,5))
axs.set_xticklabels(labels=labels)
axs.set_yscale('log')
plt.ylim(ymin=1)
plt.legend(['Total', 'Filtered'])

#save figure to pdf
plt.savefig('processed/trimmed/reports/' + filename + '_read_length_histrograms.pdf')
plt.close()

# plot histogram of position of the motif in the reads
fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
fig.suptitle('section ' + filename, fontsize=16)
axs.hist(cpg_positions, bins=max(cpg_positions)-min(cpg_positions), color='blue', edgecolor='white', linewidth=1, align='left')
axs.set_title("CpG motif distribution")
axs.set_xlabel("Position (bp)")
axs.set_ylabel("Count (log10)")
axs.set_xticks(ticks=np.arange(0,max(cpg_positions),5)) #ticks on positions
axs.set_xticklabels(labels=np.arange(0,max(cpg_positions),5)) #labels have +1 because python uses 0 as position 1
axs.set_yscale('log')
plt.ylim(ymin=1)

#save figure to pdf
plt.savefig('processed/trimmed/reports/' + filename + '_CpG_motif_histrogram.pdf')
plt.close()

# plot histogram of position of the motif in the reads
fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
fig.suptitle('section ' + filename, fontsize=16)
axs.hist(dcm_positions, bins=max(dcm_positions)-min(dcm_positions), color='blue', edgecolor='white', linewidth=1, align='left')
axs.set_title("DCM motif distribution")
axs.set_xlabel("Position (bp)")
axs.set_ylabel("Count (log10)")
axs.set_xticks(ticks=np.arange(0,max(dcm_positions),5))
axs.set_xticklabels(labels=np.arange(0,max(dcm_positions),5))
axs.set_yscale('log')
plt.ylim(ymin=1)

#save figure to pdf
plt.savefig('processed/trimmed/reports/' + filename + '_DCM_motif_histrogram.pdf')
plt.close()
