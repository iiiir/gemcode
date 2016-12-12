sj pipeline for analyzing genmcode/ chromium data

#STEPS TO RUN GEMCODE PIPELINE

##step 1: fastq to ubam
fqtoUbam_batch_cmd.sh /path/to/demux/ SAMPLE_ID SAMPLE_SM SAMPLE_LB > fqtoUbam_batch.myjobs

Required scripts:
* fqtoUbam_batch_cmd.sh
* fqtoUbam_gemcode.py

runtime: 
5 hours

##step 2: bwa-aln
run_bwa_aln_gemcode.py -b /path/to/ubam/*.bam -o b37.bam -O `pwd` --tmp /scratch_space/ -j b37.sjm    

sjm b37.sjm

Required scripts:
* run_bwa_aln_gemcode.py
* picard_sortUbam.sh
* bwa_aln_pe_qn.sh
* picard_mergeBam.sh
* barbam.sh
* barbam.py
* picard_msam.sh
* picard_mdup.sh
* samtools (flagstat)

runtime:

