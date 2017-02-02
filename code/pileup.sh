#!/bin/sh
#$ -S /bin/bash

# Print the parameters
echo Input file: $INFILE
echo Output file: $OUTFILE

cd /home/dak33/pivus/

# Stop on error
set -e

module load samtools

# INFILE=/srv/gsfs0/projects/montgomery/bballiu/PIVUS/RNAseq/Filtered_bam/PIVUS03.Aligned.out_mapq255_sorted_dedup.bam
# OUTFILE=test.gz

# mpileup: -B prevents error calculation, 
# -d1000000 sets the maximum depth to record
# -l specifies the regions (SNPs) to quantify
# pileup_filter 2 2 is a little prog I wrote to process the mpileup output (should prob use vcf/bcf output option from mpileup instead)
samtools mpileup -B -d1000000 -l hg19_snps.txt -f /srv/gsfs0/shared_data/RefGenomes/H_sapiens/hg19/hg19.fa $INFILE | ./pileup_filter 2 2 | gzip > $OUTFILE
