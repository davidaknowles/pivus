#!/bin/bash

INDIR=/srv/gsfs0/projects/montgomery/bballiu/PIVUS/DATA/RNAseq/Filtered_bam/
#INDIR=/srv/gsfs0/projects/montgomery/bballiu/PIVUS/RNAseq/STAR_G-T/

OUTDIR=/srv/gsfs0/projects/montgomery/bballiu/PIVUS/DATA/RNAseq/Filtered_bam/allelic_counts/
#OUTDIR=/srv/gsfs0/projects/montgomery/bballiu/PIVUS/RNAseq/STAR_G-T/allelic_counts/

for f in ${INDIR}*.bam; do
  BASE=`basename ${f%.bam}`
  
  outfile=$OUTDIR/$BASE.txt.gz

  if [ ! -f $outfile ]; then
      echo $f $outfile
      qsub -cwd -v INFILE=$f,OUTFILE=$outfile pileup.sh
  fi
done
