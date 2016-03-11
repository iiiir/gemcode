#!/bin/bash

echo "*** Sorting BAM by read name, output to SAME folder ***"

[[ $# -lt 1 ]] && echo "Usage: $0 <bam>" && exit 1

f=`cd \`dirname $1\`; pwd`/`basename $1`
o=${f/.bam/.qn}
osam=$o.sam

echo ">>> Sorting on BAM $f"
samtools sort -n $f $o
samtools view -h $o.bam > $osam
rm $o.bam

echo ">>> Process BAM"
sam_clean_unpaired.py $osam | samtools view -bS - > $o.pe.bam


echo "*** Finished sorting BAM by position ***"
