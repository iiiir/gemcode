#!/bin/bash

# fqtoUbam_batch_cmd.sh ../demux/ SJCOLO829_D1 SJCOLO829 SJPILOT > fqtoUbam.myjobs
# echo "Assumed file naming format: read-I1_si-TGTGCGGG_lane-001-chunk-001.fastq.gz"

[[ $# -lt 4 ]] && echo "$0 <fq_path> <id> <sm> <lb>" && exit 0
[[ ! -d $1 ]] && echo "Not a valide path: $1" && exit 0

fq_path=`cd $1;pwd`

SAMPLE_ID=$2
SAMPLE_SM=$3
SAMPLE_LB=$4

#echo $fq_path

for RAFILE in $fq_path/read-RA_si-*.fastq.gz; do
	I1FILE=`echo $RAFILE | sed s/read-RA_si-/read-I1_si-/`
	#echo $RAFILE
	#echo $I1FILE
	PRENAME=`basename $RAFILE | sed s/read-RA_si-// | sed s/\.fastq\.gz// | sed s/-/_/g`
	[[ ! -f $RAFILE ]] && echo "Not a valide filepath: $RAFILE" && exit 0
	[[ ! -f $I1FILE ]] && echo "Not a valide filepath: $I1FILE" && exit 0
	cmd="fqtoUbam_gemcode.py --fastq1 $RAFILE --fastq2 $I1FILE \
			--lb $SAMPLE_LB --sm $SAMPLE_SM --id $SAMPLE_ID --outprefix $PRENAME"
	echo $cmd
done
