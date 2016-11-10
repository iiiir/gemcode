#!/bin/bash -eu

>&2 echo "*** Merge multiple BAMs ***"

if [ $# -lt 3 ]
then
    >&2 echo "Usage: $0 <out.bam> <1.bam> <2.bam> [3.bam]"
    exit 1
fi

o=`cd \`dirname $1\`; pwd`/`basename $1`
shift

BAMs=""
for bam in $@; do
	BAMs="$BAMs I=$bam"
done

# by default, tags starts with XYZ must be specified
# all bwa tags were kept
cmd="java -Xms10g -Xmx10g -jar /hpcf/apps/picard/install/2.0.1/picard.jar MergeSamFiles \
	O=$o \
	$BAMs \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT"
echo $cmd
eval $cmd
