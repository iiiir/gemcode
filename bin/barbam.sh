#!/bin/bash -eu

echo "*** Converting bam to barcode sorted cam ***"

[[ $# -lt 2 ]] && echo "$0 <in.bam> <out_prefix>"

bam=`cd \`dirname $1\`; pwd`/`basename $1`; shift
o=`cd \`dirname $1\`; pwd`/`basename $1`; shift

[[ $o = *.cam ]] && o=${o%.cam}

cmd="barbam.py -b $bam | samtools sort - $o"
eval $cmd

mv $o.bam $o.cam && sam_index.sh $o.cam

echo "*** Finished convert bam to cam ***"
