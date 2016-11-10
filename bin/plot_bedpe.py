#!/bin/env python
import sys
import pysam
import tenkit as tk
import re
import csv
import uuid
import time
import os

# plot_bedpe.py in.bam in.bedpe

fnbam  = os.path.abspath( sys.argv[1] )
fnbed  = os.path.abspath( sys.argv[2] )

flank = 50000

SAMfile = pysam.AlignmentFile(fnbam)
fbed  = open(fnbed)
next(fbed)                                  # remove first line of header
fns   = next(fbed).lstrip("#").split("\t")  # read header seperately
BEDfile = csv.DictReader( fbed, fieldnames=fns, delimiter="\t" )

tmpfile = os.path.join("/scratch_space/swang/test/",
                       "plot_bedpe_" + time.strftime("%Y%M%d%H%M") + "_"
                       + str(uuid.uuid4()).split('-')[-1] + ".txt")
ftmp   = open(tmpfile, 'w')

for i,sv_line in enumerate(BEDfile):
    sam1 = SAMfile.fetch( sv_line['chrom1'], int(sv_line['start1']) - flank, int(sv_line['stop1']) + flank )
    sam2 = SAMfile.fetch( sv_line['chrom2'], int(sv_line['start2']) - flank, int(sv_line['stop2']) + flank )

    bc1,bc12,bc21,bc2 = tk.get_bcs_share(sam1,sam2, bed=True)


    out_beda = 'fusion%s_%s_%s_%s.bed'%(i, sv_line['chrom1'], sv_line['start1'], sv_line['stop1'])
    out_bedb = 'fusion%s_%s_%s_%s.bed'%(i, sv_line['chrom2'], sv_line['start2'], sv_line['stop2'])
    tk.print_bed12(sv_line['chrom1'], bc12, out_beda, 'left')
    tk.print_bed12(sv_line['chrom2'], bc21, out_bedb, 'right')

    print >>ftmp, "plot_fusion.R -a %s -b %s -s %s -e %s -S %s -E %s -o %s" % (out_beda, out_bedb, sv_line['start1'], sv_line['stop1'], sv_line['start2'], sv_line['stop2'], "fusion_%d.pdf" % i)
ftmp.close()
print "bsub_array_for_cmdfile.sh %s -M 5000 -q pcgp" % tmpfile
