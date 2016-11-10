#!/bin/env python
import sys
import os
import pysam
import tenkit as tk
import re
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from sklearn.neighbors import KernelDensity

# homozygous del
# run_find_deletion_support.py /nfs_exports/genomes/1/projects/WHOLEGENOME/GemcodeEval/BucketRaw/NA12878/NA12878_b37/NA12878.b37.mdup.bam 2    104186941    104187136 200
# run_find_deletion_support.py /nfs_exports/genomes/1/projects/WHOLEGENOME/GemcodeEval/BucketRaw/SJMEL829_D/SJMEL829_D/SJMEL829_D.mdup.bam 1 224782795 224786034 50000 SJMEL829_1_224782795_224786034_50000.pdf

fn    = sys.argv[1]
chrom = sys.argv[2]
pos1  = int(sys.argv[3])
pos2  = int(sys.argv[4])
flank = int(sys.argv[5])
if sys.argv[6]:
    fpdf  = sys.argv[6]
else:
    fpdf  = "DEL_%s_%s_%s_%s_%s.pdf" % (os.path.basename(fn).split('_')[0], chrom, pos1, pos2, flank)

SAMfile = pysam.AlignmentFile(fn)

sam = SAMfile.fetch(chrom, pos1-flank, pos2+flank)
bcs = tk.get_bcs_from_align(sam, outbed=True)
supp_del,ns,u = tk.mol_supporting_gap(bcs, pos1 + 50, pos2 - 50 )

# plot bed
fig = plt.figure(figsize=(8, 11))
hax = fig.add_subplot(211)
keys = supp_del.keys()
x_supp = np.array([supp_del[k][0] for k in keys])
y_supp = np.array([[i+0.5]*len(supp_del[k][0]) for i,k in enumerate(keys)])
for x,y in zip(x_supp, y_supp):
    hax.plot(x, y,'bo:', mec=None, mew=0.0, ms=3, lw=0.005)

hax = fig.add_subplot(212)
keys = ns.keys()
x_supp = np.array([ns[k][0] for k in keys])
y_supp = np.array([[i+0.5]*len(ns[k][0]) for i,k in enumerate(keys)])
for x,y in zip(x_supp, y_supp):
    hax.plot(x, y,'bo:', mec=None, mew=0.0, ms=3, lw=0.005)

#hax = fig.add_subplot(313)
#keys = u.keys()
#x_supp = np.array([u[k][0] for k in keys])
#y_supp = np.array([[i+0.5]*len(u[k][0]) for i,k in enumerate(keys)])
#for x,y in zip(x_supp, y_supp):
#    hax.plot(x, y,'bo', mec=None, mew=0.0, ms=3)

pp = PdfPages(fpdf)
plt.savefig(pp, format='pdf')
pp.close()

print >>sys.stdout, 'S|NS|U:%d|%d|%d' % ( len(supp_del), len(ns), len(u) )

## wrote three bed files ##
#bed12 = 'chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\t%d\t%s\t%s'
#prefix_name = os.path.basename( fn ).split('.')[0] + '_%s_%d_%d_flank_%d' % (chrom, pos1, pos2, flank)
#tk.print_bed12(chrom, supp_del, prefix_name+'.suppdel.bed', 'Support_'+prefix_name)
#tk.print_bed12(chrom, ns, prefix_name+'.notsupp.bed', 'NotSupp_'+prefix_name)
#tk.print_bed12(chrom, u,  prefix_name+'.uniside.bed', 'Uniside_'+prefix_name)

