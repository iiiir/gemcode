#!/bin/env python
import sys
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
# run_find_fusion_support.py /nfs_exports/genomes/1/projects/WHOLEGENOME/GemcodeEval/BucketRaw/SJMEL829_D/SJMEL829_D/SJMEL829_D.mdup.bam 18:9868617 10:7634373 10000 test.pdf

fn    = sys.argv[1]
pos1  = sys.argv[2]
pos2  = sys.argv[3]
flank = int(sys.argv[4])
fpdf  = sys.argv[5]

SAMfile = samfile1 = pysam.AlignmentFile(fn)

funs = [str,int,int]
if "-" in pos1:
    pos1 = [f(v) for f,v in zip(funs,re.split(r"[:-]",pos1))]
else:
    chrom1, pos = pos1.split(':')[0], int(pos1.split(':')[1])
    pos1 = [chrom1, pos - flank, pos + flank]
if "-" in pos2:
    pos2 = [f(v) for f,v in zip(funs,re.split(r"[:-]",pos2))]
else:
    chrom2, pos = pos2.split(':')[0], int(pos2.split(':')[1])
    pos2 = [chrom2, pos - flank, pos + flank]

sam1 = SAMfile.fetch(pos1[0], pos1[1], pos1[2])
sam2 = SAMfile.fetch(pos2[0], pos2[1], pos2[2])

bc1,bc12,bc21,bc2 = tk.get_bcs_share(sam1,sam2, bed=True)

tk.print_bed12(chrom1, bc12, '1.bed', 'left')
tk.print_bed12(chrom2, bc21, '2.bed', 'right')

if pos2[0] == pos1[0]:
    distance = str(pos2[1] - pos1[2])
else:
    distance = "INT"

print >>sys.stdout, 'S1|S2|U1|U2|D:%d|%d|%d|%d|%s' % ( len(bc12.keys()), len(bc21.keys()), len(bc1.keys()), len(bc2.keys()), distance )

# output bed plot
keys = bc12.keys()
min12 = np.array([min(bc12[k][0]) for k in keys])
max12 = np.array([max(bc12[k][0]) for k in keys])

min21 = np.array([min(bc21[k][0]) for k in keys])
max21 = np.array([max(bc21[k][0]) for k in keys])

# == plot first panel == #
fig = plt.figure(figsize=(11, 8))
hax = fig.add_subplot(411)
# all alignment (start positions)
align12x = np.array([bc12[k][0] for k in keys]) 
align12y = np.array([[i+0.5]*len(bc12[k][0]) for i,k in enumerate(keys)])
for x,y in zip(align12x, align12y):
    hax.plot(x, y,'bo:', mec=None, mew=0.0, ms=3)
#hax.hlines(range(len(min12)), min12, max12, lw=2)
hax.set_xlabel('%s'%pos1[0])

y12 = np.hstack(align12x)
hax = fig.add_subplot(412)
density12 = stats.kde.gaussian_kde(y12)
x12 = np.linspace(np.min(y12),np.max(y12) , 1000)
plt.plot(x12, density12(x12))


# == plot second panel == #
hax = fig.add_subplot(413)
align21x = np.array([bc21[k][0] for k in keys])
align21y = np.array([[i+0.5]*len(bc21[k][0]) for i,k in enumerate(keys)])
for x,y in zip(align21x, align21y):
    hax.plot(x, y,'ro:', mec=None, mew=0.0, ms=3)
#hax.hlines(range(len(min21)), min21, max21, colors='r', lw=2)
hax.set_xlabel('%s'%pos2[0])
#hax.set_title('Linked Reads at %s'%pos1[0])

# do density plot
y21 = np.hstack(align21x)
hax = fig.add_subplot(414)
density21 = stats.kde.gaussian_kde(y21)
x21 = np.linspace(np.min(y21),np.max(y21) , 1000)
plt.plot(x21, density21(x21))
# y = np.hstack(align21x)[:, np.newaxis]
#x = np.linspace(np.min(y),np.max(y) , 100)[:, np.newaxis]
#from sklearn.grid_search import GridSearchCV
#grid = GridSearchCV(KernelDensity(),
#                    {'bandwidth': np.linspace(0.1, 1.0, 30)},
#                    cv=20)
#grid.fit(y)
#print grid.best_params_
#log_dens = KernelDensity(kernel="tophat",bandwidth=1.0).fit(y).score_samples(x)
#plt.plot(x,np.exp(log_dens))

# save pdf
pp = PdfPages(fpdf)
plt.savefig(pp, format='pdf')
pp.close()

