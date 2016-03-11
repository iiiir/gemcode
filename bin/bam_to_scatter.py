import pysam
import matplotlib.pyplot as plt
import numpy as np

fn = "samtools_merged_tinybams_fusion.qn.pe.sam"
f = open(fn)
x=[]
y=[]
for l in f:
    if l.startswith('@'): continue
    lis=l.strip().split('\t')
    x.append(lis[3])
    y.append(lis[7])
