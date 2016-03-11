#!/bin/env python
import sys

# sam_clean_unpaired.py qn.sorted.sam > qn.sorted.pe.sam
# samtools view -bS test.sam > test.bam

lastread=''
lastreadname=''

f = open(sys.argv[1])

for l in f:
    if l.startswith('@'):
        print l.strip()
    else:
        readname=l.split('\t')[0]
        if readname == lastreadname:
            print lastread
            print l.strip()
            lastread=''
            lastreadname=''
        else:
            lastread = l.strip()
            lastreadname = readname
