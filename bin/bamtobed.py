#!/bin/env python

#  bamtobed.py --bam in.bam --chr 1
import os
import sys
import pysam
import argparse

p  = argparse.ArgumentParser()
p.add_argument('--bam', metavar='FILE', required=True,help='bam file')
p.add_argument('--chr', type=int, required=True, help='chromosome')
p.add_argument('--trackname', help='track anme')
p.add_argument('--cuta', type=int, help='left cutoff')
p.add_argument('--cutb', type=int, help='right cutoff')
args = p.parse_args()

def find_bc(tags):
    v = None
    for c in tags:
        if c[0] == 'BX':
            v = c[1].split('-')[0]
    return v

def cigar_to_length(cigar):
    '''
    cigar_to_length( [(0, 64), (4, 26)] )
    >>> (0, 64, 26)
    cigar_to_length( [(1, 64), (2, 123), (0, 10), (4, 26)] )
    >>> (187, 10, 26)
    '''
    # find index of first M (assume no indel in read)
    # t[0] is the cigar symbol, 0 means M
    i = next( (i for i,t in enumerate(cigar) if t[0] == 0) )
    
    # find left cut
    if i == 0: 
        left_cut = 0
    else: 
        left_cut = sum( l for k,l in cigar[0:i] )
    
    # find right cut
    if i+1 == len(cigar): 
        right_cut = 0
    else: 
        right_cut = sum( l for k,l in cigar[i+1:] )
    
    middle_piece = sum(l for k,l in cigar) - left_cut - right_cut

    return left_cut, middle_piece, right_cut

def barcode_span(barcodes, deletion_range):
    barcode_names_spans = []
    for barcode, pos_len in barcodes.items():
        positions = pos_len[0]
        # check if spans deletion
        spans = any( ( pos in deletion_range for pos in positions ) )
        if not spans:
            barcode_names_spans.append(barcode)
    span_dict = {}
    for barcode in barcode_names_spans:
        span_dict[barcode] = barcodes[barcode]
        del barcodes[barcode]
    return span_dict, barcodes

def colors():
    '''
    total 15 colors
    '''
    return ['0,0,0','255,0,0','0,255,0',
    '0,0,255','255,255,0','0,255,255',
    '255,0,255','192,192,192','128,128,128',
    '128,0,0','128,128,0','0,128,0',
    '128,0,128','0,128,128','0,0,128']

def sam_to_barcodes(fn, chrom):
    '''
    >>> barcodes['GTACGTGTGAATCA']
    [[6252069, 6252061], [len1, len2]]
    '''
    samfile = pysam.Samfile(fn,'rb')
    barcodes={}
    for read in samfile:
        k = find_bc(read.tags)
        if not k: continue
        if not read.rname == chrom: continue
        if not read.cigar: continue
        read_left, read_m, read_right = cigar_to_length(read.cigar)
        # if barcode exists, append new coordinates
        if barcodes.get(k,None):
            barcodes[k][0].append(read.pos - read_left)
            barcodes[k][1].append(read_m)
        else:
            barcodes[k] = [[],[]]
            barcodes[k][0] = [read.pos - read_left]
            barcodes[k][1] = [read_m]
    return barcodes

def barcodes_to_bed(barcodes,fnbed,mychrom):
    bedfile=open(fnbed,'w')
    mycolors = colors()
    if args.trackname:
        tname = '%s_%s'%(args.trackname,mychrom)
    else:
        tname = '%s'%(os.path.basename(args.bam.rstrip('.bam')))
    print >> bedfile, "track name=%s" % tname
    print >> bedfile, 'itemRgb="On"'
    mybed = 'chr%d\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\t%d\t%s\t%s'
    for i,b in enumerate(barcodes.keys()):
        color = mycolors[i%15]           ## only iterate between 15 colors
        allstarts = barcodes[b][0]
        start,end,name,tstart,tend = min(allstarts), max(allstarts)+150,b,min(allstarts),max(allstarts)+150
        block_counts = len(allstarts)
        bstarts = ','.join([str(i-start) for i in sorted(allstarts)])
        blocks= ','.join([str(i) for i in barcodes[b][1]])
        print >>bedfile, mybed% (mychrom, start,end,name,tstart,tend,color,block_counts,blocks,bstarts)
    bedfile.close()

def main():
    fn = args.bam
    mychrom = args.chr
    fnbed=fn.rstrip('bam') + '%d' % mychrom
    barcodes = sam_to_barcodes(fn, mychrom-1)      # samformat is 0 based even for chrom
    bed = xrange(args.cuta, args.cutb)
    if args.cuta:
        fnbed_span = fnbed + '.span.bed'
        fnbed_nospan = fnbed + '.nospan.bed'
        barcode_spans, barcode_nospans = barcode_span( barcodes, bed )
        barcodes_to_bed( barcode_spans, fnbed_span, mychrom)
        barcodes_to_bed( barcode_nospans, fnbed_nospan, mychrom)
    else:
        fnbed = fnbed + '.bed'
        barcodes_to_bed(barcodes,fnbed, mychrom)

if __name__ == '__main__':
    main()
