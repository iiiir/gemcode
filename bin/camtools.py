#!/bin/env python
import sys
from string import maketrans
import argparse
import tenkit as tk
import pysam

p = argparse.ArgumentParser('camtools.py NA12878.cam AAAAAAAAAAAAAAAA | samtools sort - a.sort')
p.add_argument('--barcodelist', '-B', action='store_true',help="barcode to query")
p.add_argument('bam_barcode', nargs='+', help="BAM or bgzipped, must be indexed")
args = p.parse_args()

fn, bc = args.bam_barcode

dna_to_decimal_dict = maketrans('ACGTacgt', '01230123')
decimal_to_dna_dict = maketrans('0123', 'ACGT')

def clean_header(original_header):
    header_line = {"HD":{"VN":"1.5","SO":"unsorted"}}
    header_line["RG"] = original_header['RG']
    header_line["PG"] = original_header['PG']
    header_line["SQ"] = [sq_line for sq_line in original_header['SQ'] if sq_line.get('M5', None)]
    return header_line

def cam_to_bam(sam_line):
    refline = sam_line.get_tag('YA').split('.')
    altline = sam_line.get_tag('YB').split('.')
    sam_line.cigarstring = sam_line.get_tag('YC')
    sam_line.reference_id, sam_line.reference_start           = int(refline[0]), int(refline[1])
    sam_line.next_reference_id, sam_line.next_reference_start = int(altline[0]), int(altline[1])
    # clean tags
    sam_line.set_tag('YA',None)
    sam_line.set_tag('YB',None)
    sam_line.set_tag('YC',None)
    return sam_line

if __name__ == "__main__":
    if args.barcodelist:
        # be aware that multipel iterator might cause overhead
        samfile = pysam.AlignmentFile(fn, multiple_iterators=True)
        outhead = clean_header( samfile.header )
        outsam  = pysam.AlignmentFile('-', 'wb', header=outhead)
        for l in open(bc):
            b = l.strip()
            chrom, start = tk.barcode_to_decimals(b)
            sam = samfile.fetch("%d"%(chrom+1), start, start+1)
            for sam_line in sam:
                outsam.write(cam_to_bam(sam_line))
    else:
        samfile = pysam.AlignmentFile(fn)
        outhead = clean_header( samfile.header )
        chrom, start = tk.barcode_to_decimals(bc)
        sam = samfile.fetch("%dbc"%(chrom+1), start, start+1)
        outsam  = pysam.AlignmentFile('-', 'wb', header=outhead)
        for sam_line in sam:
            outsam.write(cam_to_bam(sam_line))
    outsam.close()
