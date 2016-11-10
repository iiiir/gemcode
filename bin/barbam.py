#!/bin/env python
import sys
import pysam
from string import maketrans
import argparse

p = argparse.ArgumentParser(description='Converting bam to cam')
p.add_argument('-b', '--inbam', metavar='FILE', required=True, help='The BAM input')
p.add_argument('-o', '--outcam', metavar='FILE', default="-", help='The output cam')
args = p.parse_args()
fn  = args.inbam
fo  = args.outcam

dna_to_decimal_dict = maketrans('ACGTacgt', '01230123')
decimal_to_dna_dict = maketrans('0123', 'ACGT')

def barcode_to_decimals(barcode):
    try:
        chrom = int( barcode[0:3].translate( dna_to_decimal_dict ), 4)
        pos   = int( barcode[3:].translate( dna_to_decimal_dict ), 4)
    except ValueError:
        chrom = -1
        pos   = -1
    return chrom,pos

def barbam_header(original_header):
    header_line = {"HD":{"VN":"1.5","SO":"unsorted"}}
    SQ_lines = [{"SN":"%sbc"%(i+1),"LN":67108863,"AS":"%s"%(i+1),"SP":"Barcodes"} for i in range(64)]
    header_line["SQ"] = SQ_lines + original_header['SQ']
    header_line["RG"] = original_header['RG']
    header_line["PG"] = original_header['PG']
    return header_line

def convert_to_barbam(f,fo):
    outsam = pysam.AlignmentFile(fo,'wb', header = barbam_header(f.header))
    for sam_line in f:
        sam_line.set_tag('YA',"%s.%s" % (sam_line.reference_id, sam_line.reference_start), "Z")
        sam_line.set_tag('YB',"%s.%s" % (sam_line.next_reference_id, sam_line.next_reference_start), "Z")
        sam_line.set_tag('YC',"%s" % sam_line.cigarstring, "Z")
        sam_line.cigartuples = None
        barcode = sam_line.get_tag('RX')
        barcode_id, barcode_start = barcode_to_decimals(barcode)
        sam_line.reference_id, sam_line.reference_start = barcode_id, barcode_start
        sam_line.next_reference_id, sam_line.next_reference_start = barcode_id, barcode_start

        outsam.write(sam_line)
    outsam.close()

def main(fn):
    f   = pysam.AlignmentFile(fn,'rb')
    convert_to_barbam(f,fo)

if __name__ == "__main__":
    main(fn)
