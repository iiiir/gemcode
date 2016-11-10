#!/bin/env python
import pysam
import sys
import argparse
from time import gmtime, strftime
import tenkit as tk

p = argparse.ArgumentParser()
p.add_argument('--fastq1', required='true', help="Read 1 FASTQ")
p.add_argument('--fastq2', required='true', help="Read 2 FASTQ")
p.add_argument('--outprefix', required='true', help="output ubam")
p.add_argument('--lines', type=int, default=4000000, help="output ubam")
p.add_argument('--trim', type=int, help="output ubam")
p.add_argument('--suf_len', type=int, default=4, help="output ubam")
p.add_argument('--id', required='true', help="ID")
p.add_argument('--sm', required='true', help="SM")
p.add_argument('--lb', required='true', help="LB")
p.add_argument('--cn', default='STJUDE', help="LB")
p.add_argument('--pl', default='Chromium', help="LB")
p.add_argument('--debug', action='store_true', help="debug mode")
args = p.parse_args()

split_lines   = args.lines
suffix_length = args.suf_len

ctime = strftime("%Y-%m-%d", gmtime())

chromium_barcodes =  tk.Barcode()

def ubam_header():
    header = { 'HD': {'VN': '1.0',
					  'SO':'queryname'},
               'RG': [{'ID': args.id, 
					   'SM': args.sm, 
					   'LB': args.lb,
					   'CN': args.cn, 
                       'PL': args.pl, 
					   'DT': ctime }] }
    return header

def unmapped_read(fq_read, sam_flag, chromium_tags, trim=None):
    ubam_read = pysam.AlignedSegment()
    ubam_read.query_name = fq_read.name
    if trim:
        ubam_read.query_sequence = fq_read.sequence[0:trim]
        ubam_read.query_qualities = pysam.qualitystring_to_array( fq_read.quality[0:trim] )
    else:
        ubam_read.query_sequence = fq_read.sequence
        ubam_read.query_qualities = pysam.qualitystring_to_array( fq_read.quality )
    ubam_read.flag = sam_flag
    ubam_read.reference_id = -1
    ubam_read.reference_start = -1
    ubam_read.next_reference_id = -1
    ubam_read.next_reference_start = -1
    ubam_read.tags = chromium_tags 
    return ubam_read

def fq_to_ubam(fq_r1, fq_tag, opre):
    reads1  = pysam.FastxFile(fq_r1);
    reads2  = pysam.FastxFile(fq_tag);

    file_count = 0
    finished_lines = 0
    outfile = None

    for fq_r1 in reads1:
        fq_r2    = next(reads1)
        tag_read = next(reads2)

        assert  fq_r1.name == fq_r2.name, "error at line %d" % finished_lines

        fq_tag = tag_read.sequence
        fq_taq = tag_read.quality

        barcode      = fq_r1.sequence[0:16]
        barcode_qual = fq_r1.quality[0:16]
        cbarcode     = chromium_barcodes.find_barcode(barcode) 
        if cbarcode:
            tags = (('RX', cbarcode),
                    ('QX', barcode_qual),
                    ('BX', '%s-1' % cbarcode),
                    ('BC', fq_tag),
                    ('QT', fq_taq), 
                    ('RG', args.id) )
        else:
            tags = (('RX', barcode),
                    ('QX', barcode_qual),
                    ('QT', fq_taq),
                    ('RG', args.id) )
        fq_r1.sequence = fq_r1.sequence[23:]
        fq_r1.quality  = fq_r1.quality[23:]
        if args.trim:
            r1 = unmapped_read( fq_r1, 77, tags, trim=args.trim)
            r2 = unmapped_read( fq_r2, 141, tags, trim=args.trim)
        else:
            r1 = unmapped_read( fq_r1, 77, tags)
            r2 = unmapped_read( fq_r2, 141, tags)

        finished_lines += 1
        if finished_lines % split_lines == 1:
            file_count += 1
            osuff = str(file_count).zfill(suffix_length)
            ofn = '%s.s%s.%s' % (opre, osuff, 'u.bam')
            if outfile is not None:
                outfile.close()
            outfile = pysam.AlignmentFile(ofn, "wb", header = ubam_header())
            
        outfile.write(r1)
        outfile.write(r2)
    if outfile is not None:
        outfile.close()

def main(args):
    fq_to_ubam(args.fastq1, args.fastq2, args.outprefix)

if __name__ == '__main__':
    main(args)
