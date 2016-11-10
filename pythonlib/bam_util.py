import pysam
import tenkit

def print_n_reads(inbamname, obamname, n=5000, print_secondary=None):
    inbamfile = pysam.AlignmentFile(inbamname,'rb')
    obamfile  = pysam.AlignmentFile(obamname,'wb',template=inbamfile)
    sbamfile  = pysam.AlignmentFile(obamname.rstrip('.bam') + '.supp.bam','wb',template=inbamfile) # supplimentary, contain all mappings
    counter   = 0
    for read in inbamfile:
        if read.is_secondary: 
            sbamfile.write(read)
        elif read.is_supplementary:
            sbamfile.write(read)
        else:
            obamfile.write(read)
            sbamfile.write(read)
            counter = counter + 1
        if counter >= n: break
    obamfile.close()

def print_reads_by_bc(inbamfile, obamname, bc, headerfile):
    # inbamfile is just a region
    outsam = pysam.AlignmentFile(obamname,'wb', template = headerfile)
    for sam_line in inbamfile:
        if sam_line.is_secondary: continue
        if sam_line.is_supplementary: continue
        if not sam_line.has_tag('BX'): continue
        if sam_line.get_tag('BX') == bc:
            outsam.write(sam_line)
    outsam.close()

def compare_mapped(r1, r2):
    if r1.reference_name  == r2.reference_name  and \
       r1.reference_start == r2.reference_start:
        return 'M'
    elif r1.reference_name  == r2.reference_name and not \
         r1.reference_start == r2.reference_start:
        return 'C'
    else:
        return 'D'

def compare_align(r1,r2):
    assert r1.qname == r2.qname, "read name not match %s and %s" % (r1.qname, r2.qname)
    if not r1.is_unmapped and not r2.is_unmapped:
        map_report = compare_mapped(r1,r2)
    elif not r1.is_unmapped and r2.is_unmapped:
        map_report = '1'
    elif r1.is_unmapped and not r2.is_unmapped:
        map_report = '2'
    else:
        map_report = 'U'
    return map_report

def comapre_bams(bamfn1, bamfn2):
    '''
    large files might be slow for zip()
    '''
    bamfile1 = pysam.AlignmentFile(bamfn1)
    bamfile2 = pysam.AlignmentFile(bamfn2)
    for r1 in bamfile1: #, r2 in zip(bamfile1,bamfile2):
        r2 = next(bamfile2)
        # in case of secondary mapping, keep finding next primary
        while r1.is_secondary or r1.is_supplementary: r1=next(bamfile1)
        while r2.is_secondary or r2.is_supplementary: r2=next(bamfile2)
        map_report = compare_align(r1,r2)
        if map_report in ['C', 'D'] and False:
            print '%s/%s\t%s:%s(%s/%s)\t%s:%s(%s/%s)\t%s' % (r1.qname, '1' if r1.is_read1 else '2',
                                                 r1.reference_name, r1.reference_start, r1.mapq, r1.cigarstring, 
                                                 r2.reference_name, r2.reference_start, r2.mapq, r2.cigarstring, map_report)
        elif map_report == "2":
            print '%s/%s\t-\t%s:%s(%s/%s)\t%s' % (r1.qname, '1' if r1.is_read1 else '2',
                                                 r2.reference_name, r2.reference_start, r2.mapq, r2.cigarstring, map_report)
        elif map_report == "2":
            print ">%s/%s|%s.%s|%s|%s\n%s" % (r2.qname, '1' if r1.is_read1 else '2', r2.reference_name, r2.reference_start, r2.mapq, r2.cigarstring, r2.query_sequence)

def convert_to_barbam(fn,fo):
    f   = pysam.AlignmentFile(fn,'rb')
    outsam = pysam.AlignmentFile(fo,'wb', header = tenkit.barcode_bam_header())
    bc  = tenkit.Barcode()
    for sam_line in f:
        sam_line.set_tag('YA',"%s.%s" % (sam_line.reference_id, sam_line.reference_start), "Z")
        sam_line.set_tag('YB',"%s.%s" % (sam_line.next_reference_id, sam_line.next_reference_start), "Z")
        barcode = sam_line.get_tag('RX')
        if 'N' in barcode:
            barcode = bc.find_barcode( sam_line.get_tag('RX') )
        if not bc.bdict.get(barcode,None):
            print sam_line.get_tag('RX')
        barcode_id, barcode_start = tenkit.barcode_to_decimals(barcode)
        sam_line.reference_id, sam_line.reference_start = barcode_id, barcode_start
        sam_line.next_reference_id, sam_line.next_reference_start = barcode_id, barcode_start
        outsam.write(sam_line)
    f.close()
    outsam.close()

def check_barcodes(fn):
    f   = pysam.AlignmentFile(fn,'rb')
    bc  = tenkit.Barcode()
    tot = 0
    mis = []
    for sam_line in f:
        barcode = sam_line.get_tag('RX')
        tot = tot + 1
        if 'N' in barcode:
            barcode = bc.find_barcode( sam_line.get_tag('RX') )
        if not bc.bdict.get(barcode,None):
            mis.append( sam_line.get_tag('RX') )
    print "\n".join(mis)
    print "There are %d/%d not found in %s" % (len(mis), tot, fn)
    f.close()

