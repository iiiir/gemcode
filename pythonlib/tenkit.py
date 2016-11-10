from string import maketrans

## ======= constants ======= ##
dna_to_decimal_dict = maketrans('ACGTacgt', '01230123')
decimal_to_dna_dict = maketrans('0123', 'ACGT')
mycolors = ['0,0,0','255,0,0','0,255,0',
            '0,0,255','255,255,0','0,255,255',
            '255,0,255','192,192,192','128,128,128',
            '128,0,0','128,128,0','0,128,0',
            '128,0,128','0,128,128','0,0,128']

## ========   class  ======== ##
class Barcode():
    def __init__(self):
        self.barcode_file = '/home/swang/app/gemcode/barcode/4M-with-alts-february-2016.txt'
        self.bdict = {}
        self.xdict = {}
        for l in open(self.barcode_file):
            self.bdict[ l.strip() ]     = l.strip()
            self.xdict[ l.strip()[1:] ] = l.strip()

    def find_barcode(self, barcode):
        try:
            return self.bdict[barcode]
        except KeyError:
            if 'N' in barcode:
                try:
                    # first base might be N
                    return self.xdict[barcode[1:]]
                except KeyError:
                    return None
            else:
                return None

class Bed():
    def __init__(self, lin):
        lis = lin.split('\t')
        self.CHROM  = lis[0]
        self.CSTART = lis[1]
        self.CEND   = lis[2]
        self.BARC   = lis[3].split('-')[0]
        self.SCORE  = lis[4]
        self.STRAND = lis[5]
        self.TSTART = lis[6]
        self.TEND   = lis[7]
        self.COLOR  = lis[8]
        self.COUNT  = lis[9]
        self.SIZES  = lis[10].split(',')
        self.STARTS = lis[11].split(',')
        self.cstarts= [int(i)+int(lis[1]) for i in lis[11].split(',')]


## ======== functions ======= ##
def barcode_to_decimals(barcode):
    '''
    must be 16bp barcode
    Note: no "+ 1"
    '''
    chrom = int( barcode[0:3].translate( dna_to_decimal_dict ), 4)
    pos   = int( barcode[3:].translate( dna_to_decimal_dict ), 4) 
    return chrom,pos

def barcode_to_decimal_str(barcode):
    '''
    Convert 16bp barcode (ACGT) to quaternary system
    Note: the "+ 1"
    '''
    chrom = int( barcode[0:3].translate( dna_to_decimal_dict ), 4)
    pos   = int( barcode[3:].translate( dna_to_decimal_dict ), 4)
    return "%s:%s-%s" % (chrom+1,pos+1,pos+1)

def barcode_to_bed(barcodefile):
    '''
    Hangs forever, deprected
    Also see: www.biostars.org/p/49306/
    '''
    outfile = barcodefile + '.bed'
    fout    = open(outfile,'w')
    for barcode in open(barcodefile):
        barcode = barcode.rstrip()
        chrom = int( barcode[0:3].translate( dna_to_decimal_dict ), 4)
        pos   = int( barcode[3:].translate( dna_to_decimal_dict ), 4)
        print >> fout, "%s\t%s\t%s" % (chrom+1, pos, pos+1)
    fout.close()
    return outfile

def get_bcs_from_align(align_obj, min_mapq = None, outbed = False):

    bc_list = {}
    for read in align_obj:
        if not read.cigar: continue
        if read.is_duplicate: continue
        if read.is_secondary: continue
        if min_mapq and (read.mapq < min_mapq): continue

        # get barcode and make sure non-empty
        try:
            bc = read.get_tag('BX')
        except KeyError:
            continue
        if bc is None or bc == '': continue

        # add barcode or readname
        if outbed:
            if bc_list.get(bc,None):
                bc_list[bc][0].append(read.reference_start)
                bc_list[bc][1].append(read.reference_end - read.reference_start)
            else:
                bc_list[bc] = [[],[]]
                bc_list[bc][0] = [read.reference_start]
                bc_list[bc][1] = [read.reference_end - read.reference_start]

        else:
            if bc_list.get(bc,None):
                bc_list[bc].append(read.qname)
            else:    
                bc_list[bc] = [read.qname]

    # only uniq read names if not bed format
    if not outbed:
        for bc in bc_list.keys():
            bc_list[bc] = list(set(bc_list[bc]))
    return bc_list

def mol_supporting_gap(barbed, pos0, pos1):
    import copy
    bc_list = copy.deepcopy(barbed)
    barcode_names_support = []
    uniside_mols = {}
    deletion_range = range(pos0, pos1)
    for barcode, pos_len in bc_list.items():
        positions = pos_len[0]
        # skip molecules if does not cover the deletion
        if len(positions) == 1:
            uniside_mols[barcode] = bc_list.pop(barcode, None)
            continue
        all_left  = all( pos < pos0 for pos in positions )
        all_right = all( pos > pos1 for pos in positions )
        is_uniside = all_left or all_right
        if is_uniside: 
            uniside_mols[barcode] = bc_list.pop(barcode, None)
            continue
        # check if spans deletion
        support = not any( ( pos in deletion_range for pos in positions ) )
        if support:
            barcode_names_support.append(barcode)
    support_deletion     = {k:bc_list.pop(k, None) for k in barcode_names_support }
    not_support_deletion = bc_list # all supporting were removed
    return support_deletion, not_support_deletion, uniside_mols

def pms2_gap_bed(fn_bed12):
    bed = open(fn_bed12)
    next(bed) # remove header
    pms2_end = 6028676
    pms2cl_start = 6775276
    reads_in_regions = {}
    for l in bed:
        b = Bed(l)
        barcode   = b.BARC
        positions = b.cstarts
        pms2_ct   = sum([pos < pms2_end for pos in positions])
        pms2cl_ct = sum([pos > pms2cl_start for pos in positions])
        reads_in_regions[barcode] = (pms2_ct, pms2cl_ct)
    return reads_in_regions

def count_bcs_share(align_obj1, align_obj2):
    '''
    Files  (AlignmentFile obj) -> tuple of counts (barcodes_1_ct, barcodes_shared_ct, barcodes_2_ct) 
    '''
    bcs1 = get_bcs_from_align(align_obj1)
    bcs2 = get_bcs_from_align(align_obj2)
    shared_bcs = set(bcs1.keys()) & set(bcs2.keys())
    return len(bcs1), len( shared_bcs ), len(bcs2)

def get_bcs_share(align_obj1, align_obj2, bed=False):
    '''
    Files  (AlignmentFile obj) -> tuple (barcodes_1, barcodes_shared, barcodes_2) 
    '''
    bcs1  = get_bcs_from_align(align_obj1, outbed=bed)
    bcs2  = get_bcs_from_align(align_obj2, outbed=bed)
    shared_bcs = set(bcs1.keys()) & set(bcs2.keys())
    bcs12 = {k:bcs1.pop(k, None) for k in list(shared_bcs) }
    bcs21 = {k:bcs2.pop(k, None) for k in list(shared_bcs) }
    return bcs1, bcs12, bcs21, bcs2

def align_to_barbed(align_obj, chrom):
    '''
    >>> barbed['GTACGTGTGAATCA']
    [[6252069, 6252061], [128, 151]] # [[pos1,pos2],[len1,len2]]
    '''
    barbed={}
    for read in align_obj:
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

def merge_bed12(d1,d2):
    d = {}
    for k,v1 in d1.items():
        v2 = d2[k]
        d[k] = [ v1[0] + v2[0],v1[1] + v2[1] ]
    return d

def barcode_bam_header():
    '''
    >>> int('TTTTTTTTTTTTT'.translate( dna_to_decimal_dict ), 4)
    67108863
    '''
    header_line = {"HD":{"VN":"1.5","SO":"unsorted"}}
    SQ_lines = [{"SN":"%s"%(i+1),"LN":67108863,"AS":"%s"%(i+1),"SP":"Barcodes"} for i in range(64)]
    header_line["SQ"] = SQ_lines
    header_line["RG"] = [{"ID":"example","SM":"Sample","LB":"10x","PL":"Chromium"}]
    return header_line

def print_bed12(chrom, d, fn, track_name):
    f = open(fn,'w')
    bed12 = 'chr%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t%s\t%d\t%s\t%s'
    print >>f, 'track name=%s itemRgb="On"' % track_name
    for i,b in enumerate(d.keys()):
        color = mycolors[i%15]
        allstarts = d[b][0]
        start,end,name,tstart,tend = min(allstarts), max(allstarts)+150,b,min(allstarts),max(allstarts)+150
        block_counts = len(allstarts)
        bstarts = ','.join([str(i-start) for i in sorted(allstarts)])
        blocks= ','.join([str(i) for i in d[b][1]])
        print >>f, bed12% (chrom, start,end,name,tstart,tend,color,block_counts,blocks,bstarts)
    f.close()

def align_to_cam(sam_line):
    '''
    strip alignment info
    '''
    sam_line.set_tag('YA',"%s.%s" % (sam_line.reference_id, sam_line.reference_start), "Z")
    sam_line.set_tag('YB',"%s.%s" % (sam_line.next_reference_id, sam_line.next_reference_start), "Z")
    sam_line.set_tag('YC',"%s" % sam_line.cigarstring, "Z")
    barcode = sam_line.get_tag('RX')
    barcode_id, barcode_start = barcode_to_decimals(barcode)
    sam_line.reference_id, sam_line.reference_start = barcode_id, barcode_start
    sam_line.next_reference_id, sam_line.next_reference_start = barcode_id, barcode_start
    sam_line.cigartuples = None
    return sam_line
