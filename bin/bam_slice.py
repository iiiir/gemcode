#!/bin/env python

# bam_slice.py \
#	--bam phased.bam \
#	--gpos "chr1:87,337,010-chr1:87,337,012;chr10:36,119,126-chr10:36,119,128" \
#	--output tiny_bam \
#	-j tiny_bam.sjm

import re
import sys
import argparse
import sjm
import util

p = argparse.ArgumentParser('bam_slice.py --bam ... --gpos ... --output ... -j ...')
p.add_argument('--bam',help='a file contains a list of bams, one bam a line')
p.add_argument('--tempfolder', default="/scratch_space/swang/small_bams", help='Locatoin to hold temp files')
p.add_argument('--gpos', nargs='+', help='genomics position of first breakpoint')
p.add_argument('--output',default='samtools_merged_tinybams', help='output bam name (will append chr_start_end)')
p.add_argument('-A','--account', metavar='STR', help='Support for aligned and dedupped BAMs as input')
p.add_argument('-j', '--jobfile', metavar='FILE', help='The jobfile name (default: stdout)')
args = p.parse_args()

sjm.Job.name_prefix="FUSIONBAM"+"."
sjm.Job.memory="5G" # default if not provided
sjm.Job.queue="pcgp"
sjm.Job.project="CompBio"
if args.account: sjm.Job.sge_options="-A %s" % args.account

def gpos2tup(gpos):
    positions = gpos#.strip().split(";")
    pos_out   = []
    for p in positions:
        try:
            chrom, pos0, chromb,pos1 = re.split(r'[:-]', p.replace(',',''))
        except:
            chrom, pos0, pos1 = re.split(r'[:-]', p.replace(',',''))
        #chrom = chrom.lstrip('chr')
        pos_out.append((chrom,pos0,pos1))
    return pos_out

def bams2reads(bam, gpos):
    jobs = []
    pos = gpos2tup(gpos)
    for i,pos_tup in enumerate(pos):
        chrom,pos0,pos1 = pos_tup
        bamfile = util.File(bam)
        job = sjm.Job('get_fusion_bam-%s-part%d'%(bamfile.prefix,i))
        job.memory = "1G"
        job.output = util.File(bamfile.chdir(args.tempfolder).chext( '%s.%s.%s.bam' % (chrom, pos0, pos1) ))
        job.append("samtools view -hb %s %s:%s-%s > %s" % (bam, chrom, pos0, pos1, job.output) )
        jobs.append(job)
    return jobs

def merge_bams(pjobs):
    tiny_bams = [pjob.output.path for pjob in pjobs]
    gpos = gpos2tup(args.gpos)
    bam_suffix = '_'.join(['_'.join(coord) for coord in gpos])
    obam = "%s_%s.bam" % (args.output, bam_suffix)
    bamfile = util.File(obam)
    job = sjm.Job('merge_bams-%s'%(bamfile.prefix))
    job.append("samtools merge %s %s && samtools index %s" % (obam,' '.join(tiny_bams), obam))
    job.output = bamfile
    job.depend(*pjobs)
    return job

def mv_bams(pjobs):
    pjob=pjobs[0]
    gpos = gpos2tup(args.gpos)
    bam_suffix = '_'.join(['_'.join(coord) for coord in gpos])
    obam = "%s_%s.bam" % (args.output, bam_suffix)
    bamfile = util.File(obam)
    job = sjm.Job('moving-%s'%(bamfile.prefix))
    job.append("mv %s %s" % (pjob.output, obam))
    job.output = bamfile
    job.depend(pjob)
    return job

def bam_qn_pe(pjob):
    bamfile = pjob.output
    job = sjm.Job('bam_to_qnpe-%s'%(bamfile.prefix))
    job.append("sam_sort_qn.sh %s" % bamfile.path)
    job.depend(pjob)
    return job

def bam_to_bed(pjob):
    pass
    #bamtobed.py --bam DEL18_LOW_chr6_157698904_157704269.qn.pe.bam --chr 6

def main():
    if args.jobfile is None:
        jobfile=None
    else:
        jobfile=util.File(args.jobfile)

    jobs = bams2reads(args.bam, args.gpos)
    if len(args.gpos) > 1:
        job  = merge_bams(jobs)
    else:
        job  = mv_bams(jobs)
    job  = bam_qn_pe(job)
    descout = sys.stdout if jobfile is None else open(jobfile.path, "w")
    descout.write(sjm.Job().depend(job).desc())
    descout.flush()

if __name__ == '__main__':
    main()
