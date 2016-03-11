#!/bin/env python
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from math import log
import argparse
import pysam
import sys
import os

p = argparse.ArgumentParser()
p.add_argument('--samfile', help='sam file')
p.add_argument('--cuta',type=int, help='left bond')
p.add_argument('--cutb',type=int, help='right bond')
p.add_argument('--norm',action='store_true', help='if the insert size is normal distr')
p.add_argument('--debug', action='store_true',help='pop up figure instead of save figg')
p.add_argument('--flanking', type=int,help='falnking region for discoplot')
args = p.parse_args()

def sam_to_coord(fnsam):
    samfile = pysam.Samfile(fnsam,'rb')
    x=[]
    y=[]
    x_rainbow = []
    y_rainbow = []
    for r1 in samfile:
        r2 = next(samfile)
        assert r1.qname == r2.qname, "r1: %s <-> r2: %s" % (r1.qname, r2.qname)
        # remove paralogous mapping
        if r1.mapq == 0 or r2.mapq == 0: continue
        # if one and only one in the pair is reverse, orientation is proper
        proper_orientation = r1.is_reverse + r2.is_reverse == 1
        if proper_orientation: 
            x.append(r1.pos)
            y.append(r2.pos)
        else:
            x_rainbow.append(r1.pos)
            y_rainbow.append(r2.pos)            
    samfile.close()
    x=np.array(x)
    y=np.array(y)
    x_rainbow = np.array(x_rainbow)
    y_rainbow = np.array(y_rainbow)
    return x,y,x_rainbow,y_rainbow

def iqr(dist):
    return np.percentile(dist, 75) - np.percentile(dist, 25)

def plot_disco(x,y,xx,yy,bed,figname='discoplot'):
    # discodant position plot
    fig = plt.figure(figsize=(11,5))
    p = fig.add_subplot(1,2,1)
    p.scatter(x,y,marker='o', facecolors='none', edgecolor='blue',alpha=0.25)
    p.scatter(xx,yy,marker='x', color='red',alpha=0.25)
    p.set_xlim(bed)
    p.set_ylim(bed)
    p.plot(p.get_xlim(), p.get_ylim(), ls="--", c="0.3")

    vector = x-y
    if args.norm:
        shift_size = np.std(vector) * 2
    else:
        shift_size = 1.5*iqr(vector) #np.percentile(abs(vector), 75)
    p.plot(p.get_xlim() + shift_size, p.get_ylim(), ls="--", c="0.3")
    p.plot(p.get_xlim(), p.get_ylim() + shift_size, ls="--", c="0.3")

    # insertion length plot
    d = fig.add_subplot(1,2,2)
    n, bins, patches = plt.hist(vector, bins=range(-500,500,10))
    y = mlab.normpdf( bins, np.mean(vector), np.std(vector))
    l = d.plot(bins, y, 'r--', linewidth=1)

    # decide figure action
    if args.debug:
        plt.show()
    else:
        fig.savefig('%s.png' % figname,dpt=300)
    

def fn_to_bed(samfile, flanking=1000):
    '''
    note the 1000bp flanking region on each side
    '''
    lis = samfile.split('.')[0].split('_')
    return [ int(lis[-2]) - flanking ,int(lis[-1]) + flanking ]

def main():
    if not args.cuta:
        if args.flanking:
            bed = fn_to_bed( os.path.basename(args.samfile), args.flanking )
        else:
            bed = fn_to_bed( os.path.basename(args.samfile) )
    else:
        bed=[args.cuta, args.cutb]
    x,y,xx,yy = sam_to_coord(args.samfile)
    figname = os.path.basename(args.samfile).split('.')[0]
    plot_disco(x,y,xx,yy,bed,figname)

if __name__ == '__main__':
    main()
