#!/usr/bin/env python
import hmmtop
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
#import numpy.numarray as na
import numpy as na

def calculate(fasta,label,out):
    #os.environ['BIOV_DEBUG'] = 'True'
    ht = hmmtop.tools()
    ht.add_library('c',fasta)
    ht.scan_libraries()
    tms = ht.results['c']
    total = open(fasta,'r').read().count('>')
    zeros = total-len(tms.keys())
    tabs = {}
    for i in tms.values():
        count = len(i)
        tabs.setdefault(count,0)
        tabs[count]+=1
    tabs.setdefault(0,zeros)
    biggest = tabs.keys()
    biggest.sort()
    biggest=biggest[-1]
    [tabs.setdefault(i,0) for i in range(biggest)]
    xlocations = na.array(range(len(tabs.keys())))
    plt.xticks(xlocations+ 0.3/2, xlocations)
    plt.bar(xlocations,tabs.values(),width=0.3)
    plt.title(label)
    plt.savefig(out, dpi=80, format="png")
    return tabs,tms

if __name__=='__main__':
    try:
        (fasta,label,out) = sys.argv[1:4]
        tabs,tms=calculate(fasta,label,out)
        print 'TMSs\tOccurences'
        for k,v in tabs:
            print '%i\t%i'%(k,v)
    except:
        print 'Usage: tmstats.py <fasta_file> <title> <output>'
