#!/usr/bin/env python

from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
import tempfile
import hmmtop
import re
import sys

min_tms = 10
max_gap = 50
fasta_file = open(sys.argv[1],'r')
ss_results = open(sys.argv[2],'r')
subjects = open(sys.argv[3],'r')
targets = open(sys.argv[4],'r')
outfile = open(sys.argv[5],'wb+')
subjects = SeqIO.parse(subjects,'fasta')
subjects= SeqIO.to_dict(subjects)
targets= SeqIO.parse(targets,'fasta')
targets= SeqIO.to_dict(targets)
hmt = hmmtop.tools()
hmt.add_library('SX',fasta_file.name)
hmt.scan_libraries()
tmss = hmt.results['SX']

for line in ss_results:
    data = line.split("\t")
    subject = data[0]
    target = data[1]
    try:
        if len(tmss[subject]) < min_tms or len(tmss[target]) < min_tms:
            continue
    #needle.asequence='asis:'+str(subjects[subject].seq)
    #needle.bsequence='asis:'+str(targets[target].seq)
    #mygap = float(gaps.search(needle()[0]).groups()[0])
    #if mygap <= max_gap:
        outfile.write(line)
        print line,
    except:
        pass
