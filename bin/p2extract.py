#!/usr/bin/env python
from Bio import SeqIO
import protocol2
import sys
import os

usage = 'p2extract.py <minsd> <count> <outfile>'

if len(sys.argv) != 4:
    print usage
    quit()

p2 = protocol2.Compare()
p2.outdir = '.'
p2.tssearch()
p2.run_global_alignments()
minsd,count,outfile = sys.argv[1:4]
minsd,count = int(minsd),int(count)
p2.globaldata.sort(key=lambda x:x['gsat_score'],reverse=True)
subjects,targets = [],[]
for result in p2.globaldata:
    score,subject,target = result['gsat_score'],result['subject_symbol'],result['target_symbol']
    if score >= minsd:
        subjects.append(subject)
        targets.append(target)
subjects,targets = list(set(subjects)),list(set(targets))
mysubjects = subjects[0:count]
mytargets = targets[0:count]

subject_file = SeqIO.parse('subjects.faa','fasta')
target_file = SeqIO.parse('targets.faa','fasta')
subject_file = SeqIO.to_dict(subject_file)
target_file = SeqIO.to_dict(target_file)
subject_select = [subject_file[i] for i in mysubjects]
target_select = [target_file[i] for i in mytargets]
subject_select.extend(target_select)

handle = open(outfile,'wb')
for i in subject_select:
    SeqIO.write(i,handle,'fasta')
