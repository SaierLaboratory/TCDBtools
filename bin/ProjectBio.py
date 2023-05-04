#!/usr/bin/env python3
# ProjectBio is an extention to BioPython

# ParseDefline() - Will return a SeqIO Object based on just the defline.
# nr_dict() - Like SeqIO.to_dict, but will tolerate redundant keys, to generate a non-redundant dict.


from Bio import SeqIO
from Bio import AlignIO
import subprocess
import os
import tempfile,re
from numpy import sqrt,average

def ParseDefline(defline,blast_title=False):
    if blast_title:
        line= defline.split(" ")
        defline=" ".join(line[1:])
    defline = re.sub('(\w)#','\\1 #',defline)
    defline = re.sub('>','',defline)
    sequence = ">%s\nXXX" %(defline)
    fasta = tempfile.NamedTemporaryFile(mode='wt+',delete=True)
    fasta.write(sequence)
    fasta.flush()
    fasta.seek(0)
    read = SeqIO.read(fasta,'fasta')
    return read

def nr_dict(object):
    res = {}
    for i in object:
        #if i.id in list(res.keys()):
        if i.id in res.keys():
            continue
        res.setdefault(i.id,i)
    return res

def ParseTC(defline):
    findtc = re.compile(r'\d+\.\w+\.\d+\.\d+\.\d+')
    acc = re.search('(gi\||gnl\|TC-DB\|)(\w+)',defline).groups()[1]
    tcid = findtc.findall(defline)[0]
    family = '.'.join(tcid.split('.')[0:3]).upper()
    return (family,tcid,acc)

def sort_tcid(tcid):
    tcid = tcid.split('.')
    mytcid = [int(tcid[0]),str(tcid[1]),int(tcid[2]),int(tcid[3]),int(tcid[4])]
    return mytcid

def sem(mylist):
    ave,aves = average(mylist),[]
    for number in mylist:
        res = (float(number)-float(ave))
        res *= res
        aves.append(res)
    sem = sum(aves)/(float(len(aves)*(len(aves)-1)))
    return float(sqrt(sem))

class Emboss:

    def __init__(self):
        self.needle = None

    def extract_alignment(self,infile):
        aln = AlignIO.read(infile,'emboss')
        content = open(infile,'r').read()
        match = re.compile("                     ([|:. ]+)",re.DOTALL)
        mymatch = ''.join(match.findall(content))
        return (list(aln)[0].seq,mymatch,list(aln)[1].seq)

class GGSearchCommandline:

    def __init__(self):

        self.binary = 'ggsearch36'
        self.subject = None
        self.target = None
        self.outfile = tempfile.NamedTemporaryFile(delete=False).name

    def __call__(self):
        mycmd = '%s -s BL62 -m 8 -w 80 -f -8 -g -2 -b=3 -3 -k 500 %s %s'%(self.binary,self.subject,self.target)
        out = os.system(mycmd)
        print(out)


