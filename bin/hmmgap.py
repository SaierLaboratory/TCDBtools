#!/usr/bin/env python

import tempfile
import re
import subprocess
import hmmtop
from Bio import SeqIO

class annotate:

    def __init__(self):
        self.global_tms_ranges = {}
        self.local_tms_ranges = {}
        self.annotate = []
        self.hmmtop = False
        self.annotation = False

    def __call__(self,seq,gap):
        self.annotate=[]
        self.gap = str(gap)
        try:
            self.seq=SeqIO.read(seq,'fasta')
        except:
            self.seq = seq
        try:
            if self.hmmtop is False:
                hmt = hmmtop.tools()
                hmt.add_library("gap",seq)
                hmt.scan_libraries()
                self.global_tms_ranges=hmt.results['gap'][self.seq.id]
            else:
                self.global_tms_ranges=self.hmmtop[self.seq.id]
        except:
            return None
        self.parse_gap()
        self.build_label()
        return self.annotation


    def get_tms_number(self,pos):
        for k,v in list(self.global_tms_ranges.items()):
            if pos in range(v[0],v[1]):
                return k
        return False

    '''
    Within the parse_gap function, we have replaced any '*' character with a '\*',
    so that the regular expression compiles without interpretting the '*'.
    - Vasu Pranav Sai Iddamsetty
    '''

    def parse_gap(self):

        if '*' in str(self.gap):

            self.gap = str(self.gap).replace('*','\*')

        match=re.search(re.sub('-','',str(self.gap)),str(self.seq.seq))
        # Annotate
        x=match.start()-1
        for residue in self.gap:
            if residue != "-":
                self.annotate.append(self.get_tms_number(x))
                x += 1
            else:
                self.annotate.append(".")

    def build_label(self):
        string =[]
        tms = 0
        for i in self.annotate:
            if i != "." and i != False:
                if i != tms:
                    string.append(str(i))
                    ignore = len(str(i))-1
                    tms = i
                else:
                    if ignore != 0 and tms > 0:
                        ignore -= 1
                        continue
                    string.append("+")
            else:
                string.append("_")
        #print string
        self.annotation = ''.join(string)

if __name__=='__main__':
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    import sys
    if len(sys.argv) < 2:
        seq = input("Enter Sequence: ")
        gap_i= input("Enter gaps: ")
    else:
        seq = sys.argv[1]
        try:
            gap_i = sys.argv[2]
        except:
            pass
    gap = gap_i if len(gap_i)>0 else seq
    hmg = annotate()
    record = SeqRecord(Seq(seq,IUPAC.protein),id="hmg", name="hmg")
    seqfile = tempfile.NamedTemporaryFile()
    SeqIO.write(record,seqfile,'fasta')
    seqfile.seek(0)
    hmg(seqfile.name,gap)
    print(hmg.annotation)
    print(seq)
    import tkinter



