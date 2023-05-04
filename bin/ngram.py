#!/usr/bin/env python

# Generates Phylogenetic trees based on N-Grams

# [Program Flow]
# Input FASTA list, representing TC families.
# Preform BLASTP for each representative
# Generate N-Grams -> Merge into linear string of size [seek]
# Save N-Gram strings to FASTA. Preform Clustal Alignment
# Generate Phylogenetic Tree based on M.A.

# Designed by Vamsee Reddy

import sys
import re,os
import blast
#import tcdb
import shutil
from random import shuffle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class tree:

    def __init__(self):
        self.n = 3
        self.superid = int
        #(self.connection,self.cursor) = tcdb.connect()

    def ngram(self,fasta,size):
        mygram = {}
        sequence = str(fasta.seq)
        for i in range(len(sequence)-size+1):
            frame = i + size
            gram = sequence[i:frame]
            mygram.setdefault(gram,0)
            mygram[gram]+=1
        return mygram

    def define_superfamily(self):
        # First load familes within this superfamily
        q = 'SELECT family from family where super_id=%i'%self.superid
        self.cursor.execute(q)
        families = [i[0] for i in self.cursor.fetchall()]
        # Now load representatives for each family
        representatives = ["system.tcid LIKE '"+i+'.%.1'+"'" for i in families]
        select = "SELECT system.tcid,protein.acc,protein.sequence from system \
        inner join tc2acc on(system.tcid=tc2acc.tcid) \
        inner join protein on(tc2acc.acc=protein.acc) where %s" %(' OR '.join(representatives))
        self.cursor.execute(select)
        fasta = open("superfamily.faa",'wb')
        for row in self.cursor.fetchall():
            (tcid,acc,sequence) = row
            sequence = re.sub('[^A-Z]','',sequence)
            record = SeqRecord(Seq(sequence,IUPAC.protein),id=acc,name=tcid,description=tcid)
            SeqIO.write(record,fasta,'fasta')

    def get_homologs(self):
        if os.path.exists('homologs') is False:
            os.mkdir('homologs')
        b = blast.tools()
        superfamily = SeqIO.parse(open('superfamily.faa','r'),'fasta')
        for fasta in superfamily:
            b.gis = []
            b.blastp(str(fasta.seq))
            b.build_raw_fasta(fasta.description)
            shutil.copy(b.raw_fasta.name,'./homologs/'+fasta.id+'.faa')

    def build_nseq(self,file):
        mygrams = {}
        handle = SeqIO.parse(open("./homologs/"+file+'.faa','r'),'fasta')
        for i,fasta in enumerate(handle):
            ng = self.ngram(fasta,self.n).items()
            for seq,count in ng:
                mygrams.setdefault(seq,0)
                mygrams[seq]+=count
        count = i+1
        return dict([(k,float(float(v)/float(count))) for k,v in mygrams.items()])

    def calculate_distance(self,a,b):
        mya = a.items()
        myb= b.items()
        mya.sort(key=lambda x:x[1],reverse=True)
        myb.sort(key=lambda x:x[1],reverse=True)
        #besta = mya[0:100]
        #bestb = myb[0:100]
        shared = list(set(dict(mya).keys()) & set(dict(myb).keys()))
        shuffle(shared)
        shared = shared[0:5]
        if len(shared) is 0:
            return 9
        distance = float(sum([abs(a[i]-b[i]) for i in shared])/float(100))
        return distance

    def generate_matrix(self):
        saved ={}
        superfamily = SeqIO.parse("superfamily.faa",'fasta')
        superfamily = SeqIO.to_dict(superfamily)
        fastas = os.listdir('./homologs')
        myids = [i.split('.')[0] for i in fastas if bool(re.search(r'\.faa',i))]
        myids.sort()
        matrix = open('infile','ab')
        matrix.write("\t"+str(len(myids))+"\n")
        for vd in myids:
            dist = []
            a=self.build_nseq(vd)
            for hd in myids:
                b=self.build_nseq(hd)
                label = [vd,hd]
                label.sort()
                label="-".join(label)
                myscore= saved.setdefault(label,self.calculate_distance(a,b))
                score = '%0.6f'%myscore
                dist.append(score)
            line = '%s    %s\n'%(vd,"   ".join(dist))
            matrix.write(line)
        matrix.write("\n")

    def format_tree(self,file):
        supertree = SeqIO.parse('superfamily.faa','fasta')
        tree = open(file,'r').read()
        for fasta in supertree:
            tree = re.sub(fasta.id,fasta.description.split(' ')[-1],tree)
        new = open(file,'wb')
        new.write(tree)
        return



if __name__=='__main__':

    mytree = tree()
    mytree.superid = 6
    mytree.format_tree('outtree')
    quit()
    for i in range(2):
        mytree.generate_matrix()
        print i




