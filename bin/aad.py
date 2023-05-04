#!/usr/bin/env python
import re,sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import generic_protein
from Bio import SeqIO
sys.path.insert(1,'../biov/')
import hmmtop
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pylab
from math import sin,cos,radians,sqrt
from matplotlib.patches import Rectangle
from numpy import std
from math import exp
from decimal import *
from numpy import median

class Average:
    
    def __init__(self):
        self.aln_file = './staph.aln'
        self.master_seq_file = './hmmtop.in'
        self.master_seq = []
        self.hmmtop = None
        self.aaindex={'G':(-0.400,0.48),'I':(4.500,1.38),'S':(-0.800,-0.18),'Q':(-3.500,-0.85),'E':(-3.500,-0.74),\
              'A':(1.800,0.62),'M':(1.900,0.64),'T':(-0.700,-0.05),'Y':(-1.300,0.26),'H':(-3.200,-0.4),\
              'V':(4.200,1.08),'F':(2.800,1.19),'C':(2.500,0.29),'W':(-0.900,0.81),'K':(-3.900,-1.5),\
              'L':(3.800,1.06),'P':(-1.600,0.12),'N':(-3.500,-0.78),'D':(-3.500,-0.90),'R':(-4.500,-2.53),'X':(0,0)}
        self.avehydro = [] # average hydropathy for each position
        self.aveamphi = [] # average amphipathicity for each position
        self.stds = [] # Position specific standard deviation.
        self.smooth_hydro = [] # average hydropathy smoothed over window
        self.smooth_amphi = [] # average amphipathicity smoothed over window
        self.smooth_std = [] # averaged '' standard deviations ''
        self.tms_ratio = {}
        self.tms = None
        
        self.groups = {'Aliphatic':["L","I","V","A","G"],\
                       'Aromatic':['F','Y','W'],\
                       'Positive':['R','K'],\
                       'Helix-Breaking':['P','G'],\
                       'Negative':['D','Q'],\
                       'Amphipathic':['Y','H','W'],
                       'Semipolar':['Q','T','S','C','N','M'],\
                       'Branched':['I','L','V']}
        self.select_groups=[('Aliphatic','Orange'),('Helix-Breaking','Purple'),('Aromatic','Green')] # I will process these ones
        self.smooth_groups = {}
    
    def build_master_sequence(self):
        pattern = re.compile(r'(\S+)\s+([A-Z,-]+)')
        master_seq = {}
        with open(self.aln_file,'r') as h:
            contents = h.read()
            contents = re.sub(r'\*|\.|:','',contents)
            for i in pattern.findall(contents):
                if i[0].lower() == 'clustal':
                    continue
                master_seq.setdefault(i[0],[]).append(i[1])
        for key,seqs in master_seq.items():
            seqs = "".join(seqs)
            self.master_seq.append(seqs)
        return True
        
    def write_fastas(self):
        records = []
        for key,seq in enumerate(self.master_seq):
            record = SeqRecord(Seq(seq),id=str(key),description='')
            records.append(record)
        SeqIO.write(records,"hmmtop.in",'fasta')
        
    def run_hmmtop(self):
        ht = hmmtop.tools()
        ht.add_library('res',self.master_seq_file)
        ht.scan_libraries()
        self.tms = ht.results
        self.hmmtop = []
        gap_seqs = {}
        for i in SeqIO.parse('hmmtop.in','fasta'):
            gap_seqs[i.id]=str(i.seq)
        for i,j in ht.results['res'].items():
            for k in j.values():
                korrect = [self.get_pos(k[0],gap_seqs[i]),self.get_pos(k[1],gap_seqs[i])]
                point = range(korrect[0],korrect[1]+1)
                self.hmmtop.append(median(point))
        graph = {}
        for i in range(1,len(self.master_seq[0])+1):
            graph[i]=0
        for tms in self.hmmtop:
            graph[i] +=1
        for key, hits in graph.items():
            ratio = Decimal(hits)/Decimal(len(self.master_seq))
            self.tms_ratio[key]=ratio
        return True
        
    def get_pos(self,pos,gaps): # Corrects an aa residue position to fit a sequence with gaps
        findme = 0
        for i in range(1,len(gaps)+1):
            if gaps[i-1] != '-':
                #This is a letter
                findme += 1
                if findme == pos:
                    return i
        
    def average_seq(self):
        num_seq = len(self.master_seq)
        seq_len = len(self.master_seq[0])
        haves = []
        aaves = []
        # Get the average hydropathy for each position
        for i in range(seq_len): #iterating over the positions
            hsum,asum,total = 0,0,0
            getstd = []
            for seq in self.master_seq:
                aa = seq[i]
                if aa != '-':
                    getstd.append(self.aaindex[aa][0])
                    total +=1
                    hsum +=self.aaindex[aa][0]
                    asum +=self.aaindex[aa][1]
                else:
                    getstd.append(0)
            simscore = 0-std(getstd)
            self.stds.append(simscore)
            haves.append(hsum/total)
            aaves.append(asum/total)
        assert len(aaves) == seq_len
        self.avehydro,self.aveamphi = haves,aaves
        return True
        
    def calculate_hydropathy(self,window=20):
        coordinates = []
        for i in range(len(self.avehydro)-window-1):
            total = 0
            sim = 0
            for j in range(window):
                total +=self.avehydro[i+j]
                sim +=self.stds[i+j]
            coordinates.append(total/window)
            self.smooth_std.append(sim/window)
        self.smooth_hydro = coordinates
        return True
    
    def calculate_amphipathicity(self,window=20,angle=100):
        limit = len(self.aveamphi) - window-1
        amphi = []
        for i in range(limit):
          total = 0
          sumsin = 0
          sumcos = 0
          tangle = angle

          for j in range(window):
            hydro = self.aveamphi[i+j]
            sumsin += hydro * sin(radians(tangle))
            sumcos += hydro * cos(radians(tangle))
            tangle = (tangle+angle)%360

          sumsin = sumsin * sumsin
          sumcos = sumcos * sumcos
          total = (sqrt(sumsin + sumcos)/window)*(5/2) # Correction factor!
          amphi.append(total)
        self.smooth_amphi = amphi
        return amphi
        
    def process_select_groups(self,window=20):
        fstring = {}
        for key in [x[0] for x in self.select_groups]:
            # Create a frequency string
            for i in range(len(self.master_seq[0])):
                # iterating through every position
                # analyze position if of every sequence in master_seq
                total = 0
                for j in self.master_seq:
                    residue = j[i]
                    # is residue in our key?
                    if residue in self.groups[key]:
                        total += 1
                fstring.setdefault(key,[]).append(total)
        # Smooth out these values over a window
        smoothf = {}
        tempsum = {}
        for i in range(len(self.avehydro)-window-1):
            # create a dictionary
            for key in [x[0] for x in self.select_groups]:
                tempsum[key]=0
            for j in range(window):
                # i+j will yield the correct position
                # must smooth out ech selected group
                for key in [x[0] for x in self.select_groups]:
                    tempsum[key] += fstring[key][i+j] # append this
            for key in [x[0] for x in self.select_groups]:
                precise = (Decimal(tempsum[key])/Decimal(len(self.master_seq)*window))*4+2
                smoothf.setdefault(key,[]).append(precise)
        self.smooth_groups = smoothf
                

                
        
    def generate_graph(self):
        x_data = range(0, len(self.smooth_hydro))
        diff=(len(self.master_seq[0])-len(self.smooth_hydro))/2
        x1_data = range(0,len(self.smooth_groups.items()[0][-1]))
        x2_data = range(0,len(self.master_seq[0]))
        plt.figure()
        plt.axhline(y=0, color='black')
        plt.ylim(-3, 6)
        plt.xlim(right=len(self.master_seq[0]))
        plt.plot(x_data, self.smooth_hydro, linewidth=1.0, label="hydrophobicity", color='r')
        plt.plot(x_data, self.smooth_amphi, linewidth=1.0, label="amphipathicity", color='g')
        for group,color in self.select_groups:
            plt.plot(x1_data, self.smooth_groups[group], linewidth=1.0, label=group, color=color)
        plt.plot(x_data, self.smooth_std, linewidth=1.0, label="similarity", linestyle='solid', color='grey')
        
        for pos in self.hmmtop:
            plt.axvline(x=pos-1, ymin=-2, ymax = 0.1, linewidth=1, color='black',alpha=0.2)
            
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
                      ncol=3, fancybox=True, shadow=True)
        plt.xlabel("Residue Number")
        plt.ylabel("Value")
        width = (0.0265)*len(self.master_seq[0]) if len(self.master_seq[0]) > 600 else 15
        plt.grid('on')
        plt.savefig('foo.png')
        
        


if __name__=='__main__':     
    A = Average()
    A.build_master_sequence()
    A.write_fastas()
    A.average_seq()
    A.run_hmmtop()
    A.calculate_hydropathy()
    A.calculate_amphipathicity()
    A.process_select_groups()
    A.generate_graph()
