#!/usr/bin/env python

import gsat
import sys
from Bio import SeqIO
#from progressbar import ProgressBar,Bar,Percentage

class compare:

    def __init__(self):
        self.subject = False
        self.target = False
        self.shuffle = 50
        self.gapopen = 8
        self.gapextend =2
        self.quiet = True
        self.defline = False
        self.outfile = 'scalagsat_results.txt'
        self.results = []
        #self.progress = False # Some day...

    def __call__(self):
        subjects = SeqIO.parse(open(self.subject),'fasta')
        #totalsubjects = SeqIO.parse(open(self.subject),'fasta')
        #totalsubjects = SeqIO.to_dict(totalsubjects)
        #totalsubjects = len(totalsubjects)
        targets = list(SeqIO.parse(open(self.target),'fasta'))
        #numjobs = len(targets)*totalsubjects
        if self.outfile is not None:
            results = open(self.outfile,'wb+')
        i = 0
        #if self.progress:
        #    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=numjobs).start()
        for subject in subjects:
            subject_seq = str(subject.seq)
            for target in targets:
                target_seq=str(target.seq)
                try:
                    gs = gsat.compare(subject_seq,target_seq,'global',self.shuffle,self.gapopen,self.gapextend,False)
                except:
                    print 'Failed to Run:'
                    print subject.id, target.id
                    print subject_seq
                    print target_seq        
                    quit()            
                row = [subject.id,subject.description,target.id,target.description,gs.gaps,gs.zscore]
                gs.zscore = round(gs.zscore)
                if gs.zscore > 0:
                    strow =(subject.id,target.id,str(gs.gaps),str(gs.zscore))
                    if self.quiet is False:
                        print "\t".join(strow)
                    if self.outfile is not None:
                        if self.defline is False:
                            row.pop(1)
                            row.pop(2)
                        results.write("\t".join([str(i) for i in row])+"\n")
                        results.flush()
                    self.results.append(row)
                    #if self.progress:
                    #    pbar.update(i+1)
                    i += 1
        if self.outfile is not None:            
            results.close()
        self.results.sort(key=lambda x:(100-float(x[-2]),float(x[-1])),reverse=True)
        if self.outfile is not None:
            results = open(self.outfile,'wb+')
            for row in self.results:
                results.write("\t".join([str(i) for i in row])+"\n")
        return True

if __name__=='__main__':

    subject = sys.argv[1]
    target = sys.argv[2]
    outfile = sys.argv[3]

    scala=compare()
    scala.subject=subject
    scala.target=target
    scala.outfile=outfile

    scala()
