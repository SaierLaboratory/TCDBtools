#!/usr/bin/env python
import sys,os,copy
#import random
import re
import tempfile
import optparse
import textwrap
import numpy
import settings
import logging
from Bio import AlignIO
from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import generic_protein

class compare:

    def __init__(self, querys, subjects,
                 alignment='global',
                 shuffles=500,
                 gapopen=10,
                 gapextend=2,
                 verbose=False):
       
        #In python3 .encode("utf-8") for the write function to work in python3
        #self.init_logging()
        self.AA = querys.encode('utf-8')
        self.BB = subjects.encode('utf-8')
        qfile = tempfile.NamedTemporaryFile()
        qfile.write(b">A_Sequence\n"+self.AA)
        
        #set the reference point at the beginning of the file 
        qfile.seek(0)
        query=qfile.name
        sfile =tempfile.NamedTemporaryFile()
        sfile.write(b">B_Sequence\n"+self.BB)
        sfile.seek(0)
        subject=sfile.name
        self.query = query
        self.subject = subject
        self.shuffles = shuffles
        watercmd = "needle" if 'needle' not in os.environ else os.environ['needle']
        self.align = NeedleCommandline(watercmd)
        self.align.asequence = query     
        self.align.gapopen = gapopen
        self.align.gapextend = gapextend


        self.outfile = tempfile.NamedTemporaryFile()
        self.verbose = verbose

        # initialize default values for the scores and std deviation
        self.stdev = 0
        self.zscore = 0
        self.zscorep = 0
        self.average = 0
        self.on_score = []
        self.scores = []

        # initialize strings raw alignment
        self.subject_aln = str()
        self.target_aln = str()

        # initialize tempfiles
        self.shuffled_seq = tempfile.NamedTemporaryFile(mode='w+t') 
        self.shuffled_comp = tempfile.NamedTemporaryFile(mode='w+t')
        Seq_fasta2 = SeqIO.read(self.subject, "fasta")
        self.shuffle(Seq_fasta2, self.shuffles)
        self.compare()
        self.results()
        self.cleanup()
        self.grab_alignments()
        sfile.close()
        qfile.close()

    def init_logging(self):
        logger = logging.getLogger('GSAT')
        hdlr = logging.FileHandler('/var/tmp/GSAT.log')
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr) 
        logger.setLevel(logging.WARNING)
        self.logger = logger

    def shuffle(self, subject, shuffles):
        '''
        Shuffles a given Seq object a specified number of times,
        then outputs each shuffle to a file for comparisons using
        EMBOSS needle or water.
        '''

        sequence = list(subject)
        for i in range(shuffles):
            numpy.random.shuffle(sequence)
            # creates a new sequence record and writes it in fasta
            # format to /tmp/a
            #Removed bioalphabet.generic_protein second argument
            shuffled_sequence = SeqRecord(Seq("".join(sequence)), id="tmp")
            self.shuffled_seq.write(shuffled_sequence.format("fasta"))

    def compare(self):
        '''
        Pulls the shuffled sequence tempfile and the shuffled comparison
        tempfile and does a iterative comparison of the query sequence
        against each of the shuffled sequences, adds all the scores to an
        array, and computes the standard deviation of that array. Then
        compares original query vs original subject and saves that comparison
        score and alignment to an outfile, if specified, or to the temp file.
        '''
        self.align.bsequence = self.shuffled_seq.name
        self.shuffled_seq.flush()
        self.align.outfile = self.shuffled_comp.name
        stdout, stderr = self.align()
        self.align.bsequence = self.subject
        self.align.outfile = self.outfile.name
        stdout, stderr = self.align()
        #rewind the files
        self.shuffled_comp.seek(0)

    def z_score(self):
        '''
        Computes the standard (z) score assuming on_score, scores, and
        stdev are all initiated and populated. could/should be
        rewritten to take arguments.
        '''

        zscore = (self.on_score[0]-self.scores.mean())/self.stdev
        self.zscorep = round(zscore,3) # precise
        self.zscore = round(zscore)


    def find_scores(self, targetfile, scores):
        '''
        Extracts and collects the scores from a given file with one or
        more needle (or equivalent) comparisons; appends to given
        score array and returns it. This method modifies the arguments
        passed to it; specifically it populates the score[] array
        argument.

        targetfile: open handle to a given file
        scores: array to use to collect scores !This is modified
        '''

        # Find the number after 'Score:'
        regex = re.compile('Score: ([0-9,\.]{1,5})',flags = re.ASCII)
        for line in targetfile:
            line = str(line)
            tmps = regex.findall(line)
            if tmps != []: #success! add it to the scores array
                try: # but first convert it to an int or float 
                    tmps = int(tmps[0]) 
                except ValueError: #i.e., it's not an int, so try float
                    tmps = float(tmps[0])                    
                scores.append(tmps)
        return scores

    def results(self):
        '''
        Collates and (optionally) outputs the results.
        '''

        self.find_scores(self.shuffled_comp, self.scores)
        self.find_scores(self.outfile, self.on_score)

        # Calculate the standard deviation of the scores array
        self.scores = numpy.array(self.scores)
        self.stdev = numpy.std(self.scores)

        # Only call z_scores after stdev, scores, and on_score have
        # been populated
        self.z_score()
        
        self.average=numpy.average(self.scores)
        if self.verbose:
            self.outfile.seek(0)
            
            #RegEx to extract the alignment section of the output of running
            #the needle command.
            find = re.compile("#=+(.+)#-+",re.DOTALL)
            
            #In python3 .decode("utf-8") is necessary to interprte '\n' as new lines,
            #otherwise the scaped characters "\n" are printed to screen
            out=find.search(self.outfile.read().decode("utf-8")).groups()[0]
            alignment = re.compile("#=+\n(.+)#-+",re.DOTALL)
            alnresults = alignment.search(out)
            aln = alnresults.group()
            alnfixed = aln[:]
            alnfixed= alnfixed.replace('.',' ')
            print(out.replace(aln,alnfixed))
            print("============ FINISHED =============\n")
            print("Average Quality (AQ)\t%0.2f +/- %0.2f\n"%(self.average, self.stdev))
            print("Standard score (Z):\t{0}\n".format(self.zscore))
            print("Precise score (Z):\t{0}\n".format(self.zscorep))
            print("Shuffles:\t{0}\n".format(self.shuffles))

    def cleanup(self):
        self.shuffled_seq.close()
        self.shuffled_comp.close()

    def grab_alignments(self):
        self.outfile.seek(0)
        gaps = re.compile('#\s+Gaps:\s+\d+/\d+\s+\((\s+)?(\d+\.\d+)%(\s+)?\)')
        pattern = '[0-9]{1,3}\s+([-,A-Z]{1,70})\s+[0-9]{1,3}'
        gsatout=self.outfile.read().decode("utf-8")
        self.gaps = float(gaps.search(gsatout).groups()[1])
        find=re.findall(pattern,gsatout)
        for line in enumerate(find):
            if line[0]%2 == 0:
                self.subject_aln += line[1]
            else:
                self.target_aln += line[1]
        self.subject_aln = self.subject_aln.replace("-",".")
        self.target_aln = self.target_aln.replace("-",".")
        return

class cmd:
    def __init__(self):
        ''' This is a CMD wrapper for formal invocations.
            This is built as a wrapper to avoid
            Obstruction of current implementations of GSAT
        '''
        self.asequence = str
        self.bsequence = str
        self.gapopen = 8
        self.gapextend = 2
        self.shuffles = 500
        self.verbose = False
        self.outfile = None
        self.zscore = int
        self.zscorep = float
        self.gaps = float
        self.average = float

    def __call__(self):
        gs = compare(self.asequence,self.bsequence,'global',
                     self.shuffles,self.gapopen,
                     self.gapextend,self.verbose)
        (self.outfile,self.zscore,self.zscorep,self.average,self.gaps) = (gs.outfile,gs.zscore,gs.zscorep,gs.average,gs.gaps)
        return


##### here this thing runs:
if __name__=='__main__':
    try:
        asequence = sys.argv[1]
        bsequence = sys.argv[2]
        try:
            shuffles = sys.argv[3]
        except:
            shuffles = 500
    except:
        asequence = input("Sequence A: ")
        bsequence = input("Sequence B: ")
        shuffles = input("Shuffles: ")
    
    gap = compare(asequence,bsequence,'global', int(shuffles), 8,2, True)
