#!/usr/bin/env python
#Last edited by Hunter Cheng 10/04/2022
#Changes were facilitated to allow program functionality in python3
'''
==========================================================
This program uses SSearch36 for mass shuffles,
as this happens to be much quicker than
EMBOSS::Water When calculating a z-score.
==========================================================
This Class provides a simple 
listed dictionary containing all your results.
import, execute ssearch.Compare(SUBJECT.FASTA,TARGET.FASTA).
Results are stored in 'results' variable with self explanitory 
dictionary keys. Keys in results are :
{subject_symbol, target_symbol, sw_score, z_score,
subject_start, subject_end, target_start, target_end}
Coded by Vamsee Reddy <symphony.dev@gmail.com>
==========================================================
'''

import numpy
import tempfile
import subprocess
import re, os
import textwrap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC

class Compare:

    def __init__(self,subject=None,target=None,shuffle=None,verbose=False):

        # Set some defaults
        self.subject = subject
        self.target = target
        self.shuffle = 200
        self.verbose=verbose
        self.results=[]
        return

    def __call__(self):
        # Parse subjects
        subject_file = open(self.subject,"r")
        self.subject_fastas = list(SeqIO.parse(subject_file,'fasta'))
        # Parse targets
        target_file = open(self.target,"r")
        self.target_fastas = list(SeqIO.parse(target_file,'fasta'))
        # Iterate through subjects, watering each target
        subject_file.close()
        target_file.close()
        for subject in self.subject_fastas:
            for target in self.target_fastas:
                values = self.water(subject, target)
                extract = self.extract(values)
                z_scorez = self.calculate_z_score(extract)
                (target_symbol,sw_score,z_score,subject_start,subject_end,target_start,target_end) = [i for i in extract if i[0]==target.name][0]
                row = {'target_symbol':target_symbol,'subject_symbol':subject.name,'sw_score':sw_score,'z_score':z_scorez,
                'subject_start':int(subject_start),'subject_end':int(subject_end),'target_start':int(target_start),'target_end':int(target_end)}
                if self.verbose is True:
                    print (self.print_line(row))
                self.results.append(row)
        self.results.sort(key=lambda x:float(x['z_score']), reverse=True)

    def checkvars(self):
        try:
            if self.shuffle < 1 or os.path.exists(self.subject) is False or os.path.exists(self.target) is False:
                return False
            else:
                return True
        except:
            return False

    def print_line(self,row):
        line = "%s (%i-%i)\t%s (%i-%i)\t%s\t%.2f" %(row['subject_symbol'],row['subject_start'],row['subject_end'],row['target_symbol'],row['target_start'],\
        row['target_end'],row['sw_score'],float(row['z_score']))
        return line

    def water(self, subject, target):
        # Create temporary shuffled target file
        #In python3 we have to specify write mode and make sure temp file is deleted as soon as it closes
        shuffled_target_file = tempfile.NamedTemporaryFile(mode="wt+",delete = True)
        # Write original FASTA Target
        shuffled_target_file.write('>'+target.name+"\n")
        shuffled_target_file.write('\n'.join(textwrap.wrap(str(target.seq),70))+'\n')
        # Write x shuffled targets to same file
        fastas=[]
        for i in range(self.shuffle):
            random = list(target)
            numpy.random.shuffle(random)
            record = SeqRecord(Seq(''.join(random)),name=target.name,id=target.name,description=target.name)
            SeqIO.write(record,shuffled_target_file,'fasta')
        # Write Subject FASTA to temporary file
        #In python3 we have to specify write mode and make sure temp file is deleted as soon as it closes
        subject_file = tempfile.NamedTemporaryFile(mode="wt+",delete = True)
        record = SeqRecord(Seq(''.join(subject.seq)),name=subject.name,id=subject.name,description=subject.name)
        SeqIO.write(record,subject_file,'fasta')
        # Preform SSearch36 Comparision of two new files
        subject_file.flush()
        shuffled_target_file.flush()
        # Create a temporary file to push our OB to
        command = 'ssearch36 -p -q -w 80 -W 0 -m 0 -z 11 -k 1000 -f -8 -g -2 -E 1000 -s BL62 "%s" "%s"' %(subject_file.name,shuffled_target_file.name)
        #For python3 compatability, universalnewlines = True must be added to allow the search to yield in string form
        ssearch = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, universal_newlines = True)
        stdout,stderr = ssearch.communicate()
        subject_file.close()
        shuffled_target_file.close()
        return stdout

    def extract(self, res):
        # Captures tuple of 7 entries:
        # (Symbol, SW-Score, Z-Score, Subjct-Start, Subject-End,
        #  Target-Start, Target-End), very cool.
        pattern = '>>([A-Za-z0-9_.|-]+)\s+\(\w+ aa\)\s+s-w\s+opt:\s+(\w+)\s+Z-score:\s+([.0-9]+)[^\n]+\s+[^\n]+\((\w+)-(\w+):(\w+)-(\w+)\)'
        results = re.findall(pattern,res)
        return results

    def calculate_z_score(self,results):
        if self.shuffle is False:
            return False
        scores = [int(s[1]) for s in results]
        std = numpy.std(scores[1:])
        mean = numpy.mean(scores[1:])
        z_score = (scores[0] - mean)/std
        return z_score

if __name__=='__main__':
    from optparse import OptionParser
    from os import path
    desc = '''
    SSearch will compare every sequence in a subject file to every
    sequence in a target file by using a Smith Waterman search.
    Results are returned after each pair is shuffled #R number of
    times, and a Z-Score is returned to determine the quality of
    the alignment.
    Contact: Vamsee Reddy <vamsee30@gmail.com>
    '''
    desc = " ".join(desc.split())
    opts = OptionParser(description=desc)
    opts.add_option('-s',
                    action='store',
                    type='string',
                    dest='subject',
                    help="Path to your subject fasta file")
    opts.add_option('-t',
                    action='store',
                    type='string',
                    dest='target',
                    help="Path to your target fasta file")
    opts.add_option('-r',
                    action='store',
                    type='int',
                    dest='random',
                    default=300,
                    help="Number of times to shuffle each sequence (Optional. Default: 300)")
    opts.add_option('-o',action='store',type='string',dest='output',default='ss_out.txt', help="Filepath To write output")
    (cli,args)= opts.parse_args()
    ss = Compare(cli.subject,cli.target,cli.random,True)
    if ss.checkvars() is False:
        opts.print_help()
        quit()
    ss()
    result = open(cli.output,'wb+')
    for line in ss.results:
        l = ss.print_line(line)
        result.write(l+"\n")
    print ("File Written to %s" %result.name)
