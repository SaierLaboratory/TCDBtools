#!/usr/bin/env python
#Changes 10/03/2022 made by Hunter Cheng for program compatability in python3

import tempfile
import subprocess
import re,os
import ssearch
from Bio import SeqIO

class compare:

    def __init__(self,subject=None,target=None):
        self.max = 3
        self.results = []
        self.keys={}
        self.shuffle = 300
        self.subject = subject
        self.target = target
        self.verbose = False
        self.swap = False

    def __call__(self):
        # Load our Fasta Files
        subject = open(self.subject,"r")
        target = open(self.target,"r")
        subjects = SeqIO.parse(subject,'fasta')
        targets = SeqIO.parse(target,'fasta')
        subjects = SeqIO.to_dict(subjects)
        targets = SeqIO.to_dict(targets)
        subject.close()
        target.close()
        # Assign variables, we need our subject to be the shorter one
        if len(subjects) > len(targets):
            self.subjects = targets
            self.subject_file = target.name
            self.targets = subjects
            self.target_file = subject.name
            self.swap = True
        else:
            self.subjects = subjects
            self.subject_file = subject.name
            self.targets = targets
            self.target_file = target.name
        # Preform preliminary SSearch36 Scan
        self.primary_scan()
        self.extract()
        self.secondary_scan()

    def checkvars(self):
        try:
            if ( self.max < 1 or self.shuffle < 1
                 or os.path.exists(self.subject) is False
                 or os.path.exists(self.target) is False ):
                return False
            else:
                return True
        except:
            return False

    def primary_scan(self):
        cmd = "ssearch36 -s BL62 -m 8C -w 80 -f -8 -g -2 -3 -k 500 %s %s" %(self.subject_file,self.target_file)
        #universal_newlines = True was added for compatability in python3. It allows for backwards compatability
        handle = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, universal_newlines = True)
        results = handle.communicate()[0]
        self.primary_results = results
        return

    def extract(self):
        pattern = '\r|\n([A-Za-z0-9_.|-]+)\s+([A-Za-z0-9_.|-]+)\s+\d+'
        results = re.findall(pattern,self.primary_results)
        # Generates a Dictionary, with subject as key
        #  and value is a list of targets
        [self.keys.setdefault(k,[]).append(v) for k,v in results]
        return

    def secondary_scan(self):
        total = len(self.keys.keys())
        # enumerate is used in case I wanted to include a progress bar
        for i in enumerate(self.keys.items()):
            subject, targets = i[1]
            #In python 3, we need to specify write mode and make sure temp file is deleted as soon as closes
            subject_file = tempfile.NamedTemporaryFile(mode = 'wt+', delete = True)
            target_file = tempfile.NamedTemporaryFile(mode = 'wt+', delete = True)
            subject_file.write(">"+self.subjects[subject].id+"\n"+str(self.subjects[subject].seq))
            for target in targets[0:self.max]:
                target_file.write(">"+self.targets[target].id+"\n"+str(self.targets[target].seq)+"\n")
            subject_file.flush()
            target_file.flush()
            subject_file.seek(0)
            target_file.seek(0)
            myssearch=ssearch.Compare(subject_file.name,target_file.name,self.shuffle,self.verbose)
            myssearch()
            subject_file.close()
            target_file.close()
            for score in myssearch.results:
                self.results.append(score)
            if self.verbose:
                print ("Completed %i of %i" %(i[0]+1,total))
        self.results.sort(key=lambda x:float(x['z_score']), reverse=True)
        return

if __name__ =='__main__':
    from optparse import OptionParser
    desc = "TSSearch is a heuristic Smith & Waterman search tool that will rapidly find close & distant homologs between two FASTA files."
    opts = OptionParser(description=desc)
    opts.add_option('-s',
                    action='store',
                    type='string',
                    dest='subject',
                    help="Path to your subject fasta file"
    )
    opts.add_option('-t',
                    action='store',
                    type='string',
                    dest='target',
                    help="Path to your target fasta file"
    )
    opts.add_option('-o',
                    action='store',
                    type='string',
                    dest='output',
                    help="Output file to create"
    )
    opts.add_option('-m',
                    action='store',
                    type='int',
                    default=3,
                    dest='max_targets',
                    help="Targets per subject (Optional, Default: 3)"
    )
    opts.add_option('-r',
                    action='store',
                    type='int',
                    dest='random',
                    default=300,
                    help="Number of times to shuffle each sequence (Optional, Default: 300)"
    )
    (cli,args) = opts.parse_args()
    tss = compare()
    tss.subject = cli.subject
    tss.target = cli.target
    tss.output = cli.output
    tss.max = cli.max_targets
    tss.shuffle = cli.random
    tss.verbose = True
    if tss.checkvars():
        tss()
        output = open(cli.output,'wb')
        ss = ssearch.Compare()
        for line in tss.results:
            l = ss.print_line(line)
            output.write(l+"\n")
        print ("File Written to %s" %output.name)
    else:
        opts.print_help()

