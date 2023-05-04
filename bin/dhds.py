#!/usr/bin/env python
# This tool will perform a cross comparison of two FASTA files, but only using a global alignment
# And hydrophobic regions. This corrects for the blinding effect that local alignments have.

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline
from progressbar import ProgressBar,Bar,Percentage
from ProjectBio import Emboss
import hmmtop, os
import re,gsat
import optparse
import pickle
import subprocess
import tempfile
import templates
import scalagsat,copy

#os.environ['BIOV_DEBUG'] = 'True'

class extract:

    def __init__(self):
        self.subject_file = str
        self.target_file = str
        self.mysubject_file = None
        self.mytarget_file = None
        self.subjectname = None
        self.targetname = None
        self.assignments = {}
        self.alignments = {}
        self.gsatscores = []
        self.outdir = 'dhds_out'
        self.subject_loops = {}
        self.target_loops = {}
        self.flank = 5
        self.assign = 3
        self.pad = ('','XXX') # Hydrophillic-Hydrophillic & Hydrophibc-Hydrophillic Padding
        self.template = os.environ['DHDS_TEMPLATE']

        # Some ScalaGsat settings
        self.gapopen = 8
        self.gapextend = 2
        self.shuffles = 100
        self.quiet = True

    def verify(self):
        try:
            if os.path.exists(self.subject_file) and os.path.exists(self.target_file):
                return True
            return False
        except:
            return False

    def __call__(self):
        # First, run HMMTOP of course
        ht = hmmtop.tools()
        ht.add_library('subject',self.subject_file)
        ht.add_library('target',self.target_file)
        ht.scan_libraries()
        # Now, build the TMS sequence
        self.subject_loops = self.build_flanks(self.subject_file,ht.results['subject'])
        self.target_loops = self.build_flanks(self.target_file,ht.results['target'])
        # Glue the loops together and write file
        if os.path.exists(self.outdir) is False:
            os.mkdir(self.outdir)
        self.write_fastas()
        # Give names to each unit
        if self.subjectname is None:
            self.subjectname = self.mysubject_file.name.split('/')[-1]
        if self.targetname is None:
            self.targetname = self.mytarget_file.name.split('/')[-1]
        self.build_target_matrix()
        self.cross_matrix()
        self.construct_report()

    def build_flanks(self,fastafile,tms):
        fastas = SeqIO.to_dict(SeqIO.parse(fastafile,'fasta'))
        myloops = {}
        for symbol, tmss in tms.items():
            in_use=[]
            [in_use.extend(range(i[0],i[1]+1)) for i in tmss.values()]
            for tms, residues in tmss.items():
                maxlen = len(str(fastas[symbol].seq))
                ores = residues[:]
                for extend in range(self.flank):
                    left = residues[0] - 1
                    right = residues[1] + 1
                    # check validity
                    if left >= 1 and left not in in_use:
                        in_use.append(left)
                        residues[0] = left
                    if right <= maxlen and right not in in_use:
                        in_use.append(right)
                        residues[1] = right
                tmsonly = fastas[symbol].seq[ores[0]-1:ores[1]]
                loop = str(fastas[symbol].seq[residues[0]-1:residues[1]])
                replacement = '%s%s%s'%(self.pad[1],tmsonly,self.pad[1])
                loop = loop.replace(str(tmsonly),replacement)
                myloops.setdefault(symbol,[]).append(loop)
        return myloops

    def write_fastas(self):
        self.mysubject_file = open(self.outdir+'/subject.faa','wb+')
        self.mytarget_file = open(self.outdir+'/target.faa','wb+')
        for acc,seq in self.subject_loops.items():
            sequence = self.pad[0].join(seq)
            record = SeqRecord(Seq(sequence,IUPAC.protein),id=acc, description='dhds')
            SeqIO.write(record,self.mysubject_file,'fasta')
        for acc,seq in self.target_loops.items():
            sequence = self.pad[0].join(seq)
            record = SeqRecord(Seq(sequence,IUPAC.protein),id=acc,description='dhds')
            SeqIO.write(record,self.mytarget_file,'fasta')
        # Rewind & flush
        self.mysubject_file.flush()
        self.mytarget_file.flush()
        self.mysubject_file.seek(0)
        self.mytarget_file.seek(0)

    def exhaustive_search(self):
        # This will just run scalagsat on all of our files
        scala = scalagsat.compare()
        scala.subject = self.mysubject_file.name
        scala.target = self.mytarget_file.name
        scala.gapopen = self.gapopen
        scala.gapextend = self.gapextend
        scala.shuffle = self.shuffles
        scala.quiet = self.quiet
        scala.outfile = self.outdir+'/scala_gsat.txt'
        scala()

    def build_target_matrix(self):
        if self.mysubject_file.read().count('>') > self.mytarget_file.read().count('>'):
            # Ensures the subject file is the smaller one
            (self.mysubject_file,self.mytarget_file) = (self.mytarget_file,self.mysubject_file)
            (self.subjectname,self.targetname) = (self.targetname,self.subjectname)
        self.mysubject_file.seek(0)
        matrixfile = self.outdir+'/matrix'
        if os.path.exists(matrixfile):
            self.assignments = pickle.load(open(matrixfile))
            return
        count = self.mysubject_file.read().count('>')
        results = re.compile(r'# 2: (\w+).+?# Gaps:\s+\d+/\d+ \((\d+\.\d+)%\).+?# Score: (\d+\.\d+)',re.DOTALL)
        needle = NeedleCommandline()
        needle.gapopen = self.gapopen
        needle.gapextend = self.gapextend
        needle.outfile = 'stdout'
        needle.bsequence = self.mytarget_file.name
        mytargets = SeqIO.parse(self.mytarget_file.name,'fasta')
        mytargets = SeqIO.to_dict(mytargets)
        mysubjects = SeqIO.parse(self.mysubject_file.name,'fasta')
        print "Creating alignment matrix. Please wait..."
        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=count).start()
        for i,subject in enumerate(mysubjects):
            needle.asequence = 'asis:%s'%str(subject.seq)
            (stdout,stderr) = needle()
            gaps = results.findall(stdout)
            gaps.sort(key=lambda x:float(x[2]),reverse=True)
            self.assignments.setdefault(subject.id,[]).extend(gaps[0:self.assign])
            pbar.update(i+1)
        #pbar.finish()
        '''
        outfile = tempfile.NamedTemporaryFile(delete=False)
        mycmd = 'ggsearch36 -s BL62 -m 8 -w 80 -f -8 -g -2 -b=3 -3 -k 500 %s %s'%(self.mysubject_file.name,self.mytarget_file.name)
        handle = subprocess.Popen(mycmd,shell=True,stdout=subprocess.PIPE)
        (stdout,stderr) = handle.communicate()
        results = re.compile('#.+?(^[^#].+?)#',re.MULTILINE|re.DOTALL)
        res = results.findall(stdout+"#")
        rows = [i.split('\n') for i in res]
        myrows = []
        for row in rows:
            row.pop()
            row.sort(key=lambda x:float(x.split('\t')[-1]),reverse=True)
            myrows.append(row[0:self.assign])
        for i in myrows:
            [self.assignments.setdefault(sub,[]).append(tar) for sub,tar in [x.split('\t')[0:2] for x in i if bool(re.search('\t',x))]]
        '''
        pickle.dump(self.assignments,open(matrixfile,'w'))
        return

    def cross_matrix(self):
        scorefile = self.outdir+'/gsatscores'
        alignmentfile = self.outdir+'/alignments'
        if os.path.exists(scorefile):
            self.gsatscores = pickle.load(open(scorefile,'r'))
            self.alignments = pickle.load(open(alignmentfile,'r'))
            return
        (i,total) = (0,len(self.assignments) * self.assign)
        subjects = SeqIO.parse(self.mysubject_file.name,'fasta')
        targets = SeqIO.parse(self.mytarget_file.name,'fasta')
        mysubjects = SeqIO.to_dict(subjects)
        mytargets = SeqIO.to_dict(targets)
        gs = gsat.cmd()
        gs.shuffles = self.shuffles
        gs.gapopen = self.gapopen
        gs.gapextend = self.gapextend
        gs.verbose = False
        emboss = Emboss()
        print 'Performing cross comparison of selected targets...'
        pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=total).start()
        for subject, targets in self.assignments.items():
            gs.asequence = str(mysubjects[subject].seq)
            for target,gap,bit in targets:
            #for target in targets:
                gs.bsequence = str(mytargets[target].seq)
                gs()
                self.gsatscores.append((subject,target,gs.gaps,gs.zscore))
                aln=emboss.extract_alignment(gs.outfile.name)
                self.alignments.setdefault(subject+'-'+target,aln)
                pbar.update(i+1)
                i +=1
        self.gsatscores.sort(key=lambda x:x[-1], reverse=True)
        pickle.dump(self.alignments,open(alignmentfile,'wb+'))
        pickle.dump(self.gsatscores,open(scorefile,'wb+'))

    def construct_report(self):
        print 'Generating HTML Report...'
        template = open(self.template,'r').read()
        # Load resultrow & alignment modules
        resultrow = re.search(r'<%RESULTROW%>(.+)<\/\%RESULTROW%>',template,re.DOTALL)
        alignment = re.search(r'<%ALIGNMENT%>(.+)<\/\%ALIGNMENT%>',template,re.DOTALL)
        resultrow_html,alignment_html = (resultrow.groups()[0],alignment.groups()[0])
        # Generate rows first
        template = template.replace('%SUBJECTNAME%',self.subjectname)
        template = template.replace('%TARGETNAME%',self.targetname)
        myrows,myalignments = [],[]
        for subjectid, targetid, gaps, z in self.gsatscores:
            myrow = resultrow_html[:]
            myrow = myrow.replace('%SUBJECTID%',subjectid)
            myrow = myrow.replace('%TARGETID%',targetid)
            myrow = myrow.replace('%GAPS%',str(gaps))
            myrow = myrow.replace('%ZSCORE%',str(z))
            myrows.append(myrow)

            myalignment = alignment_html[:]
            (aseq,match,bseq) = self.alignments[subjectid+'-'+targetid]
            (atms,btms) = (self.annotate(aseq),self.annotate(bseq))
            myalignment = myalignment.replace('%SUBJECTID%',subjectid)
            myalignment = myalignment.replace('%TARGETID%',targetid)
            myalignment = myalignment.replace('%GAPS%',str(gaps))
            myalignment = myalignment.replace('%ZSCORE%',str(z))
            myalignment = myalignment.replace('%ATMS%',atms)
            myalignment = myalignment.replace('%ASEQUENCE%',str(aseq))
            myalignment = myalignment.replace('%MATCH%',match)
            myalignment = myalignment.replace('%BTMS%',btms)
            myalignment = myalignment.replace('%BSEQUENCE%',str(bseq))
            myalignments.append(myalignment)

        rows = ''.join(myrows)
        alns = ''.join(myalignments)
        template = re.sub('<%RESULTROW%>.+<\/%RESULTROW%>',rows,template,flags=re.DOTALL)
        template = re.sub('<%ALIGNMENT%>.+<\/%ALIGNMENT%>',alns,template,flags=re.DOTALL)
        final = open(self.outdir+'/results.html','wb+')
        final.write(template)
        return

    def annotate(self,seq):
        split = 'XXX|XX-+X|X-+XX'
        tms = re.split(split,str(seq))
        tmss = [v for i,v in enumerate(tms) if i%2 != 0]
        ranges = []
        for tms in tmss:
            match = re.search(tms,str(seq))
            ranges.append((match.start(),match.end()))
        bar = [' ' for i in range(len(str(seq)))]
        for count,tr in enumerate(ranges):
            (start,end) = tr
            first = True
            for i in range(start,end-(len(str(count+1))-1)):
                bar[i]='-' if first is False else str(count+1)
                first = False
        return ''.join(bar)

if __name__=='__main__':
    dhds = extract()
    desc = "DHDS - Distant Homolog Detection System :: Establish homology between distantly related FASTA files using only global alignments and hydrophobic regions.\
    This simplifies sequences to TMSs with padding & corrects for the blinding effect Smith-Waterman searches have when searching for TMS alignments. By Vamsee Reddy :: Part of the Bio-V Suite."
    opts = optparse.OptionParser(description=desc,version='1.0')
    opts.add_option('-s',dest='subject',type='str',help='Path to subject FASTA file')
    opts.add_option('-t',dest='target',type='str',help='Path to subject FASTA file')
    opts.add_option('-o',dest='outdir',type='str',default='dhds_out', help='Path to output directory')
    opts.add_option('--subject',dest='sname',type='str',default=None,help='Name of subject to display on result page')
    opts.add_option('--target',dest='tname',type='str',default=None,help='Name of target to display on result page')
    opts.add_option('--assign',dest='assign',type='int',default=3,help='Number of targets to assign each subject')
    (cli,args) = opts.parse_args()
    dhds.subject_file=cli.subject
    dhds.target_file=cli.target
    dhds.outdir=cli.outdir
    dhds.assign=cli.assign
    dhds.subjectname = cli.sname
    dhds.targetname = cli.tname
    if dhds.verify() is False:
        opts.print_help()
        quit()
    dhds()
