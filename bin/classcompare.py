#!/usr/bin/env python

# ClassCompare - Compare entire TC classes to one another using P2.
# This program will skip members of the same superfamily.

import protocol2
import tcdb
import urllib
import scala
import os,re
import blast
import shutil
import pickle
import templates
#import resource
from Bio import SeqIO
from multiprocessing import cpu_count
from optparse import OptionParser
#resource.setrlimit(resource.RLIMIT_NOFILE, (10000,-1))
class Compare:

    def __init__(self):
        # Runtime Variables
        self.subject = None
        self.target = None
        self.threads = cpu_count()
        self.outdir = os.environ['HOME']+'/db/icc'

        # Internal Matricies
        self.familyrelations = []
        self.names = None
        self.pool = None
        self.template = os.environ['CLASSCOMPARE_TEMPLATE']
        self.fastahome = os.environ['HOME']+'/db/families'
        if os.path.exists(self.fastahome) is False:
            os.makedirs(self.fastahome)
        if os.path.exists(self.outdir) is False:
            os.makedirs(self.outdir)

    def __call__(self):
        # Check if this ICC has been done.
        mylock = '%s/%s-%s.lock'%(self.outdir,self.subject,self.target)
        #if os.path.exists(mylock) is True:
        #    return
        # First load all interfamily relations
        url = 'http://tcdb.org/cgi-bin/projectv/classcompare.py?a=%s&b=%s'%(self.subject,self.target)
        results = urllib.urlopen(url)
        [self.familyrelations.append(i.split('\t')) for i in results.read().split('\n')]
        results.close()
        self.familyrelations.pop()
        self.names = tcdb.Names()
        # Generate FASTA files & directories
        self.generate_fastas()
        # We start multithreading after the remote BLAST searches are done.
        self.pool = scala.ThreadPool(self.threads)
        for subject,target in self.familyrelations:
            print subject,target
            self.pool.queueTask(self.subprotocol,(subject,target))
        self.pool.joinAll()
        open(mylock,'wb')
        #self.generate_report()

    def generate_fastas(self):
        if os.path.exists(self.fastahome) is False:
            os.makedirs(self.fastahome)
        families = []
        for i in self.familyrelations:
            families.extend(i)
        families = list(set(families))
        myblast=blast.tools()
        for family in families:
            myblast.gis = []
            if os.path.exists(self.fastahome+'/'+family+'.faa'):
                continue
            print 'Blasting :: %s'%self.names.get_family_abr(family)
            fastas = tcdb.define_family(family,True)
            if fastas is False:
                continue
            for fasta in fastas:
                myblast.blastp(str(fasta.seq))
            # Done, write out this family
            if bool(len(myblast.gis)) is False:
                continue
            myblast.build_raw_fasta()
            shutil.copy(myblast.raw_fasta.name,self.fastahome+'/'+family+'.faa')

    def subprotocol(self,data):
        (subject,target) = data
        thisdir = '%s/%s-%s'%(self.outdir,subject,target)
        if os.path.exists(thisdir+'/gsat.matrix'):
            return
        # Set Protocol1 Variables
        protocol = protocol2.Compare()
        protocol.subject_file = self.fastahome+'/'+subject+'.faa'
        protocol.target_file = self.fastahome+'/'+target+'.faa'
        protocol.outdir = thisdir
        protocol.subject_name = self.names.get_family_abr(subject)
        protocol.target_name = self.names.get_family_abr(target)
        if os.path.exists(protocol.subject_file) is False or os.path.exists(protocol.target_file) is False:
            return
        if os.path.exists(thisdir) is False:
            os.makedirs(thisdir)
        print "Comparison :: %s vs. %s // STARTED"%(subject,target)
        protocol()
        print "Comparison :: %s vs. %s // DONE!"%(subject,target)

    def generate_report(self):
        template = open(self.template,'r').read()
        rows = re.search('{ROW}(.+?){/ROW}',template,re.DOTALL).groups()[0]
        # Iterate through directories
        dirs = os.listdir(self.outdir+'/p2/')
        format = re.compile('(\d+\.[A-Z]\.\d+)-(\d+\.[A-Z]\.\d+)')
        folders = [i for i in dirs if bool(format.search(i))]
        myrows = []
        for folder in folders:
            matrix = self.outdir+'/p2/'+folder+'/gsat.matrix'
            if os.path.exists(matrix):
                try:
                    h = open(matrix,'r')
                    matrix = pickle.load(h)
                    h.close()
                except:
                    print 'corrupt pickle. Running again...'
                    os.remove(matrix)
                    self.subprotocol(folder.split('-'))
                    matrix = pickle.load(open(matrix,'r'))
                if bool(len(matrix)) is False:
                    continue
                matrix.sort(key=lambda x:x['gsat_score'],reverse=True)
                (subject,target) =  folder.split('-')
                row = rows[:]
                row = row.replace('%SUBJECT%',subject)
                row = row.replace('%TARGET%',target)
                row = row.replace('%SUBJECTNAME%',self.names.get_family_abr(subject))
                row = row.replace('%TARGETNAME%',self.names.get_family_abr(target))
                row = row.replace('%SCORE%', str(matrix[0]['gsat_score']))
                myrows.append((row,matrix[0]['gsat_score']))
        myrows.sort(key=lambda x:x[1],reverse=True)
        template =  re.sub('{ROW}.+?{/ROW}',"\n".join([i[0] for i in myrows]),template,flags=re.DOTALL)
        handle = open(self.outdir+'/report.html','wb')
        handle.write(template)
        print 'Done.'

if __name__=='__main__':

    opts = OptionParser(description='Compare entire TC Classes (2 Digit TCIDS). Just select two TC classes to compare. You can compare a class to itself and I will\
    automatically ignore families that are already known to be homologous to save time.',version=1.0)
    opts.add_option('-a',action='store',type='string',dest='subject_class',default=None,help='Subject class. Ex. 1.A')
    opts.add_option('-b',action='store',type='string',dest='target_class',default=None,help='Subject class. Ex. 1.B')
    opts.add_option('-o',action='store',type='string',dest='outdir',default=os.environ['HOME']+'/db/icc',help='Output Directory')
    opts.add_option('--threads',action='store',type='int',dest='threads',default=cpu_count(),help='Threads to use')
    (cli,args) = opts.parse_args()
    (a,b,outdir)=(cli.subject_class,cli.target_class,cli.outdir)
    if a is None or b is None:
        opts.print_help()
        quit()
    cc = Compare()
    cc.subject=a
    cc.target=b
    cc.outdir=outdir
    cc.threads=cli.threads
    cc()

