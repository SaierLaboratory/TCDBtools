#!/usr/bin/env python

'''
This is protocol2 - python edition.
Perform a local alignment using TSS, cut sequences, run GSAT, run HMMTOP,
and highlight alignments. Thats all.
Requires SSearch36 From the FASTA package somewhwere in your $PATH variable.
Created by Vamsee Reddy, modified by Gabo Moreno-Hagelsieb.
'''

import hmmtop
import gsat
import tss,os
import pickle
import optparse
import tempfile
import hmmgap
import logging
import templates
import re
from ProjectBio import Emboss
from Bio import SeqIO
from progressbar import ProgressBar,Bar,Percentage
#os.environ['BIOV_DEBUG'] = 'True'
class Compare:

    def __init__(self):
        # Runtime variables
        self.subject_file = None
        self.subject_name = None
        self.target_file = None
        self.target_name = None
        self.outdir = 'protocol2'
        self.assign = 3
        self.tss_shuffle = 200
        self.gsat_shuffle = 300
        self.swap = False
        self.verbose_switch = False
        self.qrestrict = None
        self.trestrict = None
        self.minAln = 60
        self.min_gsat = 9

        # Internal Matrices
        self.tssdata = {}
        self.subjects_dict = {}
        self.targets_dict = {}
        self.tms = {}
        self.globaldata = []
        self.stats = {}
        #self.tmos = []

        # Storage variables
        self.tssout = '/tss.matrix'
        self.gsatout ='/gsat.matrix'

        # Alignment variables
        self.gapopen = 8
        self.gapextend = 2

        #self.init_logging()

    def __call__(self):
        if os.path.exists(self.outdir) is False:
            os.mkdir(self.outdir)
        self.tssearch()
        # Swap subjects & targets if needed - subject is always shorter
        if self.swap:
            (self.subject_file,self.target_file) = (self.target_file,self.subject_file)
            (self.subject_name,self.target_name) = (self.target_name,self.subject_name)
            (self.qrestrict,self.trestrict) = (self.trestrict,self.qrestrict)
        # Set subject & target names
        if self.subject_name is None:
            self.subject_name = self.subject_file.split('/')[-1]
        if self.target_name is None:
            self.target_name = self.target_file.split('/')[-1]
        # Write selected subjects & targets, then run HMMTOP
        subject_ids = self.tssdata.keys()
        target_ids = []
        for i in self.tssdata.values():
            target_ids.extend([x['target_symbol'] for x in i])
        target_ids = list(set(target_ids))
        subjects = SeqIO.parse(self.subject_file,'fasta')
        self.subjects_dict = SeqIO.to_dict(subjects)
        targets = SeqIO.parse(self.target_file,'fasta')
        self.targets_dict = SeqIO.to_dict(targets)
        with open(self.outdir+'/querys.faa','w+') as sf:
            for sid in subject_ids:
                SeqIO.write(self.subjects_dict[sid],sf,'fasta')
        with open(self.outdir+'/targets.faa','w+') as tf:
            for tid in target_ids:
                SeqIO.write(self.targets_dict[tid],tf,'fasta')
        hmt = hmmtop.tools()
        hmt.add_library('subjects',self.outdir+'/querys.faa')
        hmt.add_library('targets',self.outdir+'/targets.faa')
        hmt.scan_libraries()
        self.tms = hmt.results
        # Run GSAT on all of these
        self.run_global_alignments()
        # Sort global data
        self.globaldata.sort(key=lambda x:float(x['gsat_score']),reverse=True)
        for ii,i in enumerate(self.globaldata):
            try:
                z = int(i['gsat_score'])
                z = 0 if z<0 else z
            except:
                z = 0
            self.globaldata[ii]['gsat_score'] = z
        self.globaldata.sort(key=lambda x:int(x['gsat_score']),reverse=True)
        # Generate HTML report
        self.generate_report()

    def verbose(self,msg):
        if self.verbose_switch:
            print (msg)

    def init_logging(self):
        logger = logging.getLogger('PROTOCOL2')
        hdlr = logging.FileHandler('/var/tmp/PROTOCOL2.log')
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr) 
        logger.setLevel(logging.WARNING)
        self.logger = logger

    def tssearch(self):
        if os.path.exists(self.outdir+self.tssout):
            self.verbose('Loading existing Alignment Matrix...')
            handle = open(self.outdir+self.tssout,'r')
            try:
                (self.tssdata,self.swap) = pickle.load(handle)
                #handle.close()
                #return
            except ValueError:
                # Corrupt pickle - Trash it, make another.
                os.remove(handle.name)
            handle.close()
        self.verbose('Building Alignment Matrix...')
        local = tss.compare()
        local.subject = self.subject_file
        local.target = self.target_file
        local.shuffle = self.tss_shuffle
        local.max = self.assign
        local()
        # Process our data into a subject->data dict
        for result in local.results:
            [self.tssdata.setdefault(result['subject_symbol'],[]).append(result)]

        mymatrix = [self.tssdata,local.swap]
        self.swap = local.swap
        with open(self.outdir+self.tssout,'wb') as pdump:
            pickle.dump(mymatrix,pdump)
        return

    def run_global_alignments(self):
        if os.path.exists(self.outdir+'/'+self.gsatout):
            handle = open(self.outdir+'/'+self.gsatout,'r')
            try:
                self.globaldata = pickle.load(handle)
                #handle.close()
                #return
            except ValueError:
                os.remove(handle.name)
            handle.close()
        gsatcmd = gsat.cmd()
        gsatcmd.gapopen = self.gapopen
        gsatcmd.gapextend = self.gapextend
        emboss = Emboss()
        sub_annotate = hmmgap.annotate()
        tar_annotate = hmmgap.annotate()
        sub_annotate.hmmtop = self.tms['subjects']
        tar_annotate.hmmtop = self.tms['targets']
        total = len(self.tssdata.items()) * self.assign
        i=0
        self.verbose('Performing global alignments')
        if self.verbose_switch:
            pbar = ProgressBar(widgets=[Percentage(), Bar()],maxval=total).start()
        for subject,targets in self.tssdata.items():
            subject_seq = self.subjects_dict[subject]
            for target in targets:
                mytarget = self.targets_dict[target['target_symbol']].seq[target['target_start']-1:target['target_end']]
                mysubject = subject_seq.seq[target['subject_start']-1:target['subject_end']]
                gsatcmd.asequence = str(mysubject)
                gsatcmd.bsequence = str(mytarget)
                try:
                    gsatcmd()
                except:
                    error = 'GSAT problem with:\nSubject: %s | Target: %s // \n>subject\n%s\n\n>target\n%s'%(self.subject_file,self.target_file,mysubject,mytarget)
                    print (error)
                    #self.logger.error(error)
                    continue
                (subject_aln,match_aln,target_aln) = emboss.extract_alignment(gsatcmd.outfile.name)
                subject_annotation = sub_annotate(self.subjects_dict[subject],subject_aln)
                target_annotation = tar_annotate(self.targets_dict[target['target_symbol']],target_aln)
                target['gsat_score'] = gsatcmd.zscore
                (target['subject_aln'],target['match_aln'],target['target_aln']) = (subject_aln,match_aln,target_aln)
                (target['subject_annotation'],target['target_annotation']) = (subject_annotation,target_annotation)
                target['subject_tmc'] = len(self.tms['subjects'][subject]) if subject in self.tms['subjects'] else 0
                target['target_tmc'] = len(self.tms['targets'][target['target_symbol']]) if target['target_symbol'] in self.tms['targets'] else 0
                self.globaldata.append((target))
                gsatcmd.outfile.close()
                if self.verbose_switch:
                    pbar.update(1+i)
                    i +=1
        with open(self.outdir+self.gsatout,'wb') as pdump:
            pickle.dump(self.globaldata,pdump)

    def calculate_overlap(self,subject,target,sd,mo=0.5):
        soverlap,toverlap,score = {},{},0
        tms = re.compile(r'(\d+)[+_]+\+')
        sub = tms.finditer(subject)
        tar = tms.finditer(target)
        for s in sub:
            soverlap.setdefault(s.groups()[0],(s.start(),s.end()))
        for t in tar:
            toverlap.setdefault(t.groups()[0],(t.start(),t.end()))
        totals = [len(soverlap),len(toverlap)]
        totals.sort()
        for stms,srange in soverlap.items():
            for ttms,trange in toverlap.items():
                # Get min len
                srangee=range(srange[0],srange[1])
                trangee=range(trange[0],trange[1])
                minlen=[len(srangee),len(trangee)]
                minlen.sort()
                minlen = int(minlen[0]*mo)
                score += len(set(srangee)&set(trangee))
                if len( set(srangee)&set(trangee) ) >= minlen:
                    overlap = [int(stms),int(ttms)]
                    overlap.sort()
                    overlap = "-".join([str(i) for i in overlap])
                    self.stats.setdefault(overlap,[])
                    self.stats[overlap].append(sd)
        return float(score)/float(totals[0])


    def generate_report(self):
        thandle = open(os.environ['PROTOCOL2_TEMPLATE'],'r')
        template = thandle.read()
        thandle.close()
        row_template = re.search('{ROWS}(.+?){\/ROWS}',template,re.DOTALL).groups()[0]
        data_template = re.search('{DATA}(.+?){\/DATA}',template,re.DOTALL).groups()[0]
        (myrows,mydatas,mystats) = ([],[],[])
        ##### to produce tab separated table:
        tblComment = ( "# Query: " + self.subject_name
                       + " Target: " + self.target_name)
        tblLines  = [tblComment]
        listHead = ["Query ID","Target ID","SS Z-Score","GSAT Z-Score",
                    "Query align-length","Target align-length",
                    "Query seq","Target seq","TMO-score"]
        tblHeader = "\t".join(listHead)
        tblLines.append(tblHeader)
        for result in self.globaldata:
            if self.qrestrict is not None and self.trestrict is not None:
                if ( int(result['subject_tmc']) != self.qrestrict
                     or int(result['target_tmc']) != self.trestrict ):
                    continue
            #### The minimum alingment length
            subjectSeq = str(result['subject_aln'])
            subjectSeq = subjectSeq.replace('-','')
            targetSeq  = str(result['target_aln'])
            targetSeq  = targetSeq.replace('-','')
            alnLength  = len(str(result['subject_aln']))
            gsatScore = str(result['gsat_score'])
            if ( len(subjectSeq) >= self.minAln and
                 len(targetSeq)  >= self.minAln and 
                 int(gsatScore) >= self.min_gsat):
                myrow = row_template[:]
                mydata = data_template[:]
                result['match_aln'] = str(result['match_aln']).replace('.',' ')
                myrow = myrow.replace('{SUBJECTID}',
                                      result['subject_symbol'])
                myrow = myrow.replace('{TARGETID}',
                                      result['target_symbol'])
                myrow = myrow.replace('{SSZ}',
                                      str(round(result['z_score'],2)))
                myrow = myrow.replace('{GSATZ}',
                                      str(result['gsat_score']))
                mydata = mydata.replace('{SUBJECTID}',
                                        result['subject_symbol'])
                mydata = mydata.replace('{TARGETID}',
                                        result['target_symbol'])
                mydata = mydata.replace('{SSZ}',
                                        str(round(result['z_score'],2)))
                mydata = mydata.replace('{GSATZ}',
                                        str(result['gsat_score']))
                mydata = mydata.replace('%ATMS%',
                                        str(result['subject_annotation']))
                mydata = mydata.replace('%BTMS%',
                                        str(result['target_annotation']))
                mydata = mydata.replace('%ASEQUENCE%',
                                        str(result['subject_aln']))
                mydata = mydata.replace('%BSEQUENCE%',
                                        str(result['target_aln']))
                mydata = mydata.replace('%MATCH%',
                                        str(result['match_aln']))
                mydata = mydata.replace('%STMC%',
                                        str(result['subject_tmc']))
                mydata = mydata.replace('%TTMC%',
                                        str(result['target_tmc']))
                tmo = 0.0
                try:
                    tmo=self.calculate_overlap(result['subject_annotation'],
                                               result['target_annotation'],
                                               result['gsat_score'])
                except:
                    pass
                myrows.append(myrow)
                mydatas.append(mydata)
                if int(result['gsat_score']) > 8:
                    tblLine = "\t".join([ result['subject_symbol'],
                                          result['target_symbol'],
                                          str(round(result['z_score'],2)),
                                          str(result['gsat_score']),
                                          str(len(subjectSeq)),
                                          str(len(targetSeq)),
                                          str(result['subject_aln']),
                                          str(result['target_aln']),
                                          str(format(tmo, '.2f'))])
                    tblLines.append(tblLine)
        template = template.replace('{SUBJECTNAME}',self.subject_name)
        template = template.replace('{TARGETNAME}',self.target_name)
        rows = "".join(myrows)
        data = ''.join(mydatas)
        template = re.sub('{ROWS}.+?{\/ROWS}',rows,template,flags=re.DOTALL)
        template = re.sub('{DATA}.+?{\/DATA}',data,template,flags=re.DOTALL)
        stats = self.stats.items()
        #stats.sort(key=lambda x:len(x[-1]),reverse=True)
        sorted(stats)
        for key,count in stats:
            line="%s\t%i"%(key,len(count))
            mystats.append(line)
        stats = "\n".join(mystats)
        template = re.sub('%STATS%',stats,template)
        with open(self.outdir+'/report.html','w+') as output:
            output.write(template)
        ##### print tab-separated table
        with open(self.outdir+'/report.tbl','w+') as outputtab:
            for i in tblLines:
                outputtab.write(str(i)+"\n")
        return


### here this thing runs
if __name__=='__main__':
    desc = '''
Welcome to Protocol2. This tool will allow you to rapidly locate
homologs between two fasta files.

'Subject', 'Target', and 'Outdir' Are the only mandatory options.
--Subject & --Target are used to label items on your actual report.
Example usage:

protocol2 -s 2.A.1.faa -t 2.A.3.faa -o mydir --subject='APC' --target='MFS'

// Developed by Vamsee Reddy
'''
    desc = " ".join(desc.split())
    ver = '3.0'
    opts = optparse.OptionParser(description=desc,version=ver)
    opts.add_option('-q',
                    action='store',
                    metavar="File",
                    dest='query',
                    type='string',
                    default='query.faa',
                    help='Path to file with query sequences in fasta format')
    opts.add_option('-t',
                    action='store',
                    metavar="File", 
                    dest='target',
                    type='string',
                    default='target.faa',
                    help='Path to file with target sequences in fasta format')
    opts.add_option('-o',
                    action='store',
                    metavar="Path", 
                    dest='outdir',
                    type='string',
                    default='protocol2',
                    help='Output directory where results will be saved.')
    opts.add_option('--query',
                    action='store',
                    metavar="String", 
                    dest='qname',
                    type='string',
                    default=None,
                    help='Label to identify the query in the report.')
    opts.add_option('--target',
                    action='store',
                    metavar="String",
                    dest='tname',
                    type='string',
                    default=None,
                    help='Label to identify the target in the report.')
    opts.add_option('--assign',
                    action='store',
                    metavar="Int",
                    dest='num',
                    type='int',
                    default=3,
                    help='Maximum number of target matches per protein') 
    opts.add_option('--shuffles',
                    action='store',
                    metavar="Int",
                    dest='rand',
                    type='int',
                    default=300,
                    help='Number of GSAT shuffles to apply')
    opts.add_option(
        '--qtms',
        action='store',
        metavar="Int",
        dest='qrestrict',
        type='int',
        default=None,
        help='Report will contain queries with X TMSs. TTMS must be set to work')
    opts.add_option(
        '--ttms',
        action='store',
        metavar="Int",
        dest='trestrict',
        type='int',
        default=None,
        help='Report will contain targets with Y TMSs.  QTMS must be set to work')
    opts.add_option('--minaln',
                    action='store',
                    metavar="Int",
                    dest='minAln',
                    default=60,
                    type='int',
                    help='Minimum alignment length to keep (Default 60)')
    opts.add_option('--min_gsat',
                    metavar="Int",
                    action='store',
                    dest='min_gsat',
                    default=9,
                    type='int',
                    help='Minimum gsat score to keep in the report (Default 9)')
    (cli,args) = opts.parse_args()
    if ( os.path.exists(cli.query) is False
         or os.path.exists(cli.target) is False ):
        opts.print_help()
        exit()
    protocol = Compare()
    protocol.subject_file = cli.query
    protocol.target_file = cli.target
    protocol.subject_name = cli.qname
    protocol.target_name = cli.tname
    protocol.outdir = cli.outdir
    protocol.assign = cli.num
    protocol.gsat_shuffle = cli.rand
    protocol.verbose_switch = True
    protocol.qrestrict = cli.qrestrict
    protocol.trestrict = cli.trestrict
    protocol.minAln = cli.minAln
    protocol.min_gsat = cli.min_gsat
    protocol()
