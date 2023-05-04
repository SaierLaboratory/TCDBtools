#!/usr/bin/env python
import re,sys
import matplotlib
import subprocess
matplotlib.use('TkAgg')
#sys.path.insert(1,'../biov/')
#sys.path.insert(1,'../aad/')
from aad import Average
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
import tss,os
import pickle
from ProjectBio import Emboss
import gsat
import templates
import hmmgap
import optparse

class Ancient:
    
    def __init__(self):
        self.input_aln_file = None
        self.master_seq = {}
        self.window = 40
        self.ram = {}
        self.hits = {}
        self.out='AR2OUT/'
        self.tssdata = {}
        self.swap = False
        self.tssout = 'tss.matrix'
        self.gsatout ='gsat.matrix'
        self.subjects_dict = {}
        self.targets_dict= {}
        self.globaldata=[]
        self.subject_pos = {}
        self.target_pos = {}
        self.smooth_hydro = None
        self.smooth_amphi = None
        self.hmmtop = None
        self.smooth_groups = None
        self.minAln = 20
        self.tms = None
        self.omit_after = -1
        self.omit_before = 0
        self.min_gsat = 9
        self.shuffles = 500
        self.tsv_output = 'TSV_OUTPUT'
        return

    def initialize_output_tsv( self ):
        with open( self.tsv_output, 'w' ) as output:
            #added # for easy header recognition
            output.write('#')
            output.write( '\t'.join( ['query_symbol', 'target_symbol', 'gsat_score', 'aln_length', 'query_tm_count', 'target_tm_count', 'query_aln', 'target_aln', '\n'] ))

    #u hmmtop to get the count of predicted TMSs
    def countTMSs( self, name, seq ):
        with open( "tempfasta.fasta", "w" ) as fasta:
            fasta.write(">{}\n".format( name ))
            fasta.write( str(seq) )
        hmmoutput = subprocess.check_output( ["hmmtop", "-if=tempfasta.fasta"]).decode('utf-8')
        os.remove( "tempfasta.fasta" )

        try:
            #print( "output: " + hmmoutput.split()[3] )
            return hmmoutput.split()[4]
        except:
            return -1
            #raise ValueError( "\nERROR: Could not get number of TMS for {}".format( name ))
        
    def write_output( self, result ):
            subject_symbol = result['subject_symbol']
            target_symbol = result['target_symbol']
            gsat_score = result['gsat_score']
            subject_aln = result['subject_aln']
            target_aln = result['target_aln']
            subject_n_tms = self.countTMSs( subject_symbol, subject_aln )
            target_n_tms  = self.countTMSs( target_symbol, target_aln )
            aln_len = len( subject_aln )
            newresult = [ subject_symbol, target_symbol, gsat_score, aln_len, subject_n_tms, target_n_tms, subject_aln, target_aln, '\n']
            newresult = [ str(r) for r in newresult]
            with open( self.tsv_output, 'a' ) as output:
                #add if statement to make sure gsat_min is accounted for
                #if(gsat_score >= self.min_gsat):
                output.write( '\t'.join( newresult ))

    def build_master_sequence(self):
        pattern = re.compile(r'(\S+)\s+([A-Z,-]+)')
        master_seq = {}
        with open(self.input_aln_file,'r') as h:
            contents = h.read()
            contents = re.sub(r'\*|\.|:','',contents)
            
            #Remove header line before parsing the alignment. I could not 
            #just ignore the header because 'pattern' is compiled and will
            #break if the header has more than 2 words separated by spaces.
            contents = re.sub('^clustal.*\n','', contents, flags=re.IGNORECASE)

            for i in pattern.findall(contents):
                master_seq.setdefault(i[0],[]).append(i[1])
            
           
        for key,seqs in master_seq.items():
            seqs = "".join(seqs)
            self.master_seq[key]=seqs

        return True
        
    def split(self,pos):
        #ab+ changed to w+ python3 edit
        file1 = open(self.out+'first.faa','w+')
        file2 = open(self.out+'second.faa','w+')
        for key,seq in self.master_seq.items():
            first = re.sub('-','',seq[self.omit_before:pos+1])
            second = re.sub('-','',seq[pos+1:self.omit_after])
            
            # Process first half
            desc = '1-%i'%len(first)
            record = SeqRecord(Seq(first),id=key, name=key,description=desc)
            SeqIO.write(record,file1,'fasta')
            
            # Process second half
            desc = '%i-%i'%(len(first)+1,len(first)+len(second))
            record = SeqRecord(Seq(second),id=key, name=key,description=desc)
            SeqIO.write(record,file2,'fasta')
        return

    def make_dicts(self):
        subjects = SeqIO.parse(self.out+'first.faa','fasta')
        self.subjects_dict = SeqIO.to_dict(subjects)
        for i in SeqIO.parse(self.out+'first.faa','fasta'):
            s,e = i.description.split( )[-1].split('-')
            self.subject_pos[i.id]=(int(s),int(e))
        targets = SeqIO.parse(self.out+'second.faa','fasta')
        self.targets_dict = SeqIO.to_dict(targets)
        for i in SeqIO.parse(self.out+'second.faa','fasta'):
            s,e = i.description.split( )[-1].split('-')
            self.target_pos[i.id]=(int(s),int(e))
        
        
    def tssearch(self):
        if os.path.exists(self.out+self.tssout):
            print('Loading existing Alignment Matrix...')
            #changed 'r' to 'rb' for python3
            handle = open(self.out+self.tssout,'rb')
            try:
                (self.tssdata,self.swap) = pickle.load(handle)
                handle.close()
                return
            except ValueError:
                # Corrupt pickle - Trash it, make another.
                os.remove(handle.name)
            handle.close()
        print('Building Alignment Matrix...')
        local = tss.compare()
        local.subject = self.out+'first.faa'
        local.target = self.out+'second.faa'
        local.shuffle = 1
        local.max = 3
        local()
        # Process our data into a subject->data dict
        for result in local.results:
            [self.tssdata.setdefault(result['subject_symbol'],[]).append(result)]
            
        mymatrix = [self.tssdata,local.swap]
        self.swap = local.swap
        with open(self.out+self.tssout,'wb') as pdump:
            pickle.dump(mymatrix,pdump)
        return
    
    def extract_seq(self,symbol):
        seq = self.master_seq[symbol]
        seq = re.sub('-','',seq)
        return SeqRecord(Seq(seq),id=str(symbol),description='')
        
    def write_fastas(self):
        records = []
        for key,seq in self.master_seq.items():
            record = SeqRecord(Seq(seq),id=str(key),description='')
            records.append(record)
        SeqIO.write(records,"hmmtop.in",'fasta')
    
    def run_global_alignments(self):
        if os.path.exists(self.out+self.gsatout):
            #changed 'r' to 'rb' for python3
            handle = open(self.out+self.gsatout,'rb')
            try:
                self.globaldata = pickle.load(handle)
                handle.close()
                return
            except ValueError:
                os.remove(handle.name)
            handle.close()
        gsatcmd = gsat.cmd()
        gsatcmd.gapopen = 8
        gsatcmd.gapextend = 2
        gsatcmd.shuffles = self.shuffles
        emboss = Emboss()
        sub_annotate = hmmgap.annotate()
        tar_annotate = hmmgap.annotate()
        sub_annotate.hmmtop = self.tms['res']
        tar_annotate.hmmtop = self.tms['res']
        #total = len(self.tssdata.items()) * self.assign
        i=0
        print('Performing global alignments')
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
                    print(error)
                    #self.logger.error(error)
                    continue
                (subject_aln,match_aln,target_aln) = emboss.extract_alignment(gsatcmd.outfile.name)
                subject_annotation = sub_annotate(self.extract_seq(subject),subject_aln)
                target_annotation = tar_annotate(self.extract_seq(target['target_symbol']),target_aln)
                target['gsat_score'] = gsatcmd.zscore
                (target['subject_aln'],target['match_aln'],target['target_aln']) = (subject_aln,match_aln,target_aln)
                (target['subject_annotation'],target['target_annotation']) = (subject_annotation,target_annotation)
                target['subject_tmc'] = len(self.tms['res'][subject]) if subject in self.tms['res'] else 0
                target['target_tmc'] = len(self.tms['res'][target['target_symbol']]) if target['target_symbol'] in self.tms['res'] else 0
                self.globaldata.append((target))
                gsatcmd.outfile.close()

        with open(self.out+self.gsatout,'wb') as pdump:
            pickle.dump(self.globaldata,pdump)
        
    def get_relative_positions(self,symbol,start,end,subject=True):
        if subject == True:
            s,e = self.subject_pos[symbol]
        else:
            s,e = self.target_pos[symbol]
        abs_count = 0
        rel_count = 0
        for a in self.master_seq[symbol]:
            if rel_count == s-1+start:
                continue
            abs_count +=1
            if a == '-':
                continue
            rel_count +=1
        the_start = abs_count-1
        abs_count = 0
        rel_count = 0
        for a in self.master_seq[symbol]:
            if rel_count == s-1+end:
                continue
            abs_count +=1
            if a == '-':
                continue
            rel_count +=1
        the_end=abs_count
        return (the_start,the_end)
            
            
        
    def sample_range(self,start,end):
        # Will create a probability distribution of each amino acid within a range
        key = '%i-%i'%(start,end)
        if key in self.ram.keys():
            return self.ram[key]
        AA = 'GALMFWKQESPVICYHRNDT'
        counts = {}
        ncounts = {}
        for seq in self.master_seq:
            myseq = seq[start:end+1]
            for i in myseq:
                if i =='-':
                    continue
                if i in counts.keys():
                    counts[i] += 1
                else:
                    counts[i] = 1
        normalize = (end+1-start)*len(self.master_seq)
        for k,i in counts.items():
            ncounts[k]=(i/1)
        for a in AA:
            if a not in ncounts.keys():
                ncounts[a] = 0
        self.ram[key]=ncounts
        return ncounts
        
    def divergence(self,a,b):
        error = 0
        for k,v in a.items():
            error += abs(v-b[k])
        return error
        
    def scan_window(self):
        length = len(self.master_seq[0])
        results = {}
        for i in range(0,length-self.window+1):
            dist_a = self.sample_range(i,i+self.window)
            for j in range(i+self.window,length-self.window+1):
                dist_b = self.sample_range(j,j+self.window)
                score = self.divergence(dist_a,dist_b)
                results.setdefault(score,[]).append([i,j])
        results.keys().sort()
        self.hits = results
        return
        
    def init_graph(self):
        A=Average()
        A.master_seq=[i[1] for i in self.master_seq.items()]
        A.average_seq()
        A.calculate_hydropathy()
        A.calculate_amphipathicity()
        self.write_fastas()
        A.run_hmmtop()
        A.process_select_groups()
        self.smooth_hydro = A.smooth_hydro
        self.smooth_amphi = A.smooth_amphi
        self.hmmtop = A.hmmtop
        self.tms = A.tms
        self.smooth_groups = A.smooth_groups
        return

    def generate_graph(self,subject,target,filename):
        x_data = range(0, len(self.smooth_hydro))
        mslen = len([i[1] for i in self.master_seq.items()][0])
        diff=(mslen-len(self.smooth_hydro))/2
        #to make dict_items object accessible by location in python3, the object must first be cast to the list constructor
        x1_data = range(0,len(list(self.smooth_groups.items())[0][-1]))
        x2_data = range(0,mslen)
        plt.figure()
        plt.axhline(y=0, color='black')
        plt.ylim(-3, 3)
        plt.xlim(right=mslen)
        plt.plot(x_data, self.smooth_hydro, linewidth=1.0, label="hydrophobicity", color='r')
        plt.plot(x_data, self.smooth_amphi, linewidth=1.0, label="amphipathicity", color='g')
        for pos in self.hmmtop:
            plt.axvline(x=pos-1-diff, ymin=-2, ymax = 0.1, linewidth=1, color='black',alpha=0.2)
        
        plt.axvspan(subject[0]-diff,subject[1]-diff, facecolor="orange", alpha=0.2)
        plt.axvspan(target[0]-diff,target[1]-diff, facecolor="orange", alpha=0.2)

        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
                  ncol=3, fancybox=True, shadow=True)
        plt.xlabel("Residue Number")
        plt.ylabel("Value")
        width = (0.0265)*mslen if mslen > 600 else 15
        plt.grid('on')
        plt.savefig(self.out+'/graphs/'+filename+'.png')
        plt.clf()
        plt.cla()
        plt.close()
        
    def make_graphs(self):
        for i in self.globaldata:
            sub_pos = self.get_relative_positions(i['subject_symbol'],i['subject_start'],i['subject_end'],True)
            tar_pos = self.get_relative_positions(i['target_symbol'],i['target_start'],i['target_end'],False)
            filename = '%s-%s'%(i['subject_symbol'],i['target_symbol'])
            self.generate_graph(sub_pos,tar_pos,filename)
            print('.')
        return
        
    def generate_report(self):
        for ii,i in enumerate(self.globaldata):
            try:
                z = int(i['gsat_score'])
                z = 0 if z<0 else z
            except:
                z = 0
            self.globaldata[ii]['gsat_score'] = z
        self.globaldata.sort(key=lambda x:int(x['gsat_score']),reverse=True)
        
        thandle = open(os.environ['ANCIENT_TEMPLATE'],'r')
        template = thandle.read()
        thandle.close()
        row_template = re.search('{ROWS}(.+?){\/ROWS}',template,re.DOTALL).groups()[0]
        data_template = re.search('{DATA}(.+?){\/DATA}',template,re.DOTALL).groups()[0]
        (myrows,mydatas,mystats) = ([],[],[])
        
        for result in self.globaldata:
            result_copy = result.copy()
            #### The minimum alingment length
            subjectSeq = str(result['subject_aln'])
            subjectSeq = subjectSeq.replace('-','')
            targetSeq  = str(result['target_aln'])
            targetSeq  = targetSeq.replace('-','')
            alnLength  = len(str(result['subject_aln']))
            gsatScore = str(result['gsat_score'])
            #added report has to be minimum gsat or above
            if ( len(subjectSeq) >= self.minAln and
                 len(targetSeq)  >= self.minAln and int(gsatScore) >= self.min_gsat):
                self.write_output(result_copy)
                myrow = row_template[:]
                mydata = data_template[:]
                result['match_aln'] = str(result['match_aln']).replace('.',' ')
                myrow = myrow.replace('{SUBJECTID}',
                                      result['subject_symbol'])
                myrow = myrow.replace('{TARGETID}',
                                      result['target_symbol'])
                myrow = myrow.replace('{GSATZ}',
                                      str(result['gsat_score']))
                mydata = mydata.replace('{SUBJECTID}',
                                        result['subject_symbol'])
                mydata = mydata.replace('{TARGETID}',
                                        result['target_symbol'])
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
                myrows.append(myrow)
                mydatas.append(mydata)

        rows = "".join(myrows)
        data = ''.join(mydatas)
        template = re.sub('{ROWS}.+?{\/ROWS}',rows,template,flags=re.DOTALL)
        template = re.sub('{DATA}.+?{\/DATA}',data,template,flags=re.DOTALL)
        #encode the template back to bytes form
        template_bytes = template.encode('utf-8')
        with open(self.out+'/report.html','wb+') as output:
            #write out template byte form to stream
            output.write(template_bytes)
        print("Report Generated")
        print("Graphs are being generated now. These take a while, but you can view the report while they are generated")
        return
            
            
        
    
if __name__=='__main__':
    desc = '''
Welcome to ANCIENT. This tool will find intragenic repeats using
a multiple alignment. Please make sure your multiple alignment
contains no large gaps in the middle. Gaps at the edges are ok.

FOR DETAILED INSTRUCTIONS - READ THIS MANUAL:
<http://tcdb.org/biov/AR_INSTRUCTIONS.pdf>

// Developed by Vamsee Reddy
'''
    desc = " ".join(desc.split())
    ver = '2.0'
    opts = optparse.OptionParser(description=desc,version=ver)
    opts.add_option('-i',
                    action='store',
                    metavar = 'File',
                    dest='input',
                    type='string',
                    default='input.aln',
                    help='Path to the multiple alignment file.')
    opts.add_option('-n',
                    action='store',
                    metavar = 'Position',
                    dest='split',
                    type='int',
                    default=None,
                    help='Position in multiple alignment where the cut will be made.')
    opts.add_option('-o',
                    action='store',
                    metavar = 'Path',
                    dest='out',
                    type='string',
                    default='AROUT',
                    help='Output directory where results will be saved.')
    opts.add_option('--minaln',
                    action='store',
                    metavar = 'Int',
                    dest='minaln',
                    type='int',
                    default=20,
                    help='Minimum alignment lengths to keep in report (Default 20).')
    opts.add_option('--omit_before',
                    action='store',
                    metavar = 'Int',
                    dest='omit_before',
                    type='int',
                    default=0,
                    help='Ignore amino acids BEFORE this position.')
    opts.add_option('--omit_after',
                    action='store',
                    metavar = 'Int',
                    dest='omit_after',
                    type='int',
                    default=-1,
                    help='Ignore amino acids AFTER this position.')
    opts.add_option( '--min_gsat',
                     action='store',
                     metavar = 'Int',
                     dest='min_gsat',
                     type='int',
                     default=9,
                     help='Minimum gsat score to keep in report (Default 9).')
    opts.add_option( '--shuffles',
                     action='store',
                     metavar = 'Int',
                     dest='shuffles',
                     type='int',
                     default=500,
                     help='Number of GSAT shuffles to use (Default 500).')

    (cli,args) = opts.parse_args()

    if ( os.path.exists(cli.input) is False ):
        opts.print_help()
        exit()
    
    AR = Ancient()
    AR.omit_before = cli.omit_before
    AR.omit_after= cli.omit_after
    AR.input_aln_file = cli.input
    AR.minAln = cli.minaln
    AR.out = cli.out+'/'
    AR.min_gsat = cli.min_gsat
    AR.shuffles = cli.shuffles
    AR.tsv_output = cli.out + '/report.tsv'
    if os.path.exists(AR.out) is False:
        os.mkdir(AR.out)
    if os.path.exists(AR.out+'/graphs') is False:
        os.mkdir(AR.out+'/graphs')
    AR.initialize_output_tsv()
    AR.build_master_sequence()
    AR.split(cli.split)  
    AR.tssearch()
    AR.make_dicts()
    AR.init_graph()
    AR.run_global_alignments()
    #AR.init_graph()
    AR.generate_report()
    AR.make_graphs()
    
