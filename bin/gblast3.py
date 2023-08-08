#!/usr/local/anaconda3/bin/python3
# coding=utf-8
import os, sys, re, csv
from pprint import pprint
import resource,copy
from Bio.Blast import NCBIXML
from Bio import SeqIO
from ProjectBio import ParseDefline, nr_dict, ParseTC
from Bio.Blast.Applications import NcbiblastpCommandline
from urllib.parse import quote as urlencode
from urllib.request import urlopen
import hmmtop,shutil
from math import ceil
import sys,os,pickle
import subprocess
import matplotlib
import hmmgap
import tcdb
import cdd
from decimal import Decimal
from sys import exit
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab,re,tempfile
import requests
import bz2
import re
import tempfile
import io
from PlotBio import *
rsrc = resource.RLIMIT_DATA
soft, hard = resource.getrlimit(rsrc)
##resource.setrlimit(rsrc, (1073741824, hard)) #limit to one gig, omg..

import os
import gzip
import tarfile
import zipfile
import bz2



def extract_archive(from_path):
    cmd = f"zless {from_path}"
    handle = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    results = handle.communicate()[0]
    #print(results)
    return (results)
    

class Tools:

    def __init__(self):
        self.chart3_find={}
        self.chart2_find={}


        self.dbfile = '/db/tcdb'
        self.tcdb = os.environ['HOME']+self.dbfile
        self.query = False
        self.indir = False
        self.goodresults  = []
        self.bestresults  = []
        self.notmsresults = []
        self.substrates = {}
        self.debug=True
        self.tms = {}
        self.expect = 0.001
        self.minlen = 50
        self.mincov = 40.0
        self.myqueries = False
        self.mytcdb = False
        self.rows = []
        self.data = []
        self.globalcount = 0
        self.cdd_on = False
        self.abbreviations = {}
        self.ortho = False
        self.query_gis = False
        self.target_gis = False
        self.esort = False
        self.alabel = 'Query'
        self.blabel = 'TCDB-Hit'
        #self.seqParse=SeqIO.parse(open("/System/Volumes/Data/tcdbPegasus/ResearchData/Users/clark/gblast3_Final/proteome_test.faa"),'fasta')
        #self.seqParse=SeqIO.parse(io.StringIO(self.query),'fasta')
        self.names  = tcdb.Names()
        self.loadSubstrates()

        #Sequence dictionaries
        self.queries = {}
        self.tcdbHits = {}        

        #self.substrates = tcdb.Substrates()
        tcdb.use_local()
        tcdb.use_local_betabarrel()

    def generate_tcdb(self):
        tcdb=[]
        with bz2.open("/ResearchData/pfam/tcdb.pfam-a.hmmscan.bz2", "rt") as bz_file:
            for line in bz_file:
                tcdb.append(line.rstrip())
        index_row=tcdb[2]
        head_row =tcdb[1]
        def getList(row,index_row):
            if "# target name" in row:
                index=[]
                pre_char=" "
                for ind,char in enumerate(index_row):
                    if char!=" " and pre_char==" ":
                        index.append(ind)
                    pre_char=char
                output=[]
                for ind,i in enumerate(index[:-1]):
                    output.append(row[i:index[ind+1]].strip("#").strip())
                output.append(row[index[-1]:].strip("#").strip())
            else:
                output=re.split(r'\s+',row)
                output=output[:22] + [" ".join(output[22:])]
            return(output)
        tcdb_head=getList(row=head_row,index_row=index_row)
        
        char2_final={}
        for value_row in tcdb[3:]:
            if "#" not in value_row[0]:
                value_list=getList(value_row,index_row)
                value_dict={}
                for ind,value in enumerate(value_list):
                    if tcdb_head[ind] not in value_dict:
                        value_dict[tcdb_head[ind]]=value
                    else:
                        i=0
                        while True:
                            if f"{tcdb_head[ind]}_{i}" in value_dict:
                                i+=1
                            else:
                                value_dict[f"{tcdb_head[ind]}_{i}"]=value
                                break
                value_dict["query name"]=value_dict["query name"].split('-')[-1]
                _pre_chart2=char2_final[value_dict["query name"].strip()] if value_dict["query name"].strip() in char2_final else []
            
                #char3_final[value_dict["query name"].strip()]=_pre_chart3.append(value_dict)
                char2_final[value_dict["query name"].strip()]=_pre_chart2+[value_dict]
                #print(len(char2_final[value_dict["query name"].strip()]))
              
        # acc gather
        for d,v in char2_final.items():
            _acc=[]
            for c in v:
                _acc.append(c["accession"])
            for ind,_ in enumerate(v):
                char2_final[d][ind]["accAll"]=_acc
        self.chart2=char2_final
        return char2_final
        

    #Run hmmscan and parse Pfam output
    def generate_chart3(self):  
        outPfam = os.path.join(self.indir,'pfam.out')
        os.system(f"bash run_pfam.sh {os.path.join(self.indir,'myqueries.faa')} {outPfam}")      
        with open(outPfam,"r") as f:
            lines  = f.readlines()
            chart3 = [line.rstrip() for line in lines]
        index_row=chart3[2]
        head_row =chart3[1]
        def getList(row,index_row):
            index=[]
            pre_char=" "
            for ind,char in enumerate(index_row):
                if char!=" " and pre_char==" ":
                    index.append(ind)
                pre_char=char
            output=[]
            for ind,i in enumerate(index[:-1]):
                output.append(row[i:index[ind+1]].strip("#").strip())
            output.append(row[index[-1]:].strip("#").strip())
            return(output)
        chart3_head=getList(row=head_row,index_row=index_row)
        
        char3_final={}
        for value_row in chart3[3:]:
            if "#" not in value_row[0]:
                value_list=getList(value_row,index_row)
                value_dict={}
                for ind,value in enumerate(value_list):
                    if chart3_head[ind] not in value_dict:
                        value_dict[chart3_head[ind]]=value
                    else:
                        i=0
                        while True:
                            if f"{chart3_head[ind]}_{i}" in value_dict:
                                i+=1
                            else:
                                value_dict[f"{chart3_head[ind]}_{i}"]=value
                                break
                _pre_chart3=char3_final[value_dict["query name"].strip()] if value_dict["query name"].strip() in char3_final else []
            
                #char3_final[value_dict["query name"].strip()]=_pre_chart3.append(value_dict)
                char3_final[value_dict["query name"].strip()]=_pre_chart3+[value_dict]
                #print(len(char3_final[value_dict["query name"].strip()]))
              
        # acc gather
        for d,v in char3_final.items():
            _acc=[]
            for c in v:
                _acc.append(c["accession"])
            for ind,_ in enumerate(v):
                char3_final[d][ind]["accAll"]=_acc
        self.chart3=char3_final
        return char3_final
        
    def blast_all(self):
        try:
            os.makedirs(self.indir+"/xml")
        except:
            pass
        if self.ortho is not False:
            self.prep_orthologs()
        queries = SeqIO.parse(io.StringIO(self.query),'fasta')
        for query in queries:
            
            query_file = tempfile.NamedTemporaryFile(mode='wt+',delete=True)
            SeqIO.write(query,query_file,'fasta')
            query_file.flush()
            blast_out = "'"+self.indir+'/xml/'+query.id+".xml"+"'"
            # removes the quotes in a rly cool way :)
            if os.path.exists(blast_out[1:-1]) is True:
                continue # Blast xml already exists...
            blastp = NcbiblastpCommandline(
                query=query_file.name,
                db=self.tcdb,
                evalue=self.expect,
                out=blast_out,
                outfmt=5,
                comp_based_stats='0'
            )
            blastp()
            print(("Blasted :: %s" %query.id))
            
            
    def prep_orthologs(self):
        queries = SeqIO.parse(open(self.query),'fasta')
        targets = SeqIO.parse(open(self.ortho),'fasta')
        # Load query Gis
        query_gis = open(self.query_gis,'r')
        target_gis = open(self.target_gis,'r')
        queries = SeqIO.to_dict(queries)
        targets = SeqIO.to_dict(targets)
        # Make ortho directory
        try:
            os.makedirs(self.indir+"/orthologs")
        except:
            pass
        # Write queries & targets
        myqueries = open(self.indir+'/orthologs/myqueries.faa','wb')
        mytargets = open(self.indir+'/orthologs/mytargets.faa','wb')
        for sgi in query_gis:
            SeqIO.write(queries[sgi.strip()],myqueries,'fasta')
        for tgi in target_gis:
            SeqIO.write(targets[tgi.strip()],mytargets,'fasta')
        myqueries.flush()
        mytargets.flush()
        #os.system('makeblastdb -dbtype prot -in '+mytargets.name)
        subprocess.call(('makeblastdb -dbtype prot -in '+mytargets.name).split())
        self.query = myqueries.name
        self.tcdb = mytargets.name


    def load_good(self):
        db = self.indir+"/goodresults.db"
        if os.path.exists(db) and self.debug:
            self.goodresults = pickle.load(open(db,'r'))

            return
        xml = os.listdir(self.indir+"/xml")
        rez = []
        for res in xml:
            if os.path.exists(self.indir+'/xml/'+res) is False:
                continue
            try:
                for blast in NCBIXML.parse(open(self.indir+'/xml/'+res)):
                    query =  blast.query
                    try:
                        descriptions = [i for i in blast.descriptions]
                        alignments = [i for i in blast.alignments]
                        results = list(zip(descriptions,alignments))
                        results.sort(key=lambda x:x[0].e,reverse=False)
                        record = results[0] # Contains top description and top alignment object in tuple.
                        hsps = [i for i in record[1].hsps]
                        hsps.sort(key=lambda x:x.expect)
                        hsp = hsps[0]

                        if not blast.query_length:
                                print('No Query Length')        
                                print((vars(blast)))

                        if not record[1].length:

                                print('No hit length')
                                print((vars(record[1])))

                        q_len = blast.query_length
                        h_len = record[1].length

                        qcov = float('%.1f'%(((hsp.query_end-hsp.query_start)/float(q_len))*100))
                        hcov = float('%.1f'%(((hsp.sbjct_end-hsp.sbjct_start)/float(h_len))*100))

                        if qcov >= 100:
                
                                qcov = 100.0
                        
                        if hcov >= 100:

                                hcov = 100.0

                        
                        if record[0].e > self.expect or len(hsp.match) < self.minlen:
                                continue
                        

                        if qcov >= self.mincov or hcov >= self.mincov:

                                rez.append((query,record,hsp,q_len,h_len,qcov,hcov)) # (genome ID, hit record <e,title>, hit.hsp)i
                    except:
                        pass
            except:
                continue
        self.goodresults = rez
        pickle.dump(self.goodresults,open(db,'wb'))
        return

    def write_fastas(self):
        mytcdb = []
        myquery = []
        tcdb = SeqIO.parse(self.tcdb,'fasta')
        self.tcdbHits = nr_dict(tcdb)
        queries= SeqIO.parse(io.StringIO(self.query),'fasta')
        self.queries= SeqIO.to_dict(queries)
        ##### gabo's addition
        done = dict()
        
        for query,hit,hsp,q_len,h_len,qcov,hcov in self.goodresults:
            hit = hit[0]
            query=ParseDefline(query).id
            myquery.append(self.queries[str(query)])
            subject = ParseDefline(hit.title,True).id
            if subject not in list(done.keys()):
                done[subject] = 1
                try:
                    mytcdb.append(self.tcdbHits[subject])
                except:
                    (family,tcid,acc) = ParseTC(hit.title)
                    #print tcid,acc
                    print((ParseDefline(hit.title,True).id))
                    #print hit.title
                    quit()

        query_file=open(self.indir+"/myqueries.faa",'wt+')
        tcdb_file=open(self.indir+"/mytcdb.faa",'wt+')
        SeqIO.write(myquery,query_file,'fasta')
        SeqIO.write(mytcdb,tcdb_file,'fasta')

    def hmmtop(self):
        db = self.indir+'/hmmtop.db'
        if os.path.exists(db) and self.debug:
            self.tms = pickle.load(open(db,'r'))
            return
        ht = hmmtop.tools()
        ht.add_library('queries',self.indir+"/myqueries.faa") # Genome
        ht.add_library('tcdb',self.indir+"/mytcdb.faa")
        ht.scan_libraries()
        pickle.dump(ht.results,open(db,'wb'))
        self.tms = ht.results
        return

    def calculate_tms_scores(self):
        for genome,tcdb,hsp,q_len,h_len,qcov,hcov in self.goodresults:
            genome=ParseDefline(genome).id
            tcdb=tcdb[0]
            delta  = hsp.sbjct_start-hsp.query_start
            tcdbid = ParseDefline(tcdb.title,True).id
            try:
                genome_tms = list(self.tms['queries'][genome].values())
                tcdb_tms = list(self.tms['tcdb'][tcdbid].values())

            except KeyError:
                # Genome or TCDB hit dont have a TMS.
                # These must be manually revised later!
                self.notmsresults.append((genome,tcdb,hsp,q_len,h_len,qcov,hcov,None))
                continue
            g_tms = [[i[0]+delta,i[1]+delta] for i in genome_tms]
            overlap = self.find_overlap(g_tms,tcdb_tms)
            row= (genome,tcdb,hsp,q_len,h_len,qcov,hcov,overlap)
            self.bestresults.append(row)

        if self.esort:

                self.bestresults.sort(key=lambda x:(x[1].e,x[7],self.tcsort(x[1].title)))
                self.notmsresults.sort(key=lambda x:(x[1].e,self.tcsort(x[1].title)))

        else:

                self.bestresults.sort(key=lambda x:(self.tcsort(x[1].title),x[1].e,-1.0*max(x[6],x[5])))
                self.notmsresults.sort(key=lambda x:(self.tcsort(x[1].title),x[1].e,-1.0*max(x[6],x[5])))        



    def find_overlap(self,query,target):
        queries = [list(range(i[0],i[1]+1)) for i in query]
        targets = [list(range(i[0],i[1]+1)) for i in target]
        overlap = []
        for sub in queries:
            for tar in targets:
                overlap.extend(set(sub)&set(tar))
        return float(len(set(overlap)))/20

    def title_extract(self,string): #returns: (acc,tcid)
        string = string.split(" ")
        gi = string[0].split('|')[-1]
        tcid = string[1]
        return (gi,tcid)

    def tcsort(self,title):
        if self.ortho is not False:
            return 1
        title=ParseDefline(title,True).description
        (family,tc,acc) = ParseTC(title)
        tc = tc.split('.')
        return ( int(tc[0]), str(tc[1]), int(tc[2]), int(tc[3]), int(tc[4]) )

    def build_view(self,data):
        (genome,tcdb,hsp,q_len,h_len,qcov,hcov,overlap) = data

        try:
            os.mkdir(self.indir+"/img")
        except:
            pass
        try:
            os.mkdir(self.indir+"/htmls")
        except:
            pass
        genome=ParseDefline(genome).id
        tid = ParseDefline(tcdb.title,True).id
        if os.path.exists(self.indir+"/img/"+genome+".png") is False:

            try:
                san = self.queryhmg(self.myqueries[genome],hsp.query)
                tan = self.tcdbhmg(self.mytcdb[tid],hsp.sbjct)
                self.what(hsp.query,hsp.sbjct,self.indir+"/img/"+genome+".png",[san,tan])
            except:
                print((hsp.query))
                print((hsp.sbjct))
                print(genome)
                print('error, quit')
                quit()

        (family,tcid,acc) = ParseTC(ParseDefline(tcdb.title,True).description)
        try:
            query_tms = len(self.tms['queries'][genome])
        except:
            query_tms = 0
        try:
            hit_tms = len(self.tms['tcdb'][ParseDefline(tcdb.title,True).id])
        except:
            hit_tms = 0
        self.globalcount += 1
        if self.ortho is False:
            family = self.names.get_family_abr('.'.join(tcid.split('.')[0:3]))
        else:
            family = 'family_place_holder'


        '''
        Edits made by Vasu Pranav Sai Iddamsetty (VI)
        
        -Adding another file that is a tsv(tab seperated values) file so that the output of gblast
         may be parsed by another program easily.
         
         ****************************************************
         
                This is where the substrates are found.
                They populate the 'mysubstrate' variable.        
         
         ****************************************************
        '''
        ident = round((float(hsp.identities)/len(hsp.match))*100)

        '''
        substrate_info= self.substrates.get_tcid_substrates(tcid)
        mysubstrate_html = ''
        mysubstrate_tsv = ''
        
        if substrate_info is not None:
            
            category,gen_sub,spec_sub,chebi_id,status = substrate_info
            
            chebi_info = chebi_id.replace('"','').replace(' ','').split(',')
            chebi = ''

            for i in chebi_info:
            
                chebi += ' <a href="https://www.ebi.ac.uk/chebi/searchId.do;jsessionid=A9D16DCB24C6F74339FC28A4941EFBB6?chebiId=CHEBI:{}">{}</a>'.format(i,i)
        
            mysubstrate_html = '{},{},{},{},{}'.format(category,gen_sub,spec_sub,chebi,status)
        
            mysubstrate_tsv = ",".join(substrate_info)
        
        else:
        
            mysubstrate_html = 'None'
            mysubstrate_tsv = 'None'
        '''

        substrate = ''

        try:

            substrate = self.substrates[tcid]

        except:

            substrate = 'None'
            
        ##pfamDomain = list(filter(lambda x: (x["query name"].strip()==genome.strip()), self.chart3))[0]["accession"]
        try:
            if genome.strip() in self.chart3:
                _find_chart3=self.chart3[genome.strip()]
                pfamDomain = ",".join([x.split('.')[0] for x in _find_chart3[0]["accAll"]])
                self.chart3_find[genome.strip()]=_find_chart3
                del self.chart3[genome.strip()]
            elif genome.strip() in self.chart3_find:
                _find_chart3=self.chart3_find[genome.strip()]
            pfamDomain = ",".join([x.split('.')[0] for x in _find_chart3[0]["accAll"]])
        except Exception as e:
            print(genome,e)
            self.chart3_find[genome.strip()]={}
            _find_chart3=False
            pfamDomain=f"NA"
            
        
        try:
            if acc.strip() in self.chart2:
                _find_chart2=self.chart2[acc.strip()]
                SbjpfamDomain = ",".join([x.split('.')[0] for x in _find_chart2[0]["accAll"]])
                self.chart2_find[acc.strip()]=_find_chart2
                del self.chart2[acc.strip()]
            elif acc.strip() in self.chart2_find:
                _find_chart2=self.chart2_find[acc.strip()]
            SbjpfamDomain = ",".join([x.split('.')[0] for x in _find_chart2[0]["accAll"]])
        except Exception as e:
            print(acc,e)
            #import traceback
            #print(traceback.format_exc())
            self.chart2_find[acc.strip()]={}
            _find_chart2=False
            SbjpfamDomain=f"NA"
        '''
        html_results (VI)
        '''

        #glink = '<a href="content.html#%s">%s</a>'%(genome,genome)
        #tclink = '<a href="http://tcdb.org/search/result.php?tc={}">{}</a>'.format(tcid,tcid)
        glink = '<a href="htmls/%s.html">%s</a>'%(genome,genome)
        tclink = '<a href="http://tcdb.org/search/result.php?tc={}">{}</a>'.format(tcid,tcid)
        #h_row = (glink,acc,tclink,tcdb.title,len(hsp.match),tcdb.e,ident,q_len,h_len,qcov,hcov,query_tms,\
        #hit_tms,overlap,family,substrate,self.globalcount,pfamDomain,SbjpfamDomain)
        h_row = (glink,acc,tclink,tcdb.title,len(hsp.match),tcdb.e,ident,q_len,h_len,hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,qcov,hcov,query_tms,\
        hit_tms,overlap,family,substrate,pfamDomain,SbjpfamDomain)
        h_row =['<td>'+str(i)+'</td>' for i in h_row]
        htmlrow = "<tr>\n%s\n</tr>"%("\n\t".join(h_row))


        '''
        text results (VI)
        '''

        #row = (genome,acc,tcid,tcdb.title,len(hsp.match),tcdb.e,ident,q_len,h_len,qcov,hcov,query_tms,hit_tms,overlap,family,substrate,self.globalcount,pfamDomain,SbjpfamDomain)
        row = (genome,acc,tcid,tcdb.title,len(hsp.match),tcdb.e,ident,q_len,h_len,hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,qcov,hcov,query_tms,hit_tms,overlap,family,substrate,pfamDomain,SbjpfamDomain)
        row =[str(i) for i in row]
        ##print(row)
        txtrow = "%s\n"%("\t".join(row))


        '''
        #Fig result
        '''
        import re
        
        #hsp.query= hsp.query.replace("-","")
        query=ParseDefline(genome).id
        hitID = ParseDefline(tcdb.title,True).id
        #hsp.sbj= hsp.sbjct.replace("-","")
        hsp.sbj= hsp.sbjct
        seq=self.queries[str(query)].seq
        sbj=self.tcdbHits[str(hitID)].seq


        """figFullProtein(seq=seq,seqname=genome,savefig=f"out_test/img/{genome}_qry.png",chart3=_find_chart3)
        figFullProteinSbj(seq=sbj,seqname=acc,savefig=f"out_test/img/{acc}_sbj.png",chart2=_find_chart2)
        figHydroAlignedRegions(seq=seq,seqfrag=hsp.query,sbj=sbj,sbjfrag=hsp.sbj,seqname=genome,sbjname=acc,savefig=f"out_test/img/{acc}-{genome}_hydro.png")"""
        figFullProtein(seq=seq,seqname=genome,savefig=os.path.join(self.indir,f"img/{genome}-{acc}_qry.png"),chart3=_find_chart3,startend=[hsp.query_start,hsp.query_end])
        figFullProteinSbj(seq=sbj,seqname=acc,savefig=os.path.join(self.indir,f"img/{genome}-{acc}_sbj.png"),chart2=_find_chart2,startendsbj=[hsp.sbjct_start,hsp.sbjct_end])
        figHydroAlignedRegions(seq=seq,seqfrag=hsp.query,sbj=sbj,sbjfrag=hsp.sbj,seqname=genome,sbjname=acc,savefig=os.path.join(self.indir,f"img/{genome}-{acc}_hydro.png"))
        

        '''
        content results (VI)
        '''

        #self.plotHydro(genome,tcdb,acc,hsp)

        #self.rows.append(htmlrow)
        #adjusting the hsp
        import re

        hspall=str(hsp)

        hspall=hspall.split("\n")
        hspall[0]= hspall[0].replace("expectation","E-value:").replace("alignment length","Alignment length:").replace("Score","Score:").replace("(",", Bit scores: ").replace(")", "")
        hspall[0]= hspall[0]+ f", Percentage identity: {ident}"
        hspall[1]=re.sub(r"[A-Z\-\.]{5,}",str(hsp.query),hspall[1])
        index=hspall[1].find(str(hsp.query))
        hspall[2]=hspall[2][:index]+ str(hsp.match)
        hspall[3]=re.sub(r"[A-Z\-\.]{5,}",str(hsp.sbjct),hspall[3])
        hspall="\n".join(hspall)
        if self.cdd_on is True:
            mycdd = 'Query & TC-Hit Conserved Domains:<br><img src=\'cdd/%s.1.png\'><br><img src=\'cdd/%s.2.png\'><br>'%(genome,genome)
        else:
            mycdd = ''
        ol = float(overlap) if overlap is not None else float(0)
        content = '''<div class='result' id='%s'> <h3><a name='%s'>%s</a></h3>  <p>Hit Accession: %s<br>   Hit TCID: %s</p> <p>Hit Description: %s<br>
        <br>   Mach Len: %i<br>   e:%f</p> <p>Query TMS Count : %i<br>   Hit TMS Count: %i     <br>     TMS-Overlap Score: %f<br>
        Predicted Substrates:%s <br><br>     BLAST Alignment:<br>     <pre style= "border: 2px solid black;height: 100px;width:  100%%;overflow-x: auto            ;overflow-y: hidden;margin: 1em 0;background: gray;color: white;">    %s     </pre> <br>  <table><tr><th>Protein Hydropathy Plots:</th></tr> 
        <tr><td><img src='../img/%s_qry.png'></td> <td><img src='../img/%s_sbj.png'></td></tr><br>
        <tr><th><br> Pairwise Alignment-Hydropathy Plot:<br></th></tr>
        <tr><td colspan="2" style="text-align: center;"><img src='../img/%s-%s_hydro.png'></td></tr></table><br>%s </p> </div>\n''' %(genome,genome,genome,acc,tcid,tcdb.title,\
        len(hsp.match),tcdb.e,query_tms,hit_tms,ol,substrate,str(hspall),f'{genome}-{acc}',f'{genome}-{acc}',genome,acc,mycdd)
        #print (str(hspall))
        #self.data.append(content)     71 IKLAIWLIGFTFIFDAVNHEIKINLLKF-APLIRNI---SVIAIV...VVF 351
        #print("MATCH",hsp.match)
        #print("SBJCT",hsp.sbjct)
        return htmlrow,content,txtrow

    def write_results(self):
        bestresults = copy.deepcopy(self.bestresults)
        bestresults.extend(self.notmsresults)
        
        #create and open the necessary files for results (VI)
        html_results = open(self.indir+'/results.html','wt+')
        content = open(self.indir+'/content.html','wt+')
        tsv_results = open(self.indir+'/results.tsv', 'wt+')
        
        '''
            
        We are unsure of why the local number variable( the one that is supposed to count which result it is) is the last number in the rows.  --Vasu Pranav Sai Iddamsetty
        '''
        
        
        #write to results.html and content.html (VI)
        html = '''<html><table width="100%%" border="1"> <tr> <td>Query ID</td> <td>Hit ID</td>
        <td>Hit TCID</td> <td>Hit Description</td> <td>Match Len</td> <td>e-Val</td> <td>% Identity</td> <td>Query Length</td> <td>Hit Length</td> <td>Query Start</td><td>Query End</td><td>Subject Start</td><td>Subject End</td> <td>Query Coverage</td> <td>Hit Coverage</td> <td>Query TMS#</td>
        <td>Hit TMS#</td> <td>TM-Overlap Score</td> <td>Family Abrv.</td><td>Predicted Substrate</td> <td>Query Pfam</td><td>Subject Pfam</td></tr><br>\n\n'''
        html_results.write(html)
        content.write("<html>")
        
        
        #write the column descriptions to results.tsv (VI)
        columnDescription = '''#Query_id\tHit_xid\tHit_tcid\tHit_desc\tMatch_length\te-value\t%_identity\tQuery_Length\tHit_Length\tQ_start\tQ_end\tS_start\tS_end\tQuery_Coverage\tHit_Coverage\tQuery_n_TMS\tHit_n_TMS\tTM_Overlap_Score\tFamily_Abrv\tPredicted_Substrate\tQuery_Pfam\tSubject_Pfam\n'''
        tsv_results.write(columnDescription)
        
        
        
        self.myqueries = SeqIO.parse(self.indir+'/myqueries.faa','fasta')
        self.mytcdb    = SeqIO.parse(self.indir+'/mytcdb.faa','fasta')
        self.myqueries = SeqIO.to_dict(self.myqueries)
        self.mytcdb    = SeqIO.to_dict(self.mytcdb)
        if self.cdd_on:
            self.cdd_extract()
        self.queryhmg = hmmgap.annotate()
        self.tcdbhmg  = hmmgap.annotate()
        self.queryhmg.hmmtop = self.tms['queries']
        self.tcdbhmg.hmmtop  = self.tms['tcdb']
        for res in bestresults:
            #retrieve relevant formaatted data and write it to the files (VI)
            (row,data,txt) = self.build_view(res)
            ##print(row)
            ##print(txt)
            html_results.write(row)
            #content.write(data)
            (genome,tcdb,hsp,q_len,h_len,qcov,hcov,overlap) = res
            filepath= os.path.join(self.indir,"htmls/%s.html"%(genome))
            filedata="<html>%s</html>"%(data)
            f = open(filepath,"w")
            f.write(filedata)
            f.close()
            tsv_results.write(txt)
            print(("Generated Results for :: %s" %ParseDefline(res[0]).id))
        
        #end tags for the .html files
        html_results.write("</table></html>")
        content.write("</html>")

    def loadSubstrates(self):

        print('Loading Substrates')

        #substrateData = urlopen('http://tcdb.org/cgi-bin/projectv/getSubstrates.py')
        #substrateData = urlopen('https://tcdb.org/cgi-bin/substrates/getSubstrates.py')
        substrateData =  requests.get('https://tcdb.org/cgi-bin/substrates/getSubstrates.py', allow_redirects=True).text.strip().split("\n")
        
        for line in substrateData:

            data = line.replace('\n','').split('\t')
            #print(data)
            self.substrates[data[0]] = data[1].replace('|',', ')


        return
                         

    def cdd_extract(self):
        if os.path.exists(self.indir+'/cdd') is False:
            os.mkdir(self.indir+'/cdd')
        fastas  = dict(self.myqueries,**self.mytcdb)
        thisdir = self.indir+'/cdd/'
        for query,hit,hsp in self.goodresults:
            if os.path.exists(thisdir+query.title()+'.1.png') is False:
                cdd.fetch(str(hsp.query),thisdir+query.title()+'.1.png')
            if os.path.exists(thisdir+query.title()+'.2.png') is False:
                cdd.fetch(str(hsp.sbjct),thisdir+query.title()+'.2.png')

    '''
    Plots individual Protein hydropathies (VI)
    '''
 
    def plotHydro(self,genome,hit,acc,hsp):

        #File Paths

        #queryPath = self.indir+'/img/'+genome+'_qry.png'
        queryPath = self.indir+'/img/'+f'{genome}-{acc}_qry.png'
        
        #hitPath = self.indir+'/img/'+acc+'-'+genome+'_hydro.png'
        hitPath = self.indir+'/img/'+f'{genome}-{acc}_hydro.png'

        quod = ' -s -q -d {} --width 10 --height 3 --xticks 50'.format(self.indir+"/img/")


        #Query Hydropathy
        if not os.path.exists(queryPath):

            query=ParseDefline(genome).id
            querySeq = self.queries[str(query)].seq
            query = 'quod.py {} -o {} -c blue -w {}-{} -l {}'.format(querySeq,genome+'_hydro',hsp.query_start,hsp.query_end,genome) + quod

            #os.system(query)
            subprocess.call(query.split())

        #Hit Hydropathy
        if not os.path.exists(hitPath):
        
            hitID = ParseDefline(hit.title,True).id
            hitSeq = self.tcdbHits[str(hitID)].seq
            hit = 'quod.py {} -o {} -c red -w {}-{} -l {}'.format(hitSeq, acc+'-'+genome+'_hydro',hsp.sbjct_start,hsp.sbjct_end,acc) + quod

            #os.system(hit)
            subprocess.call(hit.split())

    '''
    Within the function hydro, we added a (0,0) value for the '*' symbols, so that 
    hydropathies can still be visualized for pseudogenes. - Vasu Pranav Sai Iddamsetty
    '''

    def hydro(self,gseq):
        seq=gseq.replace('-','')
        window = 19
        prev = 0
        index = {'G':(-0.400,0.48),
                 'I':(4.500,1.38),
                 'S':(-0.800,-0.18),
                 'Q':(-3.500,-0.85),
                 'E':(-3.500,-0.74),
                 'A':(1.800,0.62),
                 'M':(1.900,0.64),
                 'T':(-0.700,-0.05),
                 'Y':(-1.300,0.26),
                 'H':(-3.200,-0.4),
                 'V':(4.200,1.08),
                 'F':(2.800,1.19),
                 'C':(2.500,0.29),
                 'W':(-0.900,0.81),
                 'K':(-3.900,-1.5),
                 'L':(3.800,1.06),
                 'P':(-1.600,0.12),
                 'N':(-3.500,-0.78),
                 'D':(-3.500,-0.90),
                 'R':(-4.500,-2.53),
                 'X':(0,0),
                 'U':(0,0),
                 '*':(0,0)}
        midpt = (window+1)/2
        length = len(seq)
        hydro = []
        for i in range(length-window+1):
            total = 0
            for j in range(window):
                total +=index[seq[i+j]][0]
            total = total/window
            hydro.append(total)
        if len(seq) == len(gseq):
            return hydro
        replace = re.finditer('(-+)',gseq)
        inserts = {}
        for i in replace:
            inserts.setdefault(i.start(),i.end()-i.start())
        first = False
        newhydro = []
        for x, h in enumerate(hydro):
            if x in list(inserts.keys()) and first is False:
                first = True
                for y in range(inserts[x]):
                    newhydro.append(0)
                newcount = x + inserts[x]
                continue
            if first is False:
                newhydro.append(h)
                continue
            if first is True and newcount in list(inserts.keys()):
                for y in range(inserts[newcount]):
                    newhydro.append(0)
                newcount += inserts[newcount]
                continue
            else:
                newhydro.append(h)
                newcount +=1

        return newhydro


    def what(self,a,b,outfile,hmt=[]):
        #quit()
        ha = self.hydro(a)
        hb = self.hydro(b)
        omg = [len(ha),len(hb)]
        readbar = re.compile(r'\d+?[+_]+\+')
        omg.sort()
        ha =ha[0:omg[0]]
        hb =hb[0:omg[0]]
        x_data=list(range(0,len(ha)))
        plt.figure()
        plt.axhline(y=0,color='black')
        plt.ylim(-3,3)
        plt.xlim(right=len(a))
        plt.plot(x_data,ha,linewidth=1,label=self.alabel,color='blue')
        plt.plot(x_data,hb,linewidth=1,label=self.blabel,color='red')
        plt.xlabel("Residue")
        plt.ylabel("Hydropathy (kcal/mol)")
        plt.legend(loc='lower right')
        # Draw TMS bars
        if len(hmt) == 2:
            sub = readbar.finditer(str(hmt[0]))
            tar = readbar.finditer(str(hmt[1]))
            for tms in sub:
                subtract = (len(hmt[0])-len(ha))/2
                if tms.end()-subtract>len(ha):
                    subtract = tms.end()-len(ha)
                plt.axvspan(tms.start()-subtract,tms.end()-subtract, facecolor="blue", alpha=0.3)
            for tms in tar:
                subtract = (len(hmt[1])-len(ha))/2
                if tms.end()-subtract>len(hb):
                    subtract = tms.end()-len(hb)
                plt.axvspan(tms.start()-subtract,tms.end()-subtract, facecolor="red", alpha=0.3)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(10,3)
        plt.savefig(outfile, dpi=100, format="png",bbox_inches='tight', pad_inches=0.003)
        plt.clf()
        plt.close(fig)



if __name__=="__main__":
    from optparse import OptionParser,OptionGroup
    desc = "Welcome to GBlast! Easily identify transporters in entire Genomes/Proteomes - By Vamsee Reddy"
    version = "GBlast V3"
    opts = OptionParser(description=desc,version=version)
    opts.add_option('-i',
                    action='store',
                    type='string',
                    dest='input',
                    help="Path to genome/proteome file"
    )
    opts.add_option('-o',
                    action='store',
                    type='string',
                    dest='output',default='genome.fsa',
                    help="Results output directory"
    )
    opts.add_option('--evalue',
                    action='store',
                    type='float',
                    dest='evalue',
                    default=0.001,
                    help="Minimum e-Value [0.001]"
    )
    opts.add_option('--cov',
                    action='store',
                    type='float',
                    dest='mincov',
                    default=40.0,
                    help="Minimum Proten Coverage [40.0]"
    )
    opts.add_option('--esort',
                    action='store_true',
                    dest='esort',
                    default=False,
                    help="Use e-value as the preliminary criteria for sorting (otherwise sorted by TCID)"
    )
    '''    
    opts.add_option('--betabarrel',
                    action='store_true',
                    dest='bb',
                    default=False,
                    help="Find Beta Barrels instead of TMS"
    )
    opts.add_option('--cdd',
                    action='store_true',
                    default=False,
                    dest='cdd',
                    help="Include CDD analysis (Takes a while)"
    )
    opts.add_option('--orthologs',
                    action='store',
                    dest='ortho',
                    default=False,
                    help="Find orthologs with this fasta file (optional)"
    )
    opts.add_option('--query_gi',
                    action='store',
                    dest='subgi',
                    default=False,
                    help="Ortholog search using list of only these queries"
    )
    opts.add_option('--target_gi',
                    action='store',
                    dest='targi',
                    default=False,
                    help="Ortholog search using list of only these targets"
    )
    '''
    (cli,args)=opts.parse_args()
    input_path= cli.input
    cli.input=extract_archive(from_path=cli.input).strip()
    if cli.input != "" and cli.output != None:
        GB = Tools()
        #if(cli.bb):
            #GB.dbfile = '/db/betabarrel'
        GB.ortho  = False
        GB.indir  = cli.output
        GB.cdd_on = False
        GB.query  = cli.input
        GB.seqParse=SeqIO.parse(io.StringIO(GB.query),'fasta')
        GB.query_gis  = False
        GB.target_gis = False
        GB.expect = cli.evalue
        GB.mincov = cli.mincov
        GB.esort=cli.esort
        GB.blast_all()
        print("Loading BLAST Results")
        GB.load_good()
        print("Writing FASTAS")
        GB.write_fastas()
        print("Running HMMTOP")
        GB.hmmtop()
        print("Calculating TMS Overlap")
        GB.calculate_tms_scores()
        print("Comparing Against Pfam")
        char3_final=GB.generate_chart3()
        print("Generating Pfam Results")
        tcdb_final=GB.generate_tcdb()
        print("Writing Results")
        GB.write_results()
        chart3_find=GB.chart3_find
        chart2_find=GB.chart2_find
    else:
        opts.print_help()
        
#path="/System/Volumes/Data/tcdbPegasus/ResearchData/Users/clark/gblast3/out_test/htmls/"
path=os.path.join(GB.indir,"htmls")
dirs=os.listdir(path)
# export file
dirs = list(filter(lambda x: ("WP"== x[:2] and ".html" == x[-5:]), dirs))
###print(dirs)
def dict2htmltable(data):
    html = ''.join('<th>' + x + '</th>' for x in data[0].keys())
    for d in data:
        html += '<tr>' + ''.join('<td>' + x + '</td>' for x in d.values()) + '</tr>'
    return "<div><table class='dom' border='1', style='width:100%'>" + html + '</table>'+ '</div>'
html_style="""
    <style type="text/css">
.dom {
   border: 2px solid black;
   height: 100px;
   width:  100%;
   overflow-x: auto;
   overflow-y: auto;
   margin: 1em 0;
   background: gray;
   color: white;
}
</style>
"""
for q in chart3_find:
    final=[]
    filename=q +".html"
    filepath=os.path.join(path,filename)        
    if os.path.exists(filepath):
        for chart3_finatmp_dict in chart3_find[q]:
            if chart3_find[q]!={}:
                getClanBuff                     = os.popen(f"getClan.sh {chart3_finatmp_dict['accession'].split('.')[0]}")
                getClan                         = getClanBuff.buffer.read().decode('gbk').strip().split("\t")
                chart3_finatmp_dict["Clan"]     = getClan[1]
                chart3_finatmp_dict["Dom_name"] = getClan[3]
                chart3_finatmp_dict["Dom_info"] = getClan[4]

                chart3_fina_dict={}
                chart3_sort_dict                = {"query name":"Query","accession":"Domain","Clan":"Clan","qlen":"Dom_Length","E-value":"E-value","from":"Dom_start","to":"Dom_end","from_1":"Q_start",
                                  "to_1":"Q_end","Dom_name":"Dom_name","Dom_info":"Dom_info"}
                for key in chart3_sort_dict:
                    if key in chart3_finatmp_dict:
                        chart3_fina_dict[chart3_sort_dict[key]]=chart3_finatmp_dict[key]      
                chart3_fina_dict["Domain"]=chart3_fina_dict["Domain"].split('.')[0]                               
                final.append(chart3_fina_dict)
            else:
                pass
                #chart3_finatmp_dict= {"query name":q,"accession":"","Clan":"","qlen":"","E-value":"","from":"","to":"","from_1":"",
                #                  "to_1":"","Dom_name":"","Dom_info":""}

        with open(filepath,"r") as f:
            regex                        = r"Hit Accession:\s([a-zA-Z\d\_]+)"
            tcdb_query_name               = re.search(regex, f.read(), re.M|re.I).group(1)
            for tcdb_finatmp_dict in chart2_find[tcdb_query_name]:
                if tcdb_finatmp_dict !={}:
                    getClanBuff                   = os.popen(f"getClan.sh {tcdb_finatmp_dict['accession'].split('.')[0]}")
                    getClan                       = getClanBuff.buffer.read().decode('gbk').strip().split("\t")
                    tcdb_finatmp_dict["Clan"]     = getClan[1]
                    tcdb_finatmp_dict["Dom_name"] = getClan[3]
                    tcdb_finatmp_dict["Dom_info"] = getClan[4]

                    tcdb_sort_dict                = {"query name":"Query","accession":"Domain","Clan":"Clan","qlen":"Dom_Length","E-value":"E-value","from":"Dom_start","to":"Dom_end","from_1":"Q_start",
                                      "to_1":"Q_end","Dom_name":"Dom_name","Dom_info":"Dom_info"}
                    tcdb_fina_dict={}
                    for key in tcdb_sort_dict:
                        if key in tcdb_finatmp_dict:
                            tcdb_fina_dict[tcdb_sort_dict[key]]=tcdb_finatmp_dict[key]
                    tcdb_fina_dict["Domain"]=tcdb_fina_dict["Domain"].split('.')[0]
                    final.append(tcdb_fina_dict)
                else:
                    pass
            # tcdb and getClan
            #print(f"getClan.sh {tcdb_finatmp_dict['accession'].split('.')[0]}")

    if final!=[]:
        html = html_style
        html += dict2htmltable(final)
        #print(html)
        with open(filepath,"a") as f:
            f.write("\n")
            f.write(html)
            f.write("\n")
