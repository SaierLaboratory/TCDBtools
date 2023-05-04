#!/usr/bin/env python

# This tool includes PsiBLAST options because NCBIWWW does not support it.
# Uses a modified version of Protocol1.PHP

# Allows user to merge multiple blast sessions into one xml -> faa file
# Entrez is hit in groups of 100 to avoid HTTP GET maxouts
# Class will be expanded to include BlastP options eventually

# Created by Vamsee Reddy
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import subprocess
import protocol1
import tempfile
import hashlib
import shutil
import urllib2
import re,sys

class tools:

    def __init__(self):
        self.blast_in = []
        self.gis = []
        self.blast_threshold= 0.005
        self.status = {}
        self.max_gi = 500
        self.fasta_file = False
        self.email = "vreddy@ucalgary.ca"

    def psiblast(self,acc):
        psi = protocol1.psi()
        psi.query = acc
        psi.init_blast()
        print psi.gis
        unique = [i for i in psi.gis if i not in self.gis]
        self.gis.extend(unique)

    def blastp(self,acc):
        try:
            gis = []
            print 'here'
            result_handle = NCBIWWW.qblast("blastp", "nr", acc, format_type="XML", expect=self.blast_threshold)
            print 'here'
            for blast_record in NCBIXML.parse(result_handle):
                for alignment in blast_record.alignments:
                    gis.append(alignment.title.split("|")[1])
            unique = [int(i.strip()) for i in gis if int(i) not in self.gis]
            self.gis.extend(unique)
        except:
            self.status.setdefault(acc,False)
        return

    def load_fasta(self,path):
        handle = open(path,'r')
        fastas = SeqIO.parse(handle,'fasta')
        gis = []
        gi = re.compile('gi\|(\d+)')
        for fasta in fastas:
            gis.append(gi.search(fasta.description).groups()[0])
        self.gis.extend(gis)
        self.gis = list(set(self.gis))

    def build_xml(self):
        self.xml_file = tempfile.NamedTemporaryFile()
        # Break the GI into groups of 100, so we don't exceed HTTP_GET limits.
        if len(self.gis) > 100:
            gi_groups = [self.gis[i:i+50] for i in range(0,len(self.gis),50)] 
        else:
            gi_groups = [self.gis]
        Entrez.email = self.email
        self.xml_file.write("<TSeqSet>\n")
        for group in gi_groups:
            try:
                handle = Entrez.efetch("protein", id=','.join([str(i) for i in group]), rettype='fasta',retmode="xml")
            except:
                print group
                quit()
            for seq_record in Entrez.parse(handle):
                self.xml_file.write("\t<TSeq>\n")
                self.xml_file.write("\t\t<TSeq_seqtype value=\"protein\">protein</TSeq_seqtype>\n")
                if 'TSeq_gi' in seq_record:
                    self.xml_file.write("\t\t<TSeq_gi>" + seq_record['TSeq_gi'] + "</TSeq_gi>\n")
                if 'TSeq_accver' in seq_record:
                    self.xml_file.write("\t\t<TSeq_accver>" + seq_record['TSeq_accver'] + "</TSeq_accver>\n")
                if 'TSeq_taxid' in seq_record:
                    self.xml_file.write("\t\t<TSeq_taxid>" + seq_record['TSeq_taxid'] + "</TSeq_taxid>\n")
                if 'TSeq_orgname' in seq_record:
                    self.xml_file.write("\t\t<TSeq_orgname>" + seq_record['TSeq_orgname'] + "</TSeq_orgname>\n")
                if 'TSeq_defline' in seq_record:
                    self.xml_file.write("\t\t<TSeq_defline>" + seq_record['TSeq_defline'] + "</TSeq_defline>\n")
                if 'TSeq_length' in seq_record:
                    self.xml_file.write("\t\t<TSeq_length>" + seq_record['TSeq_length'] + "</TSeq_length>\n")
                if 'TSeq_sequence' in seq_record:
                    self.xml_file.write("\t\t<TSeq_sequence>" + seq_record['TSeq_sequence'] + "</TSeq_sequence>\n")
                self.xml_file.write("\t</TSeq>\n")
        self.xml_file.write("</TSeqSet>\n")
        self.status['xml']= True
        return

    def build_raw_fasta(self,desc=None):
        self.raw_fasta = tempfile.NamedTemporaryFile(delete=False)
        # Break the GI into groups of 100, so we don't exceed HTTP_GET limits.
        if len(self.gis) > 100:
            gi_groups = [self.gis[i:i+100] for i in range(0,len(self.gis),100)] 
        else:
            gi_groups = [self.gis]
        Entrez.email = self.email
        for group in gi_groups:
            handle = Entrez.efetch("protein", id=','.join([str(i) for i in group]), rettype='fasta',retmode="xml")
            for seq_record in Entrez.parse(handle):
                mydesc = seq_record['TSeq_defline'] if desc is None else desc
                record = SeqRecord(Seq(seq_record['TSeq_sequence'],IUPAC.protein),name=seq_record['TSeq_gi'],\
                id=seq_record['TSeq_gi'], description=mydesc)
                #ii +=1
                SeqIO.write(record,self.raw_fasta,'fasta')
        self.status['raw']= True
        self.raw_fasta.seek(0)
        self.raw_fasta.flush()
        return self.raw_fasta

    def maketable(self,cutoff):
        cmd = 'make_table5.pl -i "%s" -s 1 -l 999999 -c %f' %(self.xml_file.name, cutoff)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout,stderr = process.communicate()
        self.fasta_file = open(self.xml_file.name +'.faa')
        self.status['maketable'] = True
        return


if __name__=='__main__':
    b = tools()
    try:
        mode = sys.argv[1]
        acc = sys.argv[2]
        c = float(sys.argv[3])
        out = sys.argv[4]
        if mode == 'p':
            b.blastp(acc)
        if mode == 'psi':
            b.psiblast(acc)
        b.build_xml()
        b.maketable(c)
        shutil.copy(b.fasta_file.name,out)
    except:
        print "BLAST Command Line Tool Usage:"
        print "blast.py <p/psi> <acc> <cd-hit-thresh> <outfile>"

