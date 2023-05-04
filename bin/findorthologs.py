#!/usr/bin/env python

from ProjectBio import ParseDefline, nr_dict
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast import NCBIXML
import os,re,pickle

class Ortho:

    def __init__(self):
        # Get input Files
        self.genomes = []
        self.standard = None
        self.identity = 0.95
        self.genome_keys = {}
        self.ignore_keys = []
        self.orthologs = [] # temp cirlce
        self.inner_cirlce = []
        self.outer_circle = []
        self.standard_keys = []
        try:
            #os.mkdir('circle')
            os.mkdir('in')
            os.mkdir('out')
        except:
            pass

    def __call__(self):
        #self.build_db(self.standard)
        #self.blast_all()
        #self.build_in()
        #self.build_genome('out','13.faa')
        self.build_genome_key()
        #self.circle_blast('out')
        #self.circle_blast('.')
        self.build_standard()
        self.mycircle = self.relate_circles('./')
        self.tsv_results()


    def build_db(self,db):
        os.system('makeblastdb -dbtype prot -in '+db)

    def do_blast(self,subject):
        print 'Running BLAST...'
        blast = blastp()
        blast.db = self.standard
        blast.query = subject
        blast.outfmt = 5
        blast.out = subject+'.xml'
        blast.comp_based_stats = "0"
        if os.path.exists(blast.out):
            return
        blast()
        # Blast is complete, now parse the results. I hate my job

    def blast_all(self):
        for genome in self.genomes:
            self.do_blast(genome)

    def build_in(self):
        for genome in self.genomes:
            xmlfile = genome+'.xml'
            results = NCBIXML.parse(open(xmlfile))
            subjects =  SeqIO.parse(open(genome),'fasta')
            subjects = SeqIO.to_dict(subjects)
            inn = []
            for result in results:
                title = result.query
                title = ParseDefline(title,False)
                title = title.id
                for aln in result.alignments:
                    hsps = list(aln.hsps)
                    hsps.sort(key=lambda x:self.hspsort(x),reverse=True)
                    i = self.hspsort(hsps[0])
                    if i >= self.identity:
                        inn.append(title)
            out = list(set(subjects.keys())-set(inn))
            infile = open('in/'+genome,'wb')
            outfile = open('out/'+genome,'wb')
            for seq in inn:
                SeqIO.write(subjects[seq],infile,'fasta')
            for seq in out:
                SeqIO.write(subjects[seq],outfile,'fasta')


    def hspsort(self,hsp):
        i = float(hsp.identities)/len(hsp.match)
        return i

    def build_genome_key(self):
        genomes=self.genomes[:]
        genomes.append(self.standard)
        for g in genomes:
            gen = SeqIO.parse(open(g),'fasta')
            for fasta in gen:
                if fasta.id not in self.genome_keys.keys():
                    self.genome_keys[fasta.id] = g

    def build_genome(self,loc,subject):
        # build genome excluding a subject
        genomes = []
        fastadb = []
        for i in self.genomes:
            if i != subject:
                genomes.append(i)
        for genome in genomes:
            handle = open(loc+'/'+genome)
            fastas = SeqIO.parse(handle,'fasta')
            for fasta in fastas:
                fastadb.append(fasta)
        circle = open(loc+'/'+'circle','wb')
        for seq in fastadb:
            SeqIO.write(seq,circle,'fasta')
        self.build_db(loc+'/circle')

    def circle_blast(self,loc):
        print 'Running Circualr BLAST...'
        blast = blastp()
        blast.db = loc+'/circle'
        blast.outfmt = 5
        blast.comp_based_stats = "0"
        for subject in self.genomes:
            self.build_genome(loc,subject)
            blast.query = subject
            blast.out = loc+'/'+subject+'.xml'
            if os.path.exists(blast.out):
                return
            blast()

    def relate(self,loc,subject):
        res = open(loc+'/'+subject+'.xml')
        results = NCBIXML.parse(res)
        for result in results:
            genhits = []
            orthologs = []
            title = result.query
            title = ParseDefline(title)
            query = title.id
            if query in self.ignore_keys:
                continue
            self.ignore_keys.append(query)
            orthologs.append(query)
            for aln in result.alignments:
                hsps = list(aln.hsps)
                title = aln.title
                title = ParseDefline(title,'True')
                hit = title.id
                if hit in self.ignore_keys:
                    continue
                hsps.sort(key=lambda x:self.hspsort(x), reverse=True)
                if self.hspsort(hsps[0]) >= self.identity:
                    # Found an ortholog, check with genhits.
                    genome=self.genome_keys[hit]
                    if genome in genhits:
                        continue
                    genhits.append(genome)
                    orthologs.append(hit)
                    self.ignore_keys.append(hit)
            if bool(len(orthologs)):
                self.orthologs.append(orthologs)

    def relate_circles(self,loc):
        self.orthologs = []
        if os.path.exists(loc+'/circle_relations'):
            self.orthologs = pickle.load(open(loc+'/circle_relations'))
            return self.orthologs
        for genome in self.genomes:
            self.relate(loc,genome)
        handle = open(loc+'/circle_relations','wb')
        pickle.dump(self.orthologs,handle)


    def build_standard(self):
        for g in self.genomes:
            handle = open('./in/'+g)
            fastas = SeqIO.parse(handle,'fasta')
            fastas = nr_dict(fastas)
            keys = fastas.keys()
            self.standard_keys.extend(keys)

    def tsv_results(self):
        self.mycircle.sort(key=lambda x:len(x),reverse=True)
        for group in self.mycircle:
            for item in group:
                genome = self.genome_keys[item]
                standard = 'Yes' if item in self.standard_keys else 'No'
                print '%s\t%s\t%s'%(item,genome,standard)
            print '--------------------------------'



if __name__=='__main__':

    ortho = Ortho()
    ortho.genomes = ['13.faa','14.faa','15.faa','16.faa','17.faa','18.faa','19.faa','20.faa']
    ortho.standard = '12.faa'
    ortho()
