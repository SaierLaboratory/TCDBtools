#!/usr/bin/env python
# Generates orthologs from FASTA list. I don't remmber how this works -V
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio import SeqIO
from ProjectBio import ParseDefline
import os
import re
import pickle

class Ortho:

    def __init__(self):
        # Load input files
        self.genomes = []
        self.subtract = None
        # Define outputs
        self.db = 'db.faa'
        # Blast settings
        self.expect = 1e-20
        self.expect_diff = 1e-20
        # Containers
        self.orthologs = []

    def build_db(self):
        if os.path.exists(self.db):
            return
        merged = []
        handle = open(self.db,'wb')
        remove = r'^[A-Za-z0-9.]'
        for org,db in self.genomes:
            records = SeqIO.parse(db,'fasta')
            for record in records:
                record.id = '%s_%s'%(org,record.id)
                #record.name = re.sub(remove,'',record.name)
                merged.append(record)
        for record in merged:
            SeqIO.write(record,handle,'fasta')
        handle.flush()
        os.system('makeblastdb -dbtype prot -in '+handle.name)

    def primary_blast(self):
        if os.path.exists('primary.xml') is True:
            return
        print 'Running BLAST...'
        blast = blastp()
        blast.db = self.db
        blast.query = self.genomes[0][1]
        blast.outfmt = 5
        blast.out = 'primary.xml'
        blast.evalue = self.expect
        blast.comp_based_stats = "0"
        blast()
        # Blast is complete, now parse the results

    def secondary_blast(self):
        if os.path.exists('orthologs.faa') is False or os.path.exists('secondary.xml') is True:
            return
        os.system('makeblastdb -dbtype prot -in orthologs.faa')
        blast = blastp()
        blast.db = 'orthologs.faa'
        blast.query = self.substract
        blast.outfmt = 5
        blast.out = 'secondary.xml'
        blast.evalue = self.expect_diff
        blast.comp_based_stats = "0"
        blast()

    def parse_blast(self):
        print 'Locating Orthologs...'
        if os.path.exists('parsed'):
            self.orthologs = pickle.load(open('parsed'))
            return
        results = NCBIXML.parse(open('primary.xml'))
        for row in results:
            group = {}
            query = row.query
            for desc in row.descriptions:
                title = ParseDefline(desc.title,True)
                e = desc.e
                org = title.id.split('_')[0]
                if org != self.genomes[0][0] and e <= self.expect:
                    group.setdefault(org,[]).append(title.id)
            if len(group.keys()) == len(self.genomes)-1:
                self.orthologs.append(group.values())
        dump = open('parsed','wb')
        pickle.dump(self.orthologs,dump)

    def parse_secondary(self):
        records = NCBIXML.parse(open('secondary.xml'))
        positives = []
        for record in records:
            title = ParseDefline(record.query)
            for desc in record.descriptions:

                if desc.e <= self.expect_diff:
                    positives.append(title.id)
        print positives
        hipster = SeqIO.parse(open(self.substract),'fasta')
        hipster = SeqIO.to_dict(hipster)
        unique = list(set(hipster.keys())-set(positives))
        handle = open('unique.faa','wb')
        for key in unique:
            record = hipster[key]
            SeqIO.write(record,handle,'fasta')

    def write_ortho(self):
        db = SeqIO.parse(open(self.db),'fasta')
        db = SeqIO.to_dict(db)
        ortho = open('orthologs.faa','wb')
        writeme = []
        for group in self.orthologs:
            for items in group:
                for item in items:
                    writeme.append(item)
        for item in list(set(writeme)):
            SeqIO.write(db[item],ortho,'fasta')

if __name__ == '__main__':
    g = [('E13','E13.txt'),('E14','E14.txt'),('E15','E15.txt'),('E16','E16.txt'),('E17','E17.txt'),('E18','E18.txt'),('E19','E19.txt'),('E20','E20.txt')]
    O = Ortho()
    O.genomes=g
    O.substract = 'E21.txt'
    O.build_db()
    O.primary_blast()
    O.parse_blast()
    O.write_ortho()
    O.secondary_blast()
    O.parse_secondary()
    # AF1TANG@GMAIL.COM
