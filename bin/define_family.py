#!/usr/bin/env python
from Bio import Entrez
import tcdb,re,sys,blast,shutil

Entrez.email = 'vreddy@ucalgary.ca'

def verify(acc):
    error = re.compile(r'error',re.IGNORECASE)
    handle = Entrez.efetch("protein", id=acc, rettype='fasta', retmode='xml')
    if bool(error.search(handle.read())):
        return False
    return True

if __name__=='__main__':

    try:
        family = sys.argv[1]
        method = sys.argv[2].lower()
        output = sys.argv[3]
    except:
        print "Usage: define_family.py FAMILY <P/PSI> OUTPUT"
        quit()

    accs = tcdb.define_family(family)
    print "Validating all members of family"
    valid = [i for i in accs if verify(i) is True]
    invalid = list(set(accs)-set(valid))
    if len(invalid) > 0:
        print "The Following IDs are invalid:"
        for bad in invalid:
            print bad
    ncbi = blast.tools()
    for good in valid:
        print "Blasting :: %s"%good
        if method == 'psi':
            ncbi.psiblast(good)
        else:
            ncbi.blastp(good)
    print "Found %i results" %len(ncbi.gis)
    c = float(raw_input("CD-HIT Threshold % :: "))
    ncbi.build_xml()
    ncbi.maketable(c)
    shutil.copy(ncbi.fasta_file.name,output)
    print "Done, file created."
