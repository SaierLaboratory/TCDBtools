# Documentation for script: _famXpander.pl_

## Summary
This program performs _blastp_ or _psiblast_ searches of query proteins 
against the NCBI non-redundant database. The sequences of top matches
are extracted and redundant sequences filtered out. This program is 
similar to protocol1 from the [BioV suite](https://github.com/SaierLaboratory/BioVx),  
but there are three main differences: **1)** famXpander performs blast searches 
locally, although it can also run them remotely; **2)** blast searches are 
performed  first for all protein queries and the redundant sequences are removed 
at the end only once; and **3)** famXpander can extract either the full sequence 
of the top hits or just the aligned regions.

## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_blast+ 2.4.0_**  
Newer versions of blast may require minor adaptations. Visit 
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

2. **NCBI non-redundant protein database**  
Given that blast runs locally the NCBI non-redundant (NR) database
must be available locally through the environment varaible _$BLASTDB_. 
You can download NR from the NCBI FTP site:  
ftp://ftp.ncbi.nlm.nih.gov/blast/db/

3. **_cd-hit 4.6_**  
Visit the [official website](http://weizhongli-lab.org/cd-hit/) to 
download the latest version.

4. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.

## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -i  input filename in fasta format (required)  
    -o  output folder (default: faaOut)  
    -n  max number of aligned sequences to keep (default: 10000)  
    -e  evalue threshold (default: 1e-7)  
    -f  psiblast evalue threshold (default: 1e-5)  
    -t  psiblast iterations (default: 1)  
    -h  keep only aligned region [T/F] (default: T)  
    -c  minimum alignment coverage of original sequence 
        (default: 0.8)  
    -s  minimal subject seq length relative to query seq length 
        (default: 0.8)
    -l  maximal subject seq length relative to query seq length.
        Option is ignored if **-h T** (default: 1.25)
    -r  identity redundancy threshold for cd-hit (default: 0.8)  
    -a  number of cpus to use.  
    -p  run remotely at ncbi [T/F] (default F)  

