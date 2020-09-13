# Documentation for script: _matchDomains.pl_


## Summary
This program identifies domains in protein sequences based on several programs (i.e., [_hmmscan_](http://hmmer.org/), [_rps-blast_](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) or [_mmseqs_](https://github.com/soedinglab/MMseqs2)) and domain databases (i.e., [Pfam](https://pfam.xfam.org/), [CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml), [COG](https://www.ncbi.nlm.nih.gov/COG/), [TIGRFAMs](http://tigrfams.jcvi.org/cgi-bin/index.cgi)). This program facilitates comparisons of domain contents between genomes and the identification of characteristic domains within families and superfamilies.


## Contributor
Gabriel Moreno-Hagelsieb


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.  

2. **_Blast+ 2.6.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

3. **_HHMER 3.2.1_**  
Package suite used for searching sequence databases for sequence homologs, and for making 
sequence alignments based on probabilistic models called profile hidden Markov models 
(profile HMMs). You can download this suite of programs from its [official site](http://hmmer.org/). 

4. **_MMseqs_**  
Open-source software suite for very fast, parallelized protein sequence searches
and clustering of huge protein sequence data sets. For more information, visit the
[official repository](https://github.com/soedinglab/MMseqs2).

5. **_Domain profile databases_**  
Download the following profile databases _Pfam_, _CDD_, _COG_, and _tigrfam_ as they 
need to be locally available. The program expects an environmental variable, DOMAINDB, pointing 
to a directory where the databases will be found. For example, if the base directory is named: 
DomainDBs, then the variable should poin to this directory, and it should contain the databases 
in subdirectories as follows: pfam and tigerfam in: DomainDBs/xfamDB; mmseqs databases 
in: DomainDBs/mmseqsDB; and cdd databases in: DomainDBs/cddDB


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    SYNOPSIS
       matchDomains.pl -q [fastaFile] -f [cog|cdd|pfam|tigrfam] [options]

    EXAMPLE
       matchDomains.pl -q GCF_000005845.faa.gz -f pfam -p mmseqs -o PFAM
       matchDomains.pl -q fastaFiles/*.faa.gz -f cdd -o CDD

    OPTIONS
     -q  query fasta file(s), required

     -f  domain database: file or [cog|cdd|pfam|tigrfam], required

     -p  program [hmmscan|rpsblast|mmseqs], except for mmseqs, can be guessed
         from database:
         -   rpsblast for cog and cdd
         -   hmmscan for Pfam and TIGRFAM
         -   mmseqs has to be specified in command line

     -o  output folder, default: scanDomains

     -e  e-value threshold, default 0.001 (NCBI uses 0.01), scientific notation
         acceptable (e.g. 1e-3)

     -x  number of CPUs to use, default: 1 (max: 4)

     -c  running in computer cluster [T|F], default 'F'
     
     -h  display this help message.
     
