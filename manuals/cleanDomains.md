# Documentation for script: _cleanDomains.pl_


## Summary
This program extracts results from[_hmmscan_](http://hmmer.org/), 
[_rpsblast_](https://www.perl.org/), or [_mmseqs_](https://github.com/soedinglab/MMseqs2)
normally obtained using [matchDomains.pl](matchDomains.md)


## Contributor
Gabriel Moreno-Hagelsieb


## Dependencies
The following dependencies need to be available in your enviroment for this 
program to run properly:

1. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:  

    SYNOPSIS
      cleanDomains.pl -q [matchDomainFile] [options]

    OPTIONS
    
     -q  file with domain matches from matchDomains.pl [e.g.
         GCF_000005845.cdd.mmseqs.bz2], required

     -f  domain family, normally [cog|cd|cdd|pfam|tigrfam], make it explicit if
         not part of the name of the file with match results [e.g.
         GCF_000005845.cdd.mmseqs.bz2]

     -o  output folder, default: cleanDomains

     -c  minimum coverage of domain model, default 0.60

     -v  maximum overlap between domains, default 0.15

     -r  reference file with all scanned sequences [e.g. GCF_000005845.faa.gz]

     -a  append annotations (T|F). If 'T' the program will produce a file with
         annotations appended (only works with cdd|cog|cd), default F
