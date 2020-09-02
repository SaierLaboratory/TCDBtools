# Documentation for script: _organizeSeqsByTMS.pl_

## Summary
This script takes a file with multiple protein sequences in FASTA fromat
and generates outfile files where each file has proteins with a given
number of HMMTOP-predicted TMSs.


## Contributor
Arturo Medrano-Soto


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl. The following modules
are required:  
  a) [Bioperl](https://bioperl.org/)  
  b) TCDB::CheckDependencies (included in the [TCDBtools distribution](https://github.com/SaierLaboratory/TCDBtools)  

2. **_hmmtop_**  
TMS are predicted with HMMTOP. Visit the [download site](http://www.enzim.hu/hmmtop/html/download.html).


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:


    -s, --seqs-file={string}
      Fasta file with the proteins sequences that will be organized.
      Argument is mandatory.

    -o, --outdir={string}
      Directory where the output files will be located. If the directory
      does not exist, it will be created.
      
    -h, --help
      Print this help. This argument takes precedence over any other
       argument.
