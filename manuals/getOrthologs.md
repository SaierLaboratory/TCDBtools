# Documentation for script: _getOrthologs.pl_


## Summary
Compares the proteomes of selected genomes and infers pairs of orthoglogs 
based on reciprobal best hits.

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

3. **_Diamond_**  
Visit the [official website](http://www.diamondsearch.org/) for details on this program.

4. **_MMseqs_**  
Open-source software suite for very fast, parallelized protein sequence searches
and clustering of huge protein sequence data sets. For more information, visit the
[official repository](https://github.com/soedinglab/MMseqs2).

5. **_Lastal_**  
THis program finds similar regions between sequences and aligns them. Visit the
[Official website](http://last.cbrc.jp/) for more details.


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    SYNOPSIS
      getOrthologs.pl -q [filename1.faa] -t [filename2.faa] [options]

    EXAMPLES
      getOrthologs.pl -q Genomes/GCF_000005845.faa.gz -t
      Genomes/GCF_000009045.faa.gz

     getOrthologs.pl -d Genomes -o AllvsAllRBH


    OPTIONS
    -q  query fasta file(s) or directory with fasta files, required. files can
        be compressed with gzip or bzip2.

    -t  target fasta file(s) or directory with fasta files, required. files can
        be compressed with gzip or bzip2.

    -d  directory with fasta files for all-vs-all comparisons, no default. If
        set -q and -t will be ignored.

    -o  directory for reciprocal best hits, default: RBH.

    -p  program for pairwise comparisons [blastp|diamond|lastal|mmseqs],
        default: diamond

    -s  sensitivity (only works for diamond and mmseqs):
          F:  low sensitivity: diamond: fast, mmseqs: -s 1
          S:  diamond: sensitive, mmseqs: -s 2
          M:  diamond: more-sensitive, mmseqs: -s 4
          V:  diamond: very-sensitive, mmseqs: -s 5.7
          U:  highest sensitivity: diamond: ultra-sensitive, mmseqs: -s 7
         
          for mmseqs the option will also accept numbers between 1 and 7,
          defaults: diamond: V; mmseqs: 5.7

    -m  directory for blastp|diamond|lastal|mmseqs results, default compRuns.

    -c  minimum coverage of shortest sequence [30 - 99], default 60.

    -f  maximum overlap between fused sequences [0 - 0.5], default 0.1 (10%) of
        shortest sequence.

    -e  maximum e-value, default 1e-06

    -a  include aligned sequences in comparison results [T/F], default F. (not
         available in lastal)

    -r  find RBH [T|F], default 'T', if 'F' the program will only run the
        blastp|diamond|lastal|mmseqs comparisons

    -k  keep prior results [T|R|F], default 'T':
         F:  repeat both alignments and gathering of RBH
         R:  only repeat gathering of RBHs
         T:  keep prior alignments and prior RBHs

    -z  maximum number of targets to find with blastp|diamond|lastal|mmseqs
        minimum of 50. Defaults to 1/7th of the sequences in the target file

    -x  number of CPUs to use, default: 4 (max: 4)
