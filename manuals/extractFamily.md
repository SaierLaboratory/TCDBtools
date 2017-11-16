# Documentation for script: _extractFamily.pl_

## Summary
This program extracts the protein sequences of any transport family, etc. in TCDB; 
including the complete TCDB database. The sequences can downloaded in several formats.


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_blast+ 2.4.0 to 2.6.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

2. **_wget_**  
Sequences are downloaded from the TCDB server using [wget](https://www.gnu.org/software/wget).

3. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.

## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -i  TCDB family ID (required)
        "-i tcdb", "-i all", or "-i full" will bring the
        complete TCDB database
    -o  output directory where data will be saved. (default: Families)  
    -f  output format: fasta|column|blast (default: fasta)
        column:  presents data in 2 tab-delimited columns,
                 TCDB ID and the sequence in one string.
        blast:   Generates a BLAST database with the downloaded
                 sequences.
    -d  path to a fasta file with all sequences in TCDB,
        as in previous runs of 'extractFamily.pl -i tcdb -f fasta'
        Default: tcdb  (online database).

