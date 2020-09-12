# Documentation for script: _extractFamily.pl_

## Summary
This program allows easy retrieval of protein sequences from TCDB. Sequences can be downloaded based on TCIDs at the subclass, family, subfamily and system levels, as well as the whole protein content in TCDB. The program can present the sequences in FASTA and 2-column formats (i.e., accession and sequence). In addition, it is possible to generate a [_BLAST database_](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) with the specified sequences. Furthermore, sequences can be extracted from a local version of TCDB, which is useful when working with projects that freeze databases in a specific version.

## Contributor
Gabriel Morneo-Hagelsieb


## How to cite this program
If you find this program useful, please cite the paper:  

  * Medrano-Soto A, Moreno-Hagelsieb G, McLaughlin D, Ye ZS, Hendargo KJ, Saier MH Jr. _Bioinformatic characterization of the Anoctamin Superfamily of Ca2+-activated ion channels and lipid scramblases._  2018. PLoS One. **13**(3):e0192851  **PMID:** [29579047](https://www.ncbi.nlm.nih.gov/pubmed/?term=29579047)  


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_blast+ 2.6.0 to 2.10.0_**  
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
    
    -F  output format: fasta|column|blast (default: fasta)
        column:  presents data in 2 tab-delimited columns,
                 TCDB ID and the sequence in one string.
        blast:   Generates a BLAST database with the downloaded
                 sequences.
                 
    -d  path to a fasta file with all sequences in TCDB,
        as in previous runs of 'extractFamily.pl -i tcdb -f fasta'
        This allows to work with frozen version of TCDB.
        Default: tcdb  (online database).

