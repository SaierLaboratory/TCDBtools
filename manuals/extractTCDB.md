# Documentation for script: _extractTCDB.pl_

## Summary
This program allows easy retrieval of protein sequences from [TCDB](http://tcdb.org). Sequences can be downloaded based on TCIDs at the class, subclass, family, subfamily and system levels, as well as the whole protein content in TCDB. The program can present the sequences several formats (e.g., fasta and 2-column). In addition, it is possible to generate a [_BLAST database_](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) with the specified sequences. Furthermore, sequences can be extracted from a local version of TCDB, which is useful when working with projects that freeze databases in a specific version.

## Contributor
Gabriel Morneo-Hagelsieb


## How to cite this program
If you find this program useful, please cite the paper:  

  * Medrano-Soto A, Moreno-Hagelsieb G, McLaughlin D, Ye ZS, Hendargo KJ, Saier MH Jr. _Bioinformatic characterization of the Anoctamin Superfamily of Ca2+-activated ion channels and lipid scramblases._  2018. PLoS One. **13**(3):e0192851  **PMID:** [29579047](https://www.ncbi.nlm.nih.gov/pubmed/?term=29579047)  


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_blast+ 2.10.0 or higher_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

2. **_wget_**  
Sequences are downloaded from the TCDB server using [wget](https://www.gnu.org/software/wget).

3. **_PERL 5.18 or higher_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -i  TCDB ID (required)
        A TCDB accesion representing a clase, sublclass,
        family, subfamily or system. (e.g., 1, 1.A, 1.A.1,
        (1.A.1.1 or 1.A.1.1.1).
        "-i tcdb", "-i all", or "-i full" will bring the
        complete TCDB database.
        
    -o  output directory where data will be saved. (default: Families) 
    
    -f  Output format: fasta|column|blast (default: fasta)
        fasta:   Most common format to represent sequences.
        column:  presents data in 2 tab-delimited columns,
                 TCDB ID and the sequence in one string.
        blast:   Generates a BLAST database with the downloaded
                 sequences.
        diamond: Generates a database compatible with DIAMOND.
                 
    -s  Path to a fasta file with all sequences in TCDB,
        as in previous runs of 'extractFamily.pl -i tcdb -f fasta'
        This allows to work with frozen version of TCDB.
        Default: tcdb  (online database).


