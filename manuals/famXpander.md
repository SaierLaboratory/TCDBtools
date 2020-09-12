# Documentation for script: _famXpander.pl_

## Summary
This program substitutes and extends the capabilities of [_Protocol1_](https://github.com/SaierLaboratory/BioVx/blob/master/manuals/BioV_manual.pdf) to retrieve homologs from NCBI. In addition to running [_PSI-BLAST_](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) against the NCBI nonredundant database, retrieving large numbers of proteins and filtering results based on E-values and degrees of sequence redundancy, _famXpander_ provides additional controls (i.e., minimal coverage of query and subject, minimal sequence lengths, retrieval of either the full sequences of subject proteins or just the aligned regions, and _BLAST_ searches can be performed either locally or remotely using NCBI servers). This provides the raw data necessary to infer relationships between families, identify repeat units within a family and investigate conserved motifs and domains.


## Contributor
Gabriel Moreno-Hagelsieb


## How to cite this program
If you find this program useful, please cite the paper:  

  * Medrano-Soto A, Moreno-Hagelsieb G, McLaughlin D, Ye ZS, Hendargo KJ, Saier MH Jr. _Bioinformatic characterization of the Anoctamin Superfamily of Ca2+-activated ion channels and lipid scramblases._  2018. PLoS One. **13**(3):e0192851  **PMID:** [29579047](https://www.ncbi.nlm.nih.gov/pubmed/?term=29579047)  

## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_blast+ 2.6.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

2. **_NCBI non-redundant protein database_**  
Given that blast runs locally, the NCBI non-redundant (NR) database
must be available locally through the environment varaible _$BLASTDB_. 
You can download NR from the NCBI FTP site: ftp://ftp.ncbi.nlm.nih.gov/blast/db/  

3. **_cd-hit 4.8_**  
Visit the [official website](http://weizhongli-lab.org/cd-hit/) to 
download the latest version.

4. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.

## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -i  input filename in fasta format (Required)
    
    -d  non-redundant database (Default: nr)
    
    -o  output folder (default: faaOut)  
    
    -n  max number of aligned sequences to keep (default: 10000)  
    
    -e  evalue threshold (default: 1e-7)  
    
    -f  psiblast evalue threshold (default: 1e-5)  
    
    -t  psiblast iterations (default: 1)  
    
    -h  keep only aligned region [T/F] (default: T)  
    
    -c  minimum alignment coverage of original sequence 
        (default: 0.8) 
        
    -x  coverage applies to either sequence (Default: F)
    
    -s  minimal subject seq length relative to query seq length 
        (default: 0.8)
        
    -l  maximal subject seq length relative to query seq length.
        Option is ignored if -h T  (default: 1.25)
        
    -r  identity redundancy threshold for cd-hit (default: 0.8)  
    
    -a  number of cpus to use. (Default: all available).
    
    -w  overwrite previous psiblast.tbl file, if it exists [T/F].
        (Default: F)
        
    -p  run remotely at ncbi [T/F] (Default F)  

