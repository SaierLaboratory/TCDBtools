# Documentation for script: _findNovelTransporters.pl_

## Summary
In transportome characterization, the main task is the identification of transporters for which there is knowledge available in [TCDB](http://tcdb.org). This program performs the other side of the analysis, that is, the identification of integral membrane proteins in the genome that are likely transporters (e.g., by having a minimal number of TMSs), but show no significant similarity to any protein in TCDB. The user can specify the minimal number of TMS in the expected transporters, the alignment coverage and the E-value cutoff. For each protein that passes these criteria the program reports hydropathy plots that can be used by curators to ultimately decide whether proteins are good candidates for novel families in TCDB.


## Contributor
Arturo Medrano-Soto


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.

2. **_blast+ 2.8.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

3. **_extractFamily.pl_**  
This program is included in the [TCDBtools distribution](https://github.com/SaierLaboratory/TCDBtools). 

4. **_hmmtop_**  
TMS are predicted with HMMTOP. Visit the [download site](http://www.enzim.hu/hmmtop/html/download.html).

5. **quod.py**  
This program in included in the [BioVx distribution](https://github.com/SaierLaboratory/BioVx). 


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:


    -wd, --workdir {path}  (Defalut: current directory)
       Path to the directory where the ouput and temporary files will be stored.

    -bdb, --blastdb {path}  (Default: ./blastdb/tcdb)
       Full path to the TCDB blast database that will be used.

    -p, --proteome {path}  (Mandatory)
       File in fasta format with the proteome that will be analyzed.
       Fasta headers must contain the ID (accession) of the
       protein followed by a space and the functional description 
       (e.g. >AKM78775.1 AP4A hydrolase ...).

    -tms, --min-tms {integer}  (Default: 4)
       Minimum number of TMS that candidate transporters should have.

    -e, --min-evalue {float}  (default: 0.01)
       Proteins with E-value greater than this value will be considered.

    -re,  --redundancy-evalue {float}  (default: 1e-5)
       Evalue threshold to consider two proteins redundant for the
       the purpose of presenting final results.

    -h, --help
       Display this help. This argument takes precedence over any other
       argument.
