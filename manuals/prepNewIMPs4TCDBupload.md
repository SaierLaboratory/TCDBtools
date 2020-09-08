# Documentation for script: _prepNewIMPs4TCDBupload.pl_

## Summary
After running [_findNovelTrasporters_](findNovelTransporters.md)  and obtaining a list of proteins
in a reference genome with little or no similarity to any protein in TCDB, this programs extracts 
homologs from NCBI for each of the new candidate transporters that will be used to expand current families
or create new families in TCDB. This list will significantly optimize  TCDB curators' time.


## Contributor
Arturo Medrano-Soto


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl. The following module is/ required:  
  b) TCDB::CheckDependencies (included in the [TCDBtools distribution](https://github.com/SaierLaboratory/TCDBtools))  
  
2. **_hmmtop_**  
TMS are predicted with HMMTOP. Visit the [download site](http://www.enzim.hu/hmmtop/html/download.html).

3. **_blast+ 2.6.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

4. **_cd-hit 4.8_**  
Visit the [official website](http://weizhongli-lab.org/cd-hit/) to 
download the latest version.

4. **_famXpander.pl_**  
This program is available in our [Perl repository](https://github.com/SaierLaboratory/TCDBtools). 


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -i, --input {list|file}  (Mandatory)
      An Accession, comma-separated list of accessions, or a file with the
      accessions (first columns) on which the perform the analysis.

    -wd, --workdir {path}  (Defalut: ./newIMPs4tcdb)
      Path to the directory where the ouput and temporary files will be stored.

    -t, --tcblastdb {path}  (Default: ~/db/blastdb/tcdb)
      Full path to the TCDB blast database that will be used.

    -gdb, --gnm-blastdb {string} (Default: nr)
      Full path to the blast DB of the whole proteome of the reference
      genome. This is necessary when proteins are not annotated with proper
      RefSeq accessions (e.g. locus_tags).

    -e, --evalue {float}  (default: 1e-3)
      NCBI proteins with E-value greater than this cutoff will be ignored.

    -re,  --redundancy-evalue {float}  (Default: 1e-15)
      Evalue threshold to consider two proteins redundant for the
      the purpose of presenting final results.

    -c, --coverage {float} (Default: 0.5)
      Minimum alignment coverage in pairwise alignments.

    -m, --mode {L|R} (Default: R)
      Indicate whether blast will be run locally (L) or remotely
      at NCBI (R).

    -nm, --num-members {int} (Defatult: 3)
     Number of non-redundant homologs to create a family in TCDB,
     or to add them to an existing family as remote homologs.

    -nh, --ncbi-hom {int} (Default: 25)
      Minimum number of homologs in NCBI to qualify for this analysis.

    -h, --help
      Display this help. This argument takes precedence over any other
      argument.
