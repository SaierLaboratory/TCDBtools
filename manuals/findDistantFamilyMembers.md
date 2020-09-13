# Documentation for script: _findDistantFamilyMembers.pl_

## Summary
This program is used when we want to identify remote members of a family in TCDB in order to increase the coverage of sequence diversity within a family. It consists of BLASTing all members within a family against the NCBI NR database and selecting hits that show borderline similarity to other members in the family, while having no significant similarity to any other families in TCDB. Alignments should satisfy a minimal coverage and involve a minimal number of TMSs compatible with the topology of the family under analysis. This produces a list of candidates that should be further evaluated by an expert. 


## Contributor
Arturo Medrano-Soto


## Dependencies
The following dependencies need to be available in your enviroment for this 
program to run properly:

1. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl. The following modules
are required:  
  a) [Bioperl](https://bioperl.org/)  
  b) _TCDB::CheckDependencies_ (included in the [TCDBtools distribution](https://github.com/SaierLaboratory/TCDBtools))

2. **_extractFamily.pl_**  
This program is available in our [Perl repository](https://github.com/SaierLaboratory/TCDBtools). 

3. **_hmmtop_**  
TMS are predicted with HMMTOP. Visit the [download site](http://www.enzim.hu/hmmtop/html/download.html).

4. **_blast+ 2.8.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the [download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). If run locally the NCBI NR needs to be available in your $BLASTDB environment varaible.  


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -d, --indir {path}
      Path to the directory where the ouput files of famXpander are located
      (i.e. psiblast.tbl and results.faa). Argument is Mandatory.

    -o, --outfile {path}
      Ouput file where results should be saved. (Mandatory)

    -od, --outdir {path}
      Output directory for temporary and final output files.
      (default: ./)

    -f, --family {string}
      Family for which remote homologs will be found. This allows to
      identify TCDB blast matches with the right family.
      (Mandatory)

    -tcbdb, --tc-blast-db-dir {path}
      Path to the blast DB with all proteins in tcdb. This database needs 
      to be created with program extractFamily and be called 'tcdb'.
      If the blast DB is not found in this location, it will be created
      for you.
      (Default ~/db/blastdb);

    -n, --min-tms-num {integer}
      Minimum number of TMS that should be found in distant homologs.
      (Default 1).

    -m, fam-tcblast-hits {integer}
      When blasting against TCDB the reference family, passed with option
      -f, should be wihtin the specified top number hits. (Default 3).

    -gt, --greater-than {float}
      Lower E-value threshold to parse famXpander output.  Hits with E-value 
      greater than this value will be extracted. (default 1e-7).

    -lt, --lower-than {path}
      Higher threshold to parse famXpander output. Hits lower than this 
      E-value will be extracted. (default 1.0)

    -ilt --id-lower-than {float}
      Identity should be lower than this value. This parameter can
      take any positive value lower than 100. By default identity
      is not taken into consideration.

    -igt --id-greater-than {float}
      Identity should be greter than this value. This parameter can
      take any positive value lower than 100. By default identity
      is not taken into consideration.

    -llt --len-lower-than {int}
      Length of subject should be lower than this value. By default
      Length is not taken into consideration.

    -lgt --len-greater-than {int}
      Length of subject should be greater than this value. By default
      Length is not taken into consideration.

    -r  --len-ratio {float}
      Length ratio between redundant distant homologs. This helps
      to determine when to keep two redundant distant homologs if
      they have very different lengths (default 1.8)

    -re  --redundant-evalue
      Evalue threshold to consider to remote homologs redundat.
      (default 1e-5).

    -s, --sort
      Sort blast results. Results can be sorted in ascending (asc),
      descending (desc) or no order (no). 
      (default is asc).

    NOTE: It's possible to specify a range of E-values by providing
          both -gt and -lt thresholds.

