# Documentation for script: _list_top_gsat_hits.pl_

## Summary
After running _protocol2_ with the option (-gr) that runs _gsat_ across the full homology
transitivity path, this program will generate a sorted list of the protocol2 hits with top
scores. This output list can then be used to focus the analysis on the protocol2
hits with the highest possibilities of capturing the "remote homology" signal between 
the two families under study.

## Contributor
Arturo Medrano-Soto


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl. The following module is/ required:  
  b) TCDB::CheckDependencies (included in the [TCDBtools distribution](https://github.com/SaierLaboratory/TCDBtools))  
  
2. **_protocol2_**  
This program is included in the [BioVx distribution](https://github.com/SaierLaboratory/BioVx).

3. **_blast+ 2.6.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

## Command line options
The following options are available. You can also run the 
script without arguments to display the options:


After running program areFamiliesHomologous, use this script to
list all protocol2 hits with high GSAT score.

Input parameters:

-p2d, --proto2-dir {path}
   The root directory for all protocol2 results. Default value
   is the current directory.

-inc, --increasing
   List top GSAT hits in increasing order. By default
   results are listed this way (This option is ignored if
   the option -b is given).

-dec, --decreasing
   List top GSAT hits in decreasing order

-b,  --best-score
   Print only the best GSAT score that will represent these
   families. It takes into account all scores A-B,B-C,C-D.

-full --full-transitivity
  List Top GSAT hits including the corrensponding comparisons
  with the original TCDB proteins to complete the transitivity
  path between two TCDB families. NOTE: Use this option only
  if famXpander generated the candidate homologous proteins.
  (optional)

-v, --validate
  If given, this option will extract the full length sequences
  of protein pairs involved in signficant protocol2 and gsat hits. 
  Then runs protocol2 again in these two sets of full-lengthproteins. 
  This is because some times protocol2 runs with the aligned blast 
  regions reported by famXpander and the TMS numbers do not reflect 
  the order of the TMS in the entire protein, just the region that 
  was aligned by psiblast. Therefore, it is necessary to realign
  the sequences so the right TMS numbers are displayed by protocol2.
  Therefore, use this option ONLY when you ran famXpander and kept
  aligned regions instead of full proteins.

-h, --help
   Display this help. Also displayed if script is run without arguments.

