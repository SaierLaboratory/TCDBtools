# Documentation for script: _interpret_proto2_top_hits.pl_

## Summary
After running [_areFamiliesHomologous_](areFamiliesHomologous.md) with the option (-gr) and 
[_list_top_gsat_hits_](list_top_gsat_hits.md), use this script to run 
[_hvordan_](https://gitlab.com/khendarg/hvordan/blob/master/docs/hvordan.md) 
and generate hydropathy plots (including [_Pfam_](https://pfam.xfam.org/) matches) across 
the homology transitivity path. This will substantially help with the interpretation of 
[_protocol2_](https://github.com/SaierLaboratory/BioVx/blob/master/manuals/BioV_manual.pdf) 
results for each pair of families analyzed. 


## Contributor
Arturo Medrano-Soto


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl. The following module is required:  
  a) TCDB::CheckDependencies (included in the [TCDBtools distribution](https://github.com/SaierLaboratory/TCDBtools))  
  
2. **_hvordan_**  
This program is included in the [BioVx distribution](https://github.com/SaierLaboratory/BioVx).

3. **_blast+ 2.8.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -fxd, --fxpand-dir  { path }  (Mandatory)
      Directory with the results of runing famXpander on the families
      under analysis (mandatory).

    -p2r, --proto2-root-dir { path }  (Optional but see -p2d)
      This is the main directory where the the results of running protocol2
      for multiple pairs of families are found. This option can't be used
      in combination with option -p2d.

    -p2d, --proto2-dir {path}  (Optional but see -p2r)
      The directory with the results of running  protocol2 for a couple
      of families. This is the directory where file report.html is located.
      This option cannot  be used in combination with option -p2r

    -tsd, --top-scores-dir  { path }  (Optional)
      Directory where the results of running list_top_gsat_hits are stored.
      By default it takes the value passed to option -p2r.

    -n, --num-top-hits { int }  (Optional; Default: 15)
      Number of significant top hits for which hvordan.py will be run.

    -p2lt, --proto2-low-cutoff {float}
       Minimum threshold to identify significant hits from Protocol_2 and GSAT
       (default is 14.0)

    -p2ht, --proto2-high-cutoff {float}
      Maximum threshold to identify significant hits from Protocol_2 and GSAT.
      Use this option along with -p2lt to focus the analysis on a range of scores.
      (default is unlimited)

    -o, --outdir { path }  (Mandatory)
      Directory with the hydropathy and blast results will be saved.

    -h, --help
      Display this help. Also displayed if script is run without arguments.
