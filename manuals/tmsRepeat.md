# Documentation for script: _tmsRepeat.pl_

## Summary
We developed tmsRepeat as a strategy to identify repeated regions of TMSs in a 
transporter sequence that allows incorporation of prior knowledge about 
the TMS topology and expected size of the repeat unit. The program cuts TMS-bundles 
of a specified length depending on the size of the repeat unit expected by the user, 
aligns the sequences of the different TMS- bundles with the Smith-Waterman algorithm, 
prints the results in a HTML report that includes the hydropathy plots of the proteins 
in a family, the limits of the bundles aligned and shows what regions within the bundles 
were contained in the alignments.


## How to cite this program
If you find this program useful, please cite the paper:  

  * Medrano-Soto A, Moreno-Hagelsieb G, McLaughlin D, Ye ZS, Hendargo KJ, Saier MH Jr. 
  _Bioinformatic characterization of the Anoctamin Superfamily of Ca2+-activated ion 
  channels and lipid scramblases._  2018. PLoS One. **13**(3):e0192851 
  **PMID:** [29579047](https://www.ncbi.nlm.nih.gov/pubmed/?term=29579047)  


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_blast+ 2.4.0 to 2.6.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

2. **_R_**  
The [R package](https://www.r-project.org/) is used. Make sure the following 
packages are installed: cluster, MCMCpack and ape.

3. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.

## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

     -i query filename in fasta format, required
     -o output folder, default: Clusters
     -p program for pairwise comparisons
        [blastp|fasta36|ssearch36|ublast], default: blastp
     -c agglomerative clustering method
        [average|complete|single|ward|weighted],
         default: ward

