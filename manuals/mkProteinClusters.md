# Documentation for script: _mkProteinClusters.pl_

## Summary
Given a number of proteins presumed to form a family or superfamily, this program performs hierarchical clustering based on their pairwise alignment bit scores and reports the results as a tree. The scores can be estimated using [_BLAST_](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [_Fasta or SSEARCH_](https://fasta.bioch.virginia.edu/fasta_www2/fasta_class.shtml). Several methods for merging clusters are supported (i.e., Ward, single, average, complete and weighted). The clustering is performed using the [_R statistical computing environment_](https://www.r-project.org/). This program is especially useful when working with large numbers of proteins, or with highly diverse proteins (as in a superfamily) where standard phylogenies cannot be constructed due to low-quality multiple alignments.


## Contributor
Gabriel Moreno-Hagelsieb  


## How to cite this program
If you find this program useful, please cite the paper:  

  * Medrano-Soto A, Moreno-Hagelsieb G, McLaughlin D, Ye ZS, Hendargo KJ, Saier MH Jr. 
  _Bioinformatic characterization of the Anoctamin Superfamily of Ca2+-activated ion 
  channels and lipid scramblases._  2018. PLoS One. **13**(3):e0192851 
  **PMID:** [29579047](https://www.ncbi.nlm.nih.gov/pubmed/?term=29579047)  


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_blast+ 2.6.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

2. **_ssearch36 and fasta36 version: 36.3.8_**  
Other versions of ssearch may require minor adaptations. Visit the
[download site](https://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml). 

3. **_R_**  
The [R package](https://www.r-project.org/) is used. Make sure the following 
packages are installed: cluster, MCMCpack and ape.

4. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.


## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -i query filename in fasta format, required
    
    -o output folder. Default: Clusters
    
    -p program for pairwise comparisons
       [blastp|fasta36|ssearch36]. Default: ssearch36
       
    -s Amino acid substitution matrix that wil be used
       by ssearch36. Default: ssearch36 default for option -s
      
    -z Algorithm to be used by ssearch36 to calculate E-values,
       default: ssearch36 default for option -z 
       
    -k Number of shuffles to be used by ssearch36 in the
       calculation of E-values. Default: ssearch36 default for option -k
       
    -c agglomerative clustering method
       [average|complete|single|ward|weighted].
       Default: ward
       
