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

1. **_ssearch36 version: 36.3.8e_**  
Other versions of ssearch may require minor adaptations. Visit the
[download site](https://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml). 

2. **_tmsplit_**  
This program is available in our [Python repository](https://github.com/SaierLaboratory/BioVx).

3. **_extractFamily.pl_**  
This program is available in our [Perl repository](https://github.com/SaierLaboratory/TCDBtools). 

4. **_quod.py_**  
This program is available in our [Python repository](https://github.com/SaierLaboratory/BioVx).

5. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.

## Command line options
The following options are available. You can also run the 
script without arguments to display the options:

    -i, --infile {path}
      Input file with id/accession(s) of the protein(s) to analyze and the coordinates
      of the TMSs in that protein(s).
      (Argument is mandatory).

    -s, --seqs {path}
      Directory to access the sequences in FASTA formate that will be used to 
      search for repeats.
     (Argument is mandatory)

    -f, --id-format {string}
      format of identifier used:
        tc    plain tcdb identifier of a system (e.g., 2.A.1.8.1)
        tca   tcdb id and accession separated by dash (e.g. 2.A.1.8.3-Q9R6U5)
        o     other, it can be refSeq, uniprot or custom, but it is requried
              that is is a single string without spaces.
      (Argument is mandatory)

    -r, --repeat-unit {int)
      Size in TMS of the repeat unit to search in the protein.
      (Argument is mandatory)

    -t, --tail-size {int}
      Number of residues to add to the beginning and end of TMS regions before
      running comparisons. Value should be less than or equal to 15 residues.
     (Default: 5);

    -e, --evalue {float}
      Maximum evalue to consider an alignment between two TMS bundles significant.
      (Default: 0.1);

    -c, --coverage {float}
      Minimum alignment coverage of the smallest bundle to consider an alignment
      signifiant.
      (default: 0.8)

    -id, --perc-identity {float}
       Minimum identity in the alignemnt to consider alignments signficant.
       (defatul 0.3);

    -gs, --gsat-shuffles {int}
       Number of shuffles that will be used to run GSAT on good matches.

