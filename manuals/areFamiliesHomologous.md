# Documentation for script: _areFamiliesHomologous_

## Summary
Our software famXpander, Protocol2 and GSAT were integrated into a pipeline 
that significantly speeds up the analysis of distant evolutionary relationships between 
families using the transitivity property of homology, thus eliminating the possibility 
of human errors (PMID: 29579047). In addition, BLASTs for all proteins in TCDB were 
pre-computed to allow rapid comparisons of large superfamilies and make analyses at the 
level of class or subclass. Users have the option to use precomputed alignments or run 
new BLASTs for their sequences.

## Description
First, famXpander is run in order to extract an expanded list of candidate homologous 
proteins for each family. Second, protocol2 is run to determine homology and TMS 
topology. Only those pairs that have a protocol2 score above a user-specified 
threshold and minimum alignment length will be singled out to run gsat. If gsat has 
a z-score above a user-specified value it will be indicated.

## How to cite this program
If you find this program useful, please cite the paper:  

  * Medrano-Soto A, Moreno-Hagelsieb G, McLaughlin D, Ye ZS, Hendargo KJ, Saier MH Jr. _Bioinformatic characterization of the Anoctamin Superfamily of Ca2+-activated ion channels and lipid scramblases._  
2018. PLoS One. **13**(3):e0192851   
PMID: [29579047](https://www.ncbi.nlm.nih.gov/pubmed/?term=29579047)  


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_extractFamily.pl_**  
This program is available in our [Perl repository](https://github.com/SaierLaboratory/TCDBtools). 

2. **_famXpander.pl_**  
This program is available in our [Perl repository](https://github.com/SaierLaboratory/TCDBtools). 

3. **_protocol2.py_**  
This program is available in our [Python repository](https://github.com/SaierLaboratory/BioVx).

4. **_gsat.py_**  
This program is available in our [Python repository](https://github.com/SaierLaboratory/BioVx).

5. **_PERL 5.18_**  
Visit the [official website](https://www.perl.org/). This program 
was not tested with more recent versions of perl.

## Command line options
The following options are available. You can also run the 
script without arguments to display the options:


    -f1, --family1 {string}
       TCDB ID string of the first family to be compared. It can be
       any string, but if the option -x is given, a TCDB ID is assumed
       (option is mandatory).

    -f2, --family2 {string}
       TCDB of the second family to be compared. It can be
       any string, but if the option -x is given, a TCDB ID is assumed
       (option is mandatory).

    -u, --prog {string}
       Specify which program will be used to find homologs:
       proto1:    protocol1 (currently broken),
       fxpand:    famXpander
       globalfxd: Use global precomputed results of famXpander for all
                  proteins in TCDB.
       If fxpand is selected, input sequences must be in fasta format. 
       If proto1 is selected, input sequences must be in 2-column format.
       (Option is mandatory)

    -d, --indir {string}
       Path to the directory with the amino acid sequences of the family
       members. If argument is not given, a directory will be created in the
       current directory based on the info provided in the -f1 and
       -f2 options (optional)

    -p1d, --proto1-dir {path}
       Directory to place the results of protocol_1. This option
       is specially useful when  more than one family will be compared.

    -gfxd, --global-fxpand-dir {path} (optional)
       Directory where precomputed results of famXpander for all families in
       TCDB are located. By default this directory is: /ResearchData/famXpander

    -p2d, --proto2-dir {path}
       Directory where results of protocol_2 and GSAT will be stored.
       If argument-p2d is not given, a directory will be created in
       the current directory based on the info provided in the
       -f1 and -f2 options (optional).

    -x, --extract (can be negated with --no-extract)
       Flag Indicating that the amino acid sequences from both families should be
       automatically extracted from TCDB. It is assumed that the
       identifiers provided by the -f1 and -f2 options are TCDB family IDs.
       If this option is not given, it is assumed that the user already
       extracted the sequences and they are available in the path provided by
       the -d option (optional)

    -ox, --only-extract-seqs
       Flag instructing only to download the sequences of the input families
       from TCDB. This is useful when prior to running protocol1 and protocol2
       the user wants to limit or edit the protein sequnces. The program
       aborts after retrieving the sequences. Option -u is mandatory in
       in order to determine in what format the sequences will be
       stored. This option is incompatible with -x.

    -n, --psiblast-it {integer}
       Number of iterations that psi-blast will perform.
       (Defauult is 1)

    -e, --evalue
       Expectation value threshold for saving psiblast hits
       Argument is optional (defult value 1e-7)

    -e2, --inc-evalue
       E-value inclusion threshold for pairwise alignments in psi-blast.
       (Default value is 1e-5)

    -a, --alignment-matches {int}
       Number of blast matches to retrieve from famXpander or protocol1.
      (Default 10000)

    -k, --keep-aln-regions {T/F}
      Option that indicates whether famXpander will feed only the psiblast
      aligned regions to protocol2. 
      (Default T)

    -c, --min-coverage {float}
      For famXpander, minimum alignment coverage of original sequence
      (Default 0.8)

    -s, --min-rseq-length {float}
      For famXpander, minimal sequence length relative to original seq length,
     (Default 0.8)

    -l, --max-seq-length {float}
      famXpander  maximal sequence length relative to original seq length.
      The higher this value, the easier it will be to detect
      protein fusions.
      (Default 5.0 (ignored if -k is given)

    -fxr, --fxpand-remote
      Run famXpander remotely at NCBI. By default famXpander is run locally.

    -r, --cdhit-cutoff {float}
      Identity redundancy threshold for cd-hit.
      (default 0.9)

    -p2lt, --proto2-low-cutoff {float}
      Minimum threshold to identify significant hits from Protocol_2 and GSAT
      (default is 14.0)

    -p2ht, --proto2-high-cutoff {float}
      Maximum threshold to identify significant hits from Protocol_2 and GSAT.
      Use this option along with -p2lt to focus the analysis on a range of scores.
     (default is unlimited)

    -tmo,  --min-tms-overlap {float}
      Minimum TMS overlap threshold to consider that a protocol2 hit is
      significant. By default it accepts any overlap grater than zero.

    -p2m, --proto2-min-aln {integer}
      Protocol2 parameter that indicates the minimum alignment length
      to keep. (Default is 80)

    -gr, --run-gsat (--no-run-gsat to negate)
      This option signals to run gsat using the parameters passed to
      options -g and -gm. By default gsat will not run in order to make
      the script finish faster and let the user look a the output to
      determine the best parameters to use for gsat.

    -gip, --gsat-ignore-proteins {string}
      List of comma-separated protein accessions that will be ignored
      when parsing Protocol2 hits and select hits on which to run GSAT.
      By default no hits are ignored.

    -gsh, --gsat-shuffles {integer >= 500}
      Number of shuffles to run GSAT (default 1000).

    -g, --gsat-cutoff
      Threshold to consider GSAT results significant.
      (default is 15)

    -gm, --gsat-mode {string}
      Indicate whether GSAT will run on the full protein sequence
      (mode full) or just the segements that aligned in protocol2
      (mode segment). Default value is: segment

    -h, --help
      Display the help of this program. Help will also be
      printed if no arguments are passed to the program.
