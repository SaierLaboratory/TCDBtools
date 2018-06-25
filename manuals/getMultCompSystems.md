# Documentation for script: _getMultCompSystems.pl_

## Summary
This program characterizes the set of multicomponent systems present in a fully sequenced genome. 
This complements gblast, which targets the best matches between the query genome and TCDB. 
All multicomponent systems available in TCDB are downloaded and BLASTed against the genome 
allowing users to specify the E-value cutoff, minimal coverage of the query and subject sequences, 
and more. Hydropathy plots are generated depicting the regions in the query and subject that are 
aligned. The program also determines the genomic context of homologs of each of the components 
within a system by identifying genes in close neighborhood, thus aiding in the process of selecting 
components that may be parts of the same operon.


## Dependencies
The following programs need to be available in your path for this 
program to run properly:

1. **_blast+ 2.4.0 to 2.6.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

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

-gdir, --genome-dir { path }  (Mandatory)
  Directory with all the files as downloaded from NCBI for the genome under analysis.

-gacc, --genome-accesssion { string } (Mandatory)
  The prefix of all the files in the genome directory. This is normally the name of
  the folder in NCBI. For example:  GCA_000995795.1_ASM99579v1

-rs, --replicon-structure {string}  (Optional; Default: circular)
  The structure of the replicon under analysis. Only the following values
  are accepted: circular or linear

-gb, --gblast { file } (Mandatory)
   The GBLAST output file in tab-separated format.

-dbn, --blastdb-name {string}  (Optional; Default: proteome)
  Name of the blast database that will be generated with the proteome
  of the genome.

-bo, --tcblast-overwrite  (Optional; by default the DB is not overwritten)
  If given, the blast DB with all the protein content in TCDB will be overwritten.
  The location of this database is ~/db/blastdb. By default the database is not
  overwritten, so be careful and see if your current database is not too old.

-of, --output-format { string } (default: tsv)
  The format of the text file that will be generated with the final report.
  By default the output is tab-separated values (tsv), but it can also be
  comma-separated values (csv).

-o, --outdir { path } (Optional; Default: MultiComponentSystemsAnalysis)
  Directory where the results of the program will be stored.

-d, --max-gene-dist { int } (Optional; Default: 15)
  Maximum distance in number of genes to consider two genes as in the
  neighborhood.

-f, --fusion-len-ratio { float } (Optional; Default: 1.8)
  Minimal length ratio beteen the query and subject proteins that specifies
  al least how much larger one protein must be must be relative to the other
  in order to consider the larger protein a fusion protein.
  (This option is experimental and may disappear

-lp, --large-protein { int } (Optional; Default: 200)
  Minimal length to identify large proteins.

-sp, --small-protein { int } (Optional; Default: 100)
  Maximal length to consider a protein small.

-hc, --high-coverage { float } (Optional; Default: 50.0);
  Assumption of good coverage between query and subject to consider
  two proteins homologous (value must be grater than for option -c).

-c, --mid-coverage { float } (Optional; Defaul: 40.0)
  Assumption of OK coverage to consider two proteins homologous
  (value must be greater than for option -mc).

-mc, --min-coverage { float } (Optional;  Default: 30.0)
  Assumption of poor coverage between two proteins but still
  it may be worth taking a look at the alignment
  (Value must be larger than 20.0)

-id, --identity { float }  (Optional;  Default: 20.0)
  Minimal alignment identity to consider a hit in the results.

-e, --evalue { float } (Optional;  Default: 1e-3)
  Worst acceptable evalue for alignments with high coverage.

-a, --min-aln-length { int } (Optional; Default: 55)
  Minimal accepted alignment length.

-h, --help
  Print this help. This option takes precedence over anyother option.