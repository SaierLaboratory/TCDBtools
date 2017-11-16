#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;

#To read command line options
use Getopt::Long;

#to check dependencies
use TCDB::CheckDependencies;
use TCDB::tcdb;



#==========================================================================
#This script runs protocol1, protocol2 and gsat to help
#determine whether two families are homologous.
#
#First, famXpander or protocol1 is run in order to extract an expanded
#list of candidate homologous proteins for each family.
#Second, protocol2 is run to determine homology and TMS topology.
#Only those pairs that have a protocol2 score above a user-specified
#threshold and minimum alignment length will be singled out to run gsat.
#If gsat has a z-score above a user-specified value it will be indicated.
#
#The script accepts the following command line options:
#
# -f1, --family1 {string}
#      TCDB ID string of the first family to be compared. It can be
#      any string, but if the option -x is given, a TCDB ID is assumed
#      (option is mandatory).
#
# -f2, --family2 {string}
#      TCDB of the second family to be compared. It can be
#      any string, but if the option -x is given, a TCDB ID is assumed
#      (option is mandatory).
#
# -u, --prog {string}
#      Specify which program will be used to find homologs, protocol1 (proto1)
#      or famXpander (fxpand). If fxpand is selected, input sequences must
#      be in fasta format. If proto1 is selected, input sequences must be in
#      2-column format (Id, protein sequence in one string).
#      (Option is mandatory)
#
# -d, --indir {string}
#      Path to the directory with the amino acid sequences of the family
#      members. If argument is not given, a directory will be created in the
#      current directory based on the info provided in the -f1 and
#      -f2 options (optional)
#
# -p1d, --proto1-dir {path}
#      Directory to place the results of protocol_1. This option
#      is specially useful when  more than one family will be compared.
#
# -p2d, --proto2-dir {path}
#      Directory where results of protocol_2 and GSAT will be stored.
#      If argument-p2d is not given, a directory will be created in
#      the current directory based on the info provided in the
#      -f1 and -f2 options (optional).
#
# -x, --extract (can be negated with --no-extract)
#      Flag Indicating that the amino acid sequences from both families should be
#      automatically extracted from TCDB. It is assumed that the
#      identifiers provided by the -f1 and -f2 options are TCDB family IDs.
#      If this option is not given, it is assumed that the user already
#      extracted the sequences and they are available in the path provided by
#      the -d option (optional)
#
# -ox, --only-extract-seqs
#      Flag instructing only to download the sequences of the input families
#      from TCDB. This is useful when prior to running protocol1 and protocol2
#      the user wants to limit or edit the protein sequnces. The program
#      aborts after retrieving the sequences. Option -u is mandatory in
#      in order to determine in what format the sequences will be
#      stored. This option is incompatible with -x.
#
# -n, --psiblast-it {integer}
#    Number of iterations that psi-blast will perform.
#    (Defauult is 1)
#
# -e, --evalue
#    Expectation value threshold for saving psiblast hits
#    Argument is optional (defult value 1e-3)
#
# -e2, --inc-evalue
#    E-value inclusion threshold for pairwise alignments in psi-blast.
#    (Default value is 1e-2)
#
# -k, --keep-aln-regions {T/F}
#    Option that indicates whether famXpander will feed only the psiblast
#    aligned regions to protocol2 (Default F).
#
# -a, --alignment-matches {int}
#    Number of blast matches to retrieve from famXpander or protocol1.
#   (Default 10000)
#
# -c, --min-coverage {float}
#    For famXpander, minimum alignment coverage of original sequence
#    (Default 0.8)
#
# -s, --min-rseq-length {float}
#    For famXpander, minimal sequence length relative to original seq length,
#    (Default 0.8)
#
# -l, --max-seq-length {float}
#    famXpander  maximal sequence length relative to original seq length.
#    The higher this value, the easier it will be to detect
#    protein fusions.
#    (Default 5.0 (ignored if -k is given)
#
# -r, --cdhit-cutoff {float}
#    Identity redundancy threshold for cd-hit.
#    (default 0.9)
#
# -p2lt, --proto2-low-cutoff {float}
#     Lower threshold to identify significant hits from Protocol_2 and GSAT
#     (default is 14.0)
#
# -p2ht,  --proto2-high-cutoff {float}
#     Max threshold to identify significant hits from protocol2 and run GSAT.
#     This option should be used to establish a range of protocol2 scores
#     to work with, which is useful in situations where there are hundreds of
#     high scoring but not relevant protocol2 hits (default value is unlimited).
#
# -tmo,  --min-tms-overlap {float}
#    Minimum TMS overlap threshold to consider that a protocol2 hit is
#    significant. By default it accepts any overlap grater than zero.
#
# -p2m, --proto2-min-aln {integer}
#     Protocol2 parameter that indicates the minimum alignment length
#     to keep. (Default is 60)
#
# -gr, --run-gsat (--no-run-gsat to negate)
#     This option signals to run gsat using the parameters passed to
#     options -g and -gm. By default gsat will not run in order to make
#     the script to finish faster and let the user look a the output and
#     determine the best parameters to use for gsat.
#
# -g, --gsat-cutoff
#     Threshold to consider GSAT results significant.
#     (default is 15)
#
# -gm, --gsat-mode {string}
#     Indicate whether GSAT will run on the full protein sequence
#     (mode full) or just the segements that aligned in protocol2
#     (mode segment). Default value is: segment
#
# -h, --help
#      Display the help of this program. Help will also be
#      printed if no arguments are passed to the program.
#
#==========================================================================



#==========================================================================
#Check dependencies


my @dependencies = ("mkdir", "grep", "wc", "extractFamily.pl",
		    "famXpander.pl",  "protocol2.py", "columnizeFaa.pl",
		    "gsat.py");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;




#==========================================================================
#Read command line options


#Identifiers for the families being compared (options -f1 and -f2)
my $fam_1 = "";
my $fam_2 = "";


#what program will be used to search for homologs: proto1 or fxpand
my $prog  = "";


#Path to the original sequences of the families being compared
#(option -d)
my $inseq_dir = "";


#If many families are going to be compared, it would be convenient to
#put all protocol_1/famXpander results in one place (option -p1d)
my $proto1_dir = "";


#Where the results of Protocol2 and GSAT will be stored (option -p2d)
my $proto2_dir = "";


#Indicate whether sequences for each family will be  automatically
#extracted from TCDB (option -x)
my $extract_seqs_from_tcdb = 0;


#Indicate whether only the sequences from the families should be
#extracted from tcdb, and then abort the program.
my $only_extract_seqs = 0;

#number of interation that psiblast will perform
my $psiblast_it = 1;


#For famXpander: E-value for considering psiblast hits (option -e)
my $evalue = 1e-7;


#For famXpander: E-value for including pairwise alignments in psiblast
#(option -e2)
my $iEvalue = 1e-5;


#Number of blast alignments to retrieve from famXpander or protocol1
my $num_aligns = 10000;


#Indicate whethere famXpander will report only align regions
#or complete proteins.
my $keep_aln_regions = 'F';


#For famXpander, minimum alignment coverage of original sequence.
my $sbj_min_cov = 0.8;


#For famXpander, minimal sequence length relative to original seq length.
my $tgt_min_length = 0.8;


#For famXpander, maximal sequence length relative to original seq length.
my $tgt_max_length = 5.0;


#For famXpander, identity threshold to determine cd-hit redundancy
my $cdhit_cutoff = 0.9;


#The z-score threshold to determine significant hits from
#protocol_2 (option -p2t)
my $proto2_low_cutoff = 14;


#The maximum z-score threshold to determine significant hits from
#protocol_2 (option -p2ht)
my $proto2_high_cutoff = 10000;


#The minimum overlap between TMS in a protocol2 hit
my $tmo_cutoff = 5;

#For protocol2, minimal alignment lent to consider a hit significant
#(option -p2m)
my $proto2_min_hit_len = 80;


#This variable indicates whether GSAT should run
my $run_gsat = 0;

my $gsat_shuffles = 1000;


#The z-score threshold to determine significant hits from
#GSAT (option -g)
my $gsat_cutoff = 15;


#Determine whether final GSAT step will run with whole proteins
#or just the segment that aligned in protocol_2 (option -gm)
my $gsat_mode = "segment";



read_command_line_arguments();


#print Data::Dumper->Dump([$fam_1, $fam_2, $prog, $inseq_dir, $proto1_dir, $proto2_dir,
#			   $extract_seqs_from_tcdb, $only_extract_seqs, $psiblast_it, $evalue, $iEvalue,
#			   $keep_aln_regions, $sbj_min_cov, $tgt_min_length, $tgt_max_length,
#			   $cdhit_cutoff, $proto2_low_cutoff, $proto2_high_cutoff, $tmo_cutoff,
#                          $proto2_min_hit_len, $gsat_cutoff, $gsat_mode],
#		    [qw(*fam_1 *fam_2 *prog *inseq_dir *proto1_dir *proto2_dir
#			*extract_seqs_from_tcdb *only_extract_seqs *psiblast_it *evalue *iEvalue
#			*keep_aln_regions *sbj_min_cov *tgt_min_length *tgt_max_length
#		        *cdhit_cutoff *proto2_low_cutoff *proto2_high_cutoff *tmo_cutoff
#                       *proto2_min_hit_len *gsat_cutoff *gsat_mode)]), "\n\n";
#
#print "Command line parsed!\n";
#exit;




#==========================================================================
#If the -x argument was given, download the sequences of each
#input family from TCDB.

TCDB::tcdb::extract_seqs_from_tcdb([$fam_1, $fam_2], $inseq_dir, $prog) if ($extract_seqs_from_tcdb || $only_extract_seqs);


#Check if user only wants to extract the sequences.
if ($only_extract_seqs) {
  print "Sequences extracted for families $fam_1 and $fam_2\n\n";
  exit;
}



#==========================================================================
#Get the number of sequences in each family


my ($seq_file_fam_1, $nseqs_fam_1) = @{ TCDB::tcdb::count_sequences_in_family($fam_1, $inseq_dir, $prog) };
my ($seq_file_fam_2, $nseqs_fam_2) = @{ TCDB::tcdb::count_sequences_in_family($fam_2, $inseq_dir, $prog) };


#Total number of pairs of proteins between both families.
#This is merely the number of family members in one family multiplied
#by the number of members in the other family, less all pairs that involve
#one protein with no blast hits.
my $total_FamPairs_cnt = $nseqs_fam_1 * $nseqs_fam_2;


#print Data::Dumper->Dump([$seq_file_fam_1, $nseqs_fam_1, $seq_file_fam_2, $nseqs_fam_2, $total_FamPairs_cnt],
#			 [qw(*seq_file_fam_1 *nseqs_fam_1 *seq_file_fam_1 *nseqs_fam_2 *total_FamPairs_cnt)]);
#exit;



#==========================================================================
#Run famXpander or protocol1 for both families


foreach my $fam ([$fam_1, $seq_file_fam_1], [$fam_2, $seq_file_fam_2]) {

  if ($prog eq "fxpand") {

    #run famXpander
    TCDB::tcdb::run_famXpander($fam->[0], $fam->[1], $proto1_dir, $evalue, $iEvalue,
			     $psiblast_it, $sbj_min_cov, $tgt_min_length,
			     $tgt_max_length, $cdhit_cutoff, $keep_aln_regions, $num_aligns);
  }
  else {

    #Run protocol1
    TCDB::tcdb::run_protocol1($fam->[1], $proto1_dir, $evalue, $cdhit_cutoff, $psiblast_it, $num_aligns);
  }
}

#print "Finished running $prog\n";
#exit;


#==========================================================================
#Run protocol_2 and GSAT for both families


print "\n\nRunning protocol_2 and GSAT...\n";

#Count the total number of pairs of proteins compared by protocol2.
my $real_FamPairs_cnt = 0;

#From Protocol2 results, count proteins pairs with high z-score
my $prot2_hsp_cnt = 0;


#Count pairs that had high enough GSAT score
my $gsat_hsp_cnt  = 0;


#Count the number of protein pairs that will be considered by protocol2
my $total_prot1_pairs = 0;


if ($prog eq 'fxpand') {

  #Run protocol2
  my $p2outdir = TCDB::tcdb::run_proto2($fam_1, $fam_2, $proto1_dir, $proto2_dir, $proto2_min_hit_len, \$total_prot1_pairs);


  #Select top protocol2 hits and run GSAT
  if ($run_gsat) {
    TCDB::tcdb::run_gsat_for_top_proto2_hits($p2outdir, $gsat_mode, $proto2_low_cutoff, $proto2_high_cutoff, $gsat_cutoff,
					   \$prot2_hsp_cnt, \$gsat_hsp_cnt, $prog, $proto1_dir, $tmo_cutoff, $gsat_shuffles);
  }
}


else {

    #Read the family sequences in column format
    my @fam1_seqs = ();
    TCDB::tcdb::read_2col_seqs($seq_file_fam_1, \@fam1_seqs);

    my @fam2_seqs = ();
    TCDB::tcdb::read_2col_seqs($seq_file_fam_2, \@fam2_seqs);


  FAM1:foreach my $f1 (@fam1_seqs) {

    FAM2:foreach my $f2 (@fam2_seqs) {

	print "Protocol2 between $f1->[0] vs $f2->[0]\n";

	#The individual protein IDs must not be the same
	next FAM2 if ($f1->[0] eq $f2->[0]);


	#Run protocol2
	my $p2outdir = TCDB::tcdb::run_proto2($f1->[0], $f2->[0], $proto1_dir, $proto2_dir, $proto2_min_hit_len, \$total_prot1_pairs);


	#Count this pair as analyzed (both had protocol_1 blast hits.
	$real_FamPairs_cnt++ unless ($p2outdir);


	#Select top protocol2 hits and run GSAT
	if ($run_gsat) {
	  TCDB::tcdb::run_gsat_for_top_proto2_hits($p2outdir, $gsat_mode, $proto2_low_cutoff, $proto2_high_cutoff, $gsat_cutoff,
						 \$prot2_hsp_cnt, \$gsat_hsp_cnt, $prog, [], $tmo_cutoff);
	}
      }

  }
}

#print Data::Dumper->Dump([$total_FamPairs_cnt, $total_prot1_pairs, $prot2_hsp_cnt, $gsat_hsp_cnt],
#			 [qw(*total_FamPairs_cnt *total_prot1_pairs *prot2_hsp_cnt *gsat_hsp_cnt)]);
#print "Protocol2 is finished\n";
#exit;



#Record final statistics in a file
open (STAT, ">", "outStats_${fam_1}_vs_${fam_2}.txt") || die $!;
print STAT "Total protein pairs between families: $total_FamPairs_cnt\n";


print STAT "Total pairs analyzed by protocol_2: $total_prot1_pairs\n";
print STAT "  Total pairs with significant protocol_2 score: $prot2_hsp_cnt\n";
print STAT "  ($prot2_hsp_cnt/$total_prot1_pairs) = ", $prot2_hsp_cnt/$total_prot1_pairs, "\n\n\n";



my $gsat_to_prot2_ratio = 0;
if ($prot2_hsp_cnt == 0) {
  $gsat_to_prot2_ratio = "Undefined";
}
else {
  $gsat_to_prot2_ratio = $gsat_hsp_cnt/$prot2_hsp_cnt;
}


print STAT "Total pairs with significant GSAT score: $gsat_hsp_cnt\n";
print STAT "   ($gsat_hsp_cnt/$total_prot1_pairs) = ", $gsat_hsp_cnt/$total_prot1_pairs, "\n\n\n";
print STAT "   ($gsat_hsp_cnt/$prot2_hsp_cnt) = ", $gsat_to_prot2_ratio, "\n\n\n";
close STAT;








#==========================================================================
#Read and validate the command line arguments


sub read_command_line_arguments {


  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }


  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "f1|family1=s"              => \&read_fam1,
      "f2|family2=s"              => \&read_fam2,
      "u|prog=s"                  => \&read_prog,
      "d|indir=s"                 => \$inseq_dir,
      "p1d|proto1-dir=s"          => \$proto1_dir,
      "p2d|proto2-dir=s"          => \$proto2_dir,
      "x|extract!"                => \$extract_seqs_from_tcdb,
      "ox|only-extract-seqs"      => \$only_extract_seqs,
      "n|psiblast-it=i"           => \$psiblast_it,
      "e|evalue=f"                => \$evalue,
      "e2|inc-evalue=f"           => \$iEvalue,
      "a|alinment-matches=i"      => \$num_aligns,
      "k|keep-aln-regions=s"      => \&read_keep_aln_regions,
      "c|min-coverage=f"          => \$sbj_min_cov,
      "s|min-rseq-length=f"       => \$tgt_min_length,
      "l|max-seq-length=f"        => \$tgt_max_length,
      "r|cdhit-cutoff=f"          => \$cdhit_cutoff,
      "p2lt|proto2-low-cutoff=f"  => \$proto2_low_cutoff,
      "p2ht|proto2-high-cutoff=f" => \$proto2_high_cutoff,
      "tmo|min-tms-overlap=f"     => \$tmo_cutoff,
      "p2m|proto2-min-aln=f"      => \$proto2_min_hit_len,
      "gr|run-gsat!"              => \$run_gsat,
      "g|gsat-cutoff=f"           => \$gsat_cutoff,
      "gm|gsat-mode=s"            => \&read_gsat_mode,
      "h|help"                    => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);



  #----------------------------------------------------------------------
  #Validate command line arguments


  #Input families
  die "Error: argument -f1 is mandatory.\n" unless ($fam_1);
  die "Error: argument -f2 is mandatory.\n" unless ($fam_2);


  #program to use for searching homologs
  die "Error: argument -u is mandatory.\n" unless ($prog);


  #The protocol2 score thresholds
  die "Error: value passed to -p2ht can't be lower than for option -p2lt" if ($proto2_high_cutoff < $proto2_low_cutoff);


  #Validate -d argument (the protein sequences)
  if($inseq_dir) {

    #If the specified directory does not exits, then the option -x or -ox must be given.
    unless (-d $inseq_dir) {
      unless ($extract_seqs_from_tcdb || $only_extract_seqs) {
	die "Error: Directory passed to -d does not exist, thus argument -x is mandatory.\n";
      }
      system "mkdir -p $inseq_dir";
    }
  }
  else {

    #If directory with proteins sequences for each family is not specified,
    #then sequences must be extracted from TCDB (options -x or -ox)
    unless ($extract_seqs_from_tcdb || $only_extract_seqs) {
      die "Error: Argument -d was not given, then argument -x is mandatory\n";
    }

    $inseq_dir = "./protein_seqs";
    system "mkdir $inseq_dir" unless (-d $inseq_dir);
  }


  if ($only_extract_seqs) {
    die "Error: options -x and -ox are not compatible\n" if ($extract_seqs_from_tcdb);
  }


  #Validate the -p1d argument (results for protocol1/famXpander)
  unless ($proto1_dir) {
    $proto1_dir = "./famXpander" if ($prog eq 'fxpand');
    $proto1_dir = "./protocol1"  if ($prog eq 'proto1');
  }
  system "mkdir -p $proto1_dir" unless (-d $proto1_dir || $only_extract_seqs);


  #Validate the -p2d argument
  $proto2_dir = "./proto2_${fam_1}_vs_${fam_2}" unless ($proto2_dir);
  system "mkdir -p $proto2_dir" unless (-d $proto2_dir || $only_extract_seqs);

}


#==========================================================================
#Read the -f1 option

sub read_fam1 {

  my ($opt, $value) = @_;

  TCDB::tcdb::validate_tcdb_id([$value]);

  $fam_1 = $value;
}


#==========================================================================
#Read the -f2 option

sub read_fam2 {

  my ($opt, $value) = @_;

  TCDB::tcdb::validate_tcdb_id([$value]);

  $fam_2 = $value;
}



#==========================================================================
#Read the -u option

sub read_prog {

  my ($opt, $value) = @_;

  unless (lc $value eq "proto1" || lc $value eq "fxpand") {
    die "Error: unknown progam passed to option $opt:  $value\n  Valid options are: proto1 or fxpand\n";
  }

  $prog = lc $value;
}




#==========================================================================
#Read the -gm option

sub read_gsat_mode {

  my ($opt, $value) = @_;

  unless (lc $value eq "segment" || lc $value eq "full") {
    die "Error: unknown mode passed to option $opt:  $value\n  Valid options are: segment or full\n";
  }

  $gsat_mode = $value;
}

#==========================================================================
#Read -k option

sub read_keep_aln_regions {
  my ($opt, $value) = @_;


  unless (lc $value eq "f" || lc $value eq "t") {
    die "Error: unknown switch passed to option $opt:  $value\n  Valid options are: T or F\n";
  }

  $keep_aln_regions = uc $value;
}


#==========================================================================
#This function will print error messages and/or the help to use this
#program. finally the program will exit the program



sub print_help {

  #
  # $errMsg: The error message to be diplayed
  #
  # $printHelp: boolean value indicating whether the help of the
  #             program will be displayed.
  #


  my $help = <<'HELP';

This script runs protocol1, protocol2 and gsat to help
determine whether two families are homologous.

First, famXpander or protocol1 is run in order to extract an expanded
list of candidate homologous proteins for each family.
Second, protocol2 is run to determine homology and TMS topology.
Only those pairs that have a protocol2 score above a user-specified
threshold and minimum alignment length will be singled out to run gsat.
If gsat has a z-score above a user-specified value it will be indicated.

The script accepts the following command line options:

 -f1, --family1 {string}
      TCDB ID string of the first family to be compared. It can be
      any string, but if the option -x is given, a TCDB ID is assumed
      (option is mandatory).

 -f2, --family2 {string}
      TCDB of the second family to be compared. It can be
      any string, but if the option -x is given, a TCDB ID is assumed
      (option is mandatory).

 -u, --prog {string}
      Specify which program will be used to find homologs, protocol1 (proto1)
      or famXpander (fxpand). If fxpand is selected, input sequences must
      be in fasta format. If proto1 is selected, input sequences must be in
      2-column format (Id, protein sequence in one string).
      (Option is mandatory)

 -d, --indir {string}
      Path to the directory with the amino acid sequences of the family
      members. If argument is not given, a directory will be created in the
      current directory based on the info provided in the -f1 and
      -f2 options (optional)

 -p1d, --proto1-dir {path}
      Directory to place the results of protocol_1. This option
      is specially useful when  more than one family will be compared.

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
      (Default F)

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

HELP


  print "$help\n";
  exit;
}



