#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

$Data::Dumper::Deepcopy = 1;

use Getopt::Long;
use LWP;
use Bio::SeqIO;

#Local libraries
use TCDB::CheckDependencies;
use TCDB::Assorted;

###########################################################################
#
# Search the famXpander results of a family for sequences with a
# specific number of TMS. The purpose is to search for evidence
# indicating which TMS were lost/gained when repeat units are not
# complete.
#
# TMSs are predicted with HMMTOP and the sequences with that number of
# TMS are compared against TCDB using GBLAST to get a nice graphical
# display of the similarities.
#
###########################################################################

#==========================================================================
#Check dependencies

my @dependencies = ("grep", "gblast3.py", "hmmtop", "blastdbcmd", "cd-hit");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;




#==========================================================================
#Read command line arguments

my $infile       = "";  #Expected to be the output of famXpander
my $family       = "";

my @proteins     = (); #Search for hits in psiblast.tbl envolving these proteins
my $fxpandir     = "";
my $psiblastFile = "";
my $fastaFile    = "";

my $outdir       = "";

#This will indicate whether sequences in $infile are complete are partial
#proteins (the fragments that aligned with blastp).
my $fullSeqs = 0;

my $nTMS     = undef;
my $evalue   = 1e-10;
my $coverage = 70;
my $covControl = "X";
my $cdhit      = 1;
my $clustID    = 0.8;  #Clustering identity value for cd-hit
my $runGBLAST  = "F";

read_command_line();

#print Data::Dumper->Dump([$infile, $family, \@proteins, $psiblastFile, $fastaFile, $outdir,
#                          $evalue, $coverage, $nTMS, $fullSeqs, $runGBLAST],
#                         [qw(*infile *family *proteins *psiblastFile *fastaFile *outdir
#                          *evalue *coverage *nTMS *fullSeqs *runGBLAST)]);
#exit;


#Create work directories
my $miscDir = "$outdir/misc";
system "mkdir -p $miscDir" unless (-d $miscDir);

my $gblastDir = "$outdir/gblast";
system "mkdir -p $gblastDir" unless (-d $gblastDir);

my $blastdbDir = "$outdir/blastdb";
system "mkdir -p $blastdbDir" unless (-d $blastdbDir);


if ($family) {}
elsif ($infile) {
  extractSequencesWithTMS($infile);
}
elsif (@proteins || $fxpandir) {
  extractHomologsForProteins(\@proteins, $psiblastFile, $fastaFile);
}
else {
  die "Something went terribly wrong!";
}




#==========================================================================
################   Subroutines definition beyond ths point   ##############
#==========================================================================


#==========================================================================
#Given a protein or list of proteins, extract their close homologs from
#file psiblast.tbl and extract their sequences.


sub extractHomologsForProteins {

  my ($prots, $psiFile, $faaFile) = @_;

  #get the query string
  my $qry = (@$prots)? join ("|", @$prots) : "";


  #Get the accessions of homologs
  my %accs = ();
  findGoodMatches($qry, $psiFile, \%accs);
  die "No homologs found for: $qry" unless (%accs);


  #Get the list of non-redundant accessions from fasta file
  my %nrAccs = ();
  readNRaccessions($faaFile, \%nrAccs);


  #Get the accessions for non-redundant homologs for the query proteins
  my %nrHom = ();
  map { $nrHom{$_} = 1 if (exists $nrAccs{$_}); } keys %accs;


  #save accessions to a file so their full sequences can be extracted.
  my $accFile1 = "$miscDir/acc1.txt";
  my @tmpAcc = keys %nrHom;
  saveAccToFile(\@tmpAcc, $accFile1);


  #Delete temporal variables
  undef %accs;
  undef %nrAccs;
  undef %nrHom;


  #Extract full sequences for accessions
  my $fullSeqsFile = "$miscDir/seqs1.faa";
  #remoteNCBIseqExtract($accFile1, $fullSeqsFile);
  localNCBIseqExtract($accFile1, $fullSeqsFile);

  filterSequencesByTMS($fullSeqsFile);
}



#==========================================================================
#Read fasta file to get the list of accessions for non-redundant
#sequences

sub readNRaccessions {

  my ($faaF, $out) = @_;

  open (my $fh, "<", $faaF) || die $!;
  while (<$fh>) {
    chomp;
    s/\s+$//;

    $out->{$1} = 1 if (/^\>(\w+)/);
  }
  close $fh;

}



#==========================================================================
#Parse PSIBLAST results and get all homologs for the query proteins

sub findGoodMatches {
  my ($str, $file, $out) = @_;

  open (my $fh, "<", $file) || die $!;
  while (<$fh>) {
    chomp;
    s/\s+$//;
    next unless ($_);
    next if (/^#/);

    if ($str) { next unless (/$str/); }

    my @cols = split (/\t/, $_);

    my $hacc = $cols[1];
    my $eval = $cols[4];
    my $qcov = ($cols[7]  - $cols[6] + 1) / $cols[8]  * 100;
    my $scov = ($cols[10] - $cols[9] + 1) / $cols[11] * 100;

    if ($eval <= $evalue && TCDB::Assorted::coverage_ok($qcov, $scov, $coverage, $covControl)) {
      $out->{$hacc} = 1;
    }
  }
  close $fh;
}




#==========================================================================
#Once famXpander is run, process results to find for sequences with the
#specified number of TMS.
#
# Input sequences must contain the correct format in the fasta header!
# (e.g. '>RefSeq  Annotations....' or just ''>RefSeq')
#

sub extractSequencesWithTMS {

  my $file = shift;

  my $fullSeqsFile = "";


  #----------------------------------------------------------------------
  #Get the full length sequences to run HMMTOP

  if ($fullSeqs) {
    $fullSeqsFile = $file;
  }
  else {

    #Extract accessions from input sequence file
    my @acc1 = ();
    extractAccFromSeqFile($file, \@acc1);


    #save accessions to a file so their full sequences can be extracted.
    my $accFile1 = "$miscDir/acc1.txt";
    saveAccToFile(\@acc1, $accFile1);


    #Extract full sequences for accessions
    $fullSeqsFile = "$miscDir/seqs1.faa";
    #remoteNCBIseqExtract($accFile1, $fullSeqsFile);
    localNCBIseqExtract($accFile1, $fullSeqsFile);
  }


  filterSequencesByTMS($fullSeqsFile);
}




sub filterSequencesByTMS {

  my $fullSeqs = shift;


  #----------------------------------------------------------------------
  #Run HMMTOP and select hits with the right number of TMS

  my $hmmtopFile = "$miscDir/tms.hmmtop";
  my $cmd1 = qq (hmmtop -if=$fullSeqs -of=$hmmtopFile -sf=FAS -pi=spred -is=pseudo);
  system $cmd1 unless (-f $hmmtopFile && !(-z $hmmtopFile));


  #parse TMS number and get the accessions of proteins with the specified number of TMSs
  my @hmmTopHits = ();
  parseHMMTOP($hmmtopFile, \@hmmTopHits);



  #----------------------------------------------------------------------
  #Extract the sequences for the proteins with the right number of TMS

  #save accessions to a file so their full sequences can be extracted.
  my $accFile2 = "$miscDir/accWithExpectedTMS.txt";
  saveAccToFile(\@hmmTopHits, $accFile2);


  #Extract full sequences for accessions
  my $fullSeqsFile2 = "$miscDir/seqsWithExpectedTMS.faa";
  localNCBIseqExtract($accFile2, $fullSeqsFile2, $fullSeqs);


  #----------------------------------------------------------------------
  #Run GBLAST with the selected sequences with the expected number of TMS

  my $cmd2 = qq(gblast3.py -i $fullSeqsFile2 -o $gblastDir --evalue=$evalue --cov=$coverage);
  system $cmd2 unless (-f "$gblastDir/results.tsv" || $runGBLAST eq 'F');
}



#==========================================================================
#Parse HMMTOP output
 
sub parseHMMTOP {

  my ($file, $out)  = @_;

  #open HMMTOP file and scan results for the desired number of TMSs.
  open (my $fh, "<", $file) || die $!;
  while (<$fh>) {
    chomp;
    if (/\>HP:\s+\d+\s+(\w+).+\s+(IN|OUT)\s+(\d+)/) {

      my $acc = $1;
      my $tms = $3;
      push (@{ $out }, $acc) if (!(defined $nTMS) || $tms == $nTMS);
    }
  }


  #Verify if sequences with the specified number of TMS were found.
  unless (scalar @{ $out } > 0) {
    print " ** No sequences found with $nTMS TMSs. **\n";
    system "touch $outdir/NO_SEQUENCES_WITH_${nTMS}_TMSs";
    exit;
  }
}



#==========================================================================
#Extract sequences from the local copy of the NR blast database

sub localNCBIseqExtract {

  my ($in, $out, $refSeqsFile) = @_;

  #The default is to use the local NR copy
  my $db = 'nr';

  #Create Blast database with the sequences
  if ($refSeqsFile) {
    my $mdb = qq(makeblastdb -dbtype prot -in $refSeqsFile -title seqsdb -parse_seqids -hash_index -out $blastdbDir/seqsdb);
    system $mdb unless (-f "$blastdbDir/seqsdb.pin");
    $db = "$blastdbDir/seqsdb";
  }

  my $randNum = int(rand(10000));
  my $tmpFile ="$miscDir/seqs${randNum}.faa";
  my $cmd = qq(blastdbcmd -db $db -entry_batch $in -target_only -outfmt '\%f' -out $tmpFile);
  system $cmd unless (-f $tmpFile && !(-z $tmpFile));


  #Remove sequences with non-standard residues
  open (my $fh, "<", $tmpFile) ||  die $!;
  my $seqs = do {local $/;  <$fh>};
  close $fh;
  die "Could not read sequences from file: $tmpFile" unless ($seqs);


  if ($cdhit) {
    my $tmpCleanSeqs = "$miscDir/clean${randNum}.faa";
    getSeqsWithoutErrors($seqs, $tmpCleanSeqs);
    die "No sequence file found:  $tmpCleanSeqs" unless (-f $tmpCleanSeqs && !(-z $tmpCleanSeqs));

    run_cdhit($tmpCleanSeqs, $out);
    system "rm $tmpCleanSeqs" if (-f $tmpCleanSeqs);
  }
  else {
    getSeqsWithoutErrors($seqs, $out);
  }
  die "No sequences found in local NCBI DB database!\n" unless (-f $out && !(-z $out));


  #Remove the temporal file
  system "rm $tmpFile" if (-f $tmpFile);

}



#==========================================================================
#Run CD-HIT on extracted sequences


sub run_cdhit {

  my ($cd_infile, $cd_outfile) = @_;

  my $cmd = qq(cd-hit -i $cd_infile -o $cd_outfile -c $clustID -g 1 -d 0);
  system $cmd unless (-f $cd_outfile);
  system "rm ${cd_outfile}.clstr" if (-f "${cd_outfile}.clstr");

}


#==========================================================================
#Remote extraction of sequences from NCBI (using the POST method)


sub remoteNCBIseqExtract {

  my ($in, $out) = @_;

  #Don't do anything if output file already exists.
  return if (-f $out && !(-z $out));

  #prepare accession list for http request
  open (my $fh, "<", $in) ||  die $!;
  chomp(my @accs = <$fh>);
  close $fh;

  #parameters
  my $accsList = join (",", @accs);
  my $db       = "protein";
  my $format   = 'fasta';

  #Example URL for the query:
  #   https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=5&rettype=fasta
  my $baseURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";


  #stantiate the user agent object
  my $browser = LWP::UserAgent->new;
  $browser->show_progress(1);


  #prepare data for query
  my %data = ( db => $db, id => $accsList, rettype => $format, retmode => 'text');


  #Query NCBI
  my $response = $browser->post($baseURL, \%data);


  #Process request output
  if ($response->is_success) {
    my $seqs = $response->content;
    getSeqsWithoutErrors($seqs, $out);
  }
  else {
    die $response->status_line;
  }
}




sub getSeqsWithoutErrors {
  my ($seqs, $out) = @_;

  my $inObj  = Bio::SeqIO->new(-string => $seqs, -format => 'fasta');
  my $outObj = Bio::SeqIO->new(-file => ">$out", -format => 'fasta');

  while (my $seq = $inObj->next_seq) {
    if (isSequenceCorrect($seq->seq)) {
      $outObj->write_seq($seq);
    }
  }
}



#Check if a protein sequence contains non standard amino acids.
sub isSequenceCorrect {
  my $seqString = shift;
  return ($seqString =~ /^[GAVLIPMFWSTNQYCDEHKR]+$/)? 1 : undef;
}



#==========================================================================
sub saveAccToFile {

  my ($accList, $oFile) = @_;

  unless (-f $oFile && !(-z $oFile)) {
    open (my $fh2, ">", $oFile) || die $!;
    print $fh2 join ("\n", @{ $accList }), "\n";
    close $fh2;

    die "Could not save Accessions to file: $oFile" unless (-f $oFile && !(-z $oFile));
  }
}



#==========================================================================
#Extract accessions from sequence file in fasta format

sub extractAccFromSeqFile {
  my ($in, $out) = @_;

  open (my $fh, "<", $in) || die $!;
  while (<$fh>) {
    chomp;
    if (/^\>([0-9a-zA-Z_-]+)/) {
      push (@{ $out }, $1);
    }
  }
  close $fh;
}



#===========================================================================
#Read command line and print help


sub read_command_line {

  print_help() unless (@ARGV);

  my $status = GetOptions(
	"i|infile=s"       => \&read_infile,
	"f|family=s"       => \&read_family,
	"p|proteins=s"     => \&read_proteins,
	"d|fxpandir=s"     => \&read_fxpandir,
	"o|outdir=s"       => \&read_outdir,
	"e|evalue=f"       => \$evalue,
	"c|coverage=f"     => \$coverage,
	"cc|cov-control=s" => \&read_covControl,
	"full!"            => \$fullSeqs,
	"n|tms=i"          => \$nTMS,
	"g|gblast=s"       => \&read_gblast,
	"h|help"           => sub { print_help(); },
	"<>"               => sub { die "Error: Unknown argument: $_[0]\n"; });
  exit unless ($status);


  #Check for incompatibilities and errors
  unless ($infile || $family || $fxpandir) {
    die "Error: only one option -i, -f, -d or -p must be given!\n";
  }
  die "Error: options -i, -f and -p are incompatible\n" if ($infile && $family);
  die "Error: options -i, -f and -p are incompatible\n" if ($infile && @proteins);
  die "Error: options -i, -f and -p are incompatible\n" if ($family && @proteins);


  #Options -p and -d must go together
  if (@proteins && !($fxpandir)) {
    die "Error: if option -p is given, then option -d is mandatory --> ";
  }



  #Number of TMS must be given
#  die "Error: option -n is mandatory!" unless (defined $nTMS);


  #Create output directory
  $outdir = "seqs_with_${nTMS}_TMS" unless ($outdir);
  system "mkdir -p $outdir" unless (-d $outdir);

}


#==========================================================================
#Option -i

sub read_infile {
  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
    die "Error in option -i: File with sequences does not exist or is empty!\n";
  }

  $infile = $value;
}

#==========================================================================
#Option -f

sub read_family {
    my ($opt, $value) = @_;

    TCDB::Assorted::validate_tcdb_id($value);
    $family = $value;
}


#==========================================================================
#Option -p

sub read_proteins {
    my ($opt, $value) = @_;

    @proteins = split(/,/, uc $value);
}


#==========================================================================
#Option -d

sub read_fxpandir {
  my ($opt, $value) = @_;

  unless (-d $value) {
    die "Error in option -d: famXpander  dir does not exist -> $value!\n";
  }

  $psiblastFile = "$value/psiblast.tbl";
  $fastaFile    = "$value/results.faa";

  unless (-f $psiblastFile && -f $fastaFile) {
    die "Error: files psiblast.tbl and results.faa must exist in famXpander dir -> $value\n";
  }

  $fxpandir = $value;
}


#==========================================================================
#Option -cc

sub read_covControl {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^[XQSB]$/) {
    die "Error in option -cc: illegal charater ($value). Valid characters are Q,S,B,X\n";
  }

  $covControl = $tmp;
}


#==========================================================================
#Option -g

sub read_gblast {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^[TF]$/) {
    die "Error in option -g: illegal charater ($value). Valid characters are T,F\n";
  }

  $runGBLAST = $tmp;
}

#==========================================================================
#Option -o

sub read_outdir {
  my ($opt, $value) = @_;

  $outdir = $value;
}


#==========================================================================
#option -h


sub print_help {

    my $help = <<'HELP';

 earch the famXpander results of a family for sequences with a
 specific number of TMS. The purpose is to search for evidence
 indicating which TMS were lost/gained when repeat units are not
 complete.

 ** NOTE: one option -i, -f or -p options must be provided.

 Command line options:

 -i, --infile {file} (Optional)
   Path to outoput file from famXpander: results.faa
   (incompatible with option -f)

 -f, --family {TCID} (optional)
   TCID for the family to analyize. The program famXpander will be run
   for this family to extract all available homologous.
   (incompatible with option -i).

 -p, --proteins {accessions} (optional)
   List of protein accessions from TCDB or UniProt for which to identify
   homologous. Option -d must also be given with this option.

 -d, --fxpandir {directory} (optional)
   Directory holding the results of famXpander. This option must be given
   if option -p is given.

 -n, --tms {integer >= 0}  (Optional; Default: any number of TMS)
   Desired nubmer of TMS in the target sequences. If not given,
   The program will take all sequences into consideration.

 -o, --outdir {directory}
   Output directory where results will be organized. If directory does not
   exists it will be created.
   (Default: seqs_with_N_tms; where N is taken from the -n option)

 -e, --evalue {FLOAT >= 0.0) (Optional)
   E-value cutoff to identify sequence matches (Default: 1e-10).

 -c, --coverage {FLOAT < 100.0) (Optional; Default: 70)
   Minimum coverage for the smaller protein (query or subject) in percentage
   values

-cc, --cov-control {char} (Optional. Defaul: X)
   Controls how the coverage threshold will be applied to the alignments
   (case insensitive):
   X:  Coverage applies to at least one protein.
   B:  Coverage applies to both proteins.
   Q:  Coverage applies to the query protein only.
   S:  Coverage applies to the subject protein only.

 --full  (Optional)
   Flag that indicates that the input file (passed through the -i option) 
   contains full sequences and not the aligned fragments.
   (Default: --no-full; assume input file contains only aligned fragments)

 -g, --gblast {FLAG: T/F} (Optional; Default: F)
   Flag indicating whether GBLAST should run for the extracted sequences.

 -h, --help
   Display this help. Also displayed if script is run without arguments.

HELP

    print $help;
    exit;
}
