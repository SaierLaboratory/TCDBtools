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


#==========================================================================
#  Extract sequences from NCBI locally or remotely
#==========================================================================

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

my @dependencies = ("blastdbcmd");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;


my $accession = "";
my $accFile   = "";
my $outfile   = "";
my $blastdb   = "nr";
my $mode      = "local";


read_command_line_arguments();
#print Data::Dumper->Dump([$acc, $accFile, $outfile, $blastdb, $mode],
#                        [qw(*acc *accFile *outfile *blastdb *mode)]);
#exit;




###########################################################################
##                         SUBROUTINES                                   ##
###########################################################################


#==========================================================================
#extract sequences remotely

sub remoteNCBIseqExtract {

  my ($in, $out) = @_;

  #Don't do anything if output file already exists.
  return if (-f $out && !(-z $out));

  #prepare accession list for http request
  open (my $fh, "<", $in) ||  die $!;
  chomp(my @accs = <$fh>);
  close $fh;

  #parameters
  my $accsList = '"' . join (",", @accs) . '"';
  my $db       = "protein";
  my $format   = 'fasta';

  #Example URL for the query:
  #   https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=5&rettype=fasta
  my $baseURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";


  #instantiate the user agent object
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
#Determine if sequence has nonstandard residues

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



#==========================================================================
#Check if a protein sequence contains non standard amino acids.

sub isSequenceCorrect {
  my $seqString = shift;
  return ($seqString =~ /^[GAVLIPMFWSTNQYCDEHKR]+$/)? 1 : undef;
}



#==========================================================================
#Read command line options

sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
    "a=s"    => \&read_accession,
    "f=s"    => \&read_infile,
    "db=s"   => \$blastdb,
    "o=s"    => \$outfile,
    "m=s"    => \&readMode,
    "h|help" => sub { print_help(); },

    #For arguments that do not look like valid options
    "<>"     => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);

  die "Option -a or -f must be given." (unless ($accession || $infile);

}


#==========================================================================
#Read the -q option

sub read_accession {

  my ($opt, $value) = @_;

  my $tmp = $value;
  $tmp =~ s/\.\d*$//;

  $accession = $tmp;
}



#==========================================================================
#Read the -s option

sub read_infile {

  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
    die "Input file with accessions not found or empty."
  }

  $infile = $value;
}



#==========================================================================
#Read the -w option

sub read_mode {

  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /R|L/) {
    die "Unknown operation mode. Valid options are: R|L";
  }
  $mode = $tmp;
}



#==========================================================================

sub print_help {

    my $help = <<"HELP";

Retrieve protein sequences from NCBI/TCDB locally or remotely.


Options:

-a  {String} (optional)
    Accession of Query protein for which to extract a sequence.
    Option -a and -f are mutually exclusive, but one must be given.

-f  {FILE} (optional)
    Input file with accessions (one per line) for which sequences will
    be extracted. Option -a and -f are mutually exclusive, but one must be given.

-db {String} (Mandatory)
    Blast DB that will be used to extract the sequences

-o  {FILE} (Default: seqs.out)
    Output file where sequences will be saved.

-m  {string} (Optional; Default: L)
    Mode of operation. Specify if sequences will be retrieved locally (L)
    or remotely (R).

-h, --help
  Print this help. This option takes precedence over anyother option.

HELP

    print $help;
    exit;
}
