#!/usr/bin/env perl

no warnings;
use strict;
use Data::Dumper;

use TCDB::Repeats;
use Getopt::Long;


my $seqsDir   = undef;      #'/Users/amedrano/Desktop/Mai_tmsRepeat/sequences';
my $seqsFile  = undef;
my $tmsFile   = undef;      #'/Users/amedrano/Desktop/Mai_tmsRepeat/tms.hmmtop';
my $outDir    = "Repeats";  #'/Users/amedrano/Desktop/Mai_tmsRepeat/RepeatUnits/ResultsOOP';

my $evalue    = 1e-2;
my $coverage  = 0.85;
my $identity  = 0.2;

my @tmsRanges = ();

read_command_line();


#print Data::Dumper->Dump([$seqsDir, $seqsFile, $tmsFile, $outDir],
#			 [qw(*seqsDir *seqsFile *tmsFile *outDir )]);
#exit;


#my $repObj = TCDB::Repeat->new('seqsDir' => $seqsDir,
#			       'tmsFile' => $tmsFile,
#			       'outDir'  => $outDir,
#			      'ranges2searchTMS' => \@TMSranges);


my @TMSranges = ([1, 3], [4, 6]);

my $repObj = TCDB::Repeat->new();

#$repObj->tmsFile($tmsFile);
#$repObj->seqsDir($seqsDir);
$repObj->seqsFile($seqsFile);
$repObj->outDir($outDir);
$repObj->evalueCutoff($evalue);
$repObj->identityCutoff($identity);
$repObj->coverageCutoff($coverage);
$repObj->TMSranges2search(\@TMSranges);

$repObj-> findRepeatsTMSranges();

#print Data::Dumper->Dump([$repObj ], [qw(*repObj)]);





#===========================================================================
#Read command line and print help


sub read_command_line {

  print_help() unless (@ARGV);

  my $status = GetOptions(
	"s|seqs-file=s" => \&read_seqsFile,
        "d|seqs-dir=s"  => \&read_seqsDir,
	"o|outdir=s"    => \&read_outdir,
	"t|tms=s"       => \&read_tmsFile,
	"e|evalue=f"    => \$evalue,
	"i|identity=f"  => \$identity,
	"c|coverage=f"  => \$coverage,
	"h|help"        => sub { print_help(); },
	"<>"            => sub { die "Error: Unknown argument: $_[0]\n"; });
  exit unless ($status);


  #Validadte input file option
  die "Error: no sequence file detected!" unless ($seqsFile);
}


#==========================================================================
#Option -s

sub read_seqsFile {
  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
    die "Error: file with sequences does not exist or is empty!\n";
  }

  $seqsFile = $value;
}


#==========================================================================
#Option -t

sub read_tmsFile {
  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
    die "Error in option -t: File with TMSs (hhmtop output) does not exist or is empty!\n";
  }

  $tmsFile = $value;
}


#==========================================================================
#Option -d

sub read_seqsDir {
  my ($opt, $value) = @_;

  die "Error: directory with sequences does not exist." unless (-d $value);

  $seqsDir = $value;
}


#==========================================================================
#Option -o

sub read_outdir {
  my ($opt, $value) = @_;

  $outDir = $value;
}


#==========================================================================
#option -h


sub print_help {

    my $help = <<'HELP';

This program searches for reapeats between different user-specified 
regions of proteins.

 Command line options:

 -s, --seqs-file {file} (mandatory)
   Path to file in fasta format with all the input sequences.
   THis option is incompatible with option -d. But one of the
   two option must be given.

 -d, --seqs-dir {path} (optional)
   Path to directory where the input sequences are located.
   This option is incompatible with options -s. But one of the
   two option must be given.

 -o, --oudir {papth} (optional)
  Path to the output directory.
  (Default: ./tmsRepeat)

 -t, --tms {file} (optional)
   File with the output of hmmtop for the input sequences, if available.
   (Default: run hmmtop on input seqeunces)

 -e, --evalue {float} (optional)
   Maximal evalue cutoff for the aligned seqments.
   (Default: 0.001)

 -i, --identity {float} (optional)
   Minimal identity in aligned regions.
  (Default: 0.2)

 -c, --coverage {float) (optional)
   Minimal coverage cutoff within the range: [0, 1] for the coverage of aligned regtions.
   (Default: 0.85)

 -h, --help
   Display this help. Also displayed if script is run without arguments.

HELP

    print $help;
    exit;
}
