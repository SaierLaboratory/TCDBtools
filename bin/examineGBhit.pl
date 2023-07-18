#!/usr/bin/env perl -w

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use TCDB::Assorted;


#==========================================================================
#  This program is a wrapper to alignSeqFiles.pl where the second
#  sequence is expected to be a TCDB protein, as is the case when
#  examining GBLAST hits.
#==========================================================================


my $query = "";
my $sAcc  = "";
my $sTC   = "";

my $outdir = ".";
my $owBlastDB = 0;
my $subMatrix = "BL50";
my $qblastdb = 'nr';

my $blastBin = "/usr/local/bin";

read_command_line_arguments();
#print Data::Dumper->Dump([$query, $sAcc, $sTC, $outdir, $owBlastDB], $qblastdb,
#                        [qw(*query *sAcc *sTC *outdir *owBlastDB)]);
#exit;


#==========================================================================
#First make sure that the TCDB blast file exists

my $blastDir = "$ENV{HOME}/db/blastdb";
system "mkdir -p $blastDir" unless (-d $blastDir);

my $cmd1 = "extractTCDB.pl -i tcdb -f blast -o $blastDir";
system $cmd1 if ($owBlastDB || !(-f "$blastDir/tcdb.pin"));



#==========================================================================
#First download the query sequence

#
#NOTE: If $query and $sAcc are the same accession, it will be necessary to modify
#      the accession of the query in the sequence file to prevent an error in
#      alignSeqFiles.pl
#

my $tmpAcc = ($query eq $sAcc)? "${query}_tmp" : $query;


my $qSeqFile = "$outdir/${tmpAcc}.faa";
my $cmd2 = qq($blastBin/blastdbcmd -db $qblastdb  -entry $query -target_only > $qSeqFile);
system $cmd2 unless (-f $qSeqFile);
die "Could not extract query sequence: $query" unless (-f $qSeqFile && !(-z $qSeqFile));


if ($query eq $sAcc) {
  my $cmd = qq(perl -i.bkp -pe 's/$query\\.*\\d*/$tmpAcc/;' $qSeqFile);
}


#==========================================================================
#Second download the TCDB protein sequence

my $sSeqFile = "$outdir/${sAcc}.faa";
my $sID = "${sTC}-$sAcc";
my $cmd3 = qq($blastBin/blastdbcmd -db tcdb -entry $sID -target_only > $sSeqFile);
system $cmd3 unless (-f $sSeqFile);
die "Could not extract TCDB sequence: $sID" unless (-f $sSeqFile && !(-z $sSeqFile));


#==========================================================================
#Compare the two sequences now

my $alnDir  = "$outdir/ssearch_${query}_vs_$sAcc";
my $repFile = "$alnDir/report.html";
my $cmd4 = qq(alignSeqsFiles.pl -q $qSeqFile -ql $query -s $sSeqFile -sl $sID -o $alnDir -e 1 -c 20 -cc X -m $subMatrix);
system $cmd4 unless (-f $repFile);


###########################################################################
##                         SUBROUTINES                                   ##
###########################################################################


sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "q=s"    => \&read_query,
      "s=s"    => \&read_sAcc,
      "t=s"    => \$sTC,
      "o=s"    => \$outdir,
      "bdb=s"  => \&read_qblastdb,
      "w=s"    => \&read_owBlastDB,
      "m=s"    => \&read_subMatrix,
      "h|help" => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"     => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);


  system "mkdir -p $outdir" unless (-d $outdir);
}


#==========================================================================
#Read the -q option

sub read_query {

  my ($opt, $value) = @_;

  my $tmp = $value;
  $tmp =~ s/\.\d*$//;

  $query = $tmp;
}


#==========================================================================
#Read the -s option

sub read_sAcc {

  my ($opt, $value) = @_;

  my $tmp = $value;
  $tmp =~ s/\.\d*$//;

  $sAcc = $tmp;
}



#==========================================================================
#Read the -w option

sub read_owBlastDB {

  my ($opt, $value) = @_;

  my $tmp = uc $value;
  if ($tmp eq "T") {
    $owBlastDB = 1;
  }
  elsif ($tmp eq "F") {
    $owBlastDB = 0;
  }
  else {
    die "Invalid argument ($tmp). Valid values for -w are: T or F"
  }
}


#==========================================================================
#Read the -bdb option. BlastDB for extracting query sequences

sub read_qblastdb {

  my ($opt, $value) = @_;

  my $tmp = "${value}.pin";

  unless  ($value eq 'nr') {
    unless (-f $tmp  && !(-z $tmp)) {
      die "BlastDB for query search not found: $value";
    }
  }
  $qblastdb = $value;
}





#==========================================================================
#Option -m (Any matrix supported by ssearch)

sub read_subMatrix {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^(BL50|BL62|P250|OPT5|VT200|VT160|P120|VT120|BL80|VT80|MD40|VT40|MD20|VT20|MD10|VT10)$/) {
    die "Error in option -$opt: illegal matrix ($value). Value should be any matrix supported by SSEARCH\n";
  }

  $subMatrix = $tmp;
}





#==========================================================================

sub print_help {

    my $help = <<"HELP";

Align the query and target proteins allowing for multiple hits. This will
help to spot potential gene duplications when analyzing gblast hits.

Options:

-q  {String} (Mandatory)
    Accession of Query protein.

-s  {String} (Mandatory)
    Subject accession (Uniprot or RefSeq of protein in TCDB).

-t  {String} (Mandatory)
    Subject system TCID (full TCID of the accession passed to option -s)

-o  {PATH} (Default: .)
    Output directory where results will be stored.

-w  {T/F} (Default: F)
    Overwrite the contents of the local TCDB blast database.

-m  {string} (Optional. Default: BL50)
   Substitution matrix to use in the alignments. Any matrix supported by
   ssearch36 can be used.

-bdb {string} (Optional. Dafault: nr)
   Blast database that will be used to extract the sequence of the query
   protein.

-h, --help
  Print this help. This option takes precedence over anyother option.

HELP

    print $help;
    exit;
}
