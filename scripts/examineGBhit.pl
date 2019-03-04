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

read_command_line_arguments();
#print Data::Dumper->Dump([$query, $sAcc, $sTC, $outdir, $owBlastDB],
#                        [qw(*query *sAcc *sTC *outdir *owBlastDB)]);
#exit;


#==========================================================================
#First make sure that the TCDB blast file exists

my $blastDir = "$ENV{HOME}/db/blastdb";
system "mkdir -p $blastDir" unless (-d $blastDir);

my $cmd1 = "extractFamily.pl -i tcdb -f blast -o $blastDir";
system $cmd1 if ($owBlastDB || !(-f "$blastDir/tcdb.pin"));



#==========================================================================
#First download the query sequence

my $qSeqFile = "$outdir/${query}.faa";
my $cmd2 = qq(blastdbcmd -db nr  -entry $query -target_only > $qSeqFile);
system $cmd2 unless (-f $qSeqFile);


#==========================================================================
#Second dowload the TCDB protein sequence

my $sSeqFile = "$outdir/${sAcc}.faa";
my $sID = "${sTC}-$sAcc";
my $cmd3 = qq(blastdbcmd -db tcdb -entry $sID -target_only > $sSeqFile);
system $cmd3 unless (-f $sSeqFile);


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

-h, --help
  Print this help. This option takes precedence over anyother option.

HELP

    print $help;
    exit;
}
