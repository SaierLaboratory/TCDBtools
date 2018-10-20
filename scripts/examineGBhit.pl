#!/usr/bin/env perl -w

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use TCDB::Assorted;



my $query = "";
my $sAcc  = "";
my $sTC   = "";

my $outdir = ".";
my $owBlastDB = 0;


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
my $cmd4 = qq(alignSeqsFiles.pl -q $qSeqFile -ql $query -s $sSeqFile -sl $sID -o $alnDir -e 10 -c 20 -cc X);
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
    Subject TCID (full TCID of the accession passed to option -s)

-o  {PATH} (Default: .)
    Output directory where results will be stored.

-w  {T/F} (Default: F)
    Overwrite the contents of the local TCDB blast database.

-h, --help
  Print this help. This option takes precedence over anyother option.

HELP

    print $help;
    exit;
}
