#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

use TCDB::Assorted;


use Getopt::Long;


my $refFams  = [];
my $cmpFams  = [];
my $refFile  = "";
my $cmpFile  = "";
my $rootDir  = "";


my $domain_cov = 0.7;
my $evalue     = 1e-5;
my $prop_prots_w_domain = 0.05;


my $tcdbSeqsFile = "";
my $pfamFile = "";
my $blastdb  = "";


readCommandLine();

#print Data::Dumper->Dump([$refFams,$cmpFams,$refFile,$cmpFile,$rootDir,$domain_cov,$evalue,$prop_prots_w_domain],
#			 [qw(*refFams *cmpFams *refFile *cmpFile *rootDir *domain_cov *evalue *prop_prots_w_domain)]);
#exit;


all_vs_all();



#==========================================================================
#Use each family as a reference and do pairwise comparisons against all
#other families against it


sub all_vs_all {

  my $fams = shift;

 REF:foreach my $ref (@$refFams) {
  CMP:foreach my $cmp (@$cmpFams) {

      next CMP if ($ref eq $cmp);
      my $pair1  = "$rootDir/${ref}_vs_$cmp";
      my $pair2  = "$rootDir/${cmp}_vs_$ref";
      my $params = qq( -dc $domain_cov -e $evalue -m $prop_prots_w_domain -sf);

      system "getDomainTopology.pl -f $ref,$cmp -o $pair1 $params" unless (-d $pair1);
      system "getDomainTopology.pl -f $cmp,$ref -o $pair2 $params" unless (-d $pair2);
    }
  }
}



#==========================================================================
#Read command line


sub readCommandLine {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "rf|ref-fams=s"      => \&readRefFams,        #TCIDs that will be used as reference
      "cf|cmp-fams=s"      => \&readCmpFams,        #TCIDs that will that will compared
      "rfile|ref-file=s"   => \&readRefFile,        #File with reference TCIDs
      "cfile|cmp-file=s"   => \&readCmpFile,        #File with comparison TCIDs
      "o|outdir=s"         => \&read_outdir,        #Ouput root directory
      "e|evalue=f"         => \$evalue,
      "dc|domain-cov=f"    => \$domain_cov,
      "m|prots-w-domain=f" => \$prop_prots_w_domain,
      "s|tcdb-seqs=s"      => \&read_tcdb_seqs,     #File with all sequences in TCDB
      "pfam=s"             => \&read_pfam,          #hmmscan output file for whole TCDB
      "b|blastdb=s"        => \&read_blastdb,       #Full path of blastdb to extract sequences
      "h|help"             => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);

  #----------------------------------------------------------------------
  #Validate command line arguments


  if (@$refFams && $refFile) {
    die "Error: options -rf and -rfile are mutually exclusive. please, provide only one\n";
  }

  if (@$cmpFams && $cmpFile) {
    die "Error: options -cf and -cfile are mutually exclusive. please, provide only one\n";
  }


  #At this point is guaranteed that only of ref and/or cmp families were given
  $refFams = TCDB::Assorted::readFamilyFile($refFile) if ($refFile);
  $cmpFams = TCDB::Assorted::readFamilyFile($cmpFile) if ($cmpFile);


  die "Error: at least -rf, -cf, -rFile, or -cFile must be given !\n" unless (@$refFams || @$cmpFams);

  #when only one type of families are given
  if (@$cmpFams && !(@$refFams)) {
    $refFams = $cmpFams;
  }
  elsif (@$refFams && !(@$cmpFams)) {
    $cmpFams = $refFams;
  }
}








#==========================================================================
#Read the -rf option

sub readRefFams {

  my ($opt, $value) = @_;

  $refFams = TCDB::Assorted::readStringFams($value);
}



#==========================================================================
#Read the -cf option

sub readCmpFams {

  my ($opt, $value) = @_;

  $cmpFams = TCDB::Assorted::readStringFams($value);
}



#==========================================================================
#Read the -rfile option

sub readRefFile {

  my ($opt, $value) = @_;

  $refFile = $value;
}



#==========================================================================
#Read the -cfile option

sub readCmpFile {

  my ($opt, $value) = @_;

  $cmpFile = $value;
}




#==========================================================================
#Read the -o option

sub read_outdir {

  my ($opt, $value) = @_;

  system "mkdir -p $value" unless (-d $value);

  $rootDir = $value;
}




#==========================================================================
#Read the -s option

sub read_tcdb_seqs {

  my ($opt, $value) = @_;

  die "File with TCDB sequences must exist and not be empty: $value" unless (-f $value && !(-z $value));

  $tcdbSeqsFile = $value;
}




#==========================================================================
#Read the -pfam option

sub read_pfam {

  my ($opt, $value) = @_;

  die "File with Pfam domains must exist and not be empty: $value" unless (-f $value && !(-z $value));

  $pfamFile = $value;
}

#==========================================================================
#Read the -b option

sub read_blastdb {

  my ($opt, $value) = @_;

  my $tmpFile = "${value}.phd";
  die "Input is not a path to a blast database: $value" unless (-f $tmpFile && !(-z $tmpFile));

  $blastdb = $value;
}



#==========================================================================
#Print help

sub print_help {

  my $help = <<'HELP';

Pairwise Pfam-domain comparisons of TCDB families. Comparisons can be
performed between all members of a known superfamily (i.e. positive control)
or between members of a positive control and a selected negative control set
of families.


-rf, --ref-fams {string} (optional)
  List of comma-separated tcdb families that will be used as reference.
  This options is imcompatible with -rfile.

-cf, --cmp-fams {string} (optional)
  List of comma-separated tcdb familes that will be used to compare against
  the set of families given through option -rf. This option is incompatible
  with -cfile.

-rfile, --ref-file {file} (optional)
  File with families that will be used a reference (i.e. --ref-fams). There
  must be only one family per line.

-cfile, --cmp-file {file} (optional)
  File with families that will be compared against the reference set
  (i.e. --cmp-familes). There must be only one family per line.

-o, --outdir {directory} (optional)
  Directory where the output of the analysis will be stored.
  Default: domainComparisons

-s, --tcdb-seqs {file} (optional)
  File with all sequences in TCDB. This files allows the analysis to be
  performed on a frozen version of TCDB. If this option is not given,
  sequences will be directly downloaded from TCDB.

--pfam {file} (optional)
  Path to the output files of hhmscan comparing the whole TCDB to
  Pfam, or atleast all the members in the families in the analysis.


Either -rf or rfile must be given. If only -cf or -cfile are given,
they will be treated as -rf.


HELP


  print $help;
  exit;
}
