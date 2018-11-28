#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;
#use List::Util qw(sum);

use TCDB::Assorted;
use TCDB::Domain::PfamParser;
use TCDB::Domain::Characterize;

use Getopt::Long;

#==========================================================================
#Global variables

#Query family or families
my @fams = ();

#This is an option for TCDB::Assorted::getSystemAccessions()
my $treatAsSuperfamily   = 0;

#Options for TCDB::Domain::PfamParser
my $domain_cov = 0.7;
my $prot_cov   = 0.1;
my $evalue     = 1e-5;
my $prop_prots_w_domain = 0.05;

#Options for TCDB::Domain::Characterize
my $rootDir              = ".";
my $tcdbSeqsFile         = "/Users/amedrano/Projects/tcdb/TOG/data/tcdb.faa"; #"/ResearchData/pfam/download/tcdb.faa";
my $pfamFile             = "/ResearchData/pfam/tcdb.pfam-a.hmmscan.bz2"; #"/Users/amedrano/Projects/tcdb/TOG/data/tcdb.pfam-a.hmmscan.bz2"; 
my $blastdb              = "/Users/amedrano/db/blastdb/tcdb";
my $prog                 = "ssearch36";
my $analysisLevel        = 'system';



#Read command line topology
read_command_line_arguments();


die "TCDB sequences file not found or empty --> $tcdbSeqsFile\n" unless (-f $tcdbSeqsFile && !(-z $tcdbSeqsFile));
die "TCDB hmmscan output file not found --> $pfamFile\n" unless (-f $pfamFile && !(-z $pfamFile));


#print Data::Dumper->Dump([\@fams, $treatAsSuperfamily,  $rootDir, $tcdbSeqsFile, $pfamFile, $blastdb, $prog, $domain_cov,
#			  $prot_cov, $evalue, $prop_prots_w_domain],
#			 [qw(*fams *treatAsSuperfamily *rootDir *tcdbSeqsFile *pfamFile *blastdb *prog *domain_cov
#			     *prot_cov *evalue *prop_prots_w_domain)]);
#exit;


#==========================================================================
#Split tcdb systems into single-component multi-component


if ($treatAsSuperfamily) {

    my $tcids = getSystemAccessions($tcdbSeqsFile, 'both', $analysisLevel, \@fams, $treatAsSuperfamily);

#    print Data::Dumper->Dump([$tcids ], [qw(*tcids )]);
#    exit;



    #==========================================================================
    #Setup the thresholds for parsing the PFAM output


    my $obj = new TCDB::Domain::PfamParser();
    $obj->pfamFile($pfamFile);
    $obj->analysisLevel($analysisLevel);
    $obj->domCovCutoff($domain_cov);
    $obj->tcCovCutoff($prot_cov);
    $obj->evalueCutoff($evalue);
    $obj->minProtsDom($prop_prots_w_domain);
    $obj->treatAsSuperfamily($treatAsSuperfamily);


    my %domFreq   = ();
    my %domCoords = ();
    $obj->getDomainStatsForUserFamilies(\@fams, $tcids, \%domFreq, \%domCoords);

#    print Data::Dumper->Dump([ \%domFreq, \%domCoords ], [qw( *domFreq *domCoords )]);
#    exit;





    #==========================================================================
    #Attempt to rescue the domains that were not recognized by PFAM in some
    #Family members


    my $rescueObj = new TCDB::Domain::Characterize();
    $rescueObj->rootDir($rootDir);
    $rescueObj->tcdbFaa($tcdbSeqsFile);
    $rescueObj->domCoords(\%domCoords);
    $rescueObj->domFreq(\%domFreq);
    $rescueObj->tcids($tcids);
    $rescueObj->searchWith($prog);
    $rescueObj->blastdb($blastdb);
    $rescueObj->evalue($evalue);
    $rescueObj->treatAsSuperfamily($treatAsSuperfamily);


    $rescueObj->rescueDomains(\@fams);

}
else {

  foreach my $fam (@fams) {

    my $tcids = getSystemAccessions($tcdbSeqsFile, 'both', $analysisLevel, [$fam], $treatAsSuperfamily);

    #print Data::Dumper->Dump([$tcids ], [qw(*tcids )]);
    #exit;



    #==========================================================================
    #Setup the thresholds for parsing the PFAM output


    my $obj = new TCDB::Domain::PfamParser();
    $obj->pfamFile($pfamFile);
    $obj->analysisLevel($analysisLevel);
    $obj->domCovCutoff($domain_cov);
    $obj->tcCovCutoff($prot_cov);
    $obj->evalueCutoff($evalue);
    $obj->minProtsDom($prop_prots_w_domain);


    my %domFreq   = ();
    my %domCoords = ();
    $obj->getDomainStatsForUserFamilies([], $tcids, \%domFreq, \%domCoords);

    #print Data::Dumper->Dump([ \%domFreq, \%domCoords ], [qw( *domFreq *domCoords )]);
    #exit;





    #==========================================================================
    #Attempt to rescue the domains that were not recognized by PFAM in some
    #Family members


    my $rescueObj = new TCDB::Domain::Characterize();
    $rescueObj->rootDir($rootDir);
    $rescueObj->tcdbFaa($tcdbSeqsFile);
    $rescueObj->domCoords(\%domCoords);
    $rescueObj->domFreq(\%domFreq);
    $rescueObj->tcids($tcids);
    $rescueObj->searchWith($prog);
    $rescueObj->blastdb($blastdb);
    $rescueObj->evalue($evalue);
    $rescueObj->domCovCutoff($domain_cov);
    $rescueObj->treatAsSuperfamily($treatAsSuperfamily);

    $rescueObj->rescueDomains();
  }
}




###########################################################################
##                         Functions                                     ##
###########################################################################



sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "f|family=s"           => \&read_fams,       #TCIDs of families to analyze (comma separated)

      #Options for TCDB::Domain::PfamParser
      "dc|domain-cov=f"      => \$domain_cov,
      "pc|protein-cov=f"     => \$prot_cov,
      "e|evalue=f"           => \$evalue,
      "m|prots-w-domain=f"   => \$prop_prots_w_domain,

      #Options for TCDB::Domain::Characterize
      "o|outdir=s"           => \&read_root_dir,        #Ouput root directory
      "s|tcdb-seqs=s"        => \&read_tcdb_seqs,       #File with all sequences in TCDB
      "sf|superfamily!"      => \$treatAsSuperfamily,	#File with the sequences of the reference family
      "pfam=s"               => \&read_pfam,            #hmmscan output file for whole TCDB
      "b|blastdb=s"          => \&read_blastdb,         #Full path of blastdb to extract sequences
      "p|rescue-prog=s"      => \&read_prog,            #Read the program that will be used to rescue domains (blastp|ssearch36)
      "h|help"               => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);

  #----------------------------------------------------------------------
  #Validate command line arguments


  die "Error: at least one TCID is mandatory!\n" unless (@fams);
}


#==========================================================================
#Read the -f option

sub read_fams {

  my ($opt, $value) = @_;

  @fams = split (/,/, $value);

  TCDB::Assorted::validate_tcdb_id(\@fams);
}


#==========================================================================
#Read the -d option

sub read_root_dir {

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

  die "File with TCDB sequences must exist and not be empty: $value" unless (-f $value && !(-z $value));

  $pfamFile = $value;
}


#==========================================================================
#Read the -b option

sub read_blastdb {

  my ($opt, $value) = @_;

  my $tmpFile = "${value}.phd";
  die "Input is not a path to blast database: $value" unless (-f $tmpFile && !(-z $tmpFile));

  $blastdb = $value;
}



#==========================================================================
#Read the - option

sub read_prog {

  my ($opt, $value) = @_;

  die "Unrecognized program: $value" unless ($value =~ /(ssearch36|blastp)/);

  $prog = $value;
}




#==========================================================================
#Print help

sub print_help {

  my $help = <<'HELP';

Identify the main protein domains in a family.

 -f, --family {string} (Mandatory)
  TCID of the family for which Pfam domain analysis will be carried out.
  If multiple TCIDs are given, they whould be comma-separated. The analysis
  will be performed individually for each family, unless the flag -sf is
  given, in which case all TCIDs will be treated as a superfamily.

 -dc, --domain-cov {float} (Default: 0.7)
  Minimum coverage of the Pfam domain to consider it a match. If coverage
  is less than the specified threshold, the coverage must apply to the
  the query protein to consider the domain hit significant.

 -pc, --protein-cov {float} (Default: 0.1)
  Minimum coverage of the query protein sequence per domain.

 -e, --evalue {float} (Default: 1e-5)
  Maximum evalue threshold to consider a Pfam domain hit.

 -m, --prots-w-domain {float) (Default: 0.05)
  Minimum proportion of the proteins in the input family that should
  contain a domain, in order to consider the domain as part of the
  family for the purpose of this analyis.

 -o, --outdir {path} (Default: .)
  Directory where results and intermediary files will be saved.

 -s, --tcdb-seqs {file} (Mandatory)
  FASTA file with all sequences in TCDB (as generated by program
  extractFamily.pl). This file will be used to extract sequences
  and TCIDs for all members of the input family. This allows to
  freeze TCDB contents at a specific date, for the purpose defining
  a project.

 -sf, --superfamily {flag} (default: negated as --no-superfamily)
  This indicates that all families passed to option -f will be
  treated as a superfamily (i.e. as a single family) for the purpose of
  the Pfam domain  analysis. This is useful when a superfamily is
  composed of two or more different family TCIDs.

 -pfam {file} (Mandatory)
  The output of running hmmscan against all TCDB, or at least the input
  family(-ies).

 -b, --blastdb {path} (Default: assumed available in $BLASTDB)
  Full path to the blast database containing all the sequences in
  TCDB, this is to easily extract segments of proteins that match
  Pfam domains.

 -p, --rescue-prog (Default: ssearch36)
  Program that will be used to project or rescue domains in
  proteins that did have direct hits with Pfam domains.

 -h, --help
  Display this help. If present, this option takes precedence over any
  other option.

HELP


  print $help;
  exit;
}
