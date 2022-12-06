#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

#To read command line options
use Getopt::Long;

#to check dependencies
use TCDB::CheckDependencies;
use TCDB::Assorted;
use Bio::SeqIO;



#==========================================================================
# This script  extracts famXpander results from the global repository:
#                  /ResearchData/famXpander
# This repository contains famXpander results for each system component
# in TCDB. The directories within this repository follow the format:
#            [TCID]-[UniProt/RefSeq accession]
#
# The script can either: 1) extract all the famXpander results under a
# a given TCID, or 2) recive a list of of component ids in the same format
# as the famXpander repository.
#
# Results will be placed in an user-defined output directory with the
# following famXpander-compatible files:
#   1. psiblast.tbl
#   2. results.faa
#   3. results.cd-hit.clstr
#--------------------------------------------------------------------------
#
#  Written by: Arturo Medrano 2021
#
#==========================================================================



#==========================================================================
#Check dependencies


my @dependencies = ("mkdir", "cat", "rm", "mv", "cd-hit");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;




#==========================================================================
#Read command line options


#Identifier for the class, subclass, family, subfamily or system of interest
my $tcid = "";

#File with individuale [TCID]-[Accsession] of the proteins of interest.
my $tcFile = "";

#This is the directory were the precomputed famXpander results for
#all TCDB are stored
my $global_fxpand_dir = "$ENV{RESEARCH_DATA}/famXpander";

#Output directory for global famXpander results
my $outdir = "./famXpander";

#For famXpander, identity threshold to determine cd-hit redundancy
my $cdhit_cutoff = 0.9;




read_command_line_arguments();


#print Data::Dumper->Dump([$tcid, $tcFile, $outdir, $cdhit_cutoff],
#			 [qw(*tcid *tcFile *outdir *cdhit_cutoff)]);
#exit;




#==========================================================================
#If the -x argument was given, download the sequences of each
#input family from TCDB.


if ($tcid) {
  get_famXpander_for_tcid ($outdir, $tcid);
}
else {
  get_famXpander_for_compList ($outdir, $tcFile);
}

###########################################################################
##                    Subroutine definitions                             ##
###########################################################################


#==========================================================================
#Given a list of [tcid]-[accesion] extract all their famXpander results
#and treat them as a single family.

sub get_famXpander_for_compList {

   my ($outDir, $inFile) = @_;

   system "mkdir -p $outDir" unless (-d $outDir);

   my $finalResFile   = "$outDir/results.faa";
   my $finalPsiFile   = "$outDir/psiblast.tbl";
   my $finalClstrFile = "$outDir/results.faa.cdhit.clstr";

   open (my $fh, "<", $inFile) || die $!;
   chomp (my @prots = <$fh>);
   close $fh;
   my %tcComponents = map {$_ => 1} @prots;

#   print Data::Dumper->Dump([\%tcComponents], [qw(*tcComponents )]);
#   exit;

   print "\n\nGetting fxpander results for input proteins\n";
  getHomologsForProteins($outDir, $finalResFile, $finalPsiFile, $finalClstrFile, \%tcComponents);

}


#==========================================================================
#Extract subfamilies and systems with famXpander results from the global
#repository, remove redundant sequences, and generate files with the same
#names than usual famXpander output files.


sub get_famXpander_for_tcid {

  my ($outDir, $tc) = @_;


  #Output files
  my $rootDir = "$outDir/$tc";
  system "mkdir -p $rootDir" unless (-d $rootDir);

  my $finalResFile   = "$rootDir/results.faa";
  my $finalPsiFile   = "$rootDir/psiblast.tbl";
  my $finalClstrFile = "$rootDir/results.faa.cdhit.clstr";


  #First get the number of elements in the input tcid
  my @elements1 = split (/\./, $tc);


  #----------------------------------------------------------------------
  #get the systems associated with $fam_id1 and $fam_id2

  my %tcComponents = ();
  opendir (my $dirh, $global_fxpand_dir) || die "Could not open dir: $global_fxpand_dir --> $_";

  unless (-f $finalResFile && -f $finalPsiFile && $finalClstrFile) {

    #Get the diectories that match the input tcdb family
    if (scalar @elements1 < 5) {

      #The second regex is in case the directory is equal to the input tcid
      %tcComponents = map { $_ => 1 } grep { /^$tc\./ || /^$tc$/} readdir $dirh;
    }
    else {
      %tcComponents = map {$_ => 1} grep { /^$tc$/ } readdir $dirh;
    }

    unless (%tcComponents) {
      print "No famXpander results found for family $tc\n";
      exit;
    }
  }
  closedir $dirh;


#  print Data::Dumper->Dump([\%tcComponents], [qw(*tcComponents)]);
#  print "Total components: ", scalar keys %tcComponents, "\n";
#  exit;

  print "\n\nGetting fxpander results for: $tc\n";
  getHomologsForProteins($rootDir, $finalResFile, $finalPsiFile, $finalClstrFile, \%tcComponents);
}



#==========================================================================
#Collapse the results of famXpander into one single file and put it in $outDir.

sub getHomologsForProteins {

  my ($oDir, $rFile, $psiFile, $cFile, $prots) = @_;


  #Clean the file contents before starting the concatenation
  system "cat /dev/null >  $rFile";
  system "cat /dev/null >  $psiFile";
  system "cat /dev/null >  $cFile";


  #Concatenated data for input proteins
  my $concResFile = "$oDir/global_results.faa";
  system "cat /dev/null > $concResFile";


  #Start collecting homologous sequences and blast results
  foreach my $component (keys %{ $prots }) {

    #path to global famXpander results.faa and psiblast.tbl file
    my $resFile   = "$global_fxpand_dir/$component/results.faa";
    my $blastFile = "$global_fxpand_dir/$component/psiblast.tbl";

    #Concatenate files here and add a new line for safety
    system "cat $resFile   >> $concResFile";
    system "cat $blastFile >> $psiFile\n";
  }

#  print "Populated files:\n  $concResFile\n  $finalPsiFile\n";
#  exit;

  #run cd-hit to remove redundant sequences
  die "Could not find file: $concResFile" unless (-f $concResFile && !(-z $concResFile));

  system "cd-hit -i $concResFile -o $rFile -c $cdhit_cutoff -M 0 -l 30 -d 0 -g 1";
  system "rm $concResFile";
  system "mv ${rFile}.clstr $cFile";

  die "Error: could not extract sequences for proteins!" unless (-f $rFile && !(-z $rFile));
}


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
      "i|tcid=s"                  => \&read_tcid,
      "f|infile=s"                => \&read_tcFile,
      "o|outdirr=s"               => \$outdir,
      "gfxd|global-fxpand-dir=s"  => \$global_fxpand_dir,
      "r|cdhit-cutoff=f"          => \$cdhit_cutoff,
      "h|help"                    => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);



  #----------------------------------------------------------------------
  #Validate command line arguments


  #Input families
  die "Error: argument -i or -f are  mandatory\n"     unless ($tcid || $tcFile);
  die "Error: arguments -i and -f are incompatible\n" if ($tcid && $tcFile);

  die "Error: Global famXpander dir not found: $global_fxpand_dir" unless (-d $global_fxpand_dir);

  system "mkdir -p $outdir" unless (-d $outdir)

}


#==========================================================================
#Read the -i option

sub read_tcid {

  my ($opt, $value) = @_;

  TCDB::Assorted::validate_tcdb_id([$value]);

  $tcid = $value;
}


#==========================================================================
#Read -f option

sub read_tcFile {
  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
    die "Error: File with component IDs not found or empty\n";
  }

  $tcFile = $value;
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

This script extracts custom famXpander results from the global
Repository.

Note that either option -i or -f must be given.

 -i, --tcid {string}
   TCDB ID that will be used to extract famXpander results. 
   all systems/components under this TCID will be processed.
   Not compatible with option -f

 -f, --infile {File}
   File with component IDs in the format [tcid]-[accession]
   that will be used to extract homologs. Not compatible with
   option -i

 -o, --outdir {Directory}
   Output directory where results will be safe
   (Default: ./famXpander)

 -gfxd, --global-fxpand-dir {Directory}
   Path to the directory with the results of running famXpander 
   for all proteins in TCDB. 
   (Default: $RESEARCH_DATA/famXpander).

 -r, --cdhit-cutoff {float}
   Identity redundancy threshold for cd-hit.
   (default 0.9)

 -h, --help
   Display the help of this program. Help will also be
   printed if no arguments are passed to the program.

HELP


  print "$help\n";
  exit;
}



