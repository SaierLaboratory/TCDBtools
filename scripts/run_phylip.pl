#!/usr/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;

use R2::CheckDependencies;


#
#This script performs phylogenetic analyses using several algorithms
#as implemented n phylip
#


#==========================================================================
#Check dependencies

my @dependencies = ("seqboot", "protdist", "neighbor", "fitch", "proml");
my $CheckDep_obj = new R2::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;



#==========================================================================
#Global variables

my @gapsPerc   = ();
my $program    = 'neighbor';
my $rootDir    = ".";
my $filePrefix = "";
my $bootstrap  = 100;
my $dirPrefix  = "gt";
my $alnDir     = "../..";

#==========================================================================
#read comand line arguments

read_command_line();

#print Data::Dumper->Dump([\@gapsPerc, $program, $rootDir, $filePrefix, $dirPrefix, $bootstrap],
#                       [qw(*gapsPerc  *program  *rootDir  *filePrefix  *dirPrefix *bootstrap)]);
#exit;



#==========================================================================
#Calculate the phylogeny for each alignment


foreach my $gap (@gapsPerc) {

  if ($program eq "neighbor") {
    run_neighbor($gap);
  }
  elsif ($program eq "fitch") {
    run_fitch($gap);
  }
  elsif ($program eq "proml") {
    run_proml($gap);
  }
}



###########################################################################
##                                                                       ##
##                        Subroutine definitions                         ##
##                                                                       ##
###########################################################################


#==========================================================================
#Run the fitch program

sub run_fitch {

  my $gapValue = shift;

  #get the working directory for this analysis
  my $workdir = "$rootDir/${dirPrefix}$gapValue";
  die "Work directory does not exist: $workdir" unless (-d $workdir);


  print "Working gap $gapValue in: $workdir\n";


  #Move to the current directory
  chdir($workdir) || die "Could not chdir to work directory: $workdir --> $!";


  #The original alignment file
  my $origAlnFile = ($alnDir)? "$alnDir/${filePrefix}${gapValue}.phy" : "${filePrefix}${gapValue}.phy";
  die "Alignment file was not found: $origAlnFile" unless (-f $origAlnFile);


  #-----------------------------------------------------------------
  #Run seqboot

  run_seqboot($gapValue);



  #-----------------------------------------------------------------
  #Run protdist

  my $paramsFile = "$workdir/params_protdist";

  unless (-f "$workdir/bootstrap_distances") {

    #Create file with the parameters for: protdist
    create_params_file($paramsFile, "bootstrap_aln", "M", "D", $bootstrap, "Y");


    #Run protdist
    system "protdist < params_protdist";


    #Save bootstrapped distances to the proper file name
    system "mv outfile bootstrap_distances";
  }



  #-----------------------------------------------------------------
  #Run fitch algorithm


  $paramsFile = "$workdir/params_fitch";

  unless (-f "$workdir/bootstrap_fitch_trees") {

    #Create file with the parameters for: fitch
    create_params_file($paramsFile, "bootstrap_distances", "G", "M", $bootstrap, 555, "Y");


    #Run Fitch-Margoliash
    system "fitch < params_fitch";


    #Saved bootstrapped trees to the proper file name
    system "mv outtree bootstrap_fitch_trees";
    system "rm outfile";
  }



  #-----------------------------------------------------------------
  #Generate the consensus tree

  $paramsFile = "$workdir/params_consense";

  unless (-f "$workdir/consense_tree.ph") {

    #Create file with the parameters for: neighbor
    create_params_file($paramsFile, "bootstrap_fitch_trees", "Y");


    #Run consense
    system "consense < params_consense > /dev/null";


    #Saved bootstrapped trees to the proper file name
    system "mv outtree consense_tree.ph";
    system "rm outfile";
  }

}






#==========================================================================
#Run the maximum-likelihood algorithm

sub run_proml {

  my $gapValue = shift;

  #get the working directory for this analysis
  my $workdir = "$rootDir/${dirPrefix}$gapValue";
  die "Work directory does not exist: $workdir" unless (-d $workdir);


  print "Working gap $gapValue in: $workdir\n";


  #Move to the current directory
  chdir($workdir) || die "Could not chdir to work directory: $workdir --> $!";


  #-----------------------------------------------------------------
  #Run seqboot

  if ($bootstrap) {
    run_seqboot($gapValue);
  }



  #-----------------------------------------------------------------
  #run proml

  my $paramsFile = "$workdir/params_proml";

  if ($bootstrap > 0) {

    unless (-f "$workdir/bootstrap_proml_trees") {

      #Create file with the parameters for: proml
      create_params_file($paramsFile, "bootstrap_aln", "S", "G", "M", "D", $bootstrap, 555, 10, "J",  "Y");


      #Run proml
      system "proml < params_proml > /dev/null";


      #Save bootstrapped distances to the proper file name
      system "mv outtree bootstrap_proml_trees; rm outfile";
    }
  }
  else {

    #Create file with the parameters for: proml
    create_params_file($paramsFile, "${filePrefix}${gapValue}.phy", "S", "G", "J", 555, 10, "Y");


    #Run proml
    system "proml < params_proml";


    #Save tree
    system "mv outtree proml_tree.ph; rm outfile";
  }



  #-----------------------------------------------------------------
  #Generate the consensus tree

  $paramsFile = "$workdir/params_consense";

  if ($bootstrap > 0) {
    unless (-f "$workdir/proml_consense_tree.ph") {

      #Create file with the parameters for: neighbor
      create_params_file($paramsFile, "bootstrap_proml_trees", "Y");


      #Run consense
      system "consense < params_consense > /dev/null";


      #Saved bootstrapped trees to the proper file name
      system "mv outtree proml_consense_tree.ph";
      system "rm outfile";
    }
  }
}





#==========================================================================
#Run the neighbor joining program

sub run_neighbor {

  my $gapValue = shift;

  #get the working directory for this analysis
  my $workdir = "$rootDir/${dirPrefix}$gapValue";
  die "Work directory does not exist: $workdir" unless (-d $workdir);


  print "Working gap $gapValue in: $workdir\n";


  #Move to the current directory
  chdir($workdir) || die "Could not chdir to work directory: $workdir --> $!";


  #The original alignment file
  my $origAlnFile = ($alnDir)? "$alnDir/${filePrefix}${gapValue}.phy" : "${filePrefix}${gapValue}.phy";
  die "Alignment file was not found: $origAlnFile" unless (-f $origAlnFile);


  #-----------------------------------------------------------------
  #Run seqboot

  run_seqboot($gapValue);



  #-----------------------------------------------------------------
  #Run protdist

  my $paramsFile = "$workdir/params_protdist";

  unless (-f "$workdir/bootstrap_distances") {

    #Create file with the parameters for: protdist
    create_params_file($paramsFile, "bootstrap_aln", "M", "D", 100, "Y");


    #Run protdist
    system "protdist < params_protdist";


    #Save bootstrapped distances to the proper file name
    system "mv outfile bootstrap_distances";
  }



  #-----------------------------------------------------------------
  #Run neighbor-joining algorithm


  $paramsFile = "$workdir/params_neighbor";

  unless (-f "$workdir/bootstrap_neighbor") {

    #Create file with the parameters for: neighbor
    create_params_file($paramsFile, "bootstrap_distances", "M", 100, 555, "Y");


    #Run neighbor
    system "neighbor < params_neighbor";


    #Saved bootstrapped trees to the proper file name
    system "mv outtree bootstrap_trees";
    system "rm outfile";
  }



  #-----------------------------------------------------------------
  #Generate the consensus tree

  $paramsFile = "$workdir/params_consense";

  unless (-f "$workdir/consense_tree.ph") {

    #Create file with the parameters for: neighbor
    create_params_file($paramsFile, "bootstrap_trees", "Y");


    #Run consense
    system "consense < params_consense";


    #Saved bootstrapped trees to the proper file name
    system "mv outtree consense_tree.ph";
    system "rm outfile";
  }

}



#==========================================================================
#Run bootstrap on a given alignment

sub run_seqboot {

  my $gapValue = shift;

  #get the working directory for this analysis
  my $workdir = "$rootDir/${dirPrefix}$gapValue";
  die "Work directory does not exist: $workdir" unless (-d $workdir);


  #Move to the current directory
  chdir($workdir) || die "Could not chdir to work directory: $workdir --> $!";


  #The original alignment file
  my $origAlnFile = ($alnDir)? "$alnDir/${filePrefix}${gapValue}.phy" : "${filePrefix}${gapValue}.phy";
  die "Alignment file was not found: $origAlnFile" unless (-f $origAlnFile);


  #-----------------------------------------------------------------
  #Run seqboot

  my $paramsFile = "$workdir/params_seqboot";

  unless (-f "$workdir/bootstrap_aln") {

    #Create file with the paramenters for: seqboot.
    create_params_file ($paramsFile, $origAlnFile, "Y", 555);


    #Run seqboot
    system "seqboot < params_seqboot > /dev/null";


    #Save bootstrapped alignments to a different file
    system "mv outfile bootstrap_aln";
  }
}



#==========================================================================
#Create parameters file for any of the phylip programs

sub create_params_file {

  my ($file, @params) = @_;

  open (my $booth, ">", $file) || die $!;
  foreach my $par (@params) {
    print $booth "$par\n"
  }
  close $booth;
}




#===========================================================================
#Read command line and print help


sub read_command_line {

    print_help() unless (@ARGV);

    my @usrFiles;

    my $status = GetOptions(
	"p|program=s"      => \&read_phylip_program,
	"g|perc-gaps=s"    => \&read_percentage_of_gaps,
	"r|root-dir=s"     => \&read_root_dir,
	"fp|file-prefix=s" => \$filePrefix,
	"dp|dir-prefix=s"  => \$dirPrefix,
        "b|bootstrap=i"    => \$bootstrap,
        "h|help"           => sub { print_help(); },
	"<>"               => sub { die "Error: Unknown argument: $_[0]\n"; });
    exit unless ($status);

    die "option -fp is mandatory" unless ($filePrefix);
}


sub read_root_dir {

  my ($opt, $value) = @_;

  $rootDir = ($value =~ /^\//)? $value : "$ENV{PWD}/$value";

  unless (-d $rootDir) {
    die "Option -r requires an existing directory";
  }
}



sub read_phylip_program {

  my ($opt, $value) = @_;

  unless ($value =~ /(neighbor|fitch|proml)/) {
    die "Unknown phylip program: $value";
  }

  $program = $value;
}


sub read_percentage_of_gaps {

  my ($opt, $value) = @_;

  @gapsPerc = split(/,/, $value);

  foreach my $n (@gapsPerc) {
    unless ($n =~ /^\d+$/) {
      die "Option -$opt requires a comma-separated list of integer values: $value";
    }
  }
}





sub print_help {

    my $help = <<'HELP';

This script runs phylip for any number of multiple alignments. The
programs currently supported are: neighbor, fitch and proml.  

It is assumed that all the phylogenies using a specific method will
have a root directory, inside that directory there must be other
subdirectories with a common prefix (e.g., for multiple alignments
with different degrees of trimming using the trimal program). The
sequencies with the reference multiple alignment can be in any
directory but also must have a common prefix.


Options:


-h, --help
   Display this help. Also displayed if script is run without arguments.


Examples:
run_phylip -p neighbor -g 60,80,90 -r ./neighbor -fp trimmed_clustalo_1.A.15_gt -dp gt -b 100
run_phylip -p fitch -g 50,55,60,65,70,75,80,85,90,95 -r ./fitch -fp trimmed_clustalo_1.A.15_gt -dp gt -b 100

HELP

    print $help;
    exit;
}

