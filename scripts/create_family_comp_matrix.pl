#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;

$Data::Dumper::Deepcopy = 1;

#To read command line options
use Getopt::Long;



#to check dependencies
use TCDB::CheckDependencies;
use TCDB::Assorted;


###########################################################################
#
#  Create a matrix of gsat comparisons between a group of families.
#  Comparison can be created all_vs_all or a reference family(es)
#  vs all.
#
###########################################################################


#==========================================================================
#Check dependencies

my @dependencies = ("sort", "grep", "list_top_gsat_hits");
my $CheckDep_obj = new R2::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;




#==========================================================================
#Read command line arguments

my $all_vs_all    = 0;
my $ref_vs_all    = 0;
my @ref_families  = ();
my $gsat_cutoff   = 15;
my $fam_div       = "_vs_";
my $rootDir       = undef;
my $exclude_fam   = '1.A.17.8';

read_command_line();

#print Data::Dumper->Dump([$all_vs_all, $ref_vs_all, \@ref_families, $gsat_cutoff, $rootDir],
#                         [qw(*all_vs_all *ref_vs_all *ref_families *gsat_cutoff *rootDir)]);
#exit;



#==========================================================================
#The reference families if applicable

my %refFamilies = ();

if (@ref_families) {
  %refFamilies =  map {$_ => 1} @ref_families;
}



#==========================================================================
#Identify the set of families that will be used to in the analysis

my %tmp_families = ();
opendir (my $rdirh, $rootDir) || die $!;

while (readdir $rdirh) {

  #Ignore dot directories
  next if (/\.+$/);

  #only work with the directories
  my $dpath = "$rootDir/$_";
  next unless (-d $dpath);


  my ($fam1, $fam2) = split (/$fam_div/, $_);
  $tmp_families{$fam1} = 1 unless ($fam1 eq $exclude_fam);
  $tmp_families{$fam2} = 1 unless ($fam2 eq $exclude_fam);

}
closedir $rdirh;


#Now validate all TCDB IDs in the comparison

my @all_families = sort {$a cmp $b} keys  %tmp_families;
TCDB::Assorted::validate_tcdb_id (\@all_families);

#print Data::Dumper->Dump([\@all_families ], [qw(*all_families)]);
#exit;


#==========================================================================
#Compare all versus all

#Where the scores will be stored
my %matrix = ();


if ($all_vs_all) {


  foreach my $idx1 (0 .. $#all_families - 1) {

    my $f1 = $all_families[$idx1];
    next if ($f1 eq $exclude_fam);

    foreach my $idx2 ( $idx1 + 1 .. $#all_families) {

      my $f2 =  $all_families[$idx2];
      next if ($f2 eq $exclude_fam);

      #Get the right directory to search for protocol2 results
      my $sdir = "";
      if (-d "$rootDir/${f1}${fam_div}$f2") {
	$sdir = "$rootDir/${f1}${fam_div}$f2";
      }
      elsif (-d "$rootDir/${f2}${fam_div}$f1") {
	$sdir = "$rootDir/${f2}${fam_div}$f1";
      }
      else {
	die "Coud not find directory comparing $f1 vs $f2.";
      }


      #Now get the top scores between each pair of families
      my $cmd = qq(list_top_gsat_hits -p2d $sdir -b);
      chomp (my $out = `$cmd`);

      my $score = undef;
      if ($out) {
	my @data = split (/\s+/, $out);

	$score = sprintf ("%.1f", $data[1]);
      }
      else {
	$score = "Failed";
      }

      $matrix{$f1}{$f2} = $score;
    }
  }

#  print Data::Dumper->Dump([\%matrix], [qw(*matrix )]);
#  exit;

  print_matrix(\%matrix, \@all_families);
}

else {


  foreach my $idx1 (0 .. $#ref_families) {

    my $f1 = $ref_families[$idx1];


    foreach my $idx2 ( 0 .. $#all_families) {

      my $f2 =  $all_families[$idx2];
      next if ($f2 eq $exclude_fam);
      next if ($f2 eq $f1);
      next if (exists $refFamilies{$f2});


      #Get the right directory to search for protocol2 results
      my $sdir = "";
      if (-d "$rootDir/${f1}${fam_div}$f2") {
 	$sdir = "$rootDir/${f1}${fam_div}$f2";
      }
      elsif (-d "$rootDir/${f2}${fam_div}$f1") {
 	$sdir = "$rootDir/${f2}${fam_div}$f1";
      }
      else {
 	die "Coud not find directory comparing $f1 vs $f2.";
      }


      #Now get the top scores between each pair of families
      my $cmd = qq(list_top_gsat_hits -p2d $sdir -b);
      chomp (my $out = `$cmd`);

      my $score = undef;
      if ($out) {
 	my @data = split (/\s+/, $out);

 	$score = sprintf ("%.1f", $data[1]);
      }
      else {
 	$score = "Failed";
      }

      $matrix{$f1}{$f2} = $score;
    }
  }

#  print Data::Dumper->Dump([\%matrix], [qw(*matrix )]);
#  exit;

  print_ref_vs_all_comparisons(\%matrix, \@ref_families, \@all_families);


}






#==========================================================================
################   Subroutines definition beyond ths point   ##############
#==========================================================================

#==========================================================================
#Print Ref vs all compariosn matrix

sub print_ref_vs_all_comparisons {

  my ($inMatrix, $refFams, $allFams) = @_;


  #Print header
  print "-,", join (",", @$refFams), "\n";


  #Print rows
  foreach my $fam (@$allFams) {

    next if ($fam eq $exclude_fam);
    next if (exists $refFamilies{$fam});

    print "$fam";

    foreach my $rfam (@$refFams) {

      next if ($rfam eq $exclude_fam);


      if (exists $inMatrix->{$rfam} && exists $inMatrix->{$rfam}->{$fam}) {
	  print ",", $inMatrix->{$rfam}->{$fam};
      }
      else {
	print ",-";
      }
    }
    print "\n";
  }

}




#==========================================================================
#Print the all vs all comparison matrix

sub print_matrix {

  my ($inMatrix, $fams) = @_;

  #Print header
  print "-,", join (",", @$fams), "\n";

  #Print rows
  foreach my $fam1 (@$fams) {

    next if ($fam1 eq $exclude_fam);

    print "$fam1";

    foreach my $fam2 (@$fams) {

      next if ($fam2 eq $exclude_fam);

      if (exists $inMatrix->{$fam1} && exists $inMatrix->{$fam1}->{$fam2}) {
	  print ",", $inMatrix->{$fam1}->{$fam2};
      }
      else {
	print ",-";
      }
    }
    print "\n";
  }

}




#===========================================================================
#Read command line and print help


sub read_command_line {

    print_help() unless (@ARGV);

    my $status = GetOptions(
	"aa|all-vs-all"        => \$all_vs_all,
	"ra|ref-vs-all"        => \$ref_vs_all,
	"rf|ref-families=s"    => \&read_ref_families,
	"d|root-dir=s"         => \&read_root_dir,
	"t|gsat-cutoff=f"      => \$gsat_cutoff,
	"h|help"               => sub { print_help(); },
	"<>"                   => sub { die "Error: Unknown argument: $_[0]\n"; });
    exit unless ($status);


    die "Error: option -d is mandatory and the given directory must exist.\n" unless ($rootDir && -d $rootDir);

    #Either $all_vs_all or $ref_vs_all is required
    unless ($all_vs_all || $ref_vs_all) {
      die "Error: either option -aa or -ra must be specified.\n";
    }


    #Incompatibilities (-aa,-ra and -aa,rf)
    if ($all_vs_all && $ref_vs_all) {
      die "Error: Options -aa and -ra are not compatible.\n";
    }
    elsif ($all_vs_all && @ref_families) {
      die "Error: Options -aa and -rf are not compatible.\n";
    }


    #If -ra is given, then -rf is required
    if ($ref_vs_all) {
      unless (@ref_families) {
	die "Error: both options -ra and -rf must be given.\n";
      }
    }


    unless (-d $rootDir) {
      die "Error: dir passed to -rd option does not exists";
    }
}



sub read_ref_families {
  my ($opt, $value) = @_;

  @ref_families = split (/,/, $value);
  TCDB::Assorted::validate_tcdb_id (\@ref_families);

}



sub read_root_dir {
  my ($opt, $value) = @_;

  unless (-d $value) {
    die "Error: directory passed to -$opt ($value) does not exist.\n";
  }
  $rootDir = $value;
}





sub print_help {

    my $help = <<'HELP';

Create a comparison matrix with top GSAT scores among
a group of tcdb families. Comparing a few families agaist
all other families is also possible (i.e calculate just a
subset of the rows in the all vs all matrix).

Input parameters:

-aa, --all-vs-all 
   Indicate that all families vs all families will be
   compared.

-ra, --ref-vs-all
   Indicate that only a subset of the families will be
   compared against all families. This option is incompatible
   with -aa.

   NOTE: either -aa or -ra must be given.

-rf, --ref-families {string}
   List of comma separated tcdb families that will be
   compared against all other families. This option is 
   mandatory if -ra is given, and incompatible with -aa.

-d, --root-dir {path}
   Root directory where the results of running protocol2 and
   GSAT are located (This option is mandatory).

-t, gsat-cutoff {float}
   The minimum gsat score to consider an alignment significant
   (Default = 15);

-h, --help
   Display this help. Also displayed if script is run without arguments.

HELP

    print $help;
    exit;
}
