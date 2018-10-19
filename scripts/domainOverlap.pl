#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

use Set::Scalar;
use Statistics::Lite;
use TCDB::Assorted;

use Getopt::Long;


#==========================================================================
#
# Based on a set domains shared between two families, calculate
# the average overlap and standard deviation between the domains shared.
#
#==========================================================================




#==========================================================================
#Global variables

#Query families for which the matrix will be built
my $refFams    = [];
my $cmpFams    = [];
my $refFile    = "";
my $cmpFile    = "";
my $inDir      = undef;

#This directory contains the domain content of individual families. This will
#allow me to calculate the overlap between domains matching different families.
my $indFamDomDir = "";


my %evidenceRank = ('noHit' => 0, 'Rescued2' => 1, 'Rescued1' => 2, 'DirectHit' => 3);


#Read command line topology
read_command_line_arguments();

#print Data::Dumper->Dump([$refFams, $cmpFams, $refFile, $cmpFile, $inDir, $indFamDomDir],
#			 [qw(*refFams *cmpFmas *refFile *cmpFile *inDir $indFamDomDir)]);
#exit;





#==========================================================================
#Read the domain within each individual family. This is the result of
#running the program getDomainTopology.pl on single families.

my %famDomains = ();
readFamilyDomains(\%famDomains, $refFams, $cmpFams);

#print Data::Dumper->Dump([\%famDomains ], [qw(*famDomains )]);
#exit;




#==========================================================================
#Get the list of shared domains


my %sharedDomains = ();
identifySharedDomains($refFams, $cmpFams, \%sharedDomains);

#print Data::Dumper->Dump([\%sharedDomains], [qw(*sharedDomains)]);
#exit;




#==========================================================================
#Calculate overlaps between domains

my %overlap = ();
getOverlap(\%sharedDomains, \%overlap);

#print Data::Dumper->Dump([\%overlap ], [qw(*overlap )]);
#exit;




#==========================================================================
#Generate matrix data for Numbers or Excel

saveData2File(\%overlap);





###########################################################################
##                         Functions                                     ##
###########################################################################

sub saveData2File {

  my $data = shift;

#  print "One\n";
#  print Data::Dumper->Dump([$data->{'9.A.14'}->{'9.B.45'}, $data->{'9.B.45'}->{'9.A.14'} ], [qw(*cov1 *cov2)]);
#  exit;

  #Print the header
  print "-\t", join("\t", @$cmpFams), "\n";


 F1:foreach my $f1 (@$refFams) {

    my $f1doms = $famDomains{$f1};

    print "$f1";

  F2:foreach my $f2 (@$cmpFams) {

      my $f2doms = $famDomains{$f2};

      #Detect instances where no common domains exist
      if ((exists $data->{$f1} && exists $data->{$f1}->{$f2} && ! %{ $data->{$f1}->{$f2} }) ||
	  (exists $data->{$f2} && exists $data->{$f1}->{$f1} && ! %{ $data->{$f2}->{$f1} })) {
	print "\t";
	next F2;
      }

      print "\t";


    D1:foreach my $d1 (keys %{ $f1doms }) {
      D2:foreach my $d2 (keys %{ $f2doms }) {

	  my ($n, $mean, $SD,$ ew) = undef;
	  my $params = {};

	  if (exists $data->{$f1} && exists $data->{$f1}->{$f2} &&
	      exists $data->{$f1}->{$f2}->{$d1} &&
	      exists $data->{$f1}->{$f2}->{$d1}->{$d2}) {

	    $params = $data->{$f1}->{$f2}->{$d1}->{$d2};
	  }
	  elsif (exists $data->{$f1} && exists $data->{$f1}->{$f2} &&
		 exists $data->{$f1}->{$f2}->{$d2} &&
		 exists $data->{$f1}->{$f2}->{$d2}->{$d1}) {

	    $params = $data->{$f1}->{$f2}->{$d2}->{$d1};
	  }
	  else { next D2 }


#	  if ($f1 eq "1.A.76" && $f2 eq "2.A.20") {
#	    print Data::Dumper->Dump([$params ], [qw(*params )]);
#	    <STDIN>;
#	  }

	  #Format numbers for printing
	  $n    = $params->{n};
	  $mean = ($params->{m})? sprintf("%.1f", $params->{m}) : undef;
	  $SD   = sprintf("%.1f", $params->{sd});
	  $ew   = sprintf("%.1f", $params->{e});


	  #Catch here cases of no overlap
	  if ($mean) {
	    print " ${d1}-${d2}:$n:$mean:$SD:$ew";
	  }
	  else {
	    print " ${d1}-${d2}:NOoverlap";
	  }
	} #d2
      } #d1
    } #f2

    print "\n";

  } #f1
}


#==========================================================================
#Calculate overlap per pair of families in the user input


sub getOverlap {

  my ($shDoms, $outOverlap) = @_;


  foreach my $f1 (keys %$shDoms) {

    my $f1doms = $famDomains{$f1};

    foreach my $f2 (keys %{ $shDoms->{$f1} }) {

      my $f2doms = $famDomains{$f2};
      my $listDoms = $shDoms->{$f1}->{$f2};

#      die "Empty list of domains for families $f1 and $f2" unless (@{ $listDoms });
#      print Data::Dumper->Dump([$f1, $f1doms, $f2, $f2doms, $listDoms], [qw(*f1 *f1doms *f2 *f2doms *listDoms)]);
#      <STDIN>;


      #-----------------------------------------------------------------
      #Extract coverage data

      my %coverageData = ();
      foreach my $hit (@$listDoms) {

	#
	#$hit has the form:
	# [['PF01027', '29-242','DirectHit'], ['Pf03595', '20-210', 'Rescued2'], ..], [...], ...]
	#

	#Explore each possible pair of domains
      D1:for(my $i1 = 0; $i1 <= ($#{$hit} - 1); $i1++) {
	  my $dom1   = $hit->[$i1];
	  my $domID1 = $dom1->[0];


	D2:for (my $i2 = $i1 + 1; $i2 <= $#{$hit}; $i2++) {
	    my $dom2   = $hit->[$i2];
	    my $domID2 = $dom2->[0];


	    #One domain must be in the family1 and the other in family2 to continue
	    next D2 unless ((exists $f1doms->{$domID1} && exists $f2doms->{$domID2}) ||
			    (exists $f1doms->{$domID2} && exists $f2doms->{$domID1}));

	    my @cov = calculateCoverage($dom1, $dom2);
	    push (@{ $coverageData{$domID1}{$domID2} },  @cov) if (@cov);
	  } #i2
	} #i1
      } #hit

#      print Data::Dumper->Dump([\%coverageData, $f1, $f2 ], [qw(*coverageData *fam1 *fam2 )]);
#      <STDIN>;



      #--------------------------------------------------------------------------
      #Calculate average coverage and evidence weight


      my $stats = getStatistics(\%coverageData);

      $outOverlap->{$f1}->{$f2} = $stats;
    } #f2
  } #f1
}


#==========================================================================
#Based on a list of coverages and evidence weights get the statistis

sub getStatistics {

  my $data = shift;

  my %out = ();

 D1:foreach my $d1 (keys %$data) {
  D2:foreach my $d2 (keys %{ $data->{$d1} }) {

      next D2 unless (@{ $data->{$d1}->{$d2} });

      my @coverage = map { $_->[0] } @{ $data->{$d1}->{$d2} };
      my @evidence = map { $_->[1] } @{ $data->{$d1}->{$d2} };

      my $mean = (scalar @coverage == 1)? $coverage[0] : Statistics::Lite::mean(@coverage) ;     #mean coverage
      my $sd   = (scalar @coverage == 1)? 0.0 : Statistics::Lite::stddev(@coverage);             #SD of coverage
      my $ev   = (scalar @evidence == 1)? $evidence[0] : Statistics::Lite::mean(@evidence);      #mean evidence weight
      my $n    = scalar(@coverage);                                                              #data points

      $out{$d1}{$d2} = {m=>$mean, sd=>$sd, n=>$n, e=>$ev};
    }
  }

  return \%out;

}


#==========================================================================
#Calculate the coverage between a pair of coordinates

sub calculateCoverage {

  my ($dHit1, $dHit2) = @_;

  #
  #dHit1 and dHit2 have the form:
  #   ['PF01027', '29-242','DirectHit']
  # Note that there can be multiple sets of coordinates
  #

  my @coords1 = grep { /\d+\-\d+/} @$dHit1;
  my @coords2 = grep { /\d+\-\d+/} @$dHit2;


  #The evidence rank
  my $rank1 = $evidenceRank{ $dHit1->[$#$dHit1] };
  my $rank2 = $evidenceRank{ $dHit2->[$#$dHit2] };
  my @rankSorted = sort {$a <=> $b} ($rank1, $rank2);


  my @out = ();

  #calculate coverage here
 C1:foreach my $c1 (@coords1) {

    my ($left1, $right1) = split (/-/, $c1);
    my $len1 = ($right1 - $left1) + 1;


  C2:foreach my $c2 (@coords2) {

      my ($left2, $right2) = split (/-/, $c2);
      my $len2 = ($right2 - $left2) + 1;

      my ($first, $second) = undef;


      #Determine which domain is leftmost
      if ($left1 <= $left2) {
	$first  = {left=>$left1, right=>$right1, length=>$len1, rank=>$rank1};
	$second = {left=>$left2, right=>$right2, length=>$len2, rank=>$rank2};
      }
      else {
	$first  = {left=>$left2, right=>$right2, length=>$len2, rank=>$rank2};
	$second = {left=>$left1, right=>$right1, length=>$len1, rank=>$rank1};
      }


      #Skip calculation if domains do not overlap
      if ($first->{right} <= $second->{left}) {
	push (@out, [0.0, $rankSorted[0]]);
	next C2;
      }


      #One domain contains the other
      my $coverage = undef;
      if (($first->{left} <= $second->{left}) &&
	  (($first->{right} >= $second->{right}) || ($first->{right} <= $second->{right}))) {

	$coverage = ($first->{length} >= $second->{length})?
	  ($second->{length}/$first->{length}) : ($first->{length}/$second->{length});
      }


      #normal overlap calculation
      elsif ($first->{right} < $second->{right}) {

	my $longerDom = ($first->{length} >= $second->{length})? $first->{length} : $second->{length};
	$coverage = ($first->{right} - $second->{left} + 1)/$longerDom;
      }
      else {
	print "Unknown condition:\n";
	print Data::Dumper->Dump([$dHit1, $dHit2], [qw(*dHit1 *dHit2)]);
	exit;
      }


      #Prepare the output
      push (@out, [$coverage, $rankSorted[0]]);
    } #c2
  } #C1

#  print Data::Dumper->Dump([$dHit1, $dHit2, \@out ], [qw(*dHit1 *dHit2 *out)]);
#  exit;

  return @out;
}



#==========================================================================
#Given a pair of families identify what domains are shared, do it at
#the level of Pfam IDs and Clans


sub identifySharedDomains {

  my ($rFams, $cFams, $shDoms) = @_;

  foreach my $f1 (@$rFams) {
    foreach my $f2 (@$cFams) {

      #Do not evaluate a family against itself
      next if ($f1 eq $f2);

      #Directories for this pair of families
      my $workDir1 = "$inDir/${f1}_vs_$f2";
      my $workDir2 = "$inDir/${f2}_vs_$f1";
      die "Error: directory not found -> $workDir1" unless (-d $workDir1);
      die "Error: directory not found -> $workDir2" unless (-d $workDir2);


      #Read shared domains between $f1 and $f2
      my @dom1 = findCommonDomains($f1, $f2, $workDir1);
      my @dom2 = findCommonDomains($f2, $f1, $workDir2);


      #combine both arrays
      my @all = (@dom1, @dom2);

      $shDoms->{$f1}->{$f2} = \@all;
    } #f1
  } #f2
}




#==========================================================================
#Given a pair of families in a specific order, identify shared domains
#and clans.

sub findCommonDomains {

  my ($rfam, $cfam, $dir) = @_;

  my $rescueFile = "$dir/$rfam/reports/${rfam}_rescuedDomains.tsv";
  die "Error: no rescue file --> $rescueFile" unless (-f  $rescueFile && !(-z $rescueFile));


  #-----------------------------------------------------------------
  #parse Domain inferences for this file

  #the list of domains hit by these two proteins
  my @out = ();


  #Parse rescue file
  open (my $fh, "<", $rescueFile) || die $!;
  while (<$fh>) {

    chomp;
    next if (/^#/);   #Ignore commented lines


    my ($hit, $prot, @domainData) = split (/\s+/);

#   print Data::Dumper->Dump([$hit, $prot, \@domainData], [qw(*hit *prot *domainData)]);
#   <STDIN>;


    my @commonDomains = ();
    foreach my $pfam (@domainData) {

      my @components = split(/\|/, $pfam);

      #Ignore domains without hits
      next if ($components[-1] eq 'Nohit');

      push (@commonDomains, \@components);
    }


    #Check if protein has at least least two domains
    if (scalar @commonDomains > 1) {
      push (@out, \@commonDomains);
    }

  }
  close $fh;


#  print Data::Dumper->Dump([\@out ], [qw(*out)]);
#  <STDIN>;


  return @out;
}



#==========================================================================
#Read the domains of multiple families and put that information in a
#hash table

sub readFamilyDomains {

  my ($outDoms, $rfams, $cfams) = @_;

  foreach my $fam (@$rfams, @$cfams) {

    my $doms = readDomains($fam, $indFamDomDir);
    $outDoms->{$fam} = $doms;
  }
}


#==========================================================================
#Read the domains matched in a single files and return an array
#with the domains in that family

sub readDomains {

  my ($family, $rootDir) = @_;

  my $infile = "$rootDir/$family/reports/${family}_rescuedDomains.tsv";
  die "Error: file not found --> $infile" unless (-f $infile);

  my %out = ();

  open (my $fh, "<", $infile) || die $!;
  while (<$fh>) {
    chomp;
    if (/^# (\w+):/) {
      $out{$1} = 1;
    }
  }
  close $fh;

  return \%out;
}






#==========================================================================
#Read command line arguments


sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "rf|ref-fams=s"        => \&readRefFams,      #TCIDs that will be used as reference
      "cf|cmp-fams=s"        => \&readCmpFams,      #TCIDs that will that will compared
      "rfile|ref-file=s"     => \&readRefFile,      #File with reference TCIDs
      "cfile|cmp-file=s"     => \&readCmpFile,      #File with comparison TCIDs
      "fd|family-doms=s"     => \&readFamDomDir,    #Directory with the domains per family
      "d|indir=s"            => \&read_indir,       #Directory with all pairwise domain comparisons
      "h|help"  => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);

  #----------------------------------------------------------------------
  #Validate command line arguments


  die "Error: input directory is mandatory!\n" unless ($inDir && -d $inDir);


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


  #Verify that the directory with the domain content per individual family was given
  die "Error: option -id is mandatory" unless ($indFamDomDir);


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
#Read the -id option

sub readFamDomDir {

  my ($opt, $value) = @_;

  die "Error: argument passed to -id is not a directory: $value \n" unless (-d $value);
  $indFamDomDir = $value;
}



#==========================================================================
#Read the -d option

sub read_indir {

  my ($opt, $value) = @_;

  die "Error: input directory must exist" unless (-d $value);

  $inDir = $value;
}





#==========================================================================
#Print help

sub print_help {

  my $help = <<'HELP';

Calculate average overlap and standard deviation of doamins shared between
families.


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

 -d, --indir {directory} (Mandatory)
  Directory with all the domain comparisons between pairs of families. These
  directories must contains the results of running script:
      relateFamiliesWithDomains.p (which runs getDomainTopology.pl)

-h, --help
  Display this help message and exit. This option takes precedence overy any other
  argument.

HELP


  print $help;
  exit;
}
