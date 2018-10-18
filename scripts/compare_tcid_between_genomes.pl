#!/usr/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use Set::Scalar;



#==========================================================================
#Global variables

my $level   = "system";
my @infiles = ();
my @labels  = ();
my $outfile = "stats.txt";
my $famAbbrevFile = "Family_abbreviations.txt";

#==========================================================================
#read comand line arguments

read_command_line();

#print Data::Dumper->Dump([$level, \@infiles, \@labels, $outfile ], [qw(*level *infiles *labels *outfile )]);
#exit;


#==========================================================================
#Read Family abbreviations

my %FamAbbrev = ();

open(my $fh, "<", $famAbbrevFile) || die $!;
while(<$fh>) {

  chomp;
  next if (/^family/);

  my ($fam, $abbrev) = split(/\t/, $_);

  if ($fam && $abbrev) {
    $FamAbbrev{$fam} = $abbrev;
  }
}
close $fh;

#print Data::Dumper->Dump([\%FamAbbrev], [qw(*FamAbbrev )]);
#exit;


#===========================================================================
#Read each file and create a hash with all the different resulting sets
#of TCIDs

my %sets = ();

foreach my $idx (0 .. $#infiles) {

  my $file  = $infiles[$idx];
  my $label = $labels[$idx];
#  print "$file\n";

  open (my $fileh, "<", $file) || die $!;
  chomp(my @tmpLines = <$fileh>);
  close $fileh;

  my @lines = ();

  if ($level eq "class") {
    foreach my $tcid (@tmpLines) {
      my @comp = split (/\./, $tcid);
      push (@lines, $comp[0]);
    }
  }
  elsif ($level eq "subclass") {
    foreach my $tcid (@tmpLines) {
      my @comp = split (/\./, $tcid);
      push (@lines, "$comp[0]\.$comp[1]");
    }
  }
  elsif ($level eq "family") {
    foreach my $tcid (@tmpLines) {
      my @comp = split (/\./, $tcid);
      push (@lines, "$comp[0]\.$comp[1]\.$comp[2]");
    }
  }
  elsif ($level eq "subfamily") {
    foreach my $tcid (@tmpLines) {
      my @comp = split (/\./, $tcid);
      push (@lines, "$comp[0]\.$comp[1]\.$comp[2]\.$comp[3]");
    }
  }
  else {
    @lines = @tmpLines;
  }

#  print Data::Dumper->Dump([\@lines ], [qw(*lines)]);
#  exit;


  my %uniqLines = map {$_ => 1} @lines;
  my @systems = keys %uniqLines;

  my $set = Set::Scalar->new(@systems);
  $sets{$label} = $set;

}

#print Data::Dumper->Dump([\%sets ], [qw( *sets )]);
#exit;



#==========================================================================
#The comparisons will be done at all levels of the TCID


#*** Intersections ***

my @results = ();

#The intesection among all genomes
global_intersection(\@labels, \%sets, \@results);


#Pairwise intersections
pairwise_intersections(\@labels, \%sets, \@results);


#Unique families per genome
unique_families(\@labels, \%sets, \@results);


#print Data::Dumper->Dump([@results ], [qw(*results )]);
#exit;



#*** print counts ***

open (my $outh, ">", $outfile) || die $!;
foreach my $r (@results) {
  print $outh $r->[0], "\n";
  print $outh print_data($r->[1]), "\n\n";
}
close $outh;


#==========================================================================


###########################################################################
##                                                                       ##
##                        Subroutine definitions                         ##
##                                                                       ##
###########################################################################



#==========================================================================
#Call back function to control the display of sets (4 families per line)

sub print_data {

  my $set = shift;

  my $out = "";
  my $cnt = 1;

  foreach my $element (sort $set->elements) {

    #Extract the family abbreveation for this tcdid
    my @comp = split (/\./, $element);

    my $abbrev = "";
    if (scalar @comp >= 3) {
      my $fam = "$comp[0]\.$comp[1]\.$comp[2]";
      $abbrev = $FamAbbrev{$fam} if (exists $FamAbbrev{$fam});
    }

    if ($cnt % 3 == 0) {

      if ($abbrev) {
	$out .= "${element} ($abbrev)\n";
      }
      else {
	$out .= "${element} (N/A)\n";
      }
    }
    else {
      if ($abbrev) {
	$out .= "${element} ($abbrev)\t";
      }
      else {
	$out .= "${element} (N/A)\t";
      }
    }
    $cnt++;
  }

  return $out;
};





#==========================================================================
#Get the families unique to each family


sub unique_families {

  my ($tags, $sets, $res) = @_;


 S1:foreach my $i1 (0 .. $#{ $tags }) {

    my $tag1 = $tags->[$i1];
    my $seto = $sets->{$tag1};


    my $label = $tag1;

  S2:foreach my $i2 (0 .. $#{ $tags }) {

      my $tag2 = $tags->[$i2];

      next S2 if ($tag1 eq $tag2);

      $label .= " - $tag2";
      my $set2 = $sets->{$tag2};


      $seto = $seto - $set2;
    }

    $label .= ": (" . scalar @$seto . ")";
    push (@$res, [$label, $seto]);
  }
}


#==========================================================================
#Get the intesection between each pair of genomes

sub pairwise_intersections {

  my ($tags, $sets, $res) = @_;


 S1:foreach my $i1 (0 .. $#{ $tags } - 1) {

    my $tag1 = $tags->[$i1];
    my $set1 = $sets->{$tag1};

  S2:foreach my $i2 ($i1 + 1 .. $#{ $tags }) {

      my $tag2 = $tags->[$i2];
      my $set2 = $sets->{$tag2};

      my $label = "($tag1 * $tag2)";
      my $seto = $set1 * $set2;


    REST:foreach my $tag (@$tags) {

	next REST if ($tag eq $tag1 || $tag eq $tag2);

	$label .= " - $tag";

	my $set3 = $sets->{$tag};
	$seto = $seto - $set3;
      }

      $label .= ": (" . scalar @$seto . ")";
      push (@$res, [$label, $seto]);

    } #S2
  }#S1
}



#==========================================================================
#Get the intersections among all genomes

sub global_intersection {

  my ($tags, $sets, $res) = @_;

  my $out = $sets->{$tags->[0]};

  foreach my $i1 (1 .. $#{ $tags }) {

    my $set2 = $sets{$tags->[$i1]};
    $out = $out * $set2;
  }

  my $tag1 = join (" * ", @$tags) . ": (" . scalar @$out . ")";

  push(@$res, [$tag1, $out]);
}



#===========================================================================
#Read command line and print help


sub read_command_line {

    print_help() unless (@ARGV);

    my @usrFiles;

    my $status = GetOptions(
	"l|level=s"     => \&read_level,
	"f|infile=s"    => \@usrFiles,
	"o|outfile=s"   => \$outfile,
	"h|help"        => sub { print_help(); },
	"<>"            => sub { die "Error: Unknown argument: $_[0]\n"; });
    exit unless ($status);

    die "Error: at least two files with TCIDs are needed." unless (scalar @usrFiles >= 2);

    #All files must exist
    foreach my $file (@usrFiles) {

      my ($label, $infile) = split(/:/, $file);
      die "No label detected for file: $file" unless ($label);
      die "File not found or empty" if (-z $infile || !(-f $infile));

      push (@infiles, $infile);
      push (@labels,  $label);
    }
}



sub read_level {

  my ($opt, $value) = @_;

  unless ($value =~ /(class|subclass|family|subfamily|system)/) {
    die "Unknown analysis level: $value";
  }

  $level = $value;
}




sub print_help {

    my $help = <<'HELP';

This script identifies sets of shared or unique tcdb IDs between different files
of identifiers:

-l, --level {string}
  Level of analysis, it can be at either: class, subclass, family,
  subfamily and system level (default: system).

-f, --infiles {path}
  File with 1 TCDB ID per line. Many files can be added by repeating the -f
  option as necessary. At least two files must be provided. The format should
  be:  Label:path_to_file
  Where Labels is a string without spaces that will help identify the genome
  From which the input file (path_to_file) was derived
  (Mandatory)

-o, --outfile {path}
   Output file where the output numbers will be saved.
   (Default: stats.txt)

-h, --help
   Display this help. Also displayed if script is run without arguments.

HELP

    print $help;
    exit;
}

