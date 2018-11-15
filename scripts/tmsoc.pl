#!/usr/bin/env perl -w


#---------------------------------------------------------------------------
# Singapore, 21st October 2011
#
# ACADEMIC SOFTWARE LICENSE for TMSOC
# ***********************************
#
# Copyright:
# Wing-Cheong Wong, Sebastian Maurer-Stroh, Georg Schneider, Frank Eisenhaber
# Bioinformatics Institute (BII) A*STAR Singapore
#
# This software is subjected to license restrictions. The software is
# provided as is. In its present form and without written consent of
# the copyright holder, the software is provided in accordance with GPL
# (www.gnu.org/licenses/gpl.html). Among other issues, this excludes any
# warranty and indemnity when using this software. For a commercial license,
# please approach the authors (via email to wongwc@bii.a-star.edu.sg).
#
# When publishing results with this software, please refer to:
# Wing-Cheong Wong, Sebastian Maurer-Stroh, Frank Eisenhaber, 2011,
# "Not all transmembrane helices are born equal: Towards the extension of
# the sequence homology concept to the membrane proteins", Biology Direct
#
# When observing bugs or strange behavior of the software, the user is
# encouraged to contact the authors.
#
#---------------------------------------------------------------------------
#
# This program was modified by Arturo Medrano (UCSD) to take command line
# options and use HHMTOP predicted TMSs, thus removing the need, but keeping
# the option, for the user to provide the TMS.
#
# Date: 11/11/2018
#---------------------------------------------------------------------------


use warnings;
use strict;	# Forces variables to be declared
use Data::Dumper;
use Getopt::Long;

use generateTMclassification;


my $i;
my $j;
my $k;
my @tmp1 = ();
my @tmp2 = ();
my @seqname = ();
my @FASTAseq = ();
my @TMsegments = ();
my $sequence = "";
my $segments = "";
my $resultsRef;
my @results;
my $newFASTAseq;

#Variables added by Arturo Medrano
my $seqFile    = "";
my $tmsFile    = "";
my $outdir     = ".";
my $outfile    = "";
my $predTMS    = 1;  #Indicates whether TMS should be predicted by HMMTOP
my $tmsFormat  = "HMMTOP";
my $onlyMasked = 0;  #Indicates whether only masked sequences will be printed

my %noTMS = ();

#Read command line
read_command_line();


#print Data::Dumper->Dump([$seqFile, $tmsFile, $predTMS, $tmsFormat, $onlyMasked, $outdir, $outfile ],
#			 [qw(*seqFile *tmsFile *predTMS *tmsFormat *onlyMasked *outdir *outfile)]);
#exit;



#If TMS are going to be predicted, make sure to remove sequences with 0 TMSs, as
#it doesn't make sense to run the program in those sequences.
predictTMS() if ($predTMS);


#Read TMS coordinates
readCVSfileCoords() unless (scalar(@TMsegments) > 0);


#Read FASTA sequences
readSequences();


#Estimate simple and complex TMSs from here
if (scalar(@FASTAseq)==scalar(@TMsegments) && scalar(@FASTAseq)>0) {

  my $fh = undef;
  if ($outfile) {
    open ($fh, ">", $outfile) || die $!;
  }
  else { $fh = *STDOUT; }


  for ($i=0; $i<scalar(@FASTAseq); $i++) {

    # Needs to offset the position values by -1; assumes first position starts at 1
    @tmp1 = split(/\s+/, $TMsegments[$i]);
    my $segment = "";
    for ($j=0; $j<scalar(@tmp1); $j++) {
      @tmp2 = split(/\,/, $tmp1[$j]);
      if ($segment eq "") { $segment = ($tmp2[0]-1).",".($tmp2[1]-1); }
      else { $segment = $segment." ".($tmp2[0]-1).",".($tmp2[1]-1); }
    }

    # Classify the TM region(s)
    ($resultsRef, $newFASTAseq) = generateTMclassification($segment, $FASTAseq[$i]);
    @results = @$resultsRef;
    if (uc($FASTAseq[$i]) eq uc($newFASTAseq)) {
      $newFASTAseq = "none";
    }

    # Output the results
    if ($onlyMasked) {
      print $fh $seqname[$i]."\n".$newFASTAseq."\n\n";
    }
    else {
      print $fh "1. TM segment(s) summary:\n";
      for ($j=0; $j<scalar(@results); $j++) {
	@tmp1 = split(/\;/,$results[$j]);
	print $fh $results[$j]."\n";
      }
      print $fh "2. Masked FASTA sequence:\n".$seqname[$i]."\n".$newFASTAseq."\n\n";
    }
  }
}
else {
  print "There are ".scalar(@FASTAseq)." sequences but ".scalar(@TMsegments)." associated TM segments. Please check input files(s).\n";
}





###########################################################################
###                             FUNCTIONS                               ###
###########################################################################


#==========================================================================
#Read sequence file and avoid sequences without TMSs.


sub readSequences {

  my $cnt = 1;

  open (my $MYFILE1, "<", $seqFile) || die $!;
  while (<$MYFILE1>) {

    chomp;
    my $line = $_;
    $line =~ s/\r//g;	# remove linefeed char
    $line =~ s/\n//g;	# remove new line char

    if ($line =~ /^>/) {

      if ($sequence ne "") {

	push(@FASTAseq, $sequence) unless (exists $noTMS{$cnt});
	$sequence = "";

#	print "|$cnt|\n", Data::Dumper->Dump([\@seqname, \@FASTAseq, \%noTMS], [qw(*seqname *FASTAseq *noTMS)]);
#	<STDIN>;

	$cnt++;
      }

      @tmp1 = split(/\s+/, $line);
      push(@seqname, $tmp1[0]) unless (exists $noTMS{$cnt});
    }
    else { $sequence = $sequence.$line; }
  }
  close($MYFILE1);

  push(@FASTAseq, $sequence) unless (exists $noTMS{$cnt});

#  print Data::Dumper->Dump([\@TMsegments, \%noTMS, \@seqname, \@FASTAseq, $cnt], [qw(*TMsegments *noTMS *seqname *FASTAseq *cnt)]);
#  exit;
}



#==========================================================================
#Reast TMS coordiantes in CSV format

sub readCVSfileCoords {

  # Read TM segments associated to FASTA sequences
  open (my $MYFILE2, "<", $tmsFile) || die $!;
  while (<$MYFILE2>) {
    chomp;
    my $line = $_;
    $line =~ s/\r//g;	# remove linefeed char
    $line =~ s/\n//g;	# remove new line char

    @tmp1 = split(/\s+/, $line);
    if (scalar(@tmp1) > 0) {
      push(@TMsegments, $line);
    }
  }
  close($MYFILE2);
}




#==========================================================================
#Run HMMTOP on sequences and put them in CSV format for TMSOC
#also remove sequences with 0 TMS, as TMSOC won't run properly on
#Those sequences.

sub predictTMS {

  #Generate random integer between 10,0000 and 100,000
  my $randNumber = int(rand 90000) + 10000;


  #run HMMTOP
  my $hmmtopFile = "$outdir/hmmtop_${randNumber}.out";
  my $cmd = qq(hmmtop -if=$seqFile -of=$hmmtopFile -sf=FAS -pi=spred -is=pseudo);
  system $cmd unless (-f $hmmtopFile && !(-z $hmmtopFile));


  #Parse HMMTOP output and get the coordinates in CSV format
  open(my $fh, "<", $hmmtopFile) || die $!;

  my $cnt = 1;
  while(<$fh>) {
    chomp;

    #Remove trailing white spaces
    s/\s+$//;

    my $n   = undef;
    my $tms = undef;

    if (/\s+(IN|OUT)\s+(\d+)/) {
      $n = $2;
      unless ($n == 0) {
	if (/(IN|OUT)\s+\d+\s+(.+)$/) {
	  $tms = $2;
	}
      }
    }

    if ($n == 0) { $noTMS{$cnt} = 1; }
    else {

      die "Error, there should be TMS: n=$n tms=$tms\nline:$_\n" unless ($tms);

      #Format the hhmtop coordiantes into CSV coordinates
      my $cvs = format_hmmtop2cvs($tms);
      die "Could not format hmmtop coords to cvs: $tms" unless ($cvs);

      push (@TMsegments, $cvs);
      $cnt++;
    }
  }
  close $fh;

  system "rm $hmmtopFile" if (-f $hmmtopFile);
}


#==========================================================================
#Format HMMTOP coordinates string to CVS for TMSOC

sub format_hmmtop2cvs {

  my $hmmtopStr = shift;

  my @tms = split(/\s+/, $hmmtopStr);
  my @cvs = ();

  for(my $i=0; $i <= $#tms - 1; $i += 2) {

    my $left  = $tms[$i];
    my $right = $tms[$i + 1];

    push(@cvs, "${left},${right}");
  }

  my $str = join("\t", @cvs);
}



#==========================================================================
#Read the command line arguments

sub read_command_line {

  print_help() unless (@ARGV);

  my $status = GetOptions(
        "s|inseqs=s"           => \&read_inseqs,
	"t|tms=s"              => \&read_tms,
	"f|tms-format=s"       => \&read_tmsFormat,
	"p|predict-tms=s"      => \&read_predict,
	"m|masked-seqs-only=s" => \&read_mask,
	"d|outdir=s"           => \$outdir,
	"o|outfile=s"          => \&read_outfile,
        "h|help"               => sub { print_help(); },
        "<>"                   => sub { die "Error: Unknown argument: $_[0]\n"; });
  exit unless ($status);

  #If TMSs are given, make sure to avoid running HMMTOP, even if the
  #option is explicitly given by the user.
  $predTMS = ($tmsFile)? 0 : 1;


  #If no output file is given, do not generate output directory because
  #output will be sent to STDOUT
  if ($outfile) {
    system "mkdir -p $outdir" unless (-d $outdir);

    #Generate the full path to the outputfile
    $outfile = "$outdir/$outfile" if ($outfile);
  }
  else {

    #If there is no outfile, generate an error if the output directory was given.
    #There is no point in generating and empty directory if output is going to STDOUT.
    if ($outdir ne ".") {
      die "Error: if output dir is different from the default, results can't be sent to STDOUT!\n";
    }
  }
}


#==========================================================================
#Read option -s

sub read_inseqs {
  my ($opt, $value) = @_;

  die "Error with Option (-$opt): sequence file does not exist or is empty!\n" unless (-f $value && !(-z $value));

  $seqFile = $value;
}


#==========================================================================
#Read option -t

sub read_tms {
  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
    die "Error with option (-$opt): File with TMS does not exist or is empty!";
  }

  $tmsFile = $value;
}


#==========================================================================
#Read option -f

sub read_tmsFormat {
  my ($opt, $value) = @_;

  my $tmp = uc($value);
  if ($tmp =~ /^(HMMTOP|CSV)$/) {

    #If option is HMMTOP, variable just takes the default value
    $tmsFormat = $tmp;
  }
  else {
    die "Error with Option (-$opt): valid options are hmmtop|csv!\n";
  }
}





#==========================================================================
#Read option -p

sub read_predict {
  my ($opt, $value) = @_;

  my $tmp = uc($value);
  if ($tmp =~ /^[TF]$/) {

    #If option is T, variable just takes the default value
    $predTMS = 0 if ($tmp eq "F");
  }
  else {
    die "Error: with Option (-$opt): valid options are T|F!\n";
  }
}


#==========================================================================
#Read option -m

sub read_mask {
  my ($opt, $value) = @_;

  my $tmp = uc($value);
  if ($tmp =~ /^[TF]$/) {

    #If option is F, variable just takes the default value
    $onlyMasked = 1 if ($tmp eq "T");
  }
  else {
    die "Error with Option (-$opt): valid options are T|F!";
  }
}


#==========================================================================
#Read option -o

sub read_outfile {
  my ($opt, $value) = @_;

  $outfile = $value;
}









#==========================================================================
#Read option -h

sub print_help {

    my $help = <<'HELP';

This program runs the script TMSOC to determine whether TMSs are
simple or complex.

Input parameters:

-s, --inseqs { FILE } (Mandatory)
   Input file in FASTA format with all sequences to be analyzed.

-t, --tms { FILE } (Optional; Default: predict with HMMTOP);
   Input file with TMS coordinates either in HMMTOP format
   or TMSOC format (comma separated pairs). Each line in this
   file must have the same order the the sequeces passed through
   the -s option.

-f, --tms-format { string } (Optional; Default: hmmtop)
   Indicates the format of the input TMS file. Options are
   'hmmtop' and 'csv'. Format csv puts cordinates separated
   by commas and TMS separated by spaces (e.g.; 3,20 26,42 ...)

-p, --predict-tms { T|F } (Optional; Default: T)
   If true (T), run HMMTOP to predic TMSs. This option is forced
   to be true (T) if no file with TMSs is provided, or false (F)
   if the file with TMSs is provided.

-m, --masked-seqs-only { T|F } (Optional; Default: F)
   Indicate whether the output will contain only sequences with
   simple TMSs masked.

-d,  --outdir { DIR } (Optional; Default: '.')
   Output directory where results will be saved. Use this option
   only if you are providing an output file name, otherwise
   an error will be generated.

-o, --outfile (Optional; Default: STDOUT)
   File in the output directory where results will be saved.
   If no output file is given, output is sent to standard output.

-h, --help
   Display this help. Also displayed if script is run without arguments.
   This option takes precedence over any other option.

HELP

    print $help;
    exit;
}
