#!/usr/bin/env perl -w


use warnings;
use strict;
use Data::Dumper;

use Getopt::Long;
use Bio::SeqIO;
use List::Util qw( shuffle );

#==========================================================================
#   This script compares two files with sequences in fasta format.
#   The files can have the sequences in random order, but if the
#   sequences are the same, this program will detect it.
#
#   The program can also randomize the order of any of the files,
#   which is useful to test the robustness of programs to meaningless
#   changes in the input sequences (e.g. ssearch36, blast+, etc.)
#==========================================================================

#f102 vs s102
#2.A.102.1.1,  2.A.112.1.14
#2.A.102.2.10, 2.A.112.1.11
#2.A.102.4.3,  2.A.112.1.12
#2.A.102.2.8,  2.A.112.1.13
#2.A.102.4.6,  2.A.112.1.13
#2.A.102.2.2,  2.A.112.1.2

#f102 vs s112
#NONE

#f112 vs s102
#2.A.102.2.8,  2.A.112.1.13
#2.A.102.2.2,  2.A.112.1.2
#2.A.102.1.1,  2.A.112.1.14
#2.A.102.4.3,  2.A.112.1.12
#2.A.102.2.10, 2.A.112.1.11
#2.A.102.4.6,  2.A.112.1.13

#GOOD
#ssearch36 -k 10000 -E 0.001 -m 8C domainFragments1it_102.faa  family-2.A.102.faa > ssearch_f102_s102.out
#ssearch36 -k 10000 -E 0.001 -m 8C domainFragments1it_112.faa  family-2.A.102.faa > ssearch_f112_s102.out


#BAD
#ssearch36 -k 10000 -E 0.001 -m 8C domainFragments1it_102.faa  family-2.A.112.faa > ssearch_f102_s112.out
#ssearch36 -k 10000 -E 0.001 -m 8C domainFragments1it_112.faa  family-2.A.112.faa > ssearch_f112_s112.out





my $file1 = undef;
my $file2 = undef;
my $randomize = 0;

read_command_line_arguments();


#print Data::Dumper->Dump([$file1, $file2, $randomize ], [qw(*file1 *file2 *randomize )]);
#exit;


if ($randomize == 0) {


  my @seqs1 = ();
  read_sequence_file($file1, \@seqs1);


  my @seqs2 = ();
  read_sequence_file($file2, \@seqs2);



  #-----------------------------------------------------------------
  #Test if the number of sequences is the same

  compareNumberOfSequences(\@seqs1, \@seqs2);



  #-----------------------------------------------------------------
  #Compare the sequences

  compareProteinSequences(\@seqs1, \@seqs2);


  print "Sequence files are identical!\n   $file1\n   $file2\n"

}


#randomize first sequence
elsif ($randomize == 1) {

  my @seqs1 = ();
  read_sequence_file($file1, \@seqs1);


  my @newOrder = shuffle(@seqs1);

  #Print shuffled sequences to STDOUT
  foreach my $seq (@newOrder) {
    print $seq->[1];
  }
}


#Randomize second sequence
elsif {

  my @seqs2 = ();
  read_sequence_file($file2, \@seqs2);


  my @newOrder = shuffle(@seqs2);

  #Print shuffled sequences to STDOUT
  foreach my $seq (@newOrder) {
    print $seq->[1];
  }
}






###########################################################################
##                    Subroutine  definitions                            ##
###########################################################################


#==========================================================================
#Compare the sequences in both files

sub compareProteinSequences {

  my ($s1, $s2) = @_;

  my %seqsFile1 = map {$_->[0] => $_->[1]} @$s1;
  my %seqsFile2 = map {$_->[0] => $_->[1]} @$s2;

#  print Data::Dumper->Dump([\%seqsFile2], [qw(*seqsFile2)]);


  #-----------------------------------------------------------------
  #Check if files have the same number of sequences with unique IDs.
  #
  #NOTE:
  #  Two files could have the same number of sequences but in one
  #  file there could be repeated sequces.

  my $nseqs1 = scalar keys %seqsFile1;
  my $nseqs2 = scalar keys %seqsFile2;

  if ($nseqs1 != $nseqs2) {

    my $errMsg1 = <<ERR1;
Files have different number of sequences with unique IDs:
  $file1 --> $nseqs1
  $file2 --> $nseqs2
ERR1

    die $errMsg1 if ($nseqs1 != $nseqs2);
  }



  #-----------------------------------------------------------------
  #Now Compare the sequences

  foreach my $id1 (keys %seqsFile1) {

    my $seq1 = $seqsFile1{$id1};
    my $seq2 = "";

    if (exists $seqsFile2{$id1}) {
      $seq2 = $seqsFile2{$id1};

      my $errMsg2 = <<ERR2;

Files have different sequences with same id:

$seq1

$seq2

ERR2

      die $errMsg2 if ($seq1 ne $seq2);
    }
    else {
      die "Files are Diferent!\n   File $file2 does not have sequence $id1\n";
    }
  }
}






#==========================================================================
#Compare number of sequences

sub compareNumberOfSequences {

  my ($s1, $s2) = @_;

  my $nseqs1 = scalar @$s1;
  my $nseqs2 = scalar @$s2;

  my $errMsg = <<ERR;
  Files have different number of sequences:
  $file1 --> $nseqs1
  $file2 --> $nseqs2
ERR

  die $errMsg if ($nseqs1 != $nseqs2);

}




#==========================================================================
#Read a file in fasta format and put the sequences in a hash table

sub read_sequence_file {

  my ($file, $output) = @_;

  my $seqObj = Bio::SeqIO->new(-file   =>  $file,
			       -format => 'Fasta');

  while (my $seq = $seqObj->next_seq()) {

    my $id    = $seq->display_id;
    my $fasta = ">$id " . $seq->desc . "\n" . $seq->seq . "\n";

    push (@{ $output }, [$id, $fasta]);
  }
}






#==========================================================================
#Read commandline arguments

sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "f1|file1=s"           => \&read_file1,       #TCIDs of families to analyze (comma separated)
      "f2|file2=s"           => \&read_file2,
      "r|randomize=s"        => \&read_randomize,
      "h|help"               => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);

  #----------------------------------------------------------------------
  #Validate command line arguments


  die "Error: both files are mandatory!\n" unless ($file1 && $file2);
}


#==========================================================================
#Read the -f option

sub read_file1 {

  my ($opt, $value) = @_;

  die "Error: First file does not exist --> $value" unless (-f $value && !(-z $value));

  $file1 = $value;
}


#==========================================================================
#Read the -d option

sub read_file2 {

  my ($opt, $value) = @_;

  die "Error: Second file does not exist --> $value" unless (-f $value && !(-z $value));

  $file2 = $value;
}


#==========================================================================
#Read the -s option

sub read_randomize {

  my ($opt, $value) = @_;

  die "Error: invalid randomization code --> $value" unless ($value =~ /^[0-2]$/);

  $randomize = $value;
}



#==========================================================================
#Print help

sub print_help {

  my $help = <<'HELP';

Compare a pair of seqeunce files in fasta format. The program can also
randomize the either of either sequence. This can be useful to test
the robustness of other programs to meaningles change in the input
(e.g. ssearch36 or blast).

-f1, --file1 {file} (Mandatory)
  First file to be compared

-f2, --file2 (file} (Mandatory)
  Second file to be compared

-r, --randomize {integer} (Optional)
  Randomization code to be used:
    0:  No randomization, proceed to file comparison (default)
    1:  Randomize only first file
    2:  Randomize only second file
  If this option is used files will not be compared, and the output will
  be the shuffled sequenced file printed directly to standard output.

HELP


  print $help;
  exit;
}
