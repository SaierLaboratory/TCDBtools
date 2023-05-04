#!/usr/bin/env perl -w

use warnings;
use strict;
use Data::Dumper;
use Bio::SeqIO;

#==========================================================================
#     Given a file with sequences in fasta format, run quod in each
#     sequence.
#==========================================================================

die "Input must be a FASTA file with protein sequences as input!\n" if (!@ARGV || (@ARGV && ($ARGV[0] =~ /-h|-help|--help/ || !(-f $ARGV[0]))));

my $seqFile = $ARGV[0];


my $sObj = Bio::SeqIO->new(-file => $seqFile,  -format => 'fasta');

while (my $seq = $sObj->next_seq) {

  my $id  = $seq->primary_id;
  my $str = $seq->seq;
  my $of  = "${id}.png";
  my $cmd = qq(quod.py -q -o $of --xticks 25 -s -- $str);

  print "$cmd\n";
  system $cmd unless (-f $of);
}





#my $hmmtopFile = $ARGV[0];
#my $seqsDir    = $ARGV[1];
#my $xticks     = $ARGV[2];
#my $outdir     = ".";
#
#$xticks = 50 unless ($xticks);
#
#open (my $tmsh, "<", $hmmtopFile) || die $!;
#
#my $cnt = 1;
#while(<$tmsh>) {
#
#  chomp;
#  s/\s+$//;
#
#  next unless ($_);
#
#
#  #Get protein ID
#  my $pid = (/(\S+)\s+(IN|OUT)\s+\d+\s+[\d\s-]+/)? $1 : undef;
#  die "Could not extract pid: $_\n" unless ($pid);
#
#
#  #Save hmmtop output for this sequence in a tmp file
#  my $tmsFile = "$outdir/${pid}.hmmtop";
#  open (my $outh, ">", $tmsFile) || die $!;
#  print $outh "$_\n";
#  close $outh;
#
#
#  #Run quod on this protein and clean trash
#  system "quod.py -q --xticks $xticks -nt +0 -lt $tmsFile -- $seqsDir/${pid}.faa";
#  system "rm $tmsFile" if (-f $tmsFile);
#}
