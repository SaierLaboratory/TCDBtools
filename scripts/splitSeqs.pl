#!/usr/bin/env perl -w

use warnings;
use strict;
use Data::Dumper;

use Bio::SearchIO;
use Bio::SeqIO;


my $infile = $ARGV[0];
my $outdir = ".";


my $obj  = Bio::SeqIO->new(-file => $infile , -format => "fasta");

while(my $seqObj = $obj->next_seq) {

  my $id  = $seqObj->primary_id;
  my $seq = $seqObj->seq;

  my $outfile = "$outdir/${id}.faa";
  open (my $outh, ">", $outfile) || die $!;
  print $outh ">$id\n$seq\n";
  close $outh;
}






