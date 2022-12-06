#!/usr/bin/env perl -w

use warnings;
use strict;
use Data::Dumper;


#==========================================================================
#THis scripts substitutes the links in file results.html to point to
#the output of program xgbhit. THis is a provisional fix so that the
#user can have data to PFAM data.
#
#This will be part of GBLAST soon.
#==========================================================================


#GBlast directory
my $GBdir        = $ARGV[0];
die "GBlast dir not found: $GBdir" unless (-d $GBdir);

#Path to results.html table.
my $htmlFile     = "$GBdir/results.html";
die "File not found: $htmlFile" unless (-f $htmlFile);

my $outFile      = "$GBdir/test.html";

my $xgbFolder    = "$GBdir/xgbhit";
if ($ARGV[1]) {
  $xgbFolder  = $ARGV[1];
}
die "Plots folder not found: $xgbFolder" unless (-d $xgbFolder);


#==========================================================================
#Now replace the accession with the locus tag in the proteome file

open (my $fh2, "<", $htmlFile) || die $!;
open (my $fh3, ">", $outFile) || die $!;
while (<$fh2>) {

  if (/\<td\>\<a href="content\S+"\>(\S+)\<\/a\>\<\/td\>/) {
    #get ssearch directory
    my $locus =$1;
    my $path = getPath2plots($locus);

    print $fh3 "	<td><a href='$path' target='_blank'>$locus</a></td>\n";

  }
  else {
    print $fh3 $_;
  }

}
close $fh2;
close $fh3;



###########################################################################
##                          Functions                                    ##
###########################################################################

#==========================================================================
#Given a genome get the path to the xgbhit results for a given gene.

sub getPath2plots {

  my $prot = shift;

  my $root = "$xgbFolder/$prot";
  die "xgbhit directory not found: $root" unless (-d $root);

  #get the list of files in the directory
  my $dirContent = qx(ls $root/*.faa | xargs basename -s .faa);
  my @accessions = split (/\n/, $dirContent);

  #Remove version from $plot find the proper html page.
  my ($cleanAcc, $version) = split(/\./, $prot);

  my $subject = undef;
  foreach my $p (@accessions) { if ($p ne $cleanAcc) {$subject = $p} }
  die "Coud not find subject protein for $prot" unless ($subject);

  #Remove tag '_tmp' if this is an alignment involving the same protein
  $subject =~ s/_tmp//;

  my $path = "$root/ssearch_${cleanAcc}_vs_$subject/report.html";
  die "Report file not found: $path" unless (-f $path);

  return $path;
}
