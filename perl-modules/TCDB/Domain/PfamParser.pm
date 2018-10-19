package TCDB::Domain::PfamParser;

no warnings;
use strict;
use Data::Dumper;

use List::Util qw(sum);
use Ref::Util qw(is_arrayref);
use TCDB::Assorted;
use Class::Struct;



struct ("TCDB::Domain::PfamParser" =>
	{
	 'analysisLevel'      => '$',  #class|subclass|family|subfamily|system
	 'treatAsSuperfamily' => '$',  #Boolean that indicates treatment of the input tcids
	 'pfamFile'           => '$',  #PFAM hmmscan output file to be parsed
	 'domCovCutoff'       => '$',  #Minimum Pfam domain coverage in the alignment
	 'tcCovCutoff'        => '$',  #Minimum coverage of TCDB protein in the alignment
	 'evalueCutoff'       => '$',  #Maximum E-value cut off for domain alignments in the protein
	 'minProtsDom'        => '$',  #Minimum number of proteins tha must have a domain present to be signifcant
	}
       );



#==========================================================================
#Assign default values for object variables



sub minProtsDom {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::PfamParser::minProtsDom';
  my $default = 1;

  if ( $value ) {
    unless($value > 0  && $value <= 1.0) {
      die "Fraction of proteins with sigificant domain hits should be: (0 < fraction <= 1.0)";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}



sub evalueCutoff {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::PfamParser::evalueCutoff';
  my $default = 10;

  if ( $value ) {
    unless($value > 0) {
      die "E-value should be a positive number";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}




sub tcCovCutoff {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::PfamParser::tcCovCutoff';
  my $default = 1e-3;

  if ( $value ) {
    unless($value > 0 || $value <= 1) {
      die "Coverage of protein should be positive and less than 1.0";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}



sub domCovCutoff {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::PfamParser::domCovCutoff';
  my $default = 1e-3;

  if ( $value ) {
    unless($value > 0 || $value <= 1) {
      die "Coverage of domain should be positive and less than 1.0";
    }
    $self->{$objPath } = $value;
  }

  unless ($self->{$objPath }) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}



sub pfamFile {
  my ($self, $file) = @_;

  my $objPath = 'TCDB::Domain::PfamParser::pfamFile';
  my $default = "";

  if ( $file ) {
    die "PFAM output file not found" unless (-f $file && !(-z $file));
    $self->{$objPath} = $file;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}



sub analysisLevel {
  my ($self, $level) = @_;

  my $objPath = 'TCDB::Domain::PfamParser::analysisLevel';
  my $default = "system";

  if ( $level ) {
    unless($level =~ /(class|subclass|family|subfamily|system)/) {
      die "Analysis level should be any of: class|subclass|family|subfamily|system";
    }
    $self->{$objPath} = $level;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}




##########################################################################
##                                                                      ##
##                     Definition of Functions                          ##
##                                                                      ##
##########################################################################



#==========================================================================
#Assumes that Analysis level is 'system'


sub getDomainStatsForUserFamilies {

  my ($self, $fams, $systems, $domFreq,  $domCoords) = @_;

  die "To use this function \$analysisLevel must be: 'system'" unless ($self->analysisLevel eq "system");

#  print Data::Dumper->Dump([$domFreq, $domCoords], [qw(*domFreq *domCoords)]);
#  exit;



  #-----------------------------------------------------------------
  #Sort user reference families.

  my @inFams = keys %{ $systems };

  my @sorted = ();
  my $lenTCIDs = scalar split(/\./, $inFams[0]);

  if    ($lenTCIDs == 1) { @sorted = sort_by_class(@inFams);     }
  elsif ($lenTCIDs == 2) { @sorted = sort_by_subclass(@inFams);  }
  elsif ($lenTCIDs == 3) { @sorted = sort_by_family(@inFams);    }
  elsif ($lenTCIDs == 4) { @sorted = sort_by_subfamily(@inFams); }
  else { @sorted = sort_by_system(@inFams); }

#  print Data::Dumper->Dump([\@sorted ], [qw(*sorted )]);
#  exit;



  #-----------------------------------------------------------------
  #Get statistics per family

  if ($self->treatAsSuperfamily) {

    #calculate Domain Statistics for this family
    my @freq   = ();
    my %coords = ();
    $self->getDomainStats($fams->[0], $fams, $systems->{$fams->[0]}, \@freq, \%coords);

    $domFreq->{$fams->[0]}   = \@freq;
    $domCoords->{$fams->[0]} = \%coords;
  }
  else {
    foreach my $tcid (@sorted) {

      my $sortedSys = $systems->{$tcid};

#      print Data::Dumper->Dump([$sortedSys ], [qw(*sortedSys )]);
#      <STDIN>;
#      next;


      #calculate Domain Statistics for this family
      my @freq   = ();
      my %coords = ();
      $self->getDomainStats($tcid, undef, $sortedSys, \@freq, \%coords);

      $domFreq->{$tcid}   = \@freq;
      $domCoords->{$tcid} = \%coords;
    }
  }
}



sub getDomainStats {

  my ($self, $tcid, $fams, $systems, $freq, $coords) = @_;

  #total number of proteins in Family
  my $totalProts = scalar @{ $systems };
#  print "Total proteins: $totalProts\n";
#  exit;


  #-----------------------------------------------------------------
  #Get the regular expression for extracting pfam hits of the TCDB ID

  die "PFAM file not found: $self->pfamFile" unless (-f $self->pfamFile);

  my ($perlGrep, $regex) = ("", "");

  if ($self->treatAsSuperfamily) {
    $regex = $self->get_regex_for_tcids($fams);
  }
  else {
    $regex = $self->get_regex_for_tcids($tcid);
  }

  $perlGrep = "bzcat " . $self->pfamFile . " | perl -ne 'print if (/($regex)/);'";

#  print "$perlGrep\n";
#  exit;



  #-----------------------------------------------------------------
  #Read PFAM data

  my %pfam = ();

  open (my $pfamh, "-|", $perlGrep) || die $!;
  while (<$pfamh>) {

    chomp;

    my @data = split(/\s+/);

    my ($pfamID, $kk) = split(/\./, $data[1]);
    my $pfamLen  = $data[2];
    my $tcLongID = $data[3];
    my $tcLength = $data[5];
    my $evalue   = $data[6];
    my $dstart   = $data[15];
    my $dend     = $data[16];
    my $tcstart  = $data[19];
    my $tcend    = $data[20];
    my $def      = substr($_, 181);



#    print Data::Dumper->Dump([$pfamID, $pfamLen, $tcLongID, $tcLength, $evalue, $dstart, $dend, $tcstart, $tcend],
#			     [qw(*pfamID *pfamLen *tcLongID *tcLength *evalue *dstart *dend *tcstart *tcend)]);
#    print Data::Dumper->Dump([$self->evalueCutoff, $self->domCovCutoff,  $self->tcCovCutoff],
#                             [qw(*evalCutoff *domCutoff *tcCutoff )]);
#    <STDIN>;



    #Calculate the domain coverage and the coverage of the sequence
    my $domainCov = ($dend  -  $dstart) / $pfamLen;
    my $tcSeqCov  = ($tcend - $tcstart) / $tcLength;



    #Determine if this is a good hit.
    #The domain coverage applies to either the domain or the protein sequence. This is because some times
    #A protein contains just one third of a domain, but that third covers the entire query protein (e.g. 1.A.76.1 vs 1.A.76.2)
    if ($evalue <= $self->evalueCutoff && ($domainCov >= $self->domCovCutoff || $tcSeqCov >= $self->domCovCutoff) && $tcSeqCov >= $self->tcCovCutoff) {
      $pfam{$tcLongID}{$pfamID}++;
      push (@{ $coords->{$tcLongID}->{$pfamID} },
	    {dlen=>$pfamLen,  dstart=>$dstart,  dend=>$dend,
	     plen=>$tcLength, pstart=>$tcstart, pend=>$tcend, def=>$def});
    }
  }
  close $pfamh;

#  print Data::Dumper->Dump([\%pfam ], [qw( *pfam )]);
#  exit;


  #There must be good pfam hits to continue
  return unless (%pfam);


  #-----------------------------------------------------------------
  #Get the counts and frequencies of each domain

  my $protsWithHits = scalar keys %pfam;
  push (@{ $freq }, ["Proteins", $protsWithHits, $protsWithHits / $totalProts]);

  #Get the counts of domain per family
  my %counts = ();
  foreach my $prot (keys %pfam) {
    foreach my $d (keys %{ $pfam{$prot} }) {

      #Each proteins contributes 1 count per different domain
      $counts{$d}++;
    }
  }


#  print Data::Dumper->Dump([ $protsWithHits, \%counts ], [qw( *protsWithDomains *counts )]);
#  exit;



  #get the output frequencies
  foreach my $dom (sort {$counts{$b} <=> $counts{$a}} keys %counts) {
    my $frac = $counts{$dom} / $protsWithHits;
    next if ($frac < $self->minProtsDom);
    push (@{ $freq }, [$dom, $counts{$dom}, $frac]);
  }

#  print Data::Dumper->Dump([\%counts, $freq ], [qw( *counts *freq)]);
#  print Data::Dumper->Dump([$coords ], [qw(*coords )]);
#  exit;
}



sub get_regex_for_tcids {

  my ($self, $tcids) = @_;

  my $regex = "";

  if (is_arrayref $tcids) {

    my @tmp = ();

    foreach my $tcid (@$tcids) {
      push (@tmp, $self->get_regex_for_tcid($tcid));
    }

    $regex = join("|", @tmp);
  }
  else {

    $regex = $self->get_regex_for_tcid($tcids);

  }

  return $regex;
}






sub get_regex_for_tcid {

  my ($self, $tcid) = @_;

  die "Function must receive a tcid: $tcid" unless ($tcid);


  my $regex = "";
  my $lenID = scalar split(/\./, $tcid);

  if ($lenID <= 4) {
    $regex = qq(\\s$tcid\\.);
  }
  else {
    $regex = qq(\\s$tcid\\s);
  }

  return $regex;

}





1;
