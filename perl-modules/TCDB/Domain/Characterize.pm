package TCDB::Domain::Characterize;

use strict;
no warnings;
use Data::Dumper;

use TCDB::Assorted;

use List::Util qw(sum);
use Ref::Util qw(is_hashref);
use Class::Struct;



#==========================================================================
# This module takes the proteins in a family containing PFAM domains,
# extracts the sequences in that domain and ssearch36 them
# against the proteins that did not match that domain.
#
# The purpose is to verify whether proteins without domain matches
# actually have the missing domains but the HMM failed to identify it.
#
# Some remote homologs may need two iterations in order to rescue their
# domains
#
#--------------------------------------------------------------------------
#
# Written by:  Arturo Medrano
# Date:        8/17/2017
#
#==========================================================================



struct ("TCDB::Domain::Characterize" =>
	{
	 'domFreq'            => '%',  #array with the frequency of domain hits
	 'domCoords'          => '%',  #Hash with the coordinates of the domains
	 'tcids'              => '%',  #Hash wit the TCIDs and accessions
	 'searchWith'         => '$',  #Program to use in the searches (default ssearch36)
	 'evalue'             => '$',  #Minimum evalue for domain-proteins alignments
	 'domCovCutoff'       => '$',  #minimum domain coverage for rescued domains
	 'treatAsSuperfamily' => '$',  #tread a least of families as a superfamily (0|1)
	 'refDomains'         => '$',  #query all input families for these domains instead of querying for a families domains.

	 #Variables with paths
	 'blastdb'            => '$',  #Full path to the blast DB that will be used to extract sequences
	 "tcdbFaa"            => '$',  #Full path to the fasta files with all sequences in TCDB
	 'rootDir'            => '$',  #The root directory where results will be saved
	 'protSeqs1it'        => '$',  #Fasta file with all sequences of the reference TCDB family
	 'domSeqs1it'         => '$',  #File with the sequences of the PFAM domains found in the family (direct hits)
	 'aln1stCycle'        => '$',  #File with the domain_vs_proteins comaprisons of the 1st rescue cycle
	 'idsFile2it'         => '$',  #File with the Ids of the proteins with missing domains after 1st rescue cycle
	 'protSeqs2it'        => '$',  #File with the sequences of the ids in variable $self->idsFile2it
	 'domSeqs2it'         => '$',  #File with the sequences of the PFAM domains found in the family (direct hits)
	 'aln2ndCycle'        => '$'   #File with the domain_vs_proteins comaprisons of the 2nd rescue cycle
	});




#==========================================================================
#Default values


#The frequences of the PFAM domains within a family of proteins
sub domFreq {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::Characterize::domFreq';
  my $default = {};

  if ( $value ) {
    unless (is_hashref $value) {
      die "domFreq: value should be a reference to a non-empty hash of domain coordinates -> ";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}





#The coordinates of the PFAM hits in terms of the domain and the protein sequence
sub domCoords {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::Characterize::domCoords';
  my $default = {};

  if ( $value ) {
    unless (is_hashref $value) {
      die "domCoords: value should be a reference to a non-empty hash of domain coordinates -> ";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}



#The different TCIDs and Accession for all proteins in a family;
#The format of the hash should be identical as the output generated with the function:
#     TCDB::Assorted::getSystemAccessions();
sub tcids {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::Characterize::tcids';
  my $default = {};

  if ( $value ) {
    unless (is_hashref $value) {
      die "tcids: value should be a reference to a non-empty hash of domain coordinates -> ";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}







#Program to use in the sequence alignments when rescueing domains (default ssearch36)
sub searchWith {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::Characterize::searchWith';
  my $default = "ssearch36";

  if ( $value ) {
    unless($value =~ /(blastp|ssearch36)/) {
      die "searchWith: this value can only be 'ssearch36' or 'blastp' -> ";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}






#Sequences of the refernece TCDB family
sub tcFamSeqs {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::Characterize::tcFamSeqs';
  my $default = "";

  if ( $value ) {
    unless(-f $value) {
      die "tcFamSeqs: file with family sequences not found or empty -> ";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}








#Full path to the blast DB that will be used to extract the sequence fragments
#corresponding to PFAM domain matches.
sub blastdb {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::Characterize::blastdb';
  my $default = "";

  if ( $value ) {
    unless(-f "${value}.phr" ) {
      die "BlastDB not detected -> ${value}.phr";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}



#This is the evalue that will be used to compare the sequence of domains against
#the full sequence of proteins in a given family
sub evalue {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::Characterize::evalue';
  my $default = 10;

  if ( $value ) {
    unless($value > 0) {
      die "Evalue should be >0.0 -> ";
    }
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}



#This is the main directory where results will be placed
sub rootDir {
  my ($self, $value) = @_;

  my $objPath = 'TCDB::Domain::Characterize::rootDir';
  my $default = ".";

  if ( $value ) {
    system "mkdir -p $value" unless (-d $value);
    $self->{$objPath} = $value;
  }

  unless ($self->{$objPath}) {
    $self->{$objPath} = $default;
  }

  return $self->{$objPath};
}




###########################################################################
##                       Function definitions
###########################################################################





sub rescueDomains {

  my ($self, $tcids) = @_;


  if ($self->treatAsSuperfamily) {
    $self->rescueDomForFamily($tcids->[0], $tcids);
  }
  else {
    my @fams = keys %{ $self->domFreq };

    foreach my $fam (@fams) {
      $self->rescueDomForFamily($fam, undef);
    }
  }
}



sub rescueDomForFamily  {

  my ($self, $fam, $fams) = @_;

  #-----------------------------------------------------------------
  #Relevant variables for this function

  #Create sequence directories
  my $seqDir   = $self->rootDir . "/$fam/sequences";
  system "mkdir -p $seqDir" unless (-d $seqDir);


  #-----------------------------------------------------------------
  #Get sequences of the reference family

  #If sequences were given as a parameter, copy them to the proper folder
  my $famSeqsFile = "$seqDir/family-${fam}.faa";
#  if (-f $self->protSeqs1it && !(-z $self->protSeqs1it)) {
#    my $cmd = "cp " . $self->protSeqs1it . " $seqDir";
#    system $cmd;
#  }

  #Extract the sequences of the full proteins in the family if they
  #were not provided by the user
  print "Rescue cycle 1 ($fam): Retrieving all family sequences\n";

  my $tcdbFaaParam = (-f $self->tcdbFaa)? "-d " . $self->tcdbFaa : "";
  if ($self->treatAsSuperfamily) {

    my @files = ();
    foreach my $f (@$fams) {
      my $file = "$seqDir/family-${f}.faa";
      system "extractFamily.pl -i $f -o $seqDir -f fasta $tcdbFaaParam" unless (-f $file);
      push(@files, $file);
    }

    my $str = join (" ", @files);
    system "cat $str > $seqDir/tmp.faa; mv $seqDir/tmp.faa $famSeqsFile";
    die "Coud not generate sequences for superfamily: " unless (-f $famSeqsFile && !(-z $famSeqsFile))
  }
  else {
    system "extractFamily.pl -i $fam -o $seqDir -f fasta $tcdbFaaParam" unless (-f $self->protSeqs1it);
  }


  #instatiate variables with sequence data
  $self->protSeqs1it($famSeqsFile) unless ($self->protSeqs1it);
  $self->domSeqs1it("$seqDir/domainFragments1it.faa");
  $self->idsFile2it("$seqDir/pids4SecondRescueCycle.faa");
  $self->protSeqs2it("$seqDir/seqs4SecondRescueCycle.faa");
  $self->domSeqs2it("$seqDir/domainFragments2it.faa");


  #-----------------------------------------------------------------
  #Create plots directory

  my $plotsDir = $self->rootDir . "/$fam/plots";
  system "mkdir -p $plotsDir" unless (-d $plotsDir);


  #-----------------------------------------------------------------
  #Create reports directory

  my $reportsDir = $self->rootDir . "/$fam/reports";
  system "mkdir -p $reportsDir" unless (-d $reportsDir);


  #-----------------------------------------------------------------
  #ReportFiles

  if ($self->searchWith eq 'ssearch36') {
    $self->aln1stCycle ("$reportsDir/alnDomProt_ssearch36_1stCycle.tsv");
    $self->aln2ndCycle ("$reportsDir/alnDomProt_ssearch36_2ndCycle.tsv");
  }
  elsif ($self->searchWith eq 'blastp') {
    $self->aln1stCycle ("$reportsDir/alnDomProt_blastp_1stCycle.tsv");
    $self->aln2ndCycle ("$reportsDir/alnDomProt_blastp_2ndCycle.tsv");
  }



  #-----------------------------------------------------------------
  #First write the report with the relevant direct domain hits

  print "Rescue cycle 1 ($fam): Extracting direct PFAM hits in the family.\n";
  $self->getRelevantDomains($fam);




  #-----------------------------------------------------------------
  #Now extract the sequence fragments corresponding to the PFAM
  #domains that are regarded as relevent for the definition of a family

  print "Rescue cycle 1 ($fam): Extracting sequences of PFAM domain hits in family.\n";

  my @domains = @{ $self->domFreq->{$fam} }[1 .. $#{ $self->domFreq->{$fam} }];
#  print Data::Dumper->Dump([$self->domFreq->{$fam} ], [qw( *domFreq)]);
#  exit;


  #Stop if no relevant domains were found
  unless (@domains) {
    my $errFile = $self->rootDir . "/$fam/NO_RESULTS.txt";
    my $msg = "Not enough proteins had characteristic domains.\n" .
      "* Check parameter minProtsDom of object TCDB::Domain::PfamParser\n";
    open (my $errh, ">", $errFile) || die $!;
    print $errh $msg;
    close $errh;
    return;
  }

  #Extract domain seqs
  unless (-f $self->domSeqs1it) {
    $self->extractDomainFragments1it($fam, \@domains);
  }

  #File with domain sequences must exist and not be empty
  unless (-f $self->domSeqs1it && !(-z $self->domSeqs1it)) {
    die "Missing or empty fiqle with matching PFAM domain sequences: ", $self->domSeqs1it, " -> ";
  }




  #-----------------------------------------------------------------
  #First rescue cycle: compare domain sequences with all the protein
  #sequences in the family

  print "Rescue cycle 1 ($fam): Comparing domain sequences with all proteins in family\n";

  unless (-f $self->aln1stCycle && !(-z $self->aln1stCycle)) {
    $self->alignDomainsVsProteins(1);
  }




  #-----------------------------------------------------------------
  #First rescue cycle: Parse alignment output and determine
  #how many proteins have been rescued and how many are still
  #without domain hits

  print "Rescue cycle 1 ($fam): Extracting domain hits per protein in family\n";

  my %protsWithDomains1 = ();
  $self->parseDomainAlignments(\%protsWithDomains1, 1);


#  print "Rescued Proteins: ", scalar keys %protsWithDomains1, "\n";
#  print Data::Dumper->Dump([\%protsWithDomains1 ], [qw(*protsWithDomains1 )]);
#  exit;



  #-----------------------------------------------------------------
  #Determine how many proteins have no domains detected after first
  #rescue iteraction in order to see if it is necessary running a
  #second domain rescue cycle.


  #First get all the IDs for the proteins in the analyzed family
  my @fullIDs = ();
  foreach my $arr (@{ $self->tcids->{$fam} }) {
    my $tc = $arr->[0];
    foreach my $acc (@{ $arr->[1] }) {
      push (@fullIDs, "${tc}-$acc");
    }
  }


  #Get the proteins still missing at least one domain
  my %missingHits = ();
  foreach my $tcid (@fullIDs) {

    my $missCnt = 0;
    foreach my $domArr (@{ $self->domFreq->{$fam} }[1 .. $#{ $self->domFreq->{$fam} }]) {

      unless (exists $protsWithDomains1{$tcid} && exists $protsWithDomains1{$tcid}{$domArr->[0]}) {
	$missCnt++;
  	$missingHits{$tcid} = $missCnt;
      }
    }
  }

#  print Data::Dumper->Dump([\%missingHits ], [qw(*missingHits )]);
#  exit;



  #-----------------------------------------------------------------
  #If there are missing domains, run a second doman rescue iteration

  my %protsWithDomains2 = ();
  if (%missingHits) {

    #get the sequences of the proteins with with at least one domain missing.
    print "\nRescue cycle 2 ($fam): Extracting sequences of proteins without domain hits.\n";
    $self->getProtSeqs4SecondRescueCycle(\%missingHits);


    #extract sequences of resqued domains
    print "Rescue cycle 2 ($fam): Extracting sequences rescued domains in cycle 1.\n";
    unless (-f $self->domSeqs2it) {
      $self->extractDomainFragments2it(\%protsWithDomains1);
    }


    #aligned rescued domains versus sequences without domain matches
    print "Rescue cycle 2 ($fam): Aligning rescued domain sequences to proteins without pfam hits.\n";
    $self->alignDomainsVsProteins(2);


    #Extract the coordinates of the resqued domains after the 2nd rescue cycle
    print "Rescue cycle 2 ($fam): Extracting coordenates of rescued domains.\n";
    $self->parseDomainAlignments(\%protsWithDomains2, 2);

  }




  #-----------------------------------------------------------------
  #Generate report of the recovered domains.

  my $reportName = "$reportsDir/${fam}_rescuedDomains.tsv";
  $self->writeReport($fam, $reportName, \@fullIDs, \%protsWithDomains1, \%protsWithDomains2);


  print "Reports generated for $fam!\n";

}




sub writeReport {

  my ($self, $fam, $outfile, $famIDs, $rescuedProts1, $rescuedProts2) = @_;


  my $totalIDs  = scalar @$famIDs;

#  print Data::Dumper->Dump([$totalIDs,$famIDs ], [qw(*famSize, *refIDs )]);
#  exit;

#  print Data::Dumper->Dump([ $rescuedProts1 ], [qw( *rescuedProts1)]);
#  exit;



  #-----------------------------------------------------------------
  #Get the number of direct hits per domain (e.g. Directly from PFAM)
  #and the rescue statistics.

#  print Data::Dumper->Dump([ $rescuedProts1 ], [qw( *rescuedProts1)]);
#  exit;


  my %stats = ();
  foreach my $domArr (@{ $self->domFreq->{$fam} }[1 .. $#{ $self->domFreq->{$fam} }]) {

    $stats{ $domArr->[0] }{directHits} = $domArr->[1];


    my %domCnt1 = ();  #Count first rescue cycle hits
    my %domCnt2 = ();  #Count second rescue cycle hits
    foreach my $prot (@{ $famIDs }) {
      if (exists $rescuedProts1->{$prot} && exists $rescuedProts1->{$prot}->{$domArr->[0]}) {
	$domCnt1{$prot} = 1;
      }
      elsif (exists $rescuedProts2->{$prot} && exists $rescuedProts2->{$prot}->{$domArr->[0]}) {
	$domCnt2{$prot} = 1;
      }
    }

    $stats{ $domArr->[0] }{cycle1rescue} = scalar keys %domCnt1;
    $stats{ $domArr->[0] }{cycle2rescue} = scalar keys %domCnt2;
  }

#  print Data::Dumper->Dump([\%stats ], [qw(*stats )]);
#  exit;



  #-----------------------------------------------------------------
  #Write report to file

  open (my $outh, ">", $outfile) || die $!;

  #Write global stats as comments
  foreach my $dom (keys %stats) {

    my $resc  = $stats{$dom}{cycle1rescue} + $stats{$dom}{cycle2rescue} - $stats{$dom}{directHits};
    my $total = sprintf("%.1f", ($stats{$dom}{cycle1rescue} + $stats{$dom}{cycle2rescue}) / $totalIDs * 100);
    print $outh "# $dom:  DirectHits:  $stats{$dom}{directHits}    Rescued Proteins:  $resc    ",
      "Prots with Domain in $fam:  ", $stats{$dom}{cycle1rescue} + $stats{$dom}{cycle2rescue},
      " (${total}% from a total of $totalIDs)\n";
  }


  my $cnt = 1;
  foreach my $prot (@$famIDs) {

    my $domStats = [];

    foreach my $domArr (@{ $self->domFreq->{$fam} }[1 .. $#{ $self->domFreq->{$fam} }]) {

      my $dom    = $domArr->[0];
      my $domPos = [];

      #Get the domain coords for this protein if rescued in the first cycle
      if (exists $rescuedProts1->{$prot} && exists $rescuedProts1->{$prot}->{$dom}) {

	my @matches = sort {$a->{left}<=>$b->{right}} @{ $rescuedProts1->{$prot}->{$dom} };

	foreach my $pos (@matches) {
	  push (@{ $domPos }, $pos->{left} . '-' . $pos->{right});
	}
      }


      #Get the domain coords for this protein if rescued in the second cycle
      elsif (exists $rescuedProts2->{$prot} && exists $rescuedProts2->{$prot}->{$dom}) {

	my @matches = sort {$a->{left}<=>$b->{right}} @{ $rescuedProts2->{$prot}->{$dom} };

	foreach my $pos (@matches) {
	  push (@{ $domPos }, $pos->{left} . '-' . $pos->{right});
	}
      }


      #format the rescued domain for printing
      if (@{ $domPos }) {
	if (exists $self->domCoords->{$fam} &&
	    exists $self->domCoords->{$fam}->{$prot} &&
	    exists $self->domCoords->{$fam}->{$prot}->{$dom}) {
	  push (@{ $domStats }, "${dom}|" . join("|", @{ $domPos }) . "|DirectHit");
	}
	elsif (exists $rescuedProts1->{$prot} && exists $rescuedProts1->{$prot}->{$dom}) {
	  push (@{ $domStats }, "${dom}|" . join("|", @{ $domPos }) . "|Rescued1");
	}
	elsif (exists $rescuedProts2->{$prot} && exists $rescuedProts2->{$prot}->{$dom}) {
	  push (@{ $domStats }, "${dom}|" . join("|", @{ $domPos }) . "|Rescued2");
	}
      }
      else {
	push (@{ $domStats }, "${dom}|Nohit");
      }
    }

    print $outh "$cnt\t$prot\t", join ("\t", @{ $domStats}), "\n";
    $cnt++;
  }

  close $outh;
}



#==========================================================================
#Extract the sequences of the proteins that did not have either direct
#hit or rescued PFAM domains


sub getProtSeqs4SecondRescueCycle {

  my ($self, $hr_prots) = @_;


  my $blastDB  = $self->blastdb;
  my $IDsFile  = $self->idsFile2it;
  my $outfile  = $self->protSeqs2it;


  #Save IDs in a file
  open (my $outh1, ">", $IDsFile) || die $!;
  print $outh1 join("\n", keys %{ $hr_prots }), "\n";
  close $outh1;


  #Extract the sequences here
  system qq(blastdbcmd -db $blastDB -dbtype prot -entry_batch $IDsFile -out $outfile);


  #Verify that results exist;
  my $eMsg = "Protein sequences file for second domain-rescue cycle not found file not found or empty:\n$outfile";
  die $eMsg unless (-f $outfile && !(-z $outfile));

}


#==========================================================================
#List the set of domains identified by PFAM including their postitions
#in the targer sequences.

sub getRelevantDomains {

  my ($self, $fam) = @_;

  #print Data::Dumper->Dump([ $self->domFreq->{$fam} ], [qw( *domFreq)]);
  #exit;


  my @domains = @{ $self->domFreq->{$fam} }[1 .. $#{ $self->domFreq->{$fam} }];

  #print Data::Dumper->Dump([$self->domCoords], ["*domains"]);
  #exit;

  
  my $outfile = $self->rootDir . "/$fam/reports/${fam}_relevantDomains.tsv";
  open (my $reph, ">", $outfile) || die $!;


 DOM:foreach my $domain (@domains) {

    my $dom = $domain->[0];  #The PFAM accession
    print $reph "===========================================================================\n";

  TCID:foreach my $tcid (sort_by_system keys %{ $self->domCoords->{$fam} }) {

      next TCID unless (exists $self->domCoords->{$fam}->{$tcid} && exists $self->domCoords->{$fam}->{$tcid}->{$dom});

      my ($tc, $acc) = split(/-/, $tcid);
      print $reph "$tcid";

    HIT:foreach my $hit (@{  $self->domCoords->{$fam}->{$tcid}->{$dom} }) {
	print $reph "\t${dom}:",$hit->{dlen}, ":",$hit->{dstart}, ":", $hit->{dend}, "|",
	  "${acc}:", $hit->{plen}, ":", $hit->{pstart}, ":", $hit->{pend};
      }
      print $reph "\n";
    }

    print $reph "\n";

  }
}






sub parseDomainAlignments {

  my ($self, $prots, $cycle) = @_;



  #-----------------------------------------------------------------
  #Parse sequence alignments

  my $infile = "";
  if    ($cycle == 1) { $infile =  $self->aln1stCycle; }
  elsif ($cycle == 2) { $infile =  $self->aln2ndCycle; }
  else { die "Unknown cycle: $cycle"; }



  my ($parsed, $qtcid, $qlen) = ({}, undef, undef);

  open (my $alnh, "<", $infile) || die $!;
  while (<$alnh>) {

    chomp;
    s/\s+$//;

    next if (!$_);


    #Get the length of the query for SSEARCH36 output (verify it if will work with BLASTP)
    if (/^#/) {
      if (/^#\s+Query:\s+(\S+)\s+-\s+(\d+)\s+aa/) {
	$qtcid = $1;
	$qlen  = $2;
      }
      next;
    }


    #Get line info
    my ($cachito, $protID, $id, $alnLen, $kk1, $kk2, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\s+/);

#    print "$_\n";
#    print Data::Dumper->Dump([$qtcid, $qlen, $cachito, $protID, $id, $alnLen, $qstart, $qend, $sstart, $send, $evalue, $bitscore ],
#			     [qw(*qtcid *qlen *cachito *protID *id *alnLen *qstart *qend *sstart *send *evalue *bitscore)]);
#    <STDIN>;


    next unless ($qtcid eq $cachito);


    my $domCoverage = $alnLen / $qlen;
    next if ($domCoverage < $self->domCovCutoff);


    #Get protein_ID with domain chachito and the cachito
    my ($protDomID, $domID) = split(/\|/, $cachito);


    my $res = [$protDomID, $sstart, $send, $domID, $qstart, $qend];
    push (@{ $parsed->{$protID} }, $res);

#    print Data::Dumper->Dump([$res ], [qw(*res)]);
#    <STDIN>;

  }
  close $alnh;

  if ($cycle == 1) {
#    print "\nCycle:$cycle\n$infile\nProteins with blastp|ssearch36 Domain Hits: ", scalar keys %{ $parsed }, "\n";
#    print Data::Dumper->Dump([$parsed ], [qw(*parsed )]);
#    exit;
  }



  #-----------------------------------------------------------------
  #Get the cachitos and assemble the domain hits.


  foreach my $protID (sort {$a cmp $b} keys %{ $parsed }) {

    my $fragments = $parsed->{$protID};

#    print "Processing: $protID\n";
    my $noOvp_match = 0;

    foreach my $hit (sort coords_by_left_pos @{ $fragments }) {

      my ($dom, $origCoords) = split(/:/, $hit->[3]);


#      print "\n\n===========================================================================\n";
#      print "$protID -> $dom\n", Data::Dumper->Dump([$prots ], [qw(*prots )]), "\n\n\n";
#      <STDIN>;



      if (exists $prots->{$protID} && exists $prots->{$protID}->{$dom}) {

	foreach my $coordSet (@{ $prots->{$protID}->{$dom} }) {

#	  print "-----------------------------------------------------------------\n";
#	  print Data::Dumper->Dump([$hit, $prots->{$protID}->{$dom} ], [qw(*hit *beforeCoords )]);

	  #Overlapping match
	  if (($hit->[1] >= $coordSet->{left} && $hit->[1] <= $coordSet->{right}) &&
	      ($hit->[2]  > $coordSet->{right})) {

	    #Expand the overlap of the protein
	    $coordSet->{right} = $hit->[2];
	  }

	  #Non overlapping match
	  elsif ($hit->[1] > $coordSet->{right}) {

	    #print Data::Dumper->Dump([ $noOvp_match], [qw( *noOvp_match)]);

	    if ( $noOvp_match ) {

	      $self->removeRedundantMatches($hit, $prots->{$protID}->{$dom});
	    }
	    else {
	      push (@{ $prots->{$protID}->{$dom} }, { left   => $hit->[1],
						      right  => $hit->[2] });
	      $noOvp_match = 1;
	    }
	  }

#	  print Data::Dumper->Dump([$prots->{$protID}->{$dom} ], [qw( *afterCoords )]);
#	  print "-----------------------------------------------------------------\n\n\n";
#	  <STDIN>;
	}
      }

      #First time domain is read
      else {
	push (@{ $prots->{$protID}->{$dom} }, { left   => $hit->[1],
						right  => $hit->[2] });
      }

    }
  }

#  if ($cycle == 1) {
#    print Data::Dumper->Dump([$prots ], [qw(*cachitosArmados )]);
#    exit;
#  }

}


#fix overlapsoverlaps
sub removeRedundantMatches {

  my ($self, $hit, $matches) = @_;

#  print "*** New non-overlapping match! ***\n";
#  print Data::Dumper->Dump([$hit, $matches ], [qw(*hit *domainsBefore )]), "\n\n";

  my $newMatch = 0;
  foreach my $seg (@{ $matches }) {

    #Non overlapping match
    if ($hit->[1] > $seg->{right}) {
      $newMatch = 1;
    }

    #Overlapping to the right match
    if (($hit->[1] >= $seg->{left} && $hit->[1] <= $seg->{right}) &&
	($hit->[2]  > $seg->{right})) {

      #Expand the overlap of the protein
      $seg->{right} = $hit->[2];
      $newMatch = 0;
    }

    #hit is contained within a previous match
    elsif (($hit->[1] >= $seg->{left} && $hit->[1] <= $seg->{right}) &&
	   ($hit->[2] <= $seg->{right})) {
      $newMatch = 0;
    }

    elsif (($hit->[1] < $seg->{left} && $hit->[2] >= $seg->{left}) &&
	   ($hit->[2] <= $seg->{right})) {
      die "Unexpected overlap to the left!... debug!\n", print Data::Dumper->Dump([$hit, $matches ], [qw(*hit *domains )]);;
    }

  }

  if ($newMatch) {
    push (@{ $matches }, { left  => $hit->[1], right => $hit->[2] });
  }

#  print Data::Dumper->Dump([$matches ], [qw(*domainsAfter )]);
#  print "*****************************************************************\n\n\n";
}




sub coords_by_left_pos {

  if ($a->[1] == $b->[1]) {
    $a->[2] <=> $b->[2];
  }
  else {
    $a->[1] <=> $b->[1];
  }
}





#==========================================================================
#Run blastp or ssearch36 of the sequence fragments that aligned with
#known domains against the full sequences of all proteins in the
#family


sub alignDomainsVsProteins {

  my ($self, $cycle) = @_;

  my $eval = $self->evalue;

  my $domSeqs = "";
  my $famSeqs = "";
  my $outfile = "";

  if ($cycle == 1) {
    $famSeqs = $self->protSeqs1it;
    $domSeqs = $self->domSeqs1it;
    $outfile = $self->aln1stCycle;
  }
  elsif ($cycle == 2) {
    $famSeqs = $self->protSeqs2it;
    $domSeqs = $self->domSeqs2it;
    $outfile = $self->aln2ndCycle;
  }
  else {
    die "Unknown domain rescue cycle: $cycle";
  }

#  print "Align Cycle: $cycle\n  $domSeqs\n  $famSeqs\n  $outfile\n";

  unless (-f $outfile && !(-z $outfile)) {
    if ($self->searchWith eq 'ssearch36') {

      #Defaults: -s BL50 -z 2 (but l1-l6 also works) ... sensitive, but error prone
      my $args = qq(-z 21 -k 1000 -E $eval -m 8C -s BL50 $domSeqs $famSeqs > $outfile);
      system "ssearch36 $args";
    }
    elsif ($self->searchWith eq 'blastp') {
      my $args = qq(-evalue $eval -use_sw_tback -max_hsps 4 -comp_based_stats 2 -outfmt 7 -query $domSeqs -subject $famSeqs > $outfile);
      system "blastp $args";
    }
  }
}



#==========================================================================
#Extract the sequence from the TCDB proteins that had no domain matches
#after the first PFAM domain rescue cycle.


sub extractDomainFragments2it {

  my ($self, $domAfter1it) = @_;

#  print Data::Dumper->Dump([$domAfter1it ], [qw(*domainHits )]);
#  exit;

  my $blastDB = $self->blastdb;

  open (my $outh, ">>", $self->domSeqs2it) || die $!;

  #For each tcid
  foreach my $tcid (keys %{ $domAfter1it }) {

    #For each PFAM ID
    foreach my $domID (keys %{ $domAfter1it->{$tcid} }) {

      #For each set of coordinates matching the PFAM domain
      foreach my $hit (@{ $domAfter1it->{$tcid}->{$domID} }) {

	my $pstart = $hit->{left};
	my $pend   = $hit->{right};

	my $args = qq(-db $blastDB -dbtype prot -entry $tcid -range ${pstart}-${pend});

	my $seq = qx(blastdbcmd $args);
	$seq =~ s/$tcid/$tcid:${pstart}-${pend}\|$domID/;

	print $outh $seq;
      } #hit
    } #domain
  } #tcid

  close $outh;
}




#==========================================================================
#Extract the sequence from the TCDB proteins that correspond to the
#Pfam or CDD hits, this is for the first Domain rescue cycle.


sub extractDomainFragments1it {

  my ($self, $fam, $relDomains) = @_;


#  print Data::Dumper->Dump([$self->domCoords->{$fam} ], [qw(*domCoords )]);
#  exit;


  my $blastDB = $self->blastdb;

  open (my $outh, ">>", $self->domSeqs1it) || die $!;

  foreach my $dom (@{ $relDomains }) {

    my $domID = $dom->[0];
    my ($did, $dv) = split(/\./, $domID);

    #now get all the proteins that have a match with this domain
    foreach my $prot (keys %{ $self->domCoords->{$fam} }) {

      if (exists $self->domCoords->{$fam}->{$prot}->{$domID}) {

	#Extract sequences for each individual domain hit
	foreach my $hit (@{ $self->domCoords->{$fam}->{$prot}->{$domID} }) {

	  my $pstart = $hit->{pstart};
	  my $pend   = $hit->{pend};
	  my $dstart = $hit->{dstart};
	  my $dend   = $hit->{dend};

	  my $args = qq(-db $blastDB -dbtype prot -entry $prot -range ${pstart}-${pend});
	  my $seq = qx(blastdbcmd $args);
	  $seq =~ s/$prot/${prot}:${pstart}-${pend}\|${did}:${dstart}-$dend/;

	  print $outh $seq;
#	  print "$domID\n", Data::Dumper->Dump([$args, $seq], [qw(*args *seq )]);
#	  <STDIN>;
	}
      }
    }
  }

  close $outh;
}



1;


