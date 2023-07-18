#!/usr/bin/env perl -w

use strict;
use warnings;

use Data::Dumper;

use TCDB::Assorted;
use Getopt::Long;


#==========================================================================
# This script gets the gblast hits with multicomponent systems and
# performs a genomic context to see if there are systems for which
# all or most of the hits are in the neighborhood.
#
#--------------------------------------------------------------------------
#
# Written by: Arturo Medrano
# Date: 8/21/2017
#
#==========================================================================


my $gnmAcc     = undef;
my $gnmDir     = undef;
my $tabFile    = undef;
my $gblastFile = undef;
my $outdir     = "genome_context";
my $repForm    = 'circular';
my $evalue     = 1e-3;
my $geneDis    = 15;
my $outFmt     = "tsv";

my $domain = "Bacteria";


read_command_line_arguments();
#print Data::Dumper->Dump([$gnmAcc, $gnmDir, $tabFile, $gblastFile, $outdir, $evalue],
#			 [qw(*gnmAcc *gnmDir *tabFile *gblastFile *outdir $evalue)]);
#exit;



#==========================================================================
#Get the TCIDs of multi-component systems in TCDB.

#Download the TCDB database in fasta format
my $tcdbFaa = "$outdir/tcdb/tcdb.faa";
system "extractTCDB.pl -i tcdb -o $outdir/tcdb -f fasta" unless (-f $tcdbFaa);
die "TCDB seqs not found: $tcdbFaa" unless (-f $tcdbFaa);

my $multSystems = getModeSystems($tcdbFaa, 'multi');

#print Data::Dumper->Dump([$multSystems ], [qw(*multSystems )]);
#exit;



#==========================================================================
#Get the multicomponent systems in the genome


#first get all TCIDs of multicomponent systems found by GBLAST. Note
#Don't filter by threshold to capture everything that GBLAST found.
my %gblastTCIDs = ();
getGBLASTmultiComponentTCIDs($gblastFile, 10, $multSystems, \%gblastTCIDs);

#print Data::Dumper->Dump([\%gblastTCIDs ], [qw(*gblastTCIDs)]);
#exit;




#==========================================================================
#Run blastp on the TCDB proteins that belong to the multi-component
#systems in GBlast against the entire reference genome


#First send to a file the TCIDs of the multi-component systems found in
#the genome
my $gnm_multSystems_tcdis_file = "$outdir/gnm_multSystems_tcids.txt";

unless (-f $gnm_multSystems_tcdis_file && !(-z $gnm_multSystems_tcdis_file)) {
  open(my $fh, ">", $gnm_multSystems_tcdis_file) || die $!;
  print $fh join("\n", keys %gblastTCIDs), "\n";
  close $fh;
}


#Now select top blastp hits and generate hydropathy plots
my $gnmProteomeFile = "$gnmDir/${gnmAcc}_protein.faa.gz";
die "Proteome file not found: $gnmProteomeFile -> " unless (-f $gnmProteomeFile);

my $whatOutDir = "$outdir/multSystemsHydropathy";
my $blastdb = "proteome";
my $args = "-g $gnmProteomeFile -q $gnm_multSystems_tcdis_file -d $blastdb -o $whatOutDir";
system "selectTopBlastpHitsAndPlotHydropathy $args" unless (-f "$whatOutDir/plot_index.html");



#==========================================================================
#For each multicompnent system hit by GBLAST determine if there are
#identifiable missing components in the neighborhood

my $blastFilteredFile = "$whatOutDir/blast_output/top_blast_hits.csv";
die "Filtered blast output not found: $blastFilteredFile -> " unless (-f $blastFilteredFile);

my $blastRawOutput    = "$whatOutDir/blast_output/blastOutput.csv";
die "Filtered blast output not found: $blastRawOutput -> " unless (-f $blastRawOutput);

checkGenomeContext($blastFilteredFile, $blastRawOutput, \%gblastTCIDs, $multSystems);
print "Done\n";





sub checkGenomeContext {

  my ($filteredBlastFile, $rawBlastFile, $gblastMCS, $tcdbMCS) = @_;


  #-----------------------------------------------------------------
  #Read BlastP data

  my %rawBlast =();
  readRawBlast($rawBlastFile, \%rawBlast, $gblastMCS);

  my %filteredBlast = ();
  readFilteredBlast($filteredBlastFile, \%filteredBlast, $gblastMCS);

#  print Data::Dumper->Dump([\%filteredBlast ], [qw(*filteredBlast )]);
#  exit;


  #-----------------------------------------------------------------
  #If a missing component does not have a filtered blast hit,
  #search the raw blast output and add as much blast hits as possible.

  complementFilteredBlastHits(\%filteredBlast, \%rawBlast, $gblastMCS, $tcdbMCS);

#  print Data::Dumper->Dump([\%filteredBlast ], [qw(*filteredBlast )]);
#  exit;



  #-----------------------------------------------------------------
  #Now check the genomic context of all the multicomponent systems
  #and identify those components that share genomic neighborhood


  verify_neighborhood(\%filteredBlast, $tabFile);

#  print Data::Dumper->Dump([\%filteredBlast ], [qw(*filteredBlast )]);
#  exit;



  #-----------------------------------------------------------------
  #Generate final report

  generate_report(\%filteredBlast);
}





sub generate_report {

  my $systems = shift;

  my $outfile = "$outdir/reportMulticomponentSystems.$outFmt";
  my $sep = "\t";

  if ($outFmt eq 'csv') {
    $sep = ",";
  }

  open (my $outh, ">", $outfile) || die $!;

  print $outh "tcid${sep}query_accession${sep}status${sep}subject_accession${sep}evalue${sep}Perc_idenity${sep}query_length${sep}subject_length${sep}alignment_length${sep}query_coverage${sep}subject_coverage${sep}neighbors\n";

  foreach my $tcid (sort_by_system keys %{ $systems }) {
    foreach my $tcAcc (sort_by_status1 ($systems->{$tcid})) {
      foreach my $hit (sort by_status2 @{ $systems->{$tcid}->{$tcAcc} }) {

	my $status    = $hit->{status};
	my $sacc      = "";
	my $eval      = "";
	my $ident     = "";
	my $qlen      = "";
	my $slen      = "";
	my $aln       = "";
	my $qcov      = "";
	my $scov      = "";
	my $neighbors = "";


	if ($status =~ /(GBLAST|Candidate|RawBlast)/) {
	  $sacc      = $hit->{sacc};
	  $eval      = $hit->{eval};
	  $ident     = $hit->{ident};
	  $qlen      = $hit->{qlen};
	  $slen      = $hit->{slen};
	  $aln       = $hit->{aln};
	  $qcov      = $hit->{qcov};
	  $scov      = sprintf("%.1f", $hit->{scov});
	  $neighbors = $hit->{neighbors};
	}


	print $outh "$tcid${sep}$tcAcc${sep}$status${sep}$sacc${sep}$eval${sep}$ident${sep}$qlen${sep}$slen${sep}$aln${sep}$qcov${sep}$scov${sep}$neighbors\n";

#	print Data::Dumper->Dump([$tcid, $tcAcc, $hit ], [qw(*tcid *tcAcc *hit )]);
      }
    }
  }
  close $outh;
}


#Sort TC accessions by  the number of hits in GBLAST, Candidate, RawBlast and NoHit

sub sort_by_status1 {

  my $compHash = shift;

#  print Data::Dumper->Dump([ $compHash], [qw( *compHash)]);
#  exit;


  my @counts = ();

  foreach my $tcAcc (keys %{ $compHash }) {

    my $gbCnt = 0;
    my $caCnt = 0;
    my $raCnt = 0;
    my $noCnt = 0;

    foreach my $hit (@{ $compHash->{$tcAcc} }) {
      if    ($hit->{status} eq 'GBLAST')    { $gbCnt++; }
      elsif ($hit->{status} eq 'Candidate') { $caCnt++; }
      elsif ($hit->{status} eq 'RawBlast')  { $raCnt++; }
      else { $noCnt++ }
    }

    push (@counts, [$tcAcc, {gblast=>$gbCnt, cand=>$caCnt, raw=>$raCnt, nohit=>$noCnt}]);
  }

  my @out = map {$_->[0]} sort by_best_hits @counts;

  return @out;

}


sub by_best_hits {

  if ($a->[1]->{gblast} > 0 || $b->[1]->{gblast} > 0) {
    $b->[1]->{gblast} <=> $a->[1]->{gblast}; #Sort by number of GBLAST hits;
  }
  elsif ($a->[1]->{cand} > 0 || $b->[1]->{cand} >0) {
    $b->[1]->{cand} <=> $a->[1]->{cand};  #Sort by number of candidate hits
  }
  elsif ($a->[1]->{raw} > 0 ||  $b->[1]->{raw} > 0) {
    $b->[1]->{raw} <=> $a->[1]->{raw};  #Sort by number of raw blast
  }
  else {
    $a->[1]->{nohit} <=> $b->[1]->{nohit}; #sort by fewer number of NoHits
  }
}



#Sort 2 blastphits by status and quality

sub by_status2 {

  my %code = (GBLAST=>0, Candidate=>1, RawBlast=>2, NoHit=>3);

  if ($a->{status} eq $b->{status}) {

    my $qcov1 = $a->{qcov};
    my $scov1 = $a->{scov};
    my $aln1  = $a->{aln};
    my $eval1 = $a->{eval};

    my $qcov2 = $b->{qcov};
    my $scov2 = $b->{scov};
    my $aln2  = $b->{aln};
    my $eval2 = $b->{eval};

    my $maxCov1 = $qcov1 + $scov1;
    my $maxCov2 = $qcov2 + $scov2;

    if ($maxCov1 == $maxCov2) {

      if ($eval1 == $eval2) {
	$aln2 <=> $aln1;  #Sort by alignment length
      }
      else {
	$eval1 <=> $eval2; #sort by $evalue
      }
    } #aln
    else {
      $maxCov2 <=> $maxCov1;  #sort by max Coverage
    }
  }
  else {
    $code{ $a->{status} } <=> $code{ $b->{status} };
  }


}



#==========================================================================
#Open the feature table from genome directory to get the order of genes
#in the genome.



sub verify_neighborhood {

  my ($systems, $featuresFile) = @_;


  #-----------------------------------------------------------------
  #Parse genome feature table

  open (my $fh1, "-|",  "zcat $featuresFile") || die $!;

  my @cds = ();
  my %pos2cds = ();
  my %cds2pos = ();
  while (<$fh1>) {
    chomp;

    next unless (/^CDS/);

    #Columns in file:
    #0)  feature
    #1)  class
    #2)  assembly
    #3)  assembly_unit
    #4)  seq_type
    #5)  chromosome
    #6)  genomic_accession
    #7)  start
    #8)  end
    #9) strand
    #10) product_accession
    #11) non-redundant_refseq
    #12) related_accession
    #13) name
    #14) symbol
    #15) GeneID
    #16) locus_tag
    #17) feature_interval_length
    #18) product_length
    #19) attributes

    my @d = split (/\t/, $_);
    my ($acc, $ver) = split(/\./, $d[10]);

    my $line = { acc=>$acc, start=>$d[7], end=>$d[8], strand=>$d[9], name=>$d[13] };

    push (@cds, $line);
  }
  close $fh1;



  #-----------------------------------------------------------------
  #Get the ranked position of each CDS

  my $pos = 1;
  foreach my $orf (sort { $a->{start} <=> $b->{start} } @cds) {
    $pos2cds{ $pos }         = $orf;
    $cds2pos{ $orf->{acc} } = $pos;
    $pos++;
  }

#  print Data::Dumper->Dump([\%cds2pos], [qw(*cds2pos )]);
#  exit;



  #-----------------------------------------------------------------
  #For each multicomponent system in the genome see which CDS
  #hitting different compnents are in the neighborhood


  #For calculating circular distances
  my $lastCDSpos = scalar keys %cds2pos;


#  print Data::Dumper->Dump([$systems ], [qw(*filteredSystems )]);
#  exit;

 TCID:foreach my $tcid (keys %{ $systems }) {

    my %neighbors = ();


    #For this system, get the accessions of the CDS matching each component
    #and their position rank in the genome
    my %candComp = ();
    getAccessionsForSystem($tcid, $systems, \%cds2pos, \%candComp);

#    print Data::Dumper->Dump([$tcid, \%candComp ], [qw(*tcid *matches)]);
#    exit;



    #Using each candidate component as references get how many of the
    #other proteins matching components in the same systems are in
    #the neigborhood.
    my @accs = sort {$candComp{$a}->[0] <=> $candComp{$b}->[0] } keys %candComp;
  ACC1:foreach my $acc1 (@accs) {

      my @nearby = ($acc1);

#      print "\n\n===========================================================================\n";

    ACC2:foreach my $acc2 (@accs) {

	next ACC2 if ($acc1 eq $acc2);



#	print "$acc1 vs $acc2\n";

	#determine if the genes are neighbors
	my @dist = areGenesCloseby($candComp{$acc1}, $candComp{$acc2}, $lastCDSpos, $geneDis, $repForm);

	if (@dist) {
	  my $cmp =  $candComp{$acc2}->[1];
	  push (@nearby, "${acc2}(". $dist[2] ."):$cmp");
	}

#	print Data::Dumper->Dump([$acc1, $acc2, \@dist ], [qw(*acc1 *acc2 *dist )]), "\n";
      }

#      print "\n------------------------------\n";
#      print Data::Dumper->Dump([$acc1, \@nearby], [qw(*acc *neighbors )]);
#      <STDIN>;



      if (scalar @nearby > 1) {
	$neighbors{$acc1} = join("|", sort {$a cmp $b} @nearby);
      }
      else {
	$neighbors{$acc1} = "None";
      }

#      print Data::Dumper->Dump([\%neighbors ], [qw( *neighbors )]);
#      <STDIN>;
    }



    #-----------------------------------------------------------------
    #With the neighbors fully estimated for this system, now add the
    #neighbors information to the filteredmulticomponent systems hits.

    foreach my $tcAcc (keys %{ $systems->{$tcid} }) {
      foreach my $hit (@{ $systems->{$tcid}->{$tcAcc} }) {
	if (exists $hit->{sacc} && exists $neighbors{ $hit->{sacc} }) {
	  $hit->{neighbors} = $neighbors{ $hit->{sacc} };
	}
      }
    }

#    print "$tcid\n", Data::Dumper->Dump([ $systems->{$tcid} ], [qw( *Filtered )]);
#    <STDIN>;
  }
}


#==========================================================================
#determine if two genes are neighbors and take into account the
#linearity or circularity of the genome under analysis.


sub areGenesCloseby {
  my ($pos1, $pos2, $repSize, $refGeneDist, $repStucture) = @_;


  my @pos = sort {$a <=> $b} ($pos1->[0], $pos2->[0]);
  my $dist = undef;


  #Replicon is circular
  if ($repStucture eq 'circular') {

    my $d1 = $pos[1] - $pos[0];
    my $d2 = $pos[0] + ($repSize - $pos[1]);

    $dist = (sort {$a <=> $b} ($d1, $d2))[0];

  }


  #Replicon is linear
  else {
    $dist = $pos[1] - $pos[0];
  }


#  print Data::Dumper->Dump([$pos1, $pos2, $dist,  $repSize], [qw(*pos1 *pos2 *dist *repSize)]);
#  <STDIN>;


  #Test if distance in in the acceptable range
  if ($dist <= $refGeneDist) {
    return ($pos1, $pos2, $dist);
  }
  else {
    return ();
  }

}




#==========================================================================
# Given at TCID of a multicomponent system, get the genome proteins that
# match all the components

sub getAccessionsForSystem {

  my ($tcid, $multSystems, $acc2pos, $out) = @_;

 COMP:foreach my $tcAcc (keys %{ $multSystems->{$tcid} }) {
  HIT:foreach my $hit (@{ $multSystems->{$tcid}->{$tcAcc} }) {

      next HIT if ($hit->{status} eq 'NoHit');

      my $gnmAcc = $hit->{sacc};
      if (exists $acc2pos->{ $gnmAcc }) {
	$out->{ $gnmAcc } =  [$acc2pos->{ $gnmAcc }, $tcAcc];
      }
      else {
	die "No position found for: $gnmAcc -> ";
      }
    }
  }
}



#==========================================================================
#Given the GBLAST and BlastP searches put together the components that
#had blastp matches and assign categories: gblast, rawblast and candidates


sub complementFilteredBlastHits {

  my ($filtered, $raw, $gblastMCS, $tcdbMCS) = @_;


 TCID:foreach my $tcid (keys %{ $gblastMCS }) {

    #Total tcdb components for system
    my $tcdbComp = $tcdbMCS->{$tcid};

  ACC:foreach my $acc (@{ $tcdbComp }) {

      if (exists $filtered->{$tcid} && exists $filtered->{$tcid}->{$acc}) {
	my @sortedHits = sort by_quality @{ $filtered->{$tcid}->{$acc} };
	$filtered->{$tcid}->{$acc} = \@sortedHits;
	next ACC;
      }

      if (exists $raw->{$tcid} && exists $raw->{$tcid}->{$acc}) {

	#sort raw hits by quality
	my @sortedHits = sort by_quality @{ $raw->{$tcid}->{$acc} };

	$filtered->{$tcid}->{$acc} = \@sortedHits;
      }
      else {
	$filtered->{$tcid}->{$acc} = [{status=>'NoHit'}];
      }
    }

#    print Data::Dumper->Dump([$tcid, $filtered->{$tcid}], [qw(*tcid *complemented)]);
#    <STDIN>;

  }

}



sub by_quality {

  my $qcov1 = $a->{qcov};
  my $scov1 = $a->{scov};
  my $aln1  = $a->{aln};
  my $eval1 = $a->{eval};

  my $qcov2 = $b->{qcov};
  my $scov2 = $b->{scov};
  my $aln2  = $b->{aln};
  my $eval2 = $b->{eval};

  my $maxCov1 = (sort {$b <=> $a} ($qcov1, $scov1))[0];
  my $maxCov2 = (sort {$b <=> $a} ($qcov2, $scov2))[0];

  if ($maxCov1 == $maxCov2) {

    if ($eval1 == $eval2) {
      $aln2 <=> $aln1;  #Sort by alignment length
    }
    else {
      $eval1 <=> $eval2; #sort by $evalue
    }
  } #aln
  else {
    $maxCov2 <=> $maxCov1;  #sort by max Coverage
  }
}


#==========================================================================
#Columns in filered blast Output:
# TCID qacc sacc qlen slen eval ident aln qcov scov

sub readFilteredBlast {
  my ($infile, $out, $gblast) = @_;

  open (my $fh, "<", $infile) || die $!;
  while (<$fh>) {
    chomp;
    s/ +$//;

    next if (/^TCID/ || !$_);
    my ($tcid, $qacc, $sid, $qlen, $slen, $eval, $ident, $alen, $qcov, $scov) = split (/,/);

    my ($sacc, $ver) = split(/\./, $sid);

#    print Data::Dumper->Dump([$tcid, $qacc, $sacc, $qlen, $slen, $eval, $ident, $alen, $qcov, $scov],
#			     [qw(*tcid *qacc *sacc *qlen *slen *eval *ident *alen *qcov *scov)]);
#    <STDIN>;


    my $status = undef;
    if (exists $gblast->{$tcid} && exists $gblast->{$tcid}->{$qacc} && exists $gblast->{$tcid}->{$qacc}->{$sacc}) {
      $status = "GBLAST";
    }
    else {
      $status = "Candidate";
    }

    push(@{ $out->{$tcid}->{$qacc} }, { sacc=>$sacc,
					qlen=>$qlen,
					slen=>$slen,
					eval=>$eval,
					ident=>$ident,
					aln=>$alen,
					qcov=>$qcov,
					scov=>$scov,
				        status=>$status});
  }
}




#==========================================================================
#Columns in raw blast Output:
#   qseqid sseqid qlen slen evalue pident length qcovs qseq sseq

sub readRawBlast {
  my ($infile, $out, $gblast) = @_;

  open (my $fh, "<", $infile) || die $!;
  while (<$fh>) {

    chomp;
    s/ +$//;

    next if (/^#/ || !$_);
    my ($qid, $sid, $qlen, $slen, $eval, $ident, $alen, $qcov, @kk) = split (/,/);

    my ($sacc, $kk1) = split (/\./, $sid);
    my $scov = $alen / $slen * 100;

#    print Data::Dumper->Dump([$qid, $sacc, $qlen, $slen, $eval, $ident, $alen, $qcov, $scov],
#			     [qw(*qid *acc *qlen *slen, *eval *ident *alen *qcov *scov)]);
#    <STDIN>;

    unless ($qid && $sacc && $qlen && $slen && $eval >= 0 && $ident && $alen && $qcov >= 0) {
      die "Coud not parse raw blast output line: $_ -> ";
    }

    #Ignore hits with not enough coverage
    next if ($qcov < 35 && $scov < 35);

    my ($tcid, $qrs) = split(/-/, $qid);
    my ($qacc, $kk2) = split(/\./, $qrs);

    my $status = undef;
    if (exists $gblast->{$tcid} && exists $gblast->{$tcid}->{$qacc} && exists $gblast->{$tcid}->{$qacc}->{$sacc}) {
      $status = "GBLAST";
    }
    else {
      $status = "RawBlast";
    }

    push(@{ $out->{$tcid}->{$qacc} }, { sacc=>$sacc,
					qlen=>$qlen,
					slen=>$slen,
					eval=>$eval,
					ident=>$ident,
					aln=>$alen,
					qcov=>$qcov,
					scov=>$scov,
					status=>$status});

  }
  close $fh;


}







###########################################################################
##
##                     Subroutines definitions
##
###########################################################################



sub read_command_line_arguments {


  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }


  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "gdir|genome-dir=s"       => \&read_genome_dir,
      "gacc|genome-accession=s" => \$gnmAcc,
      "gblast|gblast-tsv=s"     => \&read_gblast_file,
      "e|evalue=f"              => \$evalue,
      "o|outdir=s"              => \$outdir,
      "h|help"                  => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);

  die "Genome accession (option -gacc) is mandatory -> " unless ($gnmAcc);


  $gnmDir = "/ResearchData/Users/amedrano/CPR_Genomes/$domain/$gnmAcc" unless ($gnmDir);
  die "Genome directory was not found: $gnmDir" unless (-d $gnmDir);


  $tabFile = "$gnmDir/${gnmAcc}_feature_table.txt.gz";
  die "Gene feature file  not found: $tabFile" unless (-f $tabFile);


  die "GBLAST output file (with extension .tsv) is mandatory" unless ($gblastFile);


  die "Unknown replicon structure: $repForm -> " unless ($repForm =~ /(circular|linear)/);

  system "mkdir -p $outdir" unless (-d $outdir);
}




#==========================================================================
#Read the -gblast option

sub read_gblast_file {

  my ($opt, $value) = @_;

  unless (-f $value) {
    die "GBlast file was not found: $value -> ";
  }
  $gblastFile = $value;
}


#==========================================================================
#Read the -gd option

sub read_genome_dir {

  my ($opt, $value) = @_;

  unless (-d $value) {
    die "Directory with NCBI genome data does not exist -> ";
  }

  $gnmDir = $value;
}






