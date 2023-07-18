#!/usr/bin/env perl -w

use warnings;
use strict;
use Data::Dumper;

use Getopt::Long;

use TCDB::CheckDependencies;


#==========================================================================
#  Given a list of novel transporters or IMPs, get the homologs from NCBI
#  that meet the following requirements:
#
#  1. At least up to 30 homologs in other genomes, preferably across domains.
#     From this list the user will select at least 5.
#  2. Select homologs with E-value > 1e-15
#  3. Include remote homologs, as long as they share the same domain.
#  4. Generate a list of at least 5 homologs per query protein
#     for upload to TCDB.
#
#  The strategy will consist of extracting 30 candidates from which
#  five will be selected manually based on merit (i.e. TMS topology, alignment,
#  hydropathy overlap, and Pfam domain content).
#
#--------------------------------------------------------------------------
#  Written by: Arturo Medrano
#  Date:       March 2019
#==========================================================================



#==========================================================================
#Check dependencies

my $blastdbcmd = "/usr/local/bin/blastdbcmd";


my @dependencies = ('grep', 'hmmtop', 'blastdbcmd',  'psiblast', 'famXpander.pl', 'cd-hit', 'extractTCDB.pl');

my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;


#==========================================================================
#read command line

#Input protein accessions
my $input      = undef;
my @accessions = ();


#Output directory where the data related to this analysis will be stored
my $workDir = "newIMPs4tcdb";


#For blastDB acccess
my $tcblastdbDir = "$ENV{HOME}/db/blastdb";

#if the genome is not in NCBI yet, provide this manually
my $gnmBlastDB   = 'nr';

#Output file with the output of BLASTing query proteins against TCDB.
my $tcblastOutputFile = undef;


#Output reports
my $htmlReport = undef;
my $tsvReport  = undef;


#Blast parameters.
my $evalue = 1e-3;             #For parsing NCBI blast.
my $tceval = 1e-18;            #For parsing TCBLAST hits.
my $redundancyEvalue = 1e-15;  #For removing redundant sequences
my $coverage  = 0.5;           #Coverage for blast searches against NCBI
my $tcoverage = 0.2;           #Subject (target) coverage when blasting against TCDB
my $maxSlen   = undef;         #maximum length of subject relative to query (defined during commandline parsing).
my $blastMode = "-p T";        #Run famXpander remotely by default
my $owrite_tcblast = 0;        #Control whether or not the TCDB blast DB will be regenerated

my $membersNum = 50;           #maximum number of candidates to generate for manual inspection.
my $ncbiHom    = 25;           #Minimum number of non-redudant homologs from NCBI to analyze
my $minHoms    = 3;            #Number of homologus to create a new TCDB family

#Minimum number of TMS the homologs must have
my $minTMS = 4;


#How will the results be sorted (asc, desc, no)
my $sort = 'desc';

#Store here the blast results of the query vs the NCBI NR database
my %ncbiBlast = ();
my %ncbiMatches = ();

read_command_line_arguments();

#print Data::Dumper->Dump([$ncbiHom, $minHoms  ], [qw(*ncbiHom *minHoms )]);
#exit;

#print Data::Dumper->Dump([$input, \@accessions, $workDir, $tcblastdbDir,
#			  $htmlReport, $tsvReport, $evalue, $redundancyEvalue, $minTMS, $sort],
#			 [qw(*input *accessions *workDir *tcblastdbDir
#			  *htmlReport *tsvReport *evalue *redundantEvalue *minTMS *sort)]);
#exit;


ACC:foreach my $acc (@accessions) {


#  next ACC unless ($acc =~ /Dvar_30850/);

  my $accDir = "$workDir/$acc";
  system "mkdir -p $accDir" unless (-d $accDir);

  open (my $logh, ">", "$accDir/log.txt") || die $!;


  print "===========================================================================\n";
  print "Working with: $acc\n\n";


  print $logh "===========================================================================\n";
  print $logh "Working with: $acc\n\n";


  my $retok = process_accession($acc, $logh);

  unless ($retok) {
    print $logh " *** Problem when processing: $acc\n\n";
    print " *** Problem when processing: $acc\n\n";
  }

  print $logh "Status: $retok\n";
  print $logh "\n\n";

  close $logh;

}



###########################################################################
##                    Subroutine Defintions                              ##
###########################################################################



sub process_accession {

  my ($acc, $logh) = @_;

  my $accDir = "$workDir/$acc";


  #------------------------------------------------------------
  #Run blast and check for number of NCBI homologs

  print "  Blasting $acc against NR...\n";
  print $logh "  Blasting $acc against NR...\n";

  #Directory where sequences will be stored
  my $seqDir = "$accDir/sequences";
  system "mkdir -p $seqDir" unless (-d $seqDir);

  #Directory for miscelaneous files
  my $dataDir = "$accDir/data";
  system "mkdir -p $dataDir" unless (-d $dataDir);

  #Generate famXapnder dir
  my $fxpand_dir = "$accDir/famXpander";
  system "mkdir -p $fxpand_dir" unless (-d $fxpand_dir);

  my %acc2skingdom = ();
  my %skingdom2acc = ();
  my %nrCandidates = ();

  #Run Blast
  my $nseqs = run_ncbi_blast($acc, $seqDir, $fxpand_dir, $dataDir, \%ncbiBlast,  \%acc2skingdom, \%skingdom2acc, \%nrCandidates);

  unless ($nseqs) {
#    system "touch $accDir/LESS_THAN_${ncbiHom}_NCBI_HOMOLOGS";
    print $logh "    Total candidates from NCBI: $nseqs\n";
    return 0;
  }


  print "    Total candidates from NCBI: $nseqs\n";
  print $logh "    Total candidates from NCBI: $nseqs\n";

#  print Data::Dumper->Dump([\%nrCandidates ], [qw(*nrCandidates )]);
#  exit;

  #------------------------------------------------------------
  #Blast sequences against TCDB to determine whether or not
  #there are good hits with stablished families.

  print "    Comparing non-redundant homologs to TCDB...\n";
  print $logh "    Comparing non-redundant homologs to TCDB...\n";

  my %tcblast = ();
  blast_vs_tcdb(\%tcblast, $fxpand_dir, $dataDir);

  print "       Total TCDB hits: ", scalar keys %tcblast ,"\n";
  print $logh "       Total TCDB hits: ", scalar keys %tcblast ,"\n";

#  print Data::Dumper->Dump([\%tcblast ], [qw(*tcblast)]);
#  print "\n** Done! check the tcblast output **\n";
#  <STDIN>;



  #------------------------------------------------------------
  #Based on the famXpander and tcdblast results get the list of
  #members that should be added to the new tcdb family.

  print "    Selecting proteins to add to TCDB...\n";
  print $logh "    Selecting proteins to add to TCDB...\n";

  my %topHits = ();
  select_family_members(\%topHits, \%nrCandidates, \%acc2skingdom, \%skingdom2acc, \%tcblast,  $fxpand_dir, $dataDir, $seqDir, $acc);

#  print Data::Dumper->Dump([\%topHits, \%nrCandidates ], [qw(*topHits *nrCandidates )]);
#  exit;

  #See if there are at least {$minHoms} homologs
  my $nhits = scalar(keys %topHits);

#  print Data::Dumper->Dump([$nhits, \%topHits ], [qw(*nhits *topHits )]);
#  <STDIN>;

  unless ($nhits >= $minHoms) {
    system "touch $accDir/LESS_THAN_${ncbiHom}_CANDIDATE_HOMOLOGS";
    print "    Fewer than ${minHoms} non-redundant NCBI homologs: $nhits\n";
    print $logh "    Fewer than ${ncbiHom} non-redundant NCBI homologs: $nhits\n";
    return 0;
  }


  #------------------------------------------------------------
  #Generate reports if there are enough homologs

  print "    Generating report for: $acc ...\n";
  print $logh "    Generating report for: $acc ...\n";
  generate_reports($acc, \%topHits, $accDir, $seqDir);


  return 1;
}





sub generate_reports {

  my ($acc, $members, $outdir, $sDir) = @_;

  my $blast = $ncbiBlast{$acc};
  my $cnt   = 1;

#  print Data::Dumper->Dump([$acc, $members], [qw(*acc *members )]);

  #------------------------------------------------------------
  #Generate the text report

  my $repFile = "$outdir/${acc}_hom_report.txt";
  my @homs = ();

#  unless (-f $repFile && !(-z $repFile)) {
  open (my $fh, ">", $repFile) || die $!;

  #Print data for query protein
  my $skingdom = (exists $members->{$acc} && exists $members->{$acc}->{sk} && $members->{$acc}->{sk})?
    $members->{$acc}->{sk} : 'Not_found';

  print $fh "Query: \t$acc\t$skingdom\t", $members->{$acc}->{annot}, "\n";
  print $fh "Total NCBI homologs: ", $ncbiMatches{$acc}, "\n\n";

  print $fh "No\tHomolog\tE-value\tIdentity\tQ_cov\tS_cov\tTaxonomy\tFunction\n";

  #Print data for candidate homologs.
 HIT:foreach my $p (sort { unless ($a eq $acc || $b eq $acc) {$blast->{$a}->{eval} <=> $blast->{$b}->{eval} }} keys %{ $members }) {

    next if ($p eq $acc);

    #Store accession for later sequence retrieval
    push (@homs, $p);

    my $eval = sprintf("%.2e", $blast->{$p}->{eval});
    my $id   = $blast->{$p}->{id};
    my $qcov = sprintf("%.2f", $blast->{$p}->{qcov});
    my $scov = sprintf("%.2f", $blast->{$p}->{scov});
    my $sk   = $members->{$p}->{sk};
    my $func = $members->{$p}->{annot};

    print $fh "$cnt\t$p\t$eval\t$id\t$qcov\t$scov\t$sk\t$func\n";
    $cnt++;
  }
  close $fh;
#  }




  #------------------------------------------------------------
  #Generate the sequence and domain comparisons between the query
  #and its homologs

  #save the accessions of homologues to a file
  my $idFile = "$sDir/candAcc.txt";
  open (my $h1, ">", $idFile) || die $!;
  print $h1 join ("\n", @homs), "\n";
  close $h1;


  #Extract the sequences
  my $seqsF = "$sDir/candSeqs.faa";
  get_seqs_from_ncbi('', $idFile, $seqsF);


  #Now run alignSeqFiles to generate the final report
  my $qry = "$sDir/${acc}.faa";
  my $outDir = "$outdir/ssearch_${acc}_vs_homologs";

  my $rcmd = qq(alignSeqsFiles.pl -q $qry -ql $acc -s $seqsF -sl homologs -e 1 -c 0.3 -cc X -o $outDir);
  system $rcmd unless (-d $outDir);

}

#==========================================================================
#Select the most likely list of homologs that will be added to TCDB.
#Also include if there are significant with any protein in TCDB
#(This is to potentially assign the proteins to an exisiting family,
#even if the original query protein had no detectable similarity to
#anything in TCDB).

sub select_family_members {

  my ($bHits, $cand, $a2sk, $sk2a, $tcb, $fxDir, $dDir, $sDir, $acc) = @_;


  #------------------------------------------------------------
  #if the number of candidates is smaller than $minHoms,
  #prepare them all for manual inspection. Note that here I'm
  #trusting the superkingdom hash to contain all refseqs
  #(I should verify that).

  my $nCandidates = scalar keys %{ $cand };

  if ($nCandidates < $minHoms) {
    %$bHits = %$cand;
    return;
  }

#  print Data::Dumper->Dump([$nCandidates ], [qw(*nCandidates )]);
#  exit;

  #------------------------------------------------------------
  #Remove redundancies within the list of homologs, and make
  #sure the final list of candidates samples all
  #superkingdoms represented in the data.


  remove_redundant_sequences($bHits, $cand, $sDir, $fxDir, $dDir, $a2sk, $sk2a, $tcb, $acc);

}



#==========================================================================
#Remove redundancies within the list of homologs.
#(no two homologs should have E-value < $redundancyEvalue).
#
#Selected the final list of candidates such that all different
#superkingdoms are represented in the results.



sub remove_redundant_sequences {

  my ($res, $cand, $sDir, $fxDir, $dDir, $acc2sk, $sk2acc, $tcb, $acc) = @_;


  #-----------------------------------------------------------------
  #Get the available annotations for each protein

  my $seqFile  = "$fxDir/nrHomSeqs.faa";
  my $qSeq     = "$sDir/${acc}.faa";

  #Add the sequence of the query protein to the nr sequences in order to
  #add the function of the query to the %func hash table
  my $tmpFunc  = "$dDir/tmpFunc_seqs.faa";
  system "cat $seqFile $qSeq > $tmpFunc";

  #remove '>' and genome name from fasta headers, also remove the version number
  #from the accession (This code does not support TCIDs)
  my @annot = map {s/^\>//; s/\s+\[.+\]//; s/^([a-zA-Z0-9_]+)\S*/$1/; $_} split(/\n/, `grep '>' $tmpFunc`);

  #Generate a hash with the annotations per protein
  my %func = ();
  map {  $func{$1} = $2 if(/(\S+)\s+(.+)/) } @annot;
#  system "rm $tmpFunc" if (-f $tmpFunc);

#  print Data::Dumper->Dump([$acc, $seqFile, \@annot, \%func ], [qw(*acc *seqFile *annot *func )]);
#  print Data::Dumper->Dump([$tcb ], [qw(*tcb )]);
#  exit;



  #------------------------------------------------------------
  #If total number of sequences is less than $membersNum,
  #accept all sequences.

  my $nCandidates = scalar keys %{ $cand };

#  print "** Candidates for analysis: $nCandidates\t\tMembers allowed: $membersNum **\n";


  #Candidates are less than the minimum number, take them all.
  if ($nCandidates <= $membersNum) {

    foreach my $c (keys %{ $cand }) {
      $res->{$c}->{annot} = $func{$c};
      $res->{$c}->{sk}    = $acc2sk->{$c};

      #If there are blast hits with TCDB get the info of the best match
      my @best = ();
      if (exists $tcb->{$c}) {
	@best = (sort { $a->{eval}<=>$b->{eval} } @{ $tcb->{$c} })[0];
      }
      $res->{$c}->{blast} = \@best;
    }

#    print Data::Dumper->Dump([$res ], [qw(*res )]);
#    exit;
  }


  #------------------------------------------------------------
  #More homologs than nencessary, favor homologs with TCDB hits

  else {

    #Get the counts of genes per superkingdom
    my %skcounts = ();
    my $cnt1 = 0;
    foreach my $p (keys %{ $cand }) {
      my $superK = (exists $acc2sk->{$p})? $acc2sk->{$p} : "";

      if ($superK) {
	$skcounts{$superK}++;
	$cnt1++;
      }
    }

    die "No Superkingdoms detected in \%cand" unless ($cnt1);
    map { $skcounts{$_} = $skcounts{$_} / $cnt1 } keys %skcounts;



    #The expected number of accessions per superkingdom
    my %expCounts = ();
    foreach my $superK (keys %skcounts) {
      $expCounts{$superK} = int($skcounts{$superK} * $membersNum + 1);
    }


    #Now extract the numbers of homologs according to their expected distribution
    my $blast = $ncbiBlast{$acc};
    my $cnt   = 0;

#    print Data::Dumper->Dump([$blast ], [qw(*ncbiBlast )]);
#    exit;

  HIT:foreach my $c (sort { $blast->{$b}->{eval} <=> $blast->{$a}->{eval} unless (exists $blast->{$b} && exists $blast->{$a} )} keys %{ $cand }) {

      #    print Data::Dumper->Dump([$cnt, $c, $blast->{$c} ], [qw(*cnt, *c *bl )]);
      #    <STDIN>;

      last HIT if ($cnt > $membersNum);

      my $sking = $acc2sk->{$c};
      next HIT if ($expCounts{$sking} <= 0);

      $res->{$c}->{annot} = $func{$c};
      $res->{$c}->{sk}    = $sking;

      #if there is a blast hit with TCDB get the info of the best match
      my @best = ();
      if (exists $tcb->{$c}) {
	@best = (sort { $a->{ev}<=>$b->{ev} } @{ $tcb->{$c} })[0];
      }
      $res->{$c}->{blast} = \@best;

      $cnt++;
      $expCounts{$sking}--;
    }
  } #else



  #Add data for the query protein if not in the results already
  unless (exists $res->{$acc}) {

    unless (exists $func{$acc} && exists $acc2sk->{$acc}) {
      $res->{$acc}->{annot} = $func{$acc};
      $res->{$acc}->{sk}    = $acc2sk->{$acc};

      #if there is a blast hit with TCDB get the info of the best match
      my @best = ();
      if (exists $tcb->{$acc}) {
	@best = (sort { $a->{ev}<=>$b->{ev} } @{ $tcb->{$acc} })[0];
      }
      $res->{$acc}->{tcblast} = \@best;
    }
  }
}



#==========================================================================
#After runnign famXapnder and TCBLAST, this function will extract
#the sequences of the proteins that passed the filters for the
#purpose of removing redundant sequences


sub getSequencesForCandTransporters {
  my ($acc2sk, $dDir, $seqFile) = @_;

  my $tmpIDs = "$dDir/tmp_acc_list.txt";

  #send the list of accessions to a file
  open (my $oh, ">", $tmpIDs) || die $!;
  my $str = join ("\n", keys %$acc2sk);
  print $oh "$str\n";
  close $oh;

  die "Unavailable (or empty) file with sequences: $tmpIDs" unless (-f $tmpIDs && !(-z $tmpIDs));


  #Now extract the sequences based on the extracted list of accessions
  my $cmd2 = qq(blastdbcmd -db $gnmBlastDB -outfmt '\%f' -entry_batch $tmpIDs -target_only -out $seqFile  &> /dev/null);
  print "$cmd2\n";

  system $cmd2 unless (-f $seqFile);
  die "Sequence file Empty or not found: $seqFile" unless (-f $seqFile && !(-z $seqFile));
#  system "rm $tmpIDs" if (-f $tmpIDs);

}



#==========================================================================
#Blast homologs extracted from NCBI against TCDB and select the best
#candidates for TCDB.

sub blast_vs_tcdb {

  my ($out, $fxDir, $dataDir) = @_;

  #------------------------------------------------------------
  #Run blast against TCDB

  my $blastdb  = "$tcblastdbDir/tcdb";
  my $infile   = "$fxDir/nrHomSeqs.faa";
  my $blastOut = "$dataDir/tcblast.out";

  my $args1  = qq(-query $infile -db $blastdb -evalue 1 -max_target_seqs 100000 -use_sw_tback -out $blastOut -comp_based_stats F);
  my $args2  = qq(-outfmt '7 qacc sacc bitscore evalue pident qstart qend qlen sstart send slen');
  my $cmd1   = qq(blastp $args1 $args2);

#  print "$cmd1\n";
  system $cmd1 unless (-f $blastOut && !( -z $blastOut ));


  #------------------------------------------------------------
  #Parse blast output and select best hits

  open(my $bh, "<", $blastOut) || die "Could not open file: $blastOut --> $!";
  while(<$bh>) {
    chomp;
    next if (/^#/);

    my ($q, $s, $bs, $e, $pid, $qs, $qe, $ql, $ss, $se, $sl) = split(/\t/);

    #Check that all variables have values
    next unless ($q && $s && $bs && $pid && $qs && $qe && $ql && $ss && $se && $sl);

    my $qcov  = ($qe - $qs)/$ql;
    my $scov  = ($se - $ss)/$sl;
    my $rslen = $sl/$ql * 1.0;


    #Keep sequences with the right E-values, coverage and sequence length
    next unless (($e >= $tceval && $e <= $evalue) && ($qcov >= $tcoverage || $scov >= $tcoverage) && $rslen <= $maxSlen);


    #Remove version number from query
    my ($qacc, $qv) = split(/\./, $q);

    #Record alignment data
    my %data = (q=>$qacc, s=>$s, eval=>$e, qlen=>$ql, slen=>$sl, qcov=>$qcov, scov=>$scov);
    push (@{ $out->{$qacc} }, \%data);
  }
  close $bh
}



#==========================================================================
#Given a file with sequences in fasta format, extract the superkingdom data
#for the sequences in that file

sub retrieve_superkingdom {

  my ($acc, $outDir, $infile, $acc2sk, $sk2acc) = @_;


  #Extract accessions of nonredundant NCBI homologs
  my $tmpFile = "$outDir/tmp.txt";
  my $accFile = "$outDir/accList.txt";

  #Keep only accessions
  my $part1 = qq(grep '>' $infile | perl -pe 's/^\\>\([a-zA-Z0-9\._]+\).*\$/\$1/;');

  #Remove version numbers
  my $part2 = qq(perl -pe 's/\\.\\d\$//;' > $tmpFile);

  #Just in case, Include query accession in the list
  my $part3 = "";
  if ($gnmBlastDB eq 'nr') {
    $part3 = qq(echo $acc >> $tmpFile; );
  }
  $part3 .= qq(sort $tmpFile | uniq > $accFile);
  my $cmd1  = qq($part1 | $part2; $part3);


  system $cmd1 unless (-f $accFile);
  die "Error: empty file with list of RefSeq accessions --> $accFile" if (-z $accFile);
  system "rm $tmpFile" if (-f $tmpFile);


#  print "Check $accFile\n";
#  exit;


  #Extract the superkingdom data
  my $skingdomFile = "$outDir/superkingdom.txt";
  my $cmd2 = qq(blastdbcmd -db nr -outfmt '\%a \%K' -entry_batch $accFile -target_only -out $skingdomFile  &> /dev/null);

  system $cmd2 unless (-f $skingdomFile);
  die "Error: empty file with Superkingdom data" if (-z $skingdomFile);


  #Parse superkingdom data
  open (my $skh, "<", $skingdomFile) || die "Could not open: $skingdomFile --> $!";
  while(<$skh>) {
    chomp;
    my ($id, $skingdom) = split(/\s+/);
    my ($rs, $version)  = split(/\./, $id);
    $acc2sk->{$rs}  = $skingdom;
    push (@{ $sk2acc->{$skingdom} }, $rs);
  }
  close $skh;
}



#==========================================================================
#First run blast of the accesion against NCBI and bring all homologs with
#high coverage and below the specified evalue threshold


sub run_ncbi_blast {

  my ($query, $sDir, $fxDir, $dDir, $blastOut, $acc2sk, $sk2acc, $nrProts) = @_;


  my $queryDir = "$workDir/$query";


  #Extract sequence
  my $seqFile = "$sDir/${query}.faa";
  my $cmd1 = qq(blastdbcmd -db $gnmBlastDB -entry $query -target_only -out $seqFile  &> /dev/null);
  system $cmd1 unless (-f $seqFile);

  #No sequence extracted for: $query
  if (-z $seqFile) {
    system "touch $queryDir/NO_SEQUENCE_FOR_QUERY_$query";
    return 0;
  }


  #------------------------------------------------------------
  #Now run famXpander. Since we are looking for remote homologs,
  #I'll have to parse the psiblast.tbl file instead of working with
  #file results.faa

  my $subLen  = (1 - $coverage) + 1;
  my $cmd2    = qq(famXpander.pl -i $seqFile -e 1.0 -h F -c $coverage -s $coverage -x F -l $maxSlen -r 0.7 $blastMode -o $fxDir);
  my $seqsf   = "$fxDir/results.faa";
  system $cmd2 unless (-f $seqsf && !(-z $seqsf));


  #Aparently no sequences were extracted from NCBI
  unless (-f $seqsf && !(-z $seqsf)) {
    system "touch $queryDir/NO_NCBI_HOMOLOGS_FOR_$query";
    return 0;
  }



  #------------------------------------------------------------
  #Parse blast output to determine the number of sequences
  #that can be used for further analysis. Do not use the
  #results.faa file because we are not interested in close
  #homologs.

  my $blFile = "$fxDir/psiblast.tbl";

  $ncbiMatches{$query} = 0;

  open(my $fh, "<", $blFile) || die $!;
  while(<$fh>) {
    chomp;
    next if (/^[#\n]/);

    my ($qid, $sacc, $kk1, $kk2, $ev, $id, $qs, $qe, $ql, $ss, $se, $sl, @kk) = split(/\t/);

    #remove version number from query accession.
    my ($qacc, $qv) = split(/\./, $qid);

    next unless ($qid && $ev); #There must be an E-value
    next if ($qacc eq $sacc);  #ignore self-matches.

    my $qcov = ($qe - $qs) / $ql;
    my $scov = ($se - $ss) / $sl;


    next unless (($ev <= $evalue && $ev >= $tceval) && ($qcov >= $coverage || $scov >= $coverage));

    $ncbiMatches{$query}++;
    $blastOut->{$qacc}->{$sacc} = {qacc=>$qacc, sacc=>$sacc, eval=>$ev, id=>$id, qcov=>$qcov, scov=>$scov};
  }
  close $fh;


#  print Data::Dumper->Dump([$blastOut, \%ncbiMatches], [qw(*blastOut *ncbiMatches )]);
#  exit;


  #------------------------------------------------------------
  #Extract the accessions of all hits

  #Save accessions to a file
  my $homAccFile1 = "$fxDir/homTmpAcc1.txt";
  save_acc2file ($homAccFile1, $blastOut->{$query});



  #------------------------------------------------------------
  #extract the sequences of the proteins that passed the threshold.

  my $redHomSeqsFile1 = "$fxDir/homTmpSeqs.faa";
  get_seqs_from_ncbi($query, $homAccFile1, $redHomSeqsFile1);



  #------------------------------------------------------------
  #Remove any partial sequences (based on the fasta header line)

  #remove '>' and version numbers after discarding partial sequences
  my @annot = map {s/^\>//; s/^([a-zA-Z0-9_]+)\.\d/$1/; $_} split(/\n/, `grep '>' $redHomSeqsFile1 | grep -v partial`);

#  print Data::Dumper->Dump([\@annot ], [qw(*annotations )]);
#  exit;

  #fewer then the minimum number of accepted NCBI sequences
  my $nanno = scalar @annot;
  if ($nanno < $ncbiHom) {
    print "      Less than ${ncbiHom} homologs in NCBI for ${query}:  $nanno\n";
    system "touch $queryDir/LESS_THAN_${ncbiHom}_NCBI_HOMOLOGS_FOR_${query}_$nanno";
    return 0;
  }


  #Generate a hash with the accessions of full proteins
  my %goodAcc = ();
  map {  $goodAcc{$1} = 1 if(/(\S+)\s+.+/) } @annot;

#  system "rm $redHomSeqsFile1" if (-f $redHomSeqsFile1);

  #Delete query entry from %goodAcc if a custom genome blastdb was provided
  delete($goodAcc{$query}) if ($gnmBlastDB ne 'nr');


#  print Data::Dumper->Dump([\%goodAcc ], [qw(*goodAcc )]);
#  exit;




  #------------------------------------------------------------
  #Extract sequences excluding partial proteins

  #save accession to file
  my $homAccFile2 = "$fxDir/redHomAcc.txt";
  save_acc2file ($homAccFile2, \%goodAcc);

  #Extract sequences for redundant proteins
  my $redHomSeqsFile2 = "$fxDir/redHomSeqs.faa";
  get_seqs_from_ncbi($query, $homAccFile2, $redHomSeqsFile2);



  #------------------------------------------------------------
  #Get superkingdom data to help in the selection of nonredundant
  #proteins

  retrieve_superkingdom($query, $dDir, $redHomSeqsFile2, $acc2sk, $sk2acc);

#  my $redProts = scalar keys %$acc2sk;
#  print Data::Dumper->Dump([$redProts, $sk2acc ], [qw(*redProts *sk2acc )]);
#  exit;



  #------------------------------------------------------------
  #Blast all homologs vs all homologs in order to remove redundancies

  my $blastOutFile = "$fxDir/blast_all_vs_all.out";
  my $fmtBlastOut  = '7 qacc sacc qlen slen evalue pident qcovs qstart qend sstart send';


  unless (-f $blastOutFile && !(-z $blastOutFile)) {
    system "blastp -query $redHomSeqsFile2 -subject $redHomSeqsFile2 -evalue 1 -comp_based_stats 2 -seg no -max_target_seqs 100000 -max_hsps 1 -use_sw_tback -outfmt '$fmtBlastOut' -out $blastOutFile";
  }


  unless (-f $blastOutFile && !(-z $blastOutFile)) {
    print "All vs all Blast output file not generated or empty: $blastOutFile";
    system "touch $queryDir/PROBLEM_ALLvsALL_BLAST_OUTPUT";
    return 0;
  }




  #------------------------------------------------------------
  #Identify the redundant candidates. Candidates with very different size
  #are allowed. When 2 redundant candidates have similar size, remove the
  #the smallest one.


  my %redundantCandidates = ();
  my %prots = ();
  my $querySK = (exists $acc2sk->{$query})? $acc2sk->{$query} : "unknown";

  open (my $blasth, "<", $blastOutFile) || die $!;
 MATCH:while (<$blasth>) {

    chomp;
    next if (/^#/ || !$_);

    my ($qry, $sbj, $qlen, $slen, $ev, $pIdent, $qcovs, @kk) = split(/\t/, $_);

    #There must be data to continue
    die "Could not extract id for query: $qry" unless ($qry && $sbj && $qcovs);

    #remove version numbers
    my $qacc = ($qry =~ /(\S+)\.\d$/)? $1 : $qry;
    my $sacc = ($sbj =~ /(\S+)\.\d$/)? $1 : $sbj;

    #Ignore self comparisons
    next if ($qacc eq $sacc);


    #record all hits here
    $prots{$qacc} = $ev;
    $prots{$sacc} = $ev;


    #When there are redundant hits with very different size keep both.
    my $lenRatio = (sort {$b <=> $a} ($qlen/$slen, $slen/$qlen))[0];


    if ($ev < $redundancyEvalue && $lenRatio < 1.8) {

      my $q = $acc2sk->{$qacc};
      my $ssk = $acc2sk->{$sacc};


      #
      #When size is comparable, the smallest protein is redundant, unless one of
      #the proteins is the original blast query
      #

      if ($qacc ne $query && $sacc ne $query) {

	if ($qlen > $slen) {
	  $redundantCandidates{$qacc} = $ev; # if ($querySK eq $qsk);
	}
	elsif ($qlen < $slen) {
	  $redundantCandidates{$sacc} = $ev; # if ($querySK eq $ssk);
	}
	else {
	  #Same length: sort by accession and keep the first one
	  my $redAcc = (sort {$a cmp $b} ($qacc, $sacc))[0];
	  $redundantCandidates{$redAcc} = $ev; # if ($querySK eq $redAcc);
	}
      }


      #
      #The other protein is by default redundant if either the query
      #or subject proteins are the original blast query
      #

      elsif ($qacc eq $query) {
	$redundantCandidates{$sacc} = $ev;
      }
      elsif ($sacc eq $query) {
	$redundantCandidates{$qacc} = $ev;
      }
      else {
	die "Unknown condition; Blast line: $_";
      }

    }
  }
  close $blasth;

  my $total  = scalar keys %prots;
  my $redund = scalar keys %redundantCandidates;


#  print Data::Dumper->Dump([$prots{'Dvar_18500'}, $redundantCandidates{'Dvar_18500'} ], [qw(*prots *redCand )]);
#  exit;

  #Now get the set difference to identify the non_redundant genes.
  map { $nrProts->{$_} = [$prots{$_}, $acc2sk->{$_}] unless (exists $redundantCandidates{$_})  } keys %prots;

#  print Data::Dumper->Dump([$total, $redund, \%redundantCandidates, $nrProts], [qw(*total *redund *redundantCandidates *nrProts)]);
#  exit;



  #Save accessions of non-redundant genes to file
  my $nrAccList = "$fxDir/homTmpAcc1.txt";
  save_acc2file ($nrAccList, $nrProts);

  my $nrSeqs = "$fxDir/nrHomSeqs.faa";
  get_seqs_from_ncbi($query, $nrAccList, $nrSeqs);



  #------------------------------------------------------------
  #Count the number of non redudnant sequences returned by famXpander.
  #The number of homologs in NCBI must be at least $minHoms

  my $cmd5  = qq(grep -c '>' $nrSeqs);
  my $nseqs = qx($cmd5);
  chomp $nseqs;



  #Return value, at this point we want at least the number of proteins that could be
  #used to create the new family in TCDB, because this proteins are not redundant.
  ($nseqs < $minHoms)? return 0 : return $nseqs;

}



#==========================================================================
#Given a file with accessions, extract the sequences from the
#NCBI blast database.

sub get_seqs_from_ncbi {

  my ($qacc, $accFile, $seqsFile) = @_;

  my $cmd = qq(blastdbcmd -db nr -entry_batch $accFile -target_only -out $seqsFile &> /dev/null);
  system $cmd unless (-f $seqsFile);

  die "No sequences extracted: $seqsFile" unless (-f $seqsFile || !(-z $seqsFile));


  #Add the query protein sequence if genome blast database is not 'nr'
  if ($qacc &&  $gnmBlastDB ne 'nr') {
    my $qseq = "$workDir/$qacc/sequences/${qacc}.faa";
    die "Query sequence not found: $qseq" unless (-f $qseq);

    system qq(cat $qseq >> $seqsFile);
  }

  #Delete accessions file, it's no longer needed
#  system "rm $accFile" if (-f $accFile);

}


#==========================================================================
#Given a set of accessions as keys in a hash table,
#generate a file with the accession for later
#sequence extraction

sub save_acc2file {

  my ($outfile, $inhash) = @_;

  unless (-f $outfile && !(-z $outfile)) {
    open(my $ah, ">", $outfile) || die $!;
    print $ah join("\n", keys %{ $inhash }), "\n";
    close $ah;
  }
  die "File with accessions empty or not available: " unless (-f $outfile && !(-z $outfile));

}




#==========================================================================
#Read command line arguments

sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #Intermediate step to define the right path to $outfile
  my $outfileTmp = undef;

  my $status = GetOptions(
      "wd|workdir=s"           => \$workDir,
      "t|tcblastdb-dir=s"      => \$tcblastdbDir,
      "gdb|gnm-blastdb=s"      => \$gnmBlastDB,
      "i|input=s"              => \&read_input,
      "tms=s"                  => \$minTMS,
      "e|evalue=f"             => \$evalue,
      "c|coverage=f"           => \$coverage,
      "re|redundancy-evalue=f" => \$redundancyEvalue,
      "m|mode=s"               => \&read_blast_mode,
      "nh|ncbi-hom=n"          => \&read_ncbi_hom,
      "nm|num-members=n"       => \&read_num_members,
      "h|help"                 => sub { print_help(); },
      "<>"                     => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);


  #make sure working directory exists
  system "mkdir -p $workDir" unless (-d $workDir);

  #make sure the tcdb BlastDB dir exists
  system "mkdir -p $tcblastdbDir" unless (-d $tcblastdbDir);

  #Output file with final Unknown transporters in the input genome
  $tsvReport  = "$workDir/novelTransporters.txt";

  #Maximum length of the subject protein relative to the query.
  $maxSlen = (1 - $coverage) + 3;


  #Validate the tcblastdb
  my $testDBfile = "$tcblastdbDir/tcdb.pin";
  generate_tcdb_blastdb() if ($owrite_tcblast || !( -f $testDBfile ));
  die "TCDB blastDB not found in dir: $tcblastdbDir" unless (-f $testDBfile);
}


#==========================================================================
#Regenerate the TCDB blast DB

sub generate_tcdb_blastdb {
  my $cmd = qq(rm $tcblastdbDir/*; extractTCDB.pl -i tcdb -f blast -o $tcblastdbDir; extractTCDB.pl -i tcdb -o $tcblastdbDir);
  system $cmd;
}



#==========================================================================
#Option -i. Read the input accession as either list or input file

sub read_input {

  my ($opt, $value) = @_;

  $input = $value;

  #input file
  if (-f $input) {
    if (-z $input) {
      die "Input file is empty: $input";
    }

    read_acc_from_file($input);
  }

  #Comma-separated list of accessions is assumed
  else {
    @accessions = split(/,/, $input)
  }

}



sub read_acc_from_file {

  my $infile = shift;

  my $fileh = undef;

  #gzip compressed file
  if ($infile =~ /\.gz$/) {
    open ($fileh, "-|", "zcat $infile") || die $!;
  }

  #bzip2 compressed file
  elsif ($infile =~ /\.bz2$/) {
    open ($fileh, "-|", "bzcat $infile") || die $!;
  }

  #text file
  else {
    open ($fileh, "<", $infile) || die $!;
  }

  #extract first column whle ignoring anything else.
  #Also remove version numbers.
  while(<$fileh>) {
    chomp;

    next if (/^#/);

    my ($acc, @kk) = split(/\s+/);

    #Remove version number from accession. Note that this
    #will not work with TCIDs
    $acc =~ s/\.\d$//;

    push(@accessions, $acc);

  }
  close $fileh;
}


#==========================================================================
#Option -m: Read blast mode

sub read_blast_mode {

  my ($opt, $value) = @_;

  my $tmp = uc $value;
  die "Option -m only accepts L or R values." unless ($tmp =~ /^[RL]$/);


  $blastMode = ($tmp eq "L")? "-p F" : "-p T";

}


#==========================================================================
#Minimum number of homologus in NCB to consider a particular protein

sub read_ncbi_hom {

  my ($opt, $value) = @_;

  die "Minimum value for option -nh: 15" unless ($value >= 15);

  $ncbiHom = $value;

}



#==========================================================================
#Number of member proteins to create a new family in TCDB

sub read_num_members {

  my ($opt, $value) = @_;

  die "Minimum value for option -nm: 3" unless ($value >= 3);

  $minHoms = $value
}


#==========================================================================
#Print the help of the program

sub print_help {

  my $help = <<'HELP';

 Given a list of protein accessions with no significant sequence similarity
 to anything in TCDB, extract a list of homologs from NCBI that will be used
 to create a new family in TCDB.

 Input paramateres

 -i, --input {list|file}  (Mandatory)
    An Accession, comma-separated list of accessions, or a file with the
    accessions (first columns) on which the perform the analysis.

 -wd, --workdir {path}  (Defalut: ./newIMPs4tcdb)
    Path to the directory where the ouput and temporary files will be stored.

 -t, --tcblastdb {path}  (Default: ~/db/blastdb/tcdb)
    Full path to the TCDB blast database that will be used.

 -gdb, --gnm-blastdb {string} (Default: nr)
    Full path to the blast DB of the whole proteome of the reference
    genome. This is necessary when proteins are not annotated with proper
    RefSeq accessions (e.g. locus_tags).

 -e, --evalue {float}  (default: 1e-3)
    NCBI proteins with E-value greater than this cutoff will be ignored.

 -re,  --redundancy-evalue {float}  (Default: 1e-15)
    Evalue threshold to consider two proteins redundant for the
    the purpose of presenting final results.

 -c, --coverage {float} (Default: 0.5)
    Minimum alignment coverage in pairwise alignments.

 -m, --mode {L|R} (Default: R)
    Indicate whether blast will be run locally (L) or remotely
    at NCBI (R).

 -nm, --num-members {int} (Defatult: 3)
   Number of non-redundant homologs to create a family in TCDB,
   or to add them to an existing family as remote homologs.

 -nh, --ncbi-hom {int} (Default: 25)
    Minimum number of homologs in NCBI to qualify for this analysis.

 -h, --help
    Display this help. This argument takes precedence over any other
    argument.

HELP

  print $help;
  exit;
}
