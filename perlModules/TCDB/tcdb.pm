package TCDB::tcdb;


use warnings;
use strict;
use Data::Dumper;


#==========================================================================
# This module contains useful functions to handle tcdb-related
# data and programs.
#
# Functions in this module:
#
# run_famXpander()
# run_protocol1()
# run_proto2()
# run_gsat_for_top_proto2_hits()
# get_top_proto2_hits()
# read_prot2_output_faa_seqs()
# get_gsat_score()
# count_sequences_in_file()
# count_sequences_in_family()
# read_2col_seqs()
# extract_seqs_from_tcdb()
# validate_tcdb_id()
# validate_files()
# run_hmmtop()
#
#==========================================================================

#==========================================================================
#Run hmmtop



#==========================================================================
#Run famXpander for a sequence file in fasta format
#
#Input parameters:
#
# $id:  tcdb_id or file name (without extension) containing the sequences
#       in fasta format that will be analyzed by famXpander.
#
# $infile: Full path the file with the seuence file in fasta format.
#
# $p1d: Directory where the results of famXpander will be stored.
#
# $eval: e-value threhold for psiblast to include alignments in the 
#        results.
#
# $ieval: inclusion e-value for psiblast paiwise alignments.
#
# $it: Number of iterations that psiblast will perform.
#
# $s_min_cov: Minimal coverage of the query sequence in the alignments
#             (e.g. 0.8)
#
# $t_min_len: minimal sequence length of target relative to subject
#             length. Value expressed in 0 to 1 values (e.g. 0.85)
#
# $t_max_len: maximal sequence length of target relative to subject
#             length (e.g. 2.0 is twice the length of subject).
#
# $cdhit: cdhit threshold to consider two sequences reduntant (e.g. 0.95)
#
# $kar: (T/F) indicates whether famXpander will report only the psiblast
#       aligned regions ('T') or the complet proteins ('F')
#
# $rnaln: Retrieve this number of alignments from psiblast
#         (can be used to speed up the process).

sub run_famXpander {

  my ($id, $infile, $p1d, $eval, $ieval, $it, $s_min_cov,  $t_min_len, $t_max_len, $cdhit, $kar, $rnaln) = @_;


  #Output dir for famXpander program
  my $famXpand_res = "$p1d/$id";
  system "mkdir -p $famXpand_res" unless (-d $famXpand_res);


  #If famXpander output file already exists, move to the next sequence
  my $outfile = "$famXpand_res/results.faa";


  unless (-f $outfile) {

    print "Running famXpander for: $id\n";


    #Run famXpander
    my $params = "-e $eval -f $ieval -t $it -c $s_min_cov -x T ";
    $params .= "-s $t_min_len -l $t_max_len -r $cdhit ";
    $params .= "-h $kar -n $rnaln";

    my $cmd = "famXpander.pl -i $infile -o $famXpand_res $params";
    print "$cmd\n\n";
    system $cmd;

    die "No famXpander results for $id --> " unless (-f $outfile  && ! (-z $outfile));
  }
}




#==========================================================================
#Run protocol1 for each protein in a sequence file (in 2-column format)
#
#Input parameters:
#
# $famSeqFile: file in two-column format with all the sequences in one
#              family.
#
# $p1d: The directory where all the results of protocol1 will be placed.
#
# $eval: The evalue that will be used to run psi-blast.
#
# $cdhit_cut: The cdhit identity level (e.g. 0.9) at which psi-blast
#             output sequences will be clustered to eliminate redundant
#             sequences.
#
# $psib_it: Number of iteration that psi-blast will run. Note that
#           currently protocol1 only works for one iteration.
# $rnaln: Retrieve this number of alignments from blast
#         (can be used to speed up the process).


sub run_protocol1 {

  my ($famSeqFile, $p1d, $eval, $cdhit_cut, $psib_it, $rnaln) = @_;

  #Read the sequences in column format
  my @inSeqs = ();
  read_2col_seqs($famSeqFile, \@inSeqs);


#  print Data::Dumper->Dump([\@inSeqs], [qw(*inSeqs)]);
#  exit;


  foreach my $seq (@inSeqs) {

    #Outputs for protocol1 program
    my $proto1_res = "$p1d/$seq->[0]";
    system "mkdir -p $proto1_res" unless (-d $proto1_res);


    #If protocol1 output file already exists, move to the next sequence
    unless (-f "$proto1_res/results.faa") {

      print "Running protocol_1 for: ", $seq->[0], "\n";

      #length of que sequence
      my $len = length ($seq->[1]);


      #Limit blast hits to sequences at least 70% the size of the query protein
      #option --min=$lim
      my $lim = int (0.7 * $len);


      #Run Protocol 1
      my $tries = 5;
      my $cnt = 1;
      while ($cnt <= $tries) {

	#Currently protocol1 only runs one iteration
	my $it = $psib_it - 1;
	my $params = "-o $proto1_res -e $eval -c $cdhit_cut -n $rnaln --min=$lim -i $it";

	my $cmd = "protocol1.py $params -q $seq->[1]";
	print "Running:  $cmd\n";
	system $cmd;

	#output file exists. It can be empty because some proteins return no psi-blast hits
	if (-f "$proto1_res/results.faa") {
	  $cnt = 100;
	} else {
	  sleep(10);
	  $cnt++;
	}
      }
    } #unless output file exists
  } #foreach $seq

}





#==========================================================================
#Run protocol2. If successful, this function returns the directory where
#the results of running protocol2 were stored.
#
#Input parameters are:
#
#$f1: tcdb_id or custom_id of the first family that will be compared.
#     this ID must correspond to a directory name in the $p1dir input path.
#
#$f1: tcdb_id or custom_id of the second family that will be compared.
#     this ID must correspond to a directory name in the $p1dir input path.
#
#$p1dir: path to the directory where the results of running protocol1 are
#        stored.
#
#$p1dir: path to the directory where the results of running protocol2 will
#        be stored.
#
#$t_p1_pairs: variable to count the number of sequences that will be
#             compared by protocol2. It's the multiplication of the
#             number of sequences in both input files.



sub run_proto2 {

  my ($f1, $f2, $p1dir, $p2dir, $p2_min_hit_len, $t_p1_pairs) = @_;

  #Protocol1 output file for protein in first family
  my $seq1_name  = $f1;
  my $famID_1 = $f1;
  my $seq1_p1_faa = "$p1dir/$f1/results.faa";

  return undef unless (-f $seq1_p1_faa && ! (-z $seq1_p1_faa));

  #Count the number of sequences generated by protocol_1 for this family member
  my $nseqs1 = count_sequences_in_file($seq1_p1_faa);


  #Protocol1 output file for protein in second family
  my $seq2_name  = $f2;
  my $famID_2 = $f2;
  my $seq2_p1_faa = "$p1dir/$f2/results.faa";


  return undef unless (-f $seq2_p1_faa && ! (-z $seq2_p1_faa));


  #Count the number of sequences generated by protocol1 for this family member
  my $nseqs2 = count_sequences_in_file ($seq2_p1_faa);


  #Count the number of protein pairs that will be considered by protocol2
  ${ $t_p1_pairs } += $nseqs1 * $nseqs2;



  #Setup the output directory for protocol_2
  my $cmp_title = "${seq1_name}_vs_${seq2_name}";
  my $outdir    = "$p2dir/$cmp_title";
  system "mkdir $outdir" unless (-d $outdir);


  #If protocol2 results file already exist, do not re-run protocol2
  my $outFile = "$outdir/report.tbl";
  unless (-f $outFile) {

    print "  Running protocol2 for $cmp_title\n";

    #run Protocol_2
    my $params = "-s $seq1_p1_faa -t $seq2_p1_faa -o $outdir --subject=$seq1_name --target=$seq2_name ";
    $params .= "--shuffle 1000 --min=$p2_min_hit_len";

    my $cmd = "protocol2.py $params";
    print "Running:  $cmd\n";
    system $cmd;
  }

  return $outdir;
}



#==========================================================================
#This function parses protocol2's output to determine wheter there were
#significant hits, and then runs gsat on the selected hits. This function
#does not return any value.
#
#The input parameters are:
#
#$p2wd: directory where the results of protocol2 are located.
#
#$gmode: Indicates whether GSAT will run with full protein sequences (full)
#       or just with the segments aligned by protocol2.
#
#$p2lcutoff: The low cutoff from which protocol2 scores are considered significant.
#
#$p2hcutoff: The high cutoff from which protocol2 scores are considered significant.
#            This defines a range of protocol2 scores in combination with $p2lcutoff.
#
#$gcutoff: The cutoff from which the GSAT score is considered significant.
#
#$p2_hsp_cnt: Counts the number of protocol2 hits with significant score.
#
#$g_hsp_cnt: Counts the number of GSAT hits with significant score
#
#$p1_prog: Program that was run for protocol 1 (i.e. proto1 or fxpand)
#
#$p1d:  Directory with the results of proto1 or fxand. This will help to
#       extract the sequences of the original TCDB protein and the psiblast
#       hit to complete the transitivity rule.
#
#$tmoCutOff: minimum TMS overlap acceptable in a protocol 2 hit
#

sub run_gsat_for_top_proto2_hits {

  my ($p2wd, $gmode, $p2lcutoff, $p2hcutoff, $gcutoff, $p2_hsp_cnt, $g_hsp_cnt, $p1_prog, $p1d, $tmoCutOff, $shuffles, $ignore_proteins) = @_;


  #Select the hits that scored above the gsat_cutoff threshold and run GSAT
  my $outFile = "$p2wd/report.tbl";

  #In case no shuffles are given initialize it here
  #(This is to be compatible with areFamiliesHomologous)
  $shuffles = ($shuffles)? $shuffles : 1000;



  my $sbj_fam = "";
  my $tgt_fam = "";


  #identify if some proteins will be ignored
  my %ignoreProteins = ();
  if ($ignore_proteins) {
    my @ignore_ids = split (/,/, $ignore_proteins);
    %ignoreProteins = map {$_ => 1} @ignore_ids;
  }

#  print Data::Dumper->Dump([$ignore_proteins, \%ignoreProteins ], [qw(*ignore_proteins *ignoreProteins )]);
#  exit;

#  print "Mode: $gmode\n";
#  exit;


  my @top_prot2_hits = ();
  if (-f $outFile) {


    #Select the hits that scored above the gsat_cutoff threshold. The array
    #@top_prot2_hits has the format:
    #[ ['XP_006521872',  'XP_012536656',    19,
    #   'RQCGVCRQADEEFQVLANSWRYSSAFTNKVF---FASVDFDEGSDVFQMMNMNSAPTFIHFPYKGKPRKSDTYE...',
    #   'RQCSVVLSVDEEFEAAVRRLRNVGALTAPVVDREMRLVGLVTAADIIREMELEATDDILRFGGVESPLQRDDIE...',
    #   15.80 #This is the TMO value
    #  ],... ]


    get_top_proto2_hits(\@top_prot2_hits, \$sbj_fam, \$tgt_fam, $outFile, $gmode, $p2lcutoff, $p2hcutoff, $tmoCutOff, \%ignoreProteins);
#    print Data::Dumper->Dump([\@top_prot2_hits],["*top_prot2_hits"]);
#    exit;



    #Determine whether Protocol2 had significant hits or not
    if (@top_prot2_hits) {
      system "touch $p2wd/PROT2_PASSED";
    } else {
      system "touch $p2wd/PROT2_FAILED";
      return undef;
    }


    my $gsat_outdir = "$p2wd/gsat";


    #Count all protocol2 hits with high score
    $$p2_hsp_cnt += scalar(@top_prot2_hits);


    #Read the amino acid sequences of the complete proteins or just the/top
    #aligned segments for all protein pairs that were used to run Protocol2
    #Format:
    #  id => [id, seq]....

    my %top_hits_seq_1 = ();
    my %top_hits_seq_2 = ();
    read_prot2_output_faa_seqs($p2wd, \%top_hits_seq_1, \%top_hits_seq_2);
#    print Data::Dumper->Dump([\%top_hits_seq_1, \%top_hits_seq_2], [qw(*top_hits_seq_1 *top_hits_seq_2)]);
#    exit;


    die "Could not extract single-string sequences for proteins"
      unless (%top_hits_seq_1 && %top_hits_seq_2);


    #Now get the sequences for only the protein pairs that had top
    #protocol_2 scores and run GSAT
    #
    #NOTE: @top_prot2_hits contains arrays of top protocol_2 hits
    #      with the following format:
    #      ([hit1_prot1_id, hit1_prot2_id, z-score, aln_segment1, aln_segment2, TMO1],
    #       [hit2_prot3_is, hit2_prot4_id, z-score, aln_segment3, aln_segment4, TMO2]], .... )


    #Create GSAT output directory if it does not exist
    system "mkdir $gsat_outdir" unless (-d $gsat_outdir);


    foreach my $hit (@top_prot2_hits) {

      #Strings to identify each sequence
      my $xref1 = $hit->[0];
      my $xref2 = $hit->[1];


      my ($seq1, $seq2) = undef;


      #Extract sequences to run GSAT
      if ($gmode eq "segment") {
	$seq1 = $hit->[3];
	$seq2 = $hit->[4];
      } else {
	$seq1 = $top_hits_seq_1{$xref1}->[1];
	$seq2 = $top_hits_seq_2{$xref2}->[1];
      }
#      print Data::Dumper->Dump([$xref1, $seq1, $xref2, $seq2 ], [qw( $xref1, $seq1, $xref2, $seq2 )]);
#      exit;


      #Both sequences must be available
      die "Problem extracting faa sequences to run GSAT for $xref1 and $xref2" unless ($seq1 && $seq2);


      #run GSAT between high scoring protocol2 hits
      my $significantHit = run_gsat_and_verify_score($xref1, $seq1, $xref2, $seq2, $gsat_outdir, $gcutoff, $g_hsp_cnt, "p2hit", $shuffles);


      #If there was a signficant GSAT score, run GSAT again to close the transitivity loop. That is
      #the query TCDB protein family against its famXpander hit that generated
      if ($significantHit) {

	#GSAT score passed, now close the loop and run GSAT of the hit against
	#the query family
	if ($p1_prog eq 'fxpand') {

	  #This are the fragments of the original tcdb family sequences and their psiblast hits.
	  my $fxpand_outfile1 = "$p1d/$sbj_fam/psiblast.tbl";
	  my $fxpand_outfile2 = "$p1d/$tgt_fam/psiblast.tbl";


	  #both famXpander output files must exist
	  unless (-f $fxpand_outfile1 && -f $fxpand_outfile2) {
	    die "FamXpander output files do not exist: $fxpand_outfile1, $fxpand_outfile2";
	  }


	  #Extract the sequence of the original TCDB protein that generated the matches for this hit.
	  my $seqs4gsat_sbj = get_fam_seq_for_p2hit($xref1, $fxpand_outfile1);
	  my $seqs4gsat_tgt = get_fam_seq_for_p2hit($xref2, $fxpand_outfile2);
	  unless (@$seqs4gsat_sbj && @$seqs4gsat_tgt) {
	    die "Coud not extract p2hits from files: $fxpand_outfile1 and/or $fxpand_outfile2";
	  }

#	  print Data::Dumper->Dump([ $seqs4gsat_sbj, $seqs4gsat_tgt], [qw(*seqs4gsat_sbj *seqs4gsat_tgt)]);
#	  exit;



	  #--------------------------------------------------------------------------
	  #Run GSAT between the p2hits and their respective TCDB family queries, in order
	  #to complete the transitivity path of homology.

	  #This counter will prevent from counting these gsat runs as good protocol2 hits
	  my $dumb_cnt = 0;
	  foreach my $hitSbj (@$seqs4gsat_sbj) {

#	    print Data::Dumper->Dump([$hitSbj->[0], $hitSbj->[2], $hitSbj->[1], $hitSbj->[3], $gsat_outdir, $gcutoff, \$dumb_cnt, "verify", $shuffles],
#				     [qw(*hitSbj0 *hitSbj2 *hitSbj1 *hitSbj3 *gsat_outdir *gcutff *dumb_cnt *verify *shuffles )]);
#	    <STDIN>;
#	    next;

	    run_gsat_and_verify_score($hitSbj->[0], $hitSbj->[2], $hitSbj->[1], $hitSbj->[3], $gsat_outdir, $gcutoff, \$dumb_cnt, "verify", $shuffles);
	  }

	  foreach my $hitTgt (@$seqs4gsat_tgt) {

#	    print Data::Dumper->Dump([$hitTgt->[0], $hitTgt->[2], $hitTgt->[1], $hitTgt->[3], $gsat_outdir, $gcutoff, \$dumb_cnt, "verify", $shuffles],
#				     [qw(*hitTgt0 *hitTgt2 *hitTgt1 *hitTgt3 *gsat_outdir *gcutff *dumb_cnt *verify *shuffles )]);
#	    <STDIN>;
#	    next;

	    run_gsat_and_verify_score($hitTgt->[0], $hitTgt->[2], $hitTgt->[1], $hitTgt->[3], $gsat_outdir, $gcutoff, \$dumb_cnt, "verify", $shuffles);
	  }
	}


	#Program used to blast family against NR is protocol1.
	else {

	  #
	  #Since protocol1 does not register the fragment of the TCDB protein
	  #That had a blast hit against the NR database, it will be necessary
	  #To blast the sequence of both Protocol2 sequences that yielded a
	  #significant z-core against their original TCDB protein queries.
	  #(still needs to be implemented).

	}
      }
    } #foreach prot2 hit
  } #if outfile
}



#==========================================================================
#Run GSAT between protocol2 hits or between a protocol1 hit and the original
#TCDB query sequence


sub run_gsat_and_verify_score {
  my ($sbj, $seqSbj, $tgt, $seqTgt, $outdir, $gsatCutOFF, $hit_cnt, $mode, $shuffles) = @_;

  my $gsat_outFile = "$outdir/${sbj}.${tgt}.gsat";
  my $pass_outFile = "";
  my $fail_outFile = "";


  #Run GSAT for regular protocol2 hit
  if ($mode eq "p2hit") {
    $pass_outFile = "$outdir/PASSED.${sbj}.${tgt}.gsat";
    $fail_outFile = "$outdir/FAILED.${sbj}.${tgt}.gsat";
  }


  #Run GSAT for protocol1 or famXpander hits to verify the quality of the alignment
  #for the first and last step of the transitivity rule to infer remote homology (A-B and C-D).
  #SRRSW: Super Retro Re-SandWishificador
  elsif ($mode eq "verify") {
    $pass_outFile = "$outdir/SRRSW_PASSED.${sbj}.${tgt}.gsat";
    $fail_outFile = "$outdir/SRRSW_FAILED.${sbj}.${tgt}.gsat";
  }
  else {
    die "Unknow mode to run GSAT: $mode";
  }

  print "Running GSAT for $gsat_outFile\n";
#  print Data::Dumper->Dump([$gsat_outFile, $pass_outFile,  $fail_outFile],
#			   [qw( *gsat_outFile, *pass_outFile  *fail_outFile)]);
#  <STDIN>;



  #Now run GSAT on the top scoring sequences
  #print "Running GSAT:\n   gsat.py $seqSbj $seqTgt 20000 > $gsat_outFile\n";
  unless (-f $pass_outFile || -f $fail_outFile || -f $gsat_outFile) {
#    print "Running GSAT: gsat.py $seqSbj $seqTgt $shuffles > $gsat_outFile\n";
    system "gsat.py $seqSbj $seqTgt $shuffles > $gsat_outFile";
    die "Could not generate GSAT output file: $gsat_outFile --> " unless (-f $gsat_outFile);
  }


  #Extract GSAT score for this protein pair
  my $gsat_score = undef;
  if (-f $gsat_outFile) {
    $gsat_score = get_gsat_score ($gsat_outFile);
  }
  elsif (-f $pass_outFile) {
    $gsat_score = get_gsat_score ($pass_outFile);
  }
  elsif (-f $fail_outFile) {
    $gsat_score = get_gsat_score ($fail_outFile);
  }
  else { die "Error: no gsat file found at: $outdir"; }

  print "   GSAT score: $gsat_score\n\n";


  #Determine if the GSAT score is significant
  if (abs($gsat_score) < $gsatCutOFF) {
    system "mv $gsat_outFile $fail_outFile" unless (-f $fail_outFile);
    return 0;
  }
  else {
    system "mv $gsat_outFile $pass_outFile" unless (-f $pass_outFile);

    #Count the GSAT hit
    ${ $hit_cnt } += 1;

    return 1;

  }
}


#==========================================================================
#Given the NCBI identifier of a protocol2 significant hit, extract
#the secuences of both the query TCDB sequence and the subject blast hit
#from famXpander psiblast.tbl file.

sub get_fam_seq_for_p2hit {
  my ($id, $fxpand_outfile) = @_;

#  print "$fxpand_outfile\n";
#  exit;

  open (my $tblFile, "<", $fxpand_outfile) || die $!;
  my @hits = grep {chomp; $_ =~ /$id/} <$tblFile>;
  close $tblFile;

#  print Data::Dumper->Dump([$id, $fxpand_outfile, \@hits], [qw(*xref *id *fxpand_file *hits )]);
#  <STDIN>;


  #
  #When combining results of famXpander from the precomputed runs
  #in the TSCC, there can be duplicated lines, in terms of the same
  #NCBI subject sequence but with slightly different scores....
  #So, select the top hit to run GSAT
  #


  my %uniqHits = ();

  foreach my $hit (@hits) {

    #sometimes the protocol2 hit has exactly the same ID than TCDB, which brings the header line of the file
    #(e.g. '# Query: 2.A.71.2.3-Q55721') instead of blast matches
    next if ($hit =~ /^#\s+Query/);

    my @data = split (/\t/, $hit);

    my $tcid = $data[0];
    my $eval = $data[4];
    my $qseq = $data[$#data - 1];
    my $sseq = $data[$#data];

    if (exists $uniqHits{$tcid}) {

      #If hit is dupplicated, select the one with the lowest evalue
      if ($eval < $uniqHits{$tcid}->[2]) {
	$uniqHits{$tcid} = [$tcid, $id, $eval, $qseq, $sseq];
      }
    }
    else {
      $uniqHits{$tcid} = [$tcid, $id, $eval, $qseq, $sseq];
    }
  }


  #Get the results ready to return (remove the evalue)
  my @res = ();
  foreach my $h (values %uniqHits) {

    #Return IDs and the respective sequences
    push (@res, [$h->[0], $id, $h->[3], $h->[4]]);
  }

  return \@res;
}





#==========================================================================
#Select the top hits from protocol2. Top hits must have
#
#1) A minimum GSAT score (passed to the program in the variable
#   $gcutoff).
#
#2) The lenth of the alignment between pairs of proteins should be
#   better than a given threshold (e.g. 60). But this is already
#   taken care of by protocol2
#
# These are the input parameters:
# $ar_bestHits:  Variable where best hits will be stored
#
# $tbl_file:  Protocol2 output file in column format.
#
# $mode: Whether gsat will run with full proteins or just
#        on the sequence segments aligned by protocol2.
#
# $p2cut: Protocol2 cutoff to consider a hit significant
#
# $tmo_thr: minimum TMS overlap acceptable between protocol2 hits

sub get_top_proto2_hits {

  my ($ar_bestHits, $rsub, $rtgt, $tbl_file, $mode, $p2lcut, $p2hcut, $tmo_thr, $ignoreProteins) = @_;


#  print Data::Dumper->Dump([$ar_bestHits, $rsub, $rtgt, $tbl_file, $mode, $p2lcut, $p2hcut, $tmo_thr, $ignoreProteins],
#			   [qw(*ar_bestHits, *rsub, *rtgt, *tbl_file, *mode, *p2lcut, *p2hcut, *tmo_thr *ignoreProteins)]);
#  exit;



  open (my $TABFILE, $tbl_file) || die $!;

  while (<$TABFILE>) {

    chomp;

    #
    #Identify the Subject and Target Families as generated by protocol2.
    #NOTE:
    #  This is because because some times protocol2 swaps the subject and
    #  target families passed through the command line.
    #
    if (/^#\s+Subject\:\s+(\S+)\s+Target\:\s+(\S+)/) {
      $$rsub = $1;  $$rtgt = $2;
      next;
    }

    #Ignore header line
    next if (/^Subject/);


    #Retrieve data for each pair of proteins. Minimum length of the alignment
    #Is filtered out by Protocolo_2
    my ($id1, $id2, $ss, $gsat, $aln1, $aln2, $seq1, $seq2, $tmo) = split(/\t/, $_);


#    print Data::Dumper->Dump([$id1, $id2, $ss, $gsat, $tmo], [qw(*id1 *id2 *ss *gsat *tmo )]);
#    <STDIN>;


    #Check if any of the proteins in this protocol2 should be ignored!
    next if ($ignoreProteins->{$id1} || $ignoreProteins->{$id2});


    #Select proteins pairs with minimum GSAT score
    if ($gsat >= $p2lcut && $gsat <= $p2hcut) {

      #keep hits with minimum TMS overlap. Verifying TMO here prevents
      #early abortion of the while loop.
      if ($tmo >= $tmo_thr) {
	if ($mode eq "segment") {
	  push (@{ $ar_bestHits }, [$id1, $id2, $gsat, $seq1, $seq2, $tmo]);
	}
	else {
	  push (@{ $ar_bestHits }, [$id1, $id2, $gsat, $tmo]);
	}
      }
    }

    #Exit file only if protocol2 score is lower than the minimum threshold
    elsif($gsat < $p2lcut)  { last; }
  }

  close $TABFILE;
}



#==========================================================================
#Read the aminoacid sequences of both subjects and targets, this will be
#useful to run GSAT on the top hits of protocol_2
#
# Input parameters:
#
# $wd: Directory with the output of protocol2
#
# $hr_faa1: reference to a hash where the sequences of the subjects
#           will be placed.
#
# $hr_faa2: reference to a hash where the sequences of the targets
#           will be placed.
#

sub read_prot2_output_faa_seqs {

  my ($wd, $hr_faa1, $hr_faa2) = @_;


  #The fasta sequence files to transform to column format
  my $faa1_file = "$wd/subjects.faa";
  my $faa2_file = "$wd/targets.faa";
  die "Fasta file not found or empty: $faa1_file" unless (-f $faa1_file && ! (-z $faa1_file));
  die "Fasta file not found or empty: $faa2_file" unless (-f $faa2_file && ! (-z $faa2_file));


  #The sequence files in column format
  my $tbl1_file = "$wd/subjects.clm";
  my $tbl2_file = "$wd/targets.clm";


  #Convert fasta sequences to tabulated sequences and read
  foreach my $pair ([$faa1_file, $tbl1_file, $hr_faa1], [$faa2_file, $tbl2_file, $hr_faa2]) {

    my ($faa_file, $tbl_file, $out_hash, ) = @{ $pair };

    #tabulate fasta file
    system "columnizeFaa.pl $faa_file $tbl_file";
    die "Tab sequence file does not exist or is empty: $tbl_file"
      unless (-f $tbl_file && ! (-z $tbl_file));


    #Read column-formatted sequence file into a hash table
    open (my $INFILE, $tbl_file) || die $!;
    while (<$INFILE>) {
      chomp;
      next if (/^Id/);

      my ($id, $fheader, $seq) = split(/\t/, $_);
      $out_hash->{$id} = [$id, $seq];
    }
    close $INFILE;
  }
}


#==========================================================================
#Given the ouput of GSAT, parse the line with the z-score that will
#be used to determine whether the result is significant.
#
# Input parameters:
#
# $gsat_output: Path to GSAT output file
#

sub get_gsat_score {

  my $gsat_output = shift;

  my $line = `grep \'Precise\ score\' $gsat_output`;

  if ($line =~ /Precise\s+score\s+\(Z\):\s+([0-9\.+-]+)/) {
    return $1;
  }
  else {
    die "Coult not extract GSAT score from: $gsat_output --> ";
  }
}









#==========================================================================
#Count the number of sequences in a sequence file with either
#fasta or 2-column format. Format is determined based on the
#extension of the sequence file.
#
#Input parameters:
#
# $file: Path to the sequence in fasta or column format.
#

sub count_sequences_in_file {

  my $file = shift;

  #Get the extension of the file
  my @parts = split (/\./, $file);

  my $extension = $parts[$#parts];

  my $nLines = undef;

  #Count the number of lines in file depending on the format
  if (lc $extension eq 'faa') {
    my $line =  `grep \'>\' $file | wc -l`;
    if ($line =~ /\s+(\d+)/) {
      $nLines = $1;
    }
  }
  elsif (lc $extension eq 'clm') {
    my $line =  `wc -l $file`;
    if ($line =~ /\s+(\d+)/) {
      $nLines = $1;
    }
  }
  else {
    die "Unknown file extension: $extension --> ";
  }

  die "Could not count lines for: $file --> " unless ($nLines);

  return $nLines;
}



#==========================================================================
#Count the number of sequences in a file, they can be in either fasta
#or two-column format. Sequences can either have the format generated
#by the program extractTCDB.pl or a custom name provided by the user.
#
#This will help in the calculation in the number of protein pairs that
#will be processed by protocol2
#
#Input paramenters:
#
# $id: Either TCDB ID or custom ID that was used to sequence file.
#      Whatever the ID used, it should be the name (without extension)
#      of the sequence file.
#
# $indir: Directory where the sequence file is located
#
# $program: program that will use the sequences, either protocol1 (proto1)
#           or famXpander (fxpand). This determines the format of the
#           sequence files that will be generated.

sub count_sequences_in_family {

  my ($id, $indir, $program) = @_;

  #possible file names depending on the version of protocol1 that is
  #being run.
  my @file_names = ();


  #Determine the possible input sequence file names
  if ($program eq "fxpand") {
    @file_names = ("$indir/family-${id}.faa", "$indir/${id}.faa");
  }
  elsif ($program eq "proto1") {
    @file_names = ("$indir/family-${id}.clm", "$indir/${id}.clm");
  }
  else { die "Error: unknown program $program --> "; }


  #To indicate whether input file was found or not
  my $found = 0;

  foreach my $file (@file_names) {

    next unless (-f $file );

    if ($program eq "fxpand") {
      my $line =  `grep \'>\' $file | wc -l`;
      if ($line =~ /\s+(\d+)/) {
	return [$file, $1];
      }
    }
    elsif ($program eq "proto1") {
      my $line =  `wc -l $file`;
      if ($line =~ /\s+(\d+)/) {
	return [$file, $1];
      }
    }
  }
  die "Error: input files not found: ".join(", ", @file_names). " --> ";
}






#==========================================================================
#Read a files with protein sequences in two column format: [id, seq].
#Output is an array of arrays with the form ([id1, seq1], [id2, seq2]...]
#
#Input paramenters:
#
# $seqs_file: path to the sequence file in 2-column format (i.e. "id\tseq")
#
# $ar_out: reference to an array where sequences will be stored.
#

sub read_2col_seqs {

    my ($seqs_file, $ar_out) = @_;


    #Verify that file with protein sequences exists
    unless (-f $seqs_file && ! (-z $seqs_file)) {
      die "Sequence file not found or empty: $seqs_file";
    }


    open (my $inseq, $seqs_file) || die "Could not open file: $seqs_file --> $!";
    while (<$inseq>) {
	chomp;

	#Ignore header line if there is one
	next if (/^Id\s+Seq/);

	my ($id, $seq) = split(/\s+/, $_);

	#ignore line if there are no data
	next unless ($id && $seq);

	push (@{ $ar_out }, [$id, $seq]);
    }
    close $inseq;
}


#==========================================================================
#Given a set of TCDB IDs, extract their sequence(s) from TCDB.
#The program calling this function must valide the TCDB IDs.
#
#Input parameters:
#
# $fam_ids: Reference to an array with TCDB IDs.
#
# $outdir: Directory where the extracted sequences will be placed.
#
# $program:  program that will use the sequences, either protocol1 (proto1)
#            or famXpander (fxpand). This determines the format of the
#            sequence files that will be generated.

sub extract_seqs_from_tcdb {

  my ($fam_ids, $outdir, $program) = @_;


  foreach my $fam (@$fam_ids) {

    my ($fname, $fpath, $cmd)  = '';

    #the program is Vamsee's version of protocol1
    if ($program eq "proto1") {
      $fname = "family-${fam}.clm";
      $cmd = "extractTCDB.pl -i $fam -f column -o $outdir";
    }

    #if program is famXpander (Gabo's version of protocol1: famXpander.pl)
    elsif ($program eq "fxpand") {
      $fname = "family-${fam}.faa";
      $cmd = "extractTCDB.pl -i $fam -f fasta -o $outdir";
    }
    else {die "Error: Unknown program: $program --> "; }


    #Do not extract sequences if sequence file already exists
    $fpath  = "$outdir/$fname";
    next if (-f $fpath && ! -z $fpath);


    print "Extracting sequence for $fam\n";
    system $cmd;
  }

}





#==========================================================================
#Validate whether a string or array of strings correspond
#to TCDB ID(s). The function must receive a reference to
#an array.
#TCDB IDs have 5 components with the following format:
#               int.letter.int.int.int
#
#Input parameters:
#
# $tcdb_ids: Reference to an array with the tcdb ids that will be
#            validated
#

sub validate_tcdb_id {

  my ($tcdb_ids) = @_;


  foreach my $id ( @$tcdb_ids ) {

    my @components = split (/\./, $id);

    my $errMsg = "Error: argument is not a TCDB ID: $id\n";

    #ID cannot have more than 5 components
    my $parts = scalar @components;
    die $errMsg if ($parts > 5);


    foreach my $idx (1 .. $parts) {

      #Second TCDB ID component should be a letter.
      if ($idx == 2) {
	die $errMsg unless ($components[1] =~ /^[a-zA-Z]+$/);
      }


      #First, third, fourth and fifth components of the TCDB ID
      #must be integer numbers.
      else {
	die $errMsg unless ($components[$idx-1] =~ /^\d+$/);
      }
    }
  }


  #No error where found!
  return 1;
}




#==========================================================================
#Verify that a list of protein sequences files exist and are not empty
#in a specified directory
#
#Input parameters:
#
# $fileNames: Refernce to an array with the file names that will be
#             validated.
#
# $directory: The directory where the files are located.
#
# $format: The format "fasta" or "column" of the files that will be
#          validated
#

sub validate_files {
  my ($fileNames, $directory, $format) = @_;


  #Determine the right extension of the files
  my $extension = "";
  if ($format eq "fasta") {
    $extension = 'faa';
  }
  elsif ($format eq "column") {
    $extension = 'clm';
  }
  else { die "Unknown sequence file format: $format --> "; }

  #Validate files
  foreach my $file (@$fileNames) {

    my $path = "$directory/${file}.$extension";
    unless (-f $path && ! -z $path) {
      die "File not found: $path -->";
    }
  }
}



1;
