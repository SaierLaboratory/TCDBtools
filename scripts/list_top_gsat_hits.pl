#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;

$Data::Dumper::Deepcopy = 1;

#To read command line options
use Getopt::Long;

#to check dependencies
use R2::CheckDependencies;
use R2::tcdb;



#==========================================================================
#Check dependencies

my @dependencies = ("sort", "grep");
my $CheckDep_obj = new R2::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;




#==========================================================================
#Read command line arguments

my $sort_incrementally = 1;
my $fullTransitivity   = 0;
my $gsat_cutoff        = 15;
my $proto2_dir         = ".";
my $validate           = 0;
my $printBestScore     = 0;
my $validate_dir       = undef;

read_command_line();

#print Data::Dumper->Dump([$sort_incrementally, $fullTransitivity, $proto2_dir],
#                         [qw(*sort_incrementally *fullTransitivity *proto2_dir)]);
#exit;




#==========================================================================
#Determine the right argument for sort

my $sort_arg = "-nr";

$sort_arg = '-n' if ($sort_incrementally);


#==========================================================================
#Get the tcdb ids for both families, in order to discriminate
#proberly A and D families


#Extract the Subject and Target families from the report.tbl file

my $headCmd = qq( head -1  $proto2_dir/*/report.tbl );
chomp (my $cmdOut = `$headCmd`);

my ($family1, $family2) = ($cmdOut =~ /Subject:\s+(\S+)\s+Target:\s+(\S+)/) ? ($1, $2) : die "Could not extract families from report.tbl file";

my @elements1 = split (/\./, $family1);
my @elements2 = split (/\./, $family2);

#print Data::Dumper->Dump([$family1, $family2], [qw(*family1 *family2)]);
#exit;




#==========================================================================
#Get list of top protocol2 hits

#Command to get the list of protocol2 hits that passed the GSAT score
my $grep = qq(find $proto2_dir/*/gsat -name "PASS*" -print | xargs grep -H 'Precise score');
my $perl = q(perl -ne 'chomp; @a=split(/\s+/, $_); @b=split(/\:/,$a[0]); print "$a[3]\t\t$b[0]\n";');
my $sort = "sort $sort_arg";

my $cmd  = "$grep | $perl | $sort";
my $hits = `$cmd`;

#print Data::Dumper->Dump([ $hits], [qw( *hits)]);
#exit;



unless ($fullTransitivity) {
  $hits = "No significant GSAT hits found in dir: $proto2_dir\n" unless ($hits);
  print $hits;
  exit;
}



if ($fullTransitivity || $printBestScore) {

  #Get the resandwished results
  my $grep2 = qq(find $proto2_dir/*/gsat -name "SRRSW_*" -print | xargs grep -H 'Precise score');
  my $cmd2  = "$grep2 | $perl";
  my $srrsw = `$cmd2`;

#  print Data::Dumper->Dump([$srrsw ], [qw( *srrsw )]);
#  exit;


  die "No transitivity gsat results detected for: $proto2_dir" unless ($srrsw);


  #----------------------------------------------------------------------
  #Process the Super-Retro-Re-Sandwishified protocol2 hits that passed
  #THe GSAT threshold.

  my @srrsw_rows = split(/\n/, $srrsw);
  my %srrswGsat = ();

  foreach my $bun (@srrsw_rows) {

    my ($gsat, $file) = split (/\s+/, $bun);
    my ($path, $status, $prots) = split (/(SRRSW_PASSED\.|SRRSW_FAILED\.)/, $file);

    #my ($tcdb, $ncbi) = ($prots =~ /^(\d+\.[a-zA-Z]+\.\d+\.\d+.\d+-\w+)\.(\w+)\.gsat$/)? ($1, $2) : (undef, undef);
    my ($tcdb, $ncbi) = ($prots =~ /^(\d+\.[a-zA-Z]+[0-9\.]+-[a-zA-Z0-9\.]+)\.(\w+)\.gsat$/)? ($1, $2) : (undef, undef);
    die "Could not extract homology transitivity data from: $bun" unless ($tcdb && $ncbi);

    push (@{ $srrswGsat{$ncbi} }, [$tcdb, $ncbi, $gsat]);
  }

#  print Data::Dumper->Dump([\%srrswGsat ], [qw(*srrswGsat )]);
#  <STDIN>;




  #----------------------------------------------------------------------
  #Process results in the order of the protocol2 top hits

  my @p2hitRows = split(/\n/, $hits);
  my @allpassed = ();  #all scores in transitivity path are >= $gsat_cutoff
  my @failed    = ();  #not all scores in transitivity are >= $gsat_cutoff

#  print Data::Dumper->Dump([\@p2hitRows ], [qw(*p2hitRows )]);
#  exit;


 foreach my $p2hit (@p2hitRows) {

    my ($gsat, $file)  = split (/\s+/, $p2hit);
    my ($path, $prots) = split (/PASSED\./, $file);

    my ($query, $subject) = ($prots =~ /(\w+)\.(\w+)\.gsat/)? ($1, $2) : (undef, undef);
    die "Could not extract GSAT data for protocol2 hit: $p2hit" unless ($query && $subject);

#    next if ($query eq $subject);



    #Now get the resandwishing scores for both the query and subject
    die "No resandwish scores for query: $query"     unless (exists $srrswGsat{$query});
    die "No resandwish scores for subject: $subject" unless (exists $srrswGsat{$subject});
    my @query_srrsw = sort srrsw_by_gsat_score @{ $srrswGsat{$query} };
    my @sub_srrsw   = sort srrsw_by_gsat_score @{ $srrswGsat{$subject} };

#    print "-----\n1. ", Data::Dumper->Dump([\@query_srrsw,  \@sub_srrsw], [qw(*query_srrsw  *sub_srrsw)]);
#    print "1. $query eq $subject\n";    <STDIN>;


    #Put here the score corresponding to each family, these are the scores that will be
    #compared to generate the final outout
    my @fam1_scores = ();
    my @fam2_scores = ();



    #----------------------------------------------------------------------
    #If query and subject are the same, remove the other family from the arrays before
    #processing

    if ($query eq $subject) {

      my $continue = 1;


      #This is A, remove resandwishing hits with D
      while ($continue) {
      Q:foreach my $idx (0 .. $#query_srrsw) {

	  my $l = $query_srrsw[$idx];

	  if ($idx == $#query_srrsw) {
	    $continue = 0;
	  }

#	  print Data::Dumper->Dump([\@query_srrsw], [qw(*query_srrsw )]);
#	  print Data::Dumper->Dump([$idx, $#query_srrsw, $continue, $family1, $family2, $l], [qw(*idx *n *continue *family1 *family2 *l)]);
#	  <STDIN>;


	  if (scalar @elements2 < 5) {
	    if ($l && $l->[0] =~ /^$family2\./) {
	      splice (@query_srrsw, $idx, 1);
	      last Q;
	    }
	  }
	  else {
	    if ($l && $l->[0] =~ /^$family2$/) {
	      splice (@query_srrsw, $idx, 1);
	      last Q;
	    }
	  }
	}
      }

#      print Data::Dumper->Dump([\@query_srrsw], [qw(*query_srrsw )]);
#      print "Finished\n";
#      exit;



      #THis is D, remove resandwishing hits with A
      $continue = 1;
      while ($continue) {
      S:foreach my $idx (0 .. $#sub_srrsw) {

	  my $l = $sub_srrsw[$idx];

	  if ($idx == $#sub_srrsw) {
	    $continue = 0;
	  }


	  if (scalar @elements1 < 5) {
	    if ($l && $l->[0] =~ /^$family1\./) {
	      splice (@sub_srrsw, $idx, 1);
	      last S;
	    }
	  }
	  else {
	    if ($l && $l->[0] =~ /^$family1$/) {
	      splice (@sub_srrsw, $idx, 1);
	      last S;
	    }
	  }
	}
      }
    }


#    print "-----\n", Data::Dumper->Dump([\@query_srrsw,  \@sub_srrsw], [qw(*query_srrsw  *sub_srrsw)]);
#    print "2. $query vs $subject\n";
#    <STDIN>;


    #----------------------------------------------------------------------
    #If query is different from subject, and either familie A or D recognizes both B and C, make sure to collect
    #The top score of A-B and C-D for the analysis. Some times it could be possible to that one family has the
    #largest gsat scores(e.g. A-B and A-C > C-D, in this case A-B and C-D must be given)

    foreach my $q (@query_srrsw, @sub_srrsw) {

      if ($q->[1] eq $query  && $q->[0] =~ /^$family1/) {
	push (@fam1_scores, $q->[2]);
      }
      elsif ($q->[1] eq $subject  && $q->[0] =~ /^$family2/) {
	push (@fam2_scores, $q->[2]);
      }
    }


#    print "$query vs $subject\n";
#    print "Check Spliced!\n", Data::Dumper->Dump([\@query_srrsw,  \@sub_srrsw], [qw( *query_srrsw *sub_srrsw)]);
#    print "Check Scores!\n", Data::Dumper->Dump([\@fam1_scores,  \@fam2_scores], [qw( *fam1_scores *fam2_scores)]);
#    exit;




    #sort the transitivity scores ascendingly. This will allow to sort results based
    #upon the lower transitivity score score
    my @scores = sort {$a<=>$b} ($gsat, max(\@fam1_scores), max(\@fam2_scores));

#    print "-----\n",  Data::Dumper->Dump([\@scores ], [qw(*scores )]);
#    print "3. $query vs $subject\n";
#    <STDIN>;


    if ($scores[0] >= $gsat_cutoff && $scores[1] >= $gsat_cutoff && $scores[2] >= $gsat_cutoff) {

#      print "Entre, high enough score!\n";

      #Make sure that there are high resandwishing scores where A != D
      my $scores_ok = verify_resandwishing_scores(\@query_srrsw, \@sub_srrsw, $family1, $family2);
#      print Data::Dumper->Dump([$scores_ok,  $family1, $family2, \@query_srrsw, \@sub_srrsw],
#			       [qw(*scores_ok *family1 *family2 *query_srrsw *sub_srrsw )]);
#      <STDIN>;


      if ($scores_ok) {
	push (@allpassed, {q => $query, s => $subject, gsat=>$gsat, scores => \@scores,
			   q_srrsw => \@query_srrsw, s_srrsw => \@sub_srrsw });
      }
      else {
	push (@failed, {q => $query, s => $subject, gsat=>$gsat, scores => \@scores,
			q_srrsw => \@query_srrsw, s_srrsw => \@sub_srrsw });
      }
    }

    else {
      push (@failed, {q => $query, s => $subject, gsat=>$gsat, scores => \@scores,
		      q_srrsw => \@query_srrsw, s_srrsw => \@sub_srrsw });
    }
  }

#  print "Final!\n", Data::Dumper->Dump([\@allpassed, \@failed], [qw(*allpassed *failed)]);
#  exit;



  #----------------------------------------------------------------------
  #Sort results accoring to user preferences

  my @sorted_passed = sort results_by_user_preferences @allpassed if (@allpassed);
  my @sorted_failed = sort results_by_user_preferences @failed    if (@failed);

#  print Data::Dumper->Dump([\@sorted_passed ], [qw(*sorted_passed )]);
#  exit;



  #----------------------------------------------------------------------
  #print results

  if ($printBestScore) {
    print_best_score(\@sorted_passed, \@sorted_failed);
  }
  else {
    print_full_transitivity(\@sorted_passed, \@sorted_failed) unless ($validate);
  }




  #--------------------------------------------------------------------------
  #Validate results, run protocol 2 for the sequences that matched

  exit unless ($validate);

  #Directory for sequences
  my $seqs_dir = "$validate_dir/seqs";
  system "mkdir -p $seqs_dir" unless (-d $seqs_dir);

  my $q_seq_file = "$seqs_dir/query.faa";
  my $s_seq_file = "$seqs_dir/subject.faa";


  #Get the accesssions for queries and subjects
  my %qids = ();
  my %sids = ();
  foreach my $hit (@sorted_passed) {
    $qids{$hit->{q}} = 1;
    $sids{$hit->{s}} = 1;
  }

  #Print accessions to files
  my $qids_file = "$seqs_dir/qids.txt";
  my $sids_file = "$seqs_dir/sids.txt";
  open (my $qryh, ">", $qids_file) || die "Error in file $qids_file --> $!";
  open (my $subh, ">", $sids_file) || die "Error in file $sids_file --> $!";
  print $qryh join ("\n", sort {$a cmp $b} keys %qids), "\n";
  print $subh join ("\n", sort {$a cmp $b} keys %sids), "\n";
  close $qryh;
  close $subh;



  #extract sequences and format fasta headers for protocol2
  unless (-f $q_seq_file) {

    print "Downloading sequences for $qids_file\n";
    system "blastdbcmd -entry_batch $qids_file -db nr > $q_seq_file";

    print "formating fasta headers for $q_seq_file\n\n";
    my $cmd = q(perl -i.bkp -ne 'if (/^>\w+\|\w+\|\w+\|(\w+)/) {$id=$1;  s/\>/\>$id /; print;}else{print;}') . " $q_seq_file";;
    system $cmd;
  }
  unless (-f $s_seq_file) {

    print "Downloading sequences for $sids_file\n";
    system "blastdbcmd -entry_batch $sids_file -db nr > $s_seq_file";

     print "formating fasta headers for $s_seq_file\n\n";
    my $cmd = q(perl -i.bkp -ne 'if (/^>\w+\|\w+\|\w+\|(\w+)/) {$id=$1;  s/\>/\>$id /; print;}else{print;}' ) . " $s_seq_file";
    system $cmd;
  }


  #runt protocol2
  my ($f1, $f2) = split(/_vs_/,  $proto2_dir);
  my $val_p2d = "$validate_dir/protocol2";
  system "mkdir -p $val_p2d" unless (-d $val_p2d);

  print "Running Protocol2\n";
  my $params = "-s $q_seq_file -t $s_seq_file -o $val_p2d --subject=$f1 --target=$f2 --shuffle 1000 --min=60";
  system "protocol2.py $params";

  print "Done!\n";

}







#==========================================================================
################   Subroutines definition beyond ths point   ##############
#==========================================================================


#==========================================================================
#Get the maximum value within an array

sub max {

  my $aref = shift;

  my @sorted = sort {$b<=>$a} @$aref;
  return $sorted[0];

}


#==========================================================================
#Determine if the top scores are with the reference family, and if not
#there should significant scores with the reference family.

sub verify_resandwishing_scores {

  my ($scores1, $scores2, $fam1, $fam2) = @_;

  my $refam1 = $fam1;
  my $refam2 = $fam2;

  my @ok = (0, 0);


  my $continue = 1;
  while ($continue <= 2) {


    #Verify that there are significant hits with the right families
  F1:foreach my $q (@$scores1) {

      if ($q->[0] =~ /^$refam1/ && $q->[2] >= $gsat_cutoff) {
	$ok[0] = 1;
	last F1;
      }
    }


  F2:foreach my $s (@$scores2) {

      if ($s->[0] =~ /^$refam2/ && $s->[2] >= $gsat_cutoff) {
	$ok[1] = 1;
	last F2;
      }
    }

    if ($ok[0] && $ok[1]) {
      return 1;
    }
    else {
      $refam1 = $fam2;
      $refam2 = $fam1;
      $continue++;
    }
  }

  return 0;
}




sub print_best_score {
  my ($passed, $failed) = @_;


  if (exists $passed->[0]->{scores}) {
    my $gsat = $passed->[0]->{scores}->[0];
    print "$proto2_dir\t\t$gsat\n";
  }
  elsif (exists $failed->[0]->{scores}) {
    my $gsat = $failed->[0]->{scores}->[0];
    print "$proto2_dir\t\t$gsat\n";
  }
  else {
    die "An error occured and could not extract best score";
  }
}


#==========================================================================
#Print transitivity to screen


sub print_full_transitivity {

  my ($passed, $failed) = @_;

  #$passed and $failed are arrays of hashes with the format:
  #{ q => $query, s => $subject, gsat=>$gsat, scores => \@scores,
  #  q_srrsw => \@query_srrsw, s_srrsw => \@sub_srrsw }


  foreach my $hash (@{ $passed }, @{ $failed }) {

    my $query   = $hash->{q};
    my $subject = $hash->{s};
    my $gsat    = $hash->{gsat};

    #Now print the transitivity data
    print "===========================================================================\n";

    foreach my $qsrrsw (@{ $hash->{q_srrsw} }) {
      my ($tcid, $q, $score) = @{ $qsrrsw };
      print "    A-B: $tcid vs $q  ($score)\n";
    }
    print "\n";

    print "HIT B-C: $query vs $subject  ($gsat)\n\n";

    foreach my $ssrrsw (@{ $hash->{s_srrsw} }) {
      my ($tcid, $s, $score) = @{ $ssrrsw };
      print "    C-D: $s vs $tcid  ($score)\n";
    }
    print "\n\n";
  }
}



#===========================================================================
# Sorting functions



sub results_by_user_preferences {

  #Each element is a hash of the form:
  #{ q => $query, s => $subject, gsat=>$gsat, scores => \@scores,
  #  q_srrsw => \@query_srrsw, s_srrsw => \@sub_srrsw }

  if ($sort_incrementally) {
    $a->{scores}->[0] <=> $b->{scores}->[0];
  }
  else {
    $b->{scores}->[0] <=> $a->{scores}->[0];
  }

}



#==========================================================================
#Sort resandwishing alignments decreasingly by gsat score


sub srrsw_by_gsat_score {

  #
  #Gsat score is the third element. arrays have the structure:
  # [tcdb_id, ncbi_accession, gsat_score]).
  #

  $b->[2] <=> $a->[2];

}





#===========================================================================
#Read command line and print help


sub read_command_line {

    print_help() unless (@ARGV);

    my $status = GetOptions(
	"dec|decreasing"         => \&read_sort_dec,
	"inc|increasing"         => \&read_sort_inc,
	"p2d|proto2-dir=s"       => \$proto2_dir,
	"t|gsat-cutoff=f"        => \$gsat_cutoff,
	"full|full-transitivity" => \$fullTransitivity,
	"b|best-score"		 => \$printBestScore,
	"v|validate"             => \$validate,
	"h|help"                 => sub { print_help(); },
	"<>"                     => sub { die "Error: Unknown argument: $_[0]\n"; });
    exit unless ($status);


    die "Error: protocol2 directory must exist\n" unless (-d $proto2_dir);

    #options -b and -full print different info to the screen, thus they are incompatible
    if ($printBestScore && $fullTransitivity) {
      die "-full and -b options are not compatible";
    }
    elsif ($printBestScore) {
      $fullTransitivity   = 1;
      $sort_incrementally = 0;
    }

    if ($printBestScore && $validate) {
      die "options -b and -v are incompatible";
    }

    $validate_dir = "./validate_$proto2_dir";

}



sub read_sort_dec {
    my ($opt, $value) = @_;

    if ($value) {
	$sort_incrementally = 0;
    }
}


sub read_sort_inc {
    my ($opt, $value) = @_;

    if ($value) {
	$sort_incrementally = 1;
    }
}



sub print_help {

    my $help = <<'HELP';

After running program areFamiliesHomologous, use this script to
list all protocol2 hits with high GSAT score.

Input parameters:

-p2d, --proto2-dir {path}
ls   The root directory for all protocol2 results. Default value
   is the current directory.

-inc, --increasing
   List top GSAT hits in increasing order. By default
   results are listed this way (This option is ignored if
   the option -b is given).

-dec, --decreasing
   List top GSAT hits in decreasing order

-b,  --best-score
   Print only the best GSAT score that will represent these
   families. It takes into account all scores A-B,B-C,C-D.

-full --full-transitivity
  List Top GSAT hits including the corrensponding comparisons
  with the original TCDB proteins to complete the transitivity
  path between two TCDB families. NOTE: Use this option only
  if famXpander generated the candidate homologous proteins.
  (optional)

-v, --validate
  If given, this option will extract the full length sequences
  of protein pairs involved in signficant protocol2 and gsat hits. 
  Then runs protocol2 again in these two sets of full-lengthproteins. 
  This is because some times protocol2 runs with the aligned blast 
  regions reported by famXpander and the TMS numbers do not reflect 
  the order of the TMS in the entire protein, just the region that 
  was aligned by psiblast. Therefore, it is necessary to realign
  the sequences so the right TMS numbers are displayed by protocol2.
  Therefore, use this option ONLY when you ran famXpander and kept
  aligned regions instead of full proteins.

-h, --help
   Display this help. Also displayed if script is run without arguments.

HELP

    print $help;
    exit;
}
