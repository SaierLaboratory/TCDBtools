#!/usr/bin/env perl -w

use warnings;
use strict;
use Data::Dumper;

use Getopt::Long;

use TCDB::CheckDependencies;
use Bio::DB::SwissProt;


#==========================================================================
#
# This script reads the psiblast hits produced by script famXpander.pl
# (file psiblast.tbl) and returns those that are most likely distant members
# of a given TCDB family.
#
# NOTE: This programs still needs to be updatated to verify that the
#       caracteristic domains of the family are present in the candidate
#       distant members of the query family.
#==========================================================================


#==========================================================================
#Check dependencies

my @dependencies = ('grep', 'hmmtop', 'blastdbcmd',  'blastp');

my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;


#==========================================================================
#read command line

#Path to the directory with the output of famXpander
my $indir = undef;

#Output directory
my $outdir = ".";

#Path to file psiblast.tbl that program famXpander.pl generates
my $infile = undef;


#File where output results will be placed.
my $outfile = undef;


#Path to TCDB blast DB directory
my $tcbdb = "$ENV{HOME}/db/blastdb";

#Family for which remote homologues will be searched.
my $refFamily = undef;


#Minimal number of TMS expected in acceptable distant homologs
my $minTMS = 1;


#When running TC blast, the target family must be within this number
#of blast matches for a given blast hit.
my $famWithinMatches = 3;


#For queries that require results lower than this evalue
my $lowerThanEvalue  = 1;


#For queries that require results higher than this evalue
my $greaterThanEvalue = 1e-7;


#For queries that require results lower than this evalue
my $lowerThanIdentity  = undef;

#For queries that require results higher than this evalue
my $greaterThanIdentity = undef;


#Select hits that have a length below this threshold
my $lowerThanTgtLength = undef;


#Select hits that that have a length above this threshold
my $greaterThanTgtLength = undef;

#Length ratio between redundant distant homologs. This helps
#to determine when to keep two redundant distant homologs if
#they have very different lengths
my $lengthRatio = 1.8;

#Evalue threshold to consider to remote homologs redundat.
my $redundantEvalue = 1e-5;


#How will the results be sorted (asc, desc, no)
my $sort = 'asc';

read_command_line_arguments();

#print Data::Dumper->Dump([$indir, $infile, $outfile,  $refFamily, $minTMS, $famWithinMatches, $lowerThanEvalue,
# 			  $greaterThanEvalue,, $lowerThanIdentity, $greaterThanIdentity,
#                           $lowerThanTgtLength, $greaterThanTgtLength,  $lengthRatio, $redundantEvalue, $sort],
# 			 [qw(*indir *infile *outfile *refFamily *minTMS *famWithinMatches *lowerThanEvalue
# 			     *greaterThanEvalue *lowerThanIdentity *greaterThanIdentity
# 		             *lowerThanTgtLength *greaterThanTgtLength *lengthRatio *redundantEvalue *sort)]);
#exit;



#put non-redundant proteins in a file1
my $fxpand_seqsFile = "$indir/results.faa";
my $nrHitsIDs       = "$outdir/nr_hits.txt";

system qq(grep '>' $fxpand_seqsFile | perl -ne 's/\\>//; \@a=split(/\\s+/,\$_); print "\$a[0]\\n"' > $nrHitsIDs) unless (-f $nrHitsIDs);
die "File does not exist or is empty: $nrHitsIDs" unless (-f $nrHitsIDs && ! (-z $nrHitsIDs));



#get the psiblast hits for the non-redundant sequences
my $fxpand_nrBlastData = "$outdir/nr_blast_hits_data.txt";
system qq(grep -f $nrHitsIDs $infile > $fxpand_nrBlastData) unless (-f $fxpand_nrBlastData);
die "File does not exist or is empty: $fxpand_nrBlastData"  unless (-f $fxpand_nrBlastData  && ! (-z $fxpand_nrBlastData ));
#system "rm $nrHitsIDs";




#==========================================================================
#Get the psiblast hits: There can be repetitions, so keep only the top
#blast results.


open (my $nrBlastData, $fxpand_nrBlastData) || die $!;


my %top_nr_hits   = ();
my %discard_hits  = ();

my $countHits = 0;

HIT:while (<$nrBlastData>) {

  chomp;

  #Discard comment lines
  next if (/^#/);

  #Contents of array @hit:
  #0:  TCDB ID
  #1:  refseqs,gb,sp..etc.g (e.g. AFI28417)
  #2.  NR psiblast hit (e.g. gi|384931739|gb|AFI28417.1|)
  #3.  Bit score
  #4.  evalue
  #5.  % identity
  #6.  tcdb query start pos
  #7.  tcdb query end  pos
  #8.  tcdb query length
  #9.   NR target start pos
  #10.  NR target end  pos
  #11.  NR target length
  #12.  Genomes for target (including redundant genomes).
  #13.  Aligned fragment of TCDB family
  #14.  Aligned fragment of NR blast match

  my @hit = split(/\t/, $_);



  #There must be data to process
  next HIT unless ($hit[4]);


  #Keep only the top hits for each psiblast hits
  if (exists $top_nr_hits{$hit[1]}) {

    #If hit exists, only overwrite it if the stored hit has worse e-value.
    if ($hit[4] <= $lowerThanEvalue && $hit[4] >= $greaterThanEvalue) {
      if ($top_nr_hits{$hit[1]}->[2] > $hit[4]) {
	$top_nr_hits{$hit[1]} = [$hit[0], $hit[1], $hit[4], $hit[5], $hit[11], $hit[12]];
      }
    }

    #hit should be discarded, but make sure to keep the lowest E-value
    elsif ($hit[4] < $greaterThanEvalue) {
      if (exists $discard_hits{$hit[1]}) {
	if ($discard_hits{$hit[1]} > $hit[4]) {
	  $discard_hits{$hit[1]} = $hit[4];
	}
      }
      else {
	$discard_hits{$hit[1]} = $hit[4];
      }
    }
  }

  elsif ($hit[4] <= $lowerThanEvalue && $hit[4] >= $greaterThanEvalue) {

    #hit has not been recorded yet
    $top_nr_hits{$hit[1]} = [$hit[0], $hit[1], $hit[4], $hit[5], $hit[11], $hit[12]];
    $countHits++;
  }


  #hit should be discarded, but make sure to keep the lowest E-value. This will
  #Allow to skip hits in the next step
  elsif ($hit[4] < $greaterThanEvalue) {

    if (exists $discard_hits{$hit[1]}) {
      if ($discard_hits{$hit[1]} > $hit[4]) {
	$discard_hits{$hit[1]} = $hit[4];
      }
    }
    else {
      $discard_hits{$hit[1]} = $hit[4];
    }
  }
}
close $nrBlastData;

#print Data::Dumper->Dump([\%discard_hits],["*discard_hits"]);
#print Data::Dumper->Dump([\%top_nr_hits],["*top_nr_hits"]);
#exit;

#system "rm ./$fxpand_nrBlastDat";
die "No hits found in file: $fxpand_nrBlastData" unless (%top_nr_hits);




#==========================================================================
#Now remove hits with better matches with other tcdb proteins.
#
#Get the GIs of the hits that pass the filter to extract their sequences
#and blast them against TCDB.

my %results = ();

$countHits = 0;

foreach my $hit (values %top_nr_hits) {

  unless (exists $discard_hits{$hit->[1]}) {

    if (is_identity_in_range($hit->[3]) && is_length_in_range($hit->[4])) {

#      print Data::Dumper->Dump([$hit->[1], $hit], [qw( *id *hit )]);
#      <STDIN>;

      $results{$hit->[1]} = $hit;
      $countHits++;
    }
  }
}

#print Data::Dumper->Dump([ \%results], [qw( *results)]);
#print "famXpander candidate hits: ", scalar keys %results, "\n";
#exit;

die "Afters filtering by identity and length, no hits found in file: $fxpand_nrBlastData" unless ($countHits);




#==========================================================================
#Print the filtered blast results for the first round of candidate
#remote homologs. That is, psiblast of the target family against NR

my $step1_remote_homologs = "$outdir/step1_candidates_from_blast.txt";
open (my $step1h, ">", $step1_remote_homologs) || die $!;

print $step1h "#tcdb_id\tnr_ncbi_id\tEvalue\tSequence_identity\tlength_of_matched_protein\tGenomes_with_redundant_hit\n";
map { my $str = join ("\t", @{ $_ }); print $step1h "$str\n"; } values %results;

close $step1h;

print STDERR "STEP 1. psiblast hits were filtered: $step1_remote_homologs\n";





#==========================================================================
#print NCBI IDs to a file

my $blast_hits_ids_file = "$outdir/ids_remote_homologs.txt";

unless (-f $blast_hits_ids_file) {
  open (my $gih, ">", $blast_hits_ids_file) || die $!;
  map { print $gih "$_\n" } sort { $a cmp $b } keys %results;
  close $gih;
}

print STDERR "... NCBI IDs of remote homologs were saved: $blast_hits_ids_file \n";





#==========================================================================
#Extract the sequences of candidate remote homologs

my $blast_hits_seqs_file = "$outdir/seqs_remote_homologs.faa";
my $get_seq  = qq(blastdbcmd -db nr -entry_batch $blast_hits_ids_file -target_only -out $blast_hits_seqs_file);
system $get_seq unless (-f $blast_hits_seqs_file);

die "No sequences were retrieved for blast hits: $blast_hits_seqs_file" unless (-f $blast_hits_seqs_file);

print STDERR "... Sequences of remote homologs were extracted: $blast_hits_seqs_file\n";





#==========================================================================
#Run hmmtop in the sequences of candidate distant homologs


my $hmmtop_top_hits = "$outdir/candidate_homologs_tms.txt";

system "hmmtop -if=$blast_hits_seqs_file -of=$hmmtop_top_hits -is=pseudo -pi=spred -sf=FAS" unless (-f $hmmtop_top_hits);

die "No hmmtop output file found:  $hmmtop_top_hits" unless (-f $hmmtop_top_hits);


print STDERR "... Ran hmmtop on top blast hits: $hmmtop_top_hits\n";




#==========================================================================
#Parse hmmtop to get the number of TMS per GI

my %ids2tms = ();
open (my $tmsh, "<", $hmmtop_top_hits) || die $1;
while(<$tmsh>) {
  chomp;

  my ($len, $id, $ntms) = (/\s+(\d+)\s+(\w+).+\s(IN|OUT)\s+(\d+)/)? ($1, $2, $4) : ();

  if ($len && $id && $ntms) {
    $ids2tms{$id} = $ntms;
  }

}
close $tmsh;

#print Data::Dumper->Dump([ \%ids2tms ], [qw( *ids2tms )]);
#print "Total sequences with TMS predictions: ", scalar keys %ids2tms, "\n";
#exit;

print STDERR "... Parsed HHMTOP output\n";




#==========================================================================
#Download TCDB proteins and make a Blast database

my $tcbdb_name = "$tcbdb/tcdb";
system "extractFamily.pl -i tcdb -o  $tcbdb -f blast" unless (-f "${tcbdb_name}.phr");
die "Could not download TCDB blast DB in $tcbdb" unless (-f "${tcbdb_name}.phr");

print STDERR "\nSTEP 2: Blasting NR matches against TCDB\n";
print STDERR "... TCDB blast database created\n";




#==========================================================================
#Now blast the remote candidates homologs


my $blast_hits_clm_file = "$outdir/blast_against_tcdb.clm";
my $fmtBlastOut = '6 qseqid sseqid qlen slen evalue pident nident length qcovs qstart qend sstart send';

unless (-f $blast_hits_clm_file) {
  system "blastp -db $tcbdb_name -evalue 0.1 -comp_based_stats f -use_sw_tback -outfmt '$fmtBlastOut' -query $blast_hits_seqs_file -out $blast_hits_clm_file";
}


print STDERR "... Blasted candidate distant homologs against TCDB: $blast_hits_clm_file\n";


#==========================================================================
#Now discard all matches that involve lower evalues with families other
#than the reference family


my %top_tcdb_hits = ();
my %ignore_hits   = ();


open (my $blasth, "<", $blast_hits_clm_file) || die $!;
HIT:while(<$blasth>) {

  chomp;
  my ($ncbi, $tcdb_id, $qlen, $slen, $ev, @kk) = split(/\t/, $_);


  #There must be data to process
  next HIT unless ($tcdb_id && $ev);


  #extract the RefSeq in order to get the number of TMS
  my $ncbi_id = ($ncbi =~ /(\w+)/)? $1 : undef;
  die "Could not get refseq from: $_\n" unless ($ncbi_id);

#  print Data::Dumper->Dump([$ncbi_id, $_ ], [qw(*refseq *line )]);
#  <STDIN>;

  #Get the number of TMS in this gene and verify it has the minimum number of TMS
  next HIT unless (exists $ids2tms{ $ncbi_id } && $ids2tms{ $ncbi_id } >= $minTMS);


  #Keep only the top hits for each psiblast hit
  if (exists $top_tcdb_hits{ $ncbi_id } && exists $top_tcdb_hits{ $ncbi_id }{ $tcdb_id }) {

    #If hit exists, only overwrite it if the stored hit has worse e-value.
    if ($ev >= $greaterThanEvalue) {
      if ($top_tcdb_hits{ $ncbi_id }{ $tcdb_id }->[2] > $ev) {
	$top_tcdb_hits{ $ncbi_id }{ $tcdb_id } = [$ncbi_id, $tcdb_id, $ev, $qlen, $slen];
      }
    }

    #hit should be discarded, but make sure to keep the lowest E-value
    elsif ($ev < $greaterThanEvalue) {
      if (exists $ignore_hits{ $ncbi_id } && exists $ignore_hits{ $ncbi_id }{ $tcdb_id }) {
	if ($ignore_hits{ $ncbi_id }{ $tcdb_id } > $ev) {
	  $ignore_hits{ $ncbi_id }{ $tcdb_id } = $ev;
	}
      }
      else {
	$ignore_hits{ $ncbi_id }{ $tcdb_id } = $ev;
      }
    }
  }



  elsif ($ev >= $greaterThanEvalue  ) {

    #hit has not been recorded yet
    $top_tcdb_hits{ $ncbi_id }{ $tcdb_id } = [$ncbi_id, $tcdb_id, $ev, $qlen, $slen];
  }


  #Hit should be discarded, but make sure to keep the lowest E-value. This will
  #Allow to skip hits in the next iteration
  elsif ($ev < $greaterThanEvalue) {

    if (exists $ignore_hits{ $ncbi_id } && exists $ignore_hits{ $ncbi_id }{ $tcdb_id }) {
      if ($ignore_hits{ $ncbi_id }{ $tcdb_id } > $ev) {
	$ignore_hits{ $ncbi_id }{ $tcdb_id } = $ev;
      }
    }
    else {
      $ignore_hits{ $ncbi_id }{ $tcdb_id } = $ev;
    }
  }
}
close $blasth;

#print Data::Dumper->Dump([\%ignore_hits ], [qw(*ignore_hits)]);
#exit;

#print Data::Dumper->Dump([\%top_tcdb_hits ], [qw( *top_tcdb_hits )]);
#exit;

unless (%top_tcdb_hits) {
  print "Did not find useful tcdb blast hits in file: $blast_hits_clm_file\n";
  exit;
}

print STDERR "... Parsed top blast hits of remote homologs versus TCDB\n";






#==========================================================================
#Now make sure that only proteins with matches to the reference family are
#selected, and that the hits with the reference family have the best
#Evalues


#Put here the remote homologs for the family under consideration
my @remote_homologs = ();


ID:foreach my $id (keys %top_tcdb_hits) {

  #Flag indicates whether this ID had a blast match with the target family.
  #(ID will be ignored if it doesn't)
  my $has_blast_w_fam = 0;


  #The blast matches of candidate distant homologs
  my @blast_data = ();


 TC:foreach my $tcdb_id (keys %{ $top_tcdb_hits{$id} }) {

    #Skip protein if it has a higher $evalue with another family
    next TC if (exists $ignore_hits{ $id });

    #Determine here if there was a match with the target family
    $has_blast_w_fam = 1 if ($tcdb_id =~ /$refFamily/);

    push(@blast_data, $top_tcdb_hits{ $id }{ $tcdb_id });
  }

#  print Data::Dumper->Dump([\@blast_data,  $has_blast_w_fam, $id], [qw( *blast_data *has_blast_w_fam *id)]);
#  <STDIN>;


  #Ignore ID if it didn't had a blast match with target family.
  next ID unless ($has_blast_w_fam);


  #There must be blast matches to proceed
  next ID unless (@blast_data);


  #Sort the data by evalue and make sure that the target family is within the
  #top matches
  my @sorted_blast = sort {$a->[2] <=> $b->[2]} @blast_data;


  #Verify here that the target family is within the top matches as specified by
  #variable $famWithinMatches
  my $match_cnt = 0;

 MATCH:foreach my $match (@sorted_blast) {

    if ($match->[1] =~ /$refFamily/ && $match_cnt <= $famWithinMatches) {

      #Collect match and go to the next id.
      push(@remote_homologs, $match);
      last MATCH;
    }
    else {
      $match_cnt++;
      last MATCH if ($match_cnt > $famWithinMatches);
    }
  }
}


#print Data::Dumper->Dump([\@remote_homologs ], [qw( *remote_homologs)]);
#exit;

die "After filtering TCDB blast hits by E-value, could not find remote homologs." unless (@remote_homologs);

#print Data::Dumper->Dump([\@remote_homologs ], [qw(*remote_homologs )]);
#exit;

print STDERR "... Kept only the the with worst E-value against TCDB\n";





#==========================================================================
#Candidate remote homologs may contain redundant genes. A candidate
#is redundant with another if the evalue between both genes is better
#than 1e-5
#
#To remove redundant sequences from the final results, I tried usearch:
#
#--Clustering of sequences:
#usearch -cluster_agg best_candidates_ids.faa -treeout tree.phy -clusterout clusters.txt -id 0.3 -linkage max
#
#--Blasting sequences, and generating ID matrix
#usearch -search_local best_candidates_ids.faa -db best_candidates_ids.faa -evalue 1e-6 -userout usearch_local.out -userfields query+target+evalue+id+ql+tl




#extract the IDs of the candidate to extract their sequences
my $idsCandidatesFile = "$outdir/idsCandidates.txt";

open (my $finalIds, ">", $idsCandidatesFile) || die $!;
map {print $finalIds $_->[0], "\n";} @remote_homologs;
close $finalIds;



#Now extract the sequences of the candidate remote homologs
my $candidates_seqs_file = "$outdir/candidates_seqs.faa";
my $getCandIDsSeqsCmd  = qq(blastdbcmd -db nr -entry_batch  $idsCandidatesFile -target_only -out  $candidates_seqs_file);
system $getCandIDsSeqsCmd;

die "No sequences were retrieved for candidate homologs: $candidates_seqs_file" unless (-f $candidates_seqs_file);



#Now get the number of sequences in TCDB in order to calculate the right blast evalue
my $tcdb_size = undef;
my $strOut = `blastdbcmd -db $tcbdb_name -info`;
if ($strOut =~ /\s+[0-9,]+\s+sequences;\s+([0-9,]+)\s+total/) {
  $tcdb_size = $1;
  $tcdb_size =~ s/,//g; #remove commas from number
}
else { die "Could not extract size of TCDB database"; }

#print Data::Dumper->Dump([$tcdb_size], [qw(*tcdb_size )]);
#exit;




#Blast the candidates to remote homologs against themselves
my $cand_blast_clm_file = "$outdir/candidates_selfblast.clm";
$fmtBlastOut = '6 qseqid sseqid evalue qlen slen pident';

#print  "blastp -query $candidates_seqs_file -subject $candidates_seqs_file -dbsize $tcdb_size -evalue $redundantEvalue -comp_based_stats f -outfmt '$fmtBlastOut' -out $cand_blast_clm_file\n";
system "blastp -query $candidates_seqs_file -subject $candidates_seqs_file -dbsize $tcdb_size -evalue $redundantEvalue -comp_based_stats f -outfmt '$fmtBlastOut' -out $cand_blast_clm_file";




#Identify the redundant candidates. Candidates with very different size
#are allowed. When 2 redundant candidates have similar size, remove the
#the smallest one.
my %redundantCandidates = ();
open (my $cndBlast, "<", $cand_blast_clm_file) || die;
while (<$cndBlast>) {

  chomp;  next if (/^#/);
  my ($qry, $sbj, $eval, $ql, $sl, $id) = split(/\s+/, $_);

  #retrieve the refSeq for both query and subject
  my $qid = ($qry =~ /(\w+)/)? $1 : undef;
  my $sid = ($sbj =~ /(\w+)/)? $1 : undef;


  #THere must be IDs to continue
  die "Could not extract id for query: $qry" unless ($qid);
  die "Could not extract id for query: $sbj" unless ($sid);


  #Remove self comparisons
  next if ($qid eq $sid);


  #When there are redundant hits with very different size keep both.
  my $lenRatio1 = $ql/$sl;
  my $lenRatio2 = $sl/$ql;
  next if ($lenRatio1 >= $lengthRatio || $lenRatio2 >= $lengthRatio);


  #When size is similar, the smallest protein is redundant
  if ($ql > $sl) {
    $redundantCandidates{$sid} = 1;
  }
  else {
    $redundantCandidates{$qid} = 1;
  }

}
close $cndBlast;

#print Data::Dumper->Dump([ \%redundantCandidates], [qw(*redundantCandidates )]);
#exit;




#==========================================================================
#Generate list with non-redundant candidate remote homologs

my @final_candidates = ();
my @rsCandidates     = ();
foreach my $homolog (@remote_homologs) {
  unless (exists $redundantCandidates{$homolog->[0]}) {
    push (@final_candidates, $homolog);
    push (@rsCandidates, $homolog->[0]);
  }
}

#print Data::Dumper->Dump([ \@final_candidates, \@rsCandidates ], [qw( *final_candidates *rsCandidates)]);
#exit;


print STDERR "... Removed redundant distant homologs\n";




#==========================================================================
#Get the UniProt IDs for as many candidate homologs as possible

my $spObj = Bio::DB::SwissProt->new();
my $map1 = $spObj->id_mapper(-from => 'P_REFSEQ_AC',
			     -to   => 'ACC',
			     -ids  => \@rsCandidates);


my $map2 = $spObj->id_mapper(-from => 'EMBL',
			     -to   => 'ACC',
			     -ids  => \@rsCandidates);


#print Data::Dumper->Dump([ $map1, $map2], [qw( *rs2sp *embl2sp )]);
#exit;


#generate a single mapping hash
my %ids_map = ();
if (($map1 && %{ $map1 }) && ($map2 && %{ $map2 })) {
  %ids_map = (%$map1, %$map2);
}
elsif ($map1 && %{ $map1 }) {
  %ids_map = %$map1;
}
elsif ($map2 && %{ $map2 }) {
  %ids_map = %$map2;
}



#replace ids of homolog by uniprot if possible
foreach my $homolog (@final_candidates) {
  if (exists $ids_map{$homolog->[0]}) {
    $homolog->[0] =  $homolog->[0] . "," . join (",", @{ $ids_map{$homolog->[0]} });
  }
}

#print Data::Dumper->Dump([ \@final_candidates ], [qw( *final_candidates)]);
#exit;


print STDERR "... Retrieved as many uniprot IDs as possible\n";



#==========================================================================
#Print output file with final candidates for remote homologs


open (my $step2h, ">", $outfile) || die $!;

if (@final_candidates) {
  print $step2h "#ncbi_id\ttcdb_id\tE-value\tLength_of_query\tlength_of_subject\n";
  if ($sort eq 'no') {
    map { my $str = join ("\t", @{ $_ }); print $step2h "$str\n"; } @final_candidates;
  }
  elsif ($sort eq 'asc') {
    map { my $str = join ("\t", @{ $_ }); print $step2h "$str\n"; } sort asc_by_evalue @final_candidates;
  }
  else {
    map { my $str = join ("\t", @{ $_ }); print $step2h "$str\n"; } sort {$b->[2]<=>$a->[2] } @final_candidates;
  }
}
else {
  print $step2h "No remote homologs found four family $refFamily\n";
}
close $step2h;



print STDERR "FINISHED: List of remote homologs was generated: $outfile\n";



#==========================================================================
#sort results ascending by evalue and descending by identity and length

sub asc_by_evalue {

  #if evalue is equal and identity different
  if ($a->[2] == $b->[2] && $a->[3] != $b->[3]) {
    $b->[3]<=>$a->[3];  #Descending by identity
  }

  #if evalue is equal and identity is equal
  elsif ($a->[2] == $b->[2] && $b->[3] == $a->[3]) {
    $b->[4]<=>$a->[4];  #Descending by length
  }

  else {
    $a->[2]<=>$b->[2]; #Ascending by evalue
  }
}




#==========================================================================
#Check if the target's length is within the range defined by
#variables: $lowerThanTgtLength, $greaterThanTgtLength

sub is_length_in_range {

  my $length = shift;

  if (defined $lowerThanTgtLength && defined $greaterThanTgtLength) {

    if ($length >= $greaterThanTgtLength && $length <= $lowerThanTgtLength) {
      return 1;
    }
    else {
      return 0;
    }
  }
  elsif (defined $greaterThanTgtLength) {

    if ($length >= $greaterThanTgtLength) {
      return 1;
    }
    else {
      return 0;
    }
  }
  elsif (defined $lowerThanTgtLength) {

    if ($length <= $lowerThanTgtLength) {
      return 1;
    }
    else {
      return 0;
    }

  }
  else {

    #Do not take into account length
    return 1;
  }

}






#==========================================================================
#Determine whether the identity is withing the range defined by the 
#variables $lowerThanIdentity and  $lowerThanIdentity

sub is_identity_in_range {

  my $id = shift;

  #Both varaibles must have values to proceed, otherwise ignore identity info.
  return 1 unless (defined($greaterThanIdentity) && defined($lowerThanIdentity));

  if ($id >= $greaterThanIdentity && $id <= $lowerThanIdentity) {
    return 1;
  }
  else {
    return 0;
  }
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
      "d|indir=s"              => \$indir,
      "o|outfile=s"            => \$outfileTmp,
      "od|outdir=s"            => \$outdir,
      "tcbdb|tc-blast-db-dir"  => \$tcbdb,
      "f|family=s"             => \$refFamily,
      "n|min-tms=s"            => \$minTMS,
      "m|fam-tcblast-hits=i"   => \$famWithinMatches,
      "lt|lower-than=f"        => \$lowerThanEvalue,
      "gt|greater-than=f"      => \$greaterThanEvalue,
      "ilt|id-lower-than=f"    => \$lowerThanIdentity,
      "igt|id-greater-than=f"  => \$greaterThanIdentity,
      "llt|len-lower-than=i"   => \$lowerThanTgtLength,
      "lgt|len-greater-than=i" => \$greaterThanTgtLength,
      "r|len-ratio=f"          => \$lengthRatio,
      "re|redundant-evalue=f"  => \$redundantEvalue,
      "s|sort=s"               => \&sort_results,
      "h|help"                 => sub { print_help(); },
      "<>"                     => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);


  #Validate input directory
  die "Option -d is mandatory" unless ($indir && -d $indir);

  system "mkdir -p $outdir" unless (-d $outdir);

  die "TCDB blast DB directory does not exist" unless (-d $tcbdb);


  #Validate input file
  $infile = "$indir/psiblast.tbl";
  die "famXpander output file psiblast.tbl was not found in directory: $indir" unless (-f $infile);

  #Give correct path to output file
  die "Option -o is mandatory" unless ($outfileTmp);
  $outfile = "$outdir/$outfileTmp";

  die "Reference family is mandatory" unless ($refFamily);
}




#Read the sorting option -s
sub sort_results {

  my ($option, $value) = @_;

  my $val = lc $value;

  unless ($val eq 'asc' || $val eq 'desc' || $val eq 'no') {
    die "Unknown sorting option: $value\n  Acceptable values are: asc, desc, no\n";
  }

  $sort = $val;

}




sub print_help {

  #
  # $errMsg: The error message to be diplayed
  #
  # $printHelp: boolean value indicating whether the help of the
  #             program will be displayed.
  #

  my $help = <<'HELP';

 This script reads the psiblast hits produced by script famXpander.pl
 (file psiblast.tbl) and returns those that are most likely distant members
 of a given TCDB family.


 Input paramateres

 -d, --indir {path}
    Path to the directory where the ouput files of famXpander are located
    (i.e. psiblast.tbl and results.faa). Argument is Mandatory.

 -o, --outfile {path}
    Ouput file where results should be saved. (Mandatory)

 -od, --outdir {path}
    Output directory for temporary and final output files.
    (default: ./)

 -f, --family {string}
    Family for which remote homologs will be found. This allows to
    identify TCDB blast matches with the right family.
    (Mandatory)

 -tcbdb, --tc-blast-db-dir {path}
    Path to the blast DB with all proteins in tcdb. This database needs 
    to be created with program extractFamily and be called 'tcdb'.
    If the blast DB is not found in this location, it will be created
    for you.
    (Default ~/db/blastdb);

 -n, --min-tms-num {integer}
    Minimum number of TMS that should be found in distant homologs.
    (Default 1).

 -m, fam-tcblast-hits {integer}
    When blasting against TCDB the reference family, passed with option
    -f, should be wihtin the specified top number hits. (Default 3).

 -gt, --greater-than {float}
    Lower E-value threshold to parse famXpander output.  Hits with E-value 
    greater than this value will be extracted. (default 1e-7).

 -lt, --lower-than {path}
    Higher threshold to parse famXpander output. Hits lower than this 
    E-value will be extracted. (default 1.0)

 -ilt --id-lower-than {float}
    Identity should be lower than this value. This parameter can
    take any positive value lower than 100. By default identity
    is not taken into consideration.

 -igt --id-greater-than {float}
    Identity should be greter than this value. This parameter can
    take any positive value lower than 100. By default identity
    is not taken into consideration.

  -llt --len-lower-than {int}
    Length of subject should be lower than this value. By default
    Length is not taken into consideration.

  -lgt --len-greater-than {int}
    Length of subject should be greater than this value. By default
    Length is not taken into consideration.

  -r  --len-ratio {float}
    Length ratio between redundant distant homologs. This helps
    to determine when to keep two redundant distant homologs if
    they have very different lengths (default 1.8)

  -re  --redundant-evalue
    Evalue threshold to consider to remote homologs redundat.
    (default 1e-5).

 -s, --sort
    Sort blast results. Results can be sorted in ascending (asc),
    descending (desc) or no order (no). 
    (default is asc).

 NOTE: It's possible to specify a range of E-values by providing
       both -gt and -lt thresholds.

HELP

  print $help;
  exit;
}
