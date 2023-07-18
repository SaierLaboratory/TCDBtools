#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(max min);
use Bio::SeqIO;
use Bio::SearchIO;

use TCDB::Assorted;

use Data::Dumper;
$Data::Dumper::Deepcopy = 1;

###########################################################################
#
# In genome analysis, after running GBLAST there will be multicomponent
# systems with missing components. This is approached by running blastp
# with the TCDB proteins of the missing components against the genome
# of interest, and by checking whether proteins matching other components
# within the same system are located nearby in the genome.
#
# This script will run blastp/ssearch36 of a series of TCDB proteins in
# multicomponent systems. Top hits that will be selected must meet the
# following criteria to be called "Candidates":
#
# 1) id% >= 24%
# 2) coverage of at least 45% of either the query or the subject
# 3) If the smaller protein is less than 100aa long, the coverage muust
#    be at least 70% of the smallest protein
# 4) E-value < 1.0
#
# After the top hits are filtered, and the genome neighborhood is determined
# for each protein, then the program WHAT (quod.py) is run to generate plots
# for the alignment and the full sequence of the query (TCDB protein) and
# the subject (protein in the reference genome).
#
# BUG Report:
# When one smaller protein has multiple hits with another because the
# larger protein is a repeat of the smaller the program only shows the
# plots for the last match. This needs to be fixed so that there is a plot
# for every match.
#--------------------------------------------------------------------------
#
# Writen by: Arturo Medrano
# Created:  08/2017
# Modified: 08/2019
#
###########################################################################


#
#Improvements for next update:
#
#   * Incorporate Scaffold IDs into the results to make the script compatible
#     with genomic context searches in incomplete  metagenomes!!
#
#   * Make sure genomic context brings hits with all other components in the
#     same system. For an example check Desulfosarcina_variabilis
#     a) System 2.A.63.1.3 is missing from the output table
#     b) 2.A.63.1.4 should hit component O05228 with Dvar_54160
#     c) The genomic context for 2.A.63.1.4 is not showing all the genomic
#        proteins that match other components.
#


#==========================================================================
#Command line arguments

my $file_tcdb_queries = undef;     #The TCIDs of multicomponent systems

#input files and overall controls
my $gnmAcc          = undef;       #This is prefix to identify different files in the genome (e.g. proteome, GFF, GenBank, etc.)
my $gnmDir          = undef;       #Assumed to be the NCBI dir for the genomic data
my $repForm         = 'circular';  #or linear
my $blastdb_name    = "proteome";  #name that will be used to name the genome's blast DB
my $file_proteome   = undef;
my $gnmFeatTable    = undef;
my $gblastFile      = undef;
my $blastOverWrite  = 0;            #Recreate the blastDB of all TCDB proteins
my $outFmt          = "tsv";        #can also be "csv"
my $prog            = "ssearch36";  #Alignment program to use (blastp|ssearch36)
my $subMatrix       = 'BL50';       #Substitution matrix fro ssearch36
my $outdir          = "./MultiComponentSystemsAnalysis";
my $tcdbFaa         = undef;
my $blastBinDir     = '/usr/local/bin';

#Variables that will be used to select good blastp hits (Candidates)
my $geneDis         = 20;    #How many genes away to look for neighbors
my $fusion_lenRatio = 1.8;   #At what point a protein is considered fused to another
my $large_protein   = 200;
my $small_protein   = 100;

my $high_coverage   = 60.0;  #Assumption of good coverage for candidates 60
my $mid_coverage    = 50.0;  #Assumption of ok coverage  for candidates  50
my $min_coverage    = 40.0;  #Assumption of low coverage for candidate   40
my $minCovDiscard   = 30.0;  #don't look at  blast hits if one of the proteins has this or lower coverage
my $minEvalDiscard  = 1e-3;  #Do not look at anything worse than this evalue

my $worst_evalue    = 1e-4;  #E-value for high coverage >=$high_coverage ($minEvalDiscard/100)
my $badCov_evalue   = 1e-10; #E-value for really bad coverage >=$minCovDiscard
my $shortCov_evalue = 1e-6;  #E-value for short coverage alignments >=$min_coverage
my $midCov_evalue   = 1e-6;  #E-value for mid coveage proteins >=$mid_coverage

my $min_identity    = 18.0;  #minimum accepted identity to consider a good hit

my $compStats       = 'F';   #Use composition statistis to filter blast output
my $min_alignment   = 55;    #minimal alignment length
my $minCompMatches  = 0.25;  #At least this fraction of the components in the system must have a blast/ssearch match




read_command_line();

#print Data::Dumper->Dump([$gnmAcc, $gnmDir, $repForm, $blastdb_name, $file_proteome,
#			  $gnmFeatTable, $gblastFile, $blastOverWrite, $outFmt, $outdir,
#			  $geneDis, $fusion_lenRatio, $large_protein, $small_protein,
#			  $high_coverage, $mid_coverage, $min_coverage, $min_identity,
#			  $max_evalue, $min_alignment, $tcdbFaa],
#			 [qw(*gnmAcc *gnmDir *repForm *blastdb_name *file_proteome
#			  *gnmFeatTable *gblastFile *blastOverWrite *outFmt *outdir
#			  *geneDis *fusion_lenRatio *large_protein *small_protein
#			  *high_coverage *mid_coverage *min_coverage *min_identity
#			  *max_evalue *min_alignment *tcdbFaa)]);
#exit;





#==========================================================================
#Create working directories and declare

my $blastdb_dir = "$outdir/blastdb";
system "mkdir $blastdb_dir" unless (-d $blastdb_dir);

my $plotDir = "$outdir/plots";
system "mkdir $plotDir" unless (-d $plotDir);

my $seqDir = "$outdir/sequences";
system "mkdir $seqDir" unless (-d $seqDir);

my $htmlDir = "$outdir/html";
system "mkdir $htmlDir" unless (-d $htmlDir);

my $blastoutDir = "$outdir/files";
system "mkdir $blastoutDir" unless (-d $blastoutDir);


my $outfile = "$outdir/plot_index.html";



#==========================================================================
#Get the TCIDs of multi-component systems in TCDB. This information will
#be used to extract the TCIDs of all multicomponent systems in TCDB.

print "Retrieving TCDB sequences.\n";


#Download the TCDB database in fasta format
$tcdbFaa = "$outdir/tcdb/tcdb.faa" unless ($tcdbFaa);
system "extractTCDB.pl -i tcdb -o $outdir/tcdb -f fasta" unless (-f $tcdbFaa);
die "TCDB seqs not found: $tcdbFaa" unless (-f $tcdbFaa);

my $multSystems = getModeSystems($tcdbFaa, 'multi');

#my $kk = scalar keys %$multSystems;
#print Data::Dumper->Dump([$multSystems->{'1.A.1.2.4'}, $kk ], [qw(*multSystems *nSystems )]);
#exit;



#==========================================================================
#Get the multicomponent systems identified by GBLAST

print "Extracting multicomponent systems in TCDB\n";


#First get all TCIDs of multicomponent systems found by GBLAST.
#
#NOTE:
#  Don't filter by threshold to capture everything that GBLAST found.
#
my %gblastTCIDs = ();
getGBLASTmultiComponentTCIDs($gblastFile, 10, $multSystems, \%gblastTCIDs);

#my $kk = scalar keys %gblastTCIDs;
#print Data::Dumper->Dump([\%gblastTCIDs, $kk ], [qw(*gblastTCIDs *nSystems)]);
#exit;






#==========================================================================
#Download the sequences of the multicomponent systems and put them all
#in a single file that will be used as a query for blastp/ssearch36 searches.

print "Downloading sequences for multicomponent systems.\n";


my $file_blastp_query = "$seqDir/systems.faa";
unless (-f $file_blastp_query && !(-z $file_blastp_query)) {
  downloadSeqsForMCS(\%gblastTCIDs, $multSystems, $file_blastp_query);
}

#print Data::Dumper->Dump([$file_blastp_query], [qw(*file_blastp_query )]);
#exit;


#==========================================================================
#Download the TCDB database in fasta format in order to extract
#the sequence of full proteins for generating hydropathy plots

print "Downloading and generating TCDB blast DB\n";



my $tcdbBlastDBdir = "$ENV{HOME}/db/blastdb";
system "mkdir -p $tcdbBlastDBdir" unless (-d $tcdbBlastDBdir);


#Delete old files if they exist
my $blastTestFile = "$tcdbBlastDBdir/tcdb.pin";
if ($blastOverWrite) {
  system "rm $tcdbBlastDBdir/*" if (-f $blastTestFile);
}

#download tcdb sequences and make blast database
if ($blastOverWrite || !(-f $blastTestFile)) {
  system "extractTCDB.pl -i tcdb -o $tcdbBlastDBdir -f blast -d $tcdbFaa";
}



#==========================================================================
#Create the blast database with the reference genome. This will
#allow searching for missing components and extracting sequences
#of individual proteins in the genome

print "Creating blast database for genome!\n";


my $internalProteome = "$seqDir/proteome.faa";

if ($file_proteome =~ /.+\.gz$/) {
  system qq(gunzip -c $file_proteome > $internalProteome);
}
elsif ($file_proteome =~ /.+\.bz2$/) {
  system qq(bunzip2 -c $file_proteome > $internalProteome);
}
else {
  system qq(cp $file_proteome > $internalProteome);
}


unless (-f $internalProteome && !(-z $internalProteome)) {
  die "Could not get extract sequences for proteome:  $file_proteome -> $internalProteome";
}



#Prepare command
my $blastdb = "$blastdb_dir/$blastdb_name";
my $bdb_cmd =
  qq(makeblastdb -dbtype prot -input_type fasta -title 'Sequences: $blastdb_name' ) .
  qq(-parse_seqids -hash_index -out $blastdb);


#generate the database
print qq($bdb_cmd -in $internalProteome \n);
system qq($bdb_cmd -in $internalProteome) unless (-f "${blastdb}.phr");

#Verify that the database was successfully generated
die "Blast databse $blastdb_name was not created" unless (-f "${blastdb}.phr");

#print "Check blastdb: $blastdb\n";
#exit;





#==========================================================================
#blastp the TCDB proteins (multicomponent systems) against the genome


my %filteredBlast =();
my $blast_file = "";

#Run blastp
if ($prog eq "blastp") {

  $blast_file = "$blastoutDir/blastOutput.csv";


  #----------------------------------------------------------------------
  #Prepare the blast command

  my $outfmt = qq(10 qacc sacc qlen slen evalue pident length qcovs qstart qend sstart send qseq sseq);
  my $blast_cmd =
    qq(blastp -db $blastdb -query $file_blastp_query -evalue $minEvalDiscard -use_sw_tback ) .
    qq( -comp_based_stats $compStats -max_hsps 1 -outfmt '$outfmt' -out $blast_file );


  #----------------------------------------------------------------------
  #Run blastp

  print "Aligning multi-component systems against the genome with $prog!\n";


  system $blast_cmd unless (-f $blast_file && !(-z  $blast_file));
  die "Error: Blastp output file not found -> $blast_file " unless (-f $blast_file && !(-z  $blast_file));



  #----------------------------------------------------------------------
  #Parse the blastp output file, select top blastp hits, and get the protein
  #sequences of the aligned regions

  print "Parsing blastp output\n";

  filterBlastOutput($blast_file, \%filteredBlast, \%gblastTCIDs, $multSystems);

}

#Run ssearch36
else {

  $blast_file = "$blastoutDir/ssearch.out";

  #----------------------------------------------------------------------
  #Prepare the ssearch command

  print "Aligning multi-component systems against the genome with $prog!\n";

  my $params = qq(-z 11 -k 1000 -m 0 -W 0 -E $minEvalDiscard -s $subMatrix );
  #my $ssearch_cmd = qq(ssearch36 $params $file_blastp_query '$blastdb 12' > $blast_file);
  my $ssearch_cmd = qq(ssearch36 $params $file_blastp_query $internalProteome > $blast_file);

  print "   $ssearch_cmd\n";
  system $ssearch_cmd unless (-f $blast_file && !(-z  $blast_file));
  die "Error: ssearch36 output file not found -> $blast_file " unless (-f $blast_file && !(-z  $blast_file));

#  print "Check sssearch results: $blast_file\n";
#  exit;


  #----------------------------------------------------------------------
  #Parse ssearch output

  print "Parsing ssearch36 output\n";

  filterSSEARCHoutput($blast_file, \%filteredBlast, \%gblastTCIDs, $multSystems);

#  print Data::Dumper->Dump([\%filteredBlast], [qw(*filterdBlast)]);
#  exit;


}

#print Data::Dumper->Dump([$filteredBlast{'2.A.78.2.1'}], [qw(*katie)]);
#die "Finished parsing.\n";



#==========================================================================
#For each genome protein and its respective match in TCDB that passed the
#filtering process, infer the number of TMS. This will help the user to
#determine which of the components are membrane proteins.

print "Running HMMTOP.....\n";

my %TMS = ();
runHMMTOP(\%filteredBlast, \%TMS);


#print Data::Dumper->Dump([\%TMS], [qw(*ntms)]);
#exit;




#==========================================================================
#Now check the genomic context of all multicomponent systems and
#idenfiy those TCDB components that  share genomic neighborhood
#in the genome under analysis.

print "Analyzing genome context\n";
my %gnmAnnotations = ();
verify_neighborhood(\%filteredBlast, $gnmFeatTable, \%gnmAnnotations);

#print Data::Dumper->Dump([\%filteredBlast ], [qw(*filtered_blast)]);
#exit;






#==========================================================================
#Create text file with top hits and genome context results


print "Saving results to text and html reports\n";

generate_reports(\%filteredBlast, \%gnmAnnotations);

print "Done!\n";






###########################################################################
##                                                                       ##
##                        Subroutine definitions                         ##
##                                                                       ##
###########################################################################



#==========================================================================
#Generate the text and html output files with the results of the analysis


sub generate_reports {

  my ($systems, $annotations) = @_;

  my $outfile  = "$outdir/reportMulticomponentSystems.$outFmt";
  my $htmlFile = "$outdir/reportMulticomponentSystems.html";

  #define separator for text file
  my $sep = "\t";
  if ($outFmt eq 'csv') {
    $sep = ",";
  }


  #base URLs for the HTML table
  my $sysURL  = "http://tcdb.org/search/result.php?tc=";
  my $accURL  = "http://tcdb.org/search/result.php?acc=";
  my $ncbiURL = "https://www.ncbi.nlm.nih.gov/protein/";



  #Get the HTML header
  my @hHeader = qw (#tcid  query_accession  status  subject_accession  hydropathy
		    query_length  subject_length  query_tms subject_tms evalue  perc_idenity
		    alignment_length  query_coverage  subject_coverage
		    neighbors gnm_annotation);
  my $htmlHeader = <<HEADER;
<html>
  <head><title>Multicomponent Systems Analysis</title></head>
  <h1 style='text-align:center'>Milticomponent Systems, Genomic Context and Hydropathy Plots</h1>
  <body>
  <table border='1', style='width:100%'>
    <tr>
HEADER

  my $hColor = '#CC9966';
  $htmlHeader .= "      <th style='background-color:$hColor;'>" .
    join ("</th>\n      <th style='background-color:$hColor;'>", @hHeader) .
    "</th>\n    </tr>\n";
#  my $noHitColSpan = scalar(@hHeader) - 2;


  #Get the text header
  my @tHeader = qw (tcid  query_accession  status  subject_accession query_length  subject_length
                    query_tms subject_tms evalue  perc_idenity alignment_length  query_coverage
                    subject_coverage neighbors gnm_annotation);

  my $txtHeader = join ($sep, @tHeader) . "\n";



  open  (my $htmh, ">", $htmlFile) || die $!;
  print $htmh $htmlHeader;




  open (my $outh, ">", $outfile) || die $!;
  print $outh $txtHeader;


  my $rowCnt = 1;

 TC:foreach my $tcid (sort_by_system keys %{ $systems }) {  #Based on the entire TCDB


    #Ignore TCID if there are no detected hits.
    next TC unless (system_has_hits($tcid));


    my $rowColor = "";
    #Get the color for this system
    if ($rowCnt % 2) {
      $rowColor = "#FFFFFF";
    }
    else {
      $rowColor = "#FFEBCD";
    }

    foreach my $tcAcc (sort_by_status1 ($systems->{$tcid})) {
      foreach my $hit (sort by_status2 @{ $systems->{$tcid}->{$tcAcc} }) {

	my $status    = $hit->{status};
	my $sacc      = "";
	my $eval      = "";
	my $ident     = "";
	my $qlen      = "";
	my $slen      = "";
	my $qtms      = " ";
	my $stms      = " ";
	my $aln       = "";
	my $qcov      = "";
	my $scov      = "";
	my $neighbors = "";
	my $annot     = " ";


	my $qtcid = $tcid . '-' . $tcAcc;
	unless (exists $TMS{$qtcid}) {
	  print "No number of TMSs for query: $qtcid\n";
	  print Data::Dumper->Dump([$hit], [qw(*hit )]);
	  exit;
	}
	$qtms = $TMS{$qtcid}->[0];
	$qlen = $TMS{$qtcid}->[1];

#	print Data::Dumper->Dump([$hit, $qtcid, $qtms, $qlen], [qw(*hit *qtcid *qtms *qlen)]);
#	<STDIN>;


	if ($status =~ /(GBLAST|Candidate|RawBlast)/) {
	  $sacc      = $hit->{sacc};
	  $qlen      = $hit->{qlen};
	  $slen      = $hit->{slen};
	  $eval      = $hit->{eval};
	  $ident     = $hit->{ident};
	  $aln       = $hit->{aln};
	  $qcov      = $hit->{qcov};
	  $scov      = sprintf("%.1f", $hit->{scov});
	  $neighbors = $hit->{neighbors};
	  $annot     = (exists $annotations->{$sacc})? $annotations->{$sacc} : die "No name annotation for gene: $sacc";

	  unless (exists $TMS{$sacc}) {
	    print "No number of TMSs for subject: $sacc\n";
	    print Data::Dumper->Dump([$hit], [qw(*hit )]);
	    exit;
	  }

	  $stms = $TMS{$sacc}->[0];
	}


	#print the text version of the row
	print $outh "$tcid${sep}$tcAcc${sep}$status${sep}$sacc${sep}$qlen${sep}$slen${sep}$qtms${sep}$stms${sep}$eval${sep}$ident${sep}$aln${sep}$qcov${sep}$scov${sep}$neighbors${sep}$annot\n";


	#Instead of running quod as before, now run xgbhit.sh for each match.
	if ($status =~ /(GBLAST|Candidate|RawBlast)/) {

	  #generate hydropathy plots and print the corresponding HTML row.
	  my $plotFile = "${tcAcc}_vs_$sacc/ssearch_${sacc}_vs_${tcAcc}/report.html";

	  print "  Creating plots comparing:  ${tcid}-${tcAcc} vs ${sacc}\n";

	  my $good = run_quod($tcid, $tcAcc, $sacc);
	  unless ($good) {
	    print Data::Dumper->Dump([$tcid, $tcAcc, $sacc, $hit->{qstart}, $hit->{qend}, $hit->{sstart}, $hit->{send}, $hit->{qseq},$hit->{sseq}, $annotations->{$sacc}],
				     [qw(*tcid *tcAcc *sacc *qstart *qend *sstart *send *qseq *sseq *annotation)]);
	    die "Could not generate plot for: ${tcid}-$tcAcc vs $sacc -> ";
	  }

	  my @data = ($tcid, $tcAcc, $status, $sacc, $qlen, $slen, $qtms, $stms, $eval, $ident, $aln, $qcov, $scov, $neighbors, $annot);
	  my $hitRaw = getMatchRowHTMLstring($sysURL, $accURL, $ncbiURL, $plotFile, $rowColor, \@data);

	  print $htmh $hitRaw;
	}

	#print html row with no hit
	else {
	  my $noHitRaw = getNoHitRowHTMLstring($sysURL, $accURL, $tcid, $tcAcc, $status, $qlen, $qtms, $rowColor);
	  print $htmh $noHitRaw;
	}

#	print Data::Dumper->Dump([$tcid, $tcAcc, $hit ], [qw(*tcid *tcAcc *hit )]);
      }
    }
    $rowCnt++;
  }
  close $outh;


  my $footer = <<FOOTER;
  </table>
  </body>
</html>
FOOTER

  print $htmh $footer;

}



#==========================================================================
#Return true if system has at least one component with a hit in the
#genome

sub system_has_hits {

  my $tcid = shift;

  my $okHits = 0;

  my $nComp = scalar keys %{ $filteredBlast{$tcid} };

  #at least $minCompMatches fraction of the components in the systems
  #must have an acceptable hit
 CMP:foreach my $component (keys %{ $filteredBlast{$tcid} }) {

    my @hits = @{ $filteredBlast{$tcid}->{$component} };

  HIT:foreach my $hit (@hits) {

      if ($hit->{status} =~ /NoHit/) {
	next HIT;
      }
      else {
	$okHits++;
      }
    }
  }

  if (($okHits / $nComp) < $minCompMatches) {
    return 0;
  }
  else {
    return $okHits;
  }
}



#==========================================================================
#get a string with the HTML code for a row of data

#Print an HTML raw with the info of blastp match
sub getMatchRowHTMLstring {

  my ($sysURL, $accURL, $ncbiURL, $plotFile, $color, $data) = @_;


  my ($tcid, $tcAcc, $status, $sacc, $qlen, $slen, $qtms, $stms, $eval, $ident, $aln, $qcov, $scov, $neighbors, $annot) = @$data;

  my $identity = sprintf ("%.1f", $ident);

  $neighbors =~ s/\|/\<br \\>/g;

  my $row = <<"NOHIT";
    <tr>
      <td style='text-align:left;   background-color:$color;'><a href='$sysURL${tcid}' target='_blank'>$tcid</a></td>
      <td style='text-align:center; background-color:$color;'><a href='$accURL${tcAcc}' target='_blank'>$tcAcc</a></td>
      <td style='text-align:center; background-color:$color;'>$status</td>
      <td style='text-align:center; background-color:$color;'><a href='${ncbiURL}$sacc' target='_blank'>$sacc</a></td>
      <td style='text-align:center; background-color:$color;'><a href='plots/$plotFile' target='_blank'>plots</a></td>
      <td style='text-align:right;  background-color:$color;'>$qlen</td>
      <td style='text-align:right;  background-color:$color;'>$slen</td>
      <td style='text-align:right;  background-color:$color;'>$qtms</td>
      <td style='text-align:right;  background-color:$color;'>$stms</td>
      <td style='text-align:center; background-color:$color;'>$eval</td>
      <td style='text-align:right;  background-color:$color;'>$identity</td>
      <td style='text-align:right;  background-color:$color;'>$aln</td>
      <td style='text-align:right;  background-color:$color;'>$qcov</td>
      <td style='text-align:right;  background-color:$color;'>$scov</td>
      <td style='text-align:center; background-color:$color;'>$neighbors</td>
      <td style='text-align:center; background-color:$color;'>$annot</td>
    </tr>
NOHIT

  return $row;

}




#==========================================================================
#get the HTML string to print a row that had no hits

sub getNoHitRowHTMLstring {

  my ($tcURL, $accURL, $tcid, $tcAcc, $status, $qlen, $qtms, $color) = @_;

  my $row = <<"NOHIT";
    <tr>
      <td style='text-align:left; background-color:$color;'><a href='$tcURL${tcid}' target='_blank'>$tcid</a></td>
      <td style='text-align:center; background-color:$color;'><a href='$accURL${tcAcc}' target='_blank'>$tcAcc</a></td>
      <td style='text-align:center; background-color:$color; padding-left:10px;'>$status</td>
      <td style='text-align:center; background-color:$color;'></td>
      <td style='text-align:center; background-color:$color;'></td>
      <td style='text-align:right;  background-color:$color;'>$qlen</td>
      <td style='text-align:right;  background-color:$color;'></td>
      <td style='text-align:right;  background-color:$color;'>$qtms</td>
      <td style='text-align:right;  background-color:$color;'></td>
      <td style='text-align:center; background-color:$color;'></td>
      <td style='text-align:right;  background-color:$color;'></td>
      <td style='text-align:right;  background-color:$color;'></td>
      <td style='text-align:right;  background-color:$color;'></td>
      <td style='text-align:right;  background-color:$color;'></td>
      <td style='text-align:center; background-color:$color;'></td>
      <td style='text-align:center; background-color:$color;'></td>
    </tr>
NOHIT

  return $row;
}





#==========================================================================
#Sort TC accessions by  the number of hits in GBLAST, Candidate,
#RawBlast and NoHit

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


#Subordinate funtion to sort_by_status1
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


#==========================================================================
#Sort 2 blastp hits by status and quality

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

    my $maxCov1 = ($qcov1 + $scov1)/2;
    my $maxCov2 = ($qcov2 + $scov2)/2;


    if ($eval1 == $eval2) {
      $maxCov2 <=> $maxCov1;  #sort by max Coverage
    }
    else {
      $eval1 <=> $eval2; #sort by $evalue
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

  my ($systems, $featuresFile, $annot) = @_;

  #
  # $systems contains the filtered blast ouput
  #


  #-----------------------------------------------------------------
  #Parse genome feature table

  open (my $fh1, "-|",  "gunzip -c $featuresFile") || die $!;

  my @cds = ();
  my %pos2cds = ();
  my %cds2pos = ();
  my %protPerScaffold = ();

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
    #6)  genomic_accession (or scaffold in metagenomes)
    #7)  start
    #8)  end
    #9) strand
    #10) product_accession
    #11) non-redundant_refseq
    #12) related_accession
    #13) name (functional annotation)
    #14) symbol
    #15) GeneID
    #16) locus_tag
    #17) feature_interval_length
    #18) product_length
    #19) attributes

    #The protein accession
    my @d = split (/\t/, $_);
    next unless($d[10]);  #Ignore line if there is no protein accession
    my ($acc, $ver) = split(/\./, $d[10]);


    #The replicon or scaffold accession
    my ($scaffacc, $v2) = split(/\./, $d[6]);


    my $line = { scaffacc=>$scaffacc, acc=>$acc, start=>$d[7], end=>$d[8], strand=>$d[9], name=>$d[13], gene=>$d[14] };
    push (@{ $protPerScaffold{$scaffacc} }, $line);
    $annot->{$acc} = $line->{name};
  }
  close $fh1;


#  print Data::Dumper->Dump([\%protPerScaffold], [qw( *protPerScaffold )]);
#  exit;


  #-----------------------------------------------------------------
  #Get the ranked position of each CDS per scaffold


  foreach my $scaffold (sort {$a cmp $b} keys %protPerScaffold) {

    #gene position within the scaffold
    my $pos = 1;

    #NOTE: In theory the CDS is already sorted by scaffold and position,
    #      but it's better to sort it to prevent problems with unexpected
    #      misformatting of the reference features table.

    foreach my $orf (sort { $a->{start} <=> $b->{start} } @{ $protPerScaffold{$scaffold} }) {

#      print "$scaffold\t",$orf->{acc}, "\n";

      $pos2cds{ $scaffold }{ $pos } = $orf;
      $cds2pos{ $orf->{acc} } = [$pos, $scaffold];
      $pos++;
    }
  }


#  print Data::Dumper->Dump([\%cds2pos], [qw(*cds2pos)]);
#  exit;


  #-----------------------------------------------------------------
  #For each multicomponent system in the genome see which CDS
  #hitting different components are in the neighborhood

#  print Data::Dumper->Dump([$systems ], [qw(*filteredSystems )]);
#  exit;


 TCID:foreach my $tcid (keys %{ $systems }) {

    my %neighbors = ();

    #For debugging purpose
#    next TCID unless ($tcid eq "2.A.63.1.4");

    #For this system, get the accessions of the CDS matching each component
    #and their position rank in the replicon
    my %candComp = ();
    getAccessionsForSystem($tcid, $systems, \%cds2pos, \%candComp);

#    print Data::Dumper->Dump([$tcid, \%candComp ], [qw(*tcid *matches)]);
#    exit;



    #----------------------------------------------------------------------
    #Using each candidate component as references get how many of the
    #other proteins matching components in the same systems are in
    #the neigborhood.


    my @accs = ();
    sort_hits_by_rep_and_position (\%candComp, \@accs);

#    print Data::Dumper->Dump([$tcid, \@accs], [qw(*tcid *sorted)]);
#    exit

   ACC1:foreach my $acc1 (@accs) {

       my @nearby = ($acc1);


#       print "\n\n===========================================================================\n";
#       print "$tcid\n\n";


     ACC2:foreach my $acc2 (@accs) {

	 next ACC2 if ($acc1 eq $acc2);


#	 print "$acc1 vs $acc2\n";

	 #determine if the genes are neighbors
	 my @dist = areGenesCloseby($candComp{$acc1}, $candComp{$acc2}, $geneDis, $repForm, \%protPerScaffold);


#	 print Data::Dumper->Dump([\@dist, $candComp{$acc2}], [qw(*dist *candAcc2)]);
#	 <STDIN>;


	if (@dist) {

#	  print Data::Dumper->Dump([\@dist], [qw(*dist )]), "\nCheck Distances\n";;
#	  <STDIN>;
#	  print Data::Dumper->Dump([$candComp{$acc1}, $candComp{$acc2}], [qw(*cand1 *cand2 )]), "\nCheck candidates\n";
#	  <STDIN>;

	  my @cmpHit = ();
	  foreach my $acc2hit (@{ $dist[1] }) {
	    push (@cmpHit, $acc2hit->[1]);
	  }
	  my $cmpStr = join(",", @cmpHit);
	  push (@nearby, "${acc2}(". $dist[2] ."):$cmpStr");

	}

#	print Data::Dumper->Dump([$acc1, $acc2, \@dist, \@nearby ], [qw(*acc1 *acc2 *dist *nearby)]), "\n";
#	<STDIN>;

       } #ACC2


#       print "\n------------------------------\n";
#       print Data::Dumper->Dump([$acc1, \@nearby], [qw(*acc *neighbors )]);
#       <STDIN>;


       #register whether there were neighbors or not
       if (scalar @nearby > 1) {
	 $neighbors{$acc1} = join("|", sort {$a cmp $b} @nearby);
       }
       else {
	 $neighbors{$acc1} = "None";
       }

#       print Data::Dumper->Dump([\%candComp, \%neighbors ], [qw(*tcmatches *neighbors )]);
#       <STDIN>;
     } #ACC1

#    "Katie\n";
#    exit;

    #-----------------------------------------------------------------
    #With the neighbors fully estimated for this system, now add the
    #neighbors information to the filtered multicomponent systems hits.

    foreach my $tcAcc (keys %{ $systems->{$tcid} }) {
      foreach my $hit (@{ $systems->{$tcid}->{$tcAcc} }) {
	if (exists $hit->{sacc} && exists $neighbors{ $hit->{sacc} }) {
	  $hit->{neighbors} = $neighbors{ $hit->{sacc} };
	}
      }
    }

#    print "$tcid\n", Data::Dumper->Dump([ $systems->{$tcid} ], [qw( *Filtered )]);
#    <STDIN>;

  } #TCID
}


#==========================================================================
#Sort matches accoring to their position within their replicon/scaffold
#
# Strcuture of $matches:
#%matches = ('Dvar_61810' => [[6181,'REP001'],'P51149'],
#            'Dvar_42640' => [[4264,'REP001'],'Q9UBQ0'],
#            'Dvar_72140' => [[7213,'REP001'],'Q9ULJ7'],
#	    ...
#           );



sub sort_hits_by_rep_and_position {

  my ($matches, $out) = @_;

  @{$out} = sort
    {
      #Check whether proteins belong to the same replicon.
      if ($matches->{$a}->[0]->[1] eq $matches->{$b}->[0]->[1]) {
	$matches->{$a}->[0]->[0] <=> $matches->{$b}->[0]->[0];
      }
      else {
	$matches->{$a}->[0]->[1] cmp $matches->{$b}->[0]->[1];
      }
    }
    keys %{ $matches };
}



#==========================================================================
#determine if two genes are neighbors and take into account the
#Replicon or scaffold and linearity or circularity of the genome under analysis.


sub areGenesCloseby {

  my ($pos1, $pos2, $refGeneDist, $repStucture, $repHash) = @_;


  #First determine if genes are in the same replicon
  my $rep1 = $pos1->[0]->[0]->[1];
  my $rep2 = $pos2->[0]->[0]->[1];


  #Do not calculate distance if replicon or scaffold in different
  return () unless ($rep1 eq $rep2);


  #Number of proteins in the replicon or Scaffold
  my $repSize = scalar @{ $repHash->{$rep1} };


  #The $pos1 and $pos2 arrays are too complecated now and some data is repeated
  #unecessarily *** I can simplify them in the future ***
  my @pos = sort {$a <=> $b} ($pos1->[0]->[0]->[0], $pos2->[0]->[0]->[0]);
  my $dist = undef;


  #Replicon is circular
  if ($repStucture eq 'circular') {

    my $d1 = $pos[1] - $pos[0];
    my $d2 = $pos[0] + ($repSize - $pos[1]);

    $dist = (sort {$a <=> $b} ($d1, $d2))[0];
  }


  #Replicon is linear (this should happen for scaffolds)
  else {
    $dist = $pos[1] - $pos[0];
  }


#  print Data::Dumper->Dump([$pos1, $pos2, $dist,  $repSize], [qw(*pos1 *pos2 *dist *repSize)]);
#  <STDIN>;


  #Test if distance is in the acceptable range
  if ($dist <= $refGeneDist) {
    return ($pos1, $pos2, $dist);
  }
  else {
    return ();
  }

}




#==========================================================================
# Given a TCID of a multicomponent system, get the genome proteins that
# match all the components and their respective position rank in the genome.


sub getAccessionsForSystem {

  my ($tcid, $multSystems, $acc2pos, $out) = @_;

 COMP:foreach my $tcAcc (keys %{ $multSystems->{$tcid} }) {
  HIT:foreach my $hit (@{ $multSystems->{$tcid}->{$tcAcc} }) {

      next HIT if ($hit->{status} eq 'NoHit');

      my $gnmAcc = $hit->{sacc};

      #Need to save hits in an array because some proteins may hit
      #more than one component in the same system.
      if (exists $acc2pos->{ $gnmAcc }) {
	push (@{ $out->{ $gnmAcc }}, [$acc2pos->{ $gnmAcc }, $tcAcc]);
      }
      else {
	die "No position found for: $gnmAcc -> ";
      }
    }
  }
}




#==========================================================================
#Run HMMTOP on the proteins that had blast matches and extract the
#number of TMSs

sub runHMMTOP {

  my ($bdata, $tms) = @_;

  my %tcacc = ();
  my %gacc  = ();


  #----------------------------------------------------------------------
  #Collect the accessions for which hmmtop will run

  foreach my $tc (keys %{ $bdata }) {

    my @comp = keys %{ $bdata->{$tc} };
    foreach my $cmp (@comp) {

      #TC-ID
      my $qacc = "${tc}-$cmp";
      $tcacc{$qacc} = 1;

      foreach my $hit (@{ $bdata->{$tc}->{$cmp} }) {

	my $sacc = 'NA';
	if (exists $hit->{sacc}) {
	  $sacc = $hit->{sacc};
	  $gacc{$sacc} = 1;
	}
      }
    }
  }


  #----------------------------------------------------------------------
  #Get sequences for each set of proteins, combine them in one file,
  #run HMMTOP and parse the output putting the results within the same
  #input hash.

  #To avoid accumulating sequences remove seq file if it exists
  my $allSeqFile = "$seqDir/allAccSeqs.faa";
  system qq(rm $allSeqFile) if (-f $allSeqFile);

  #Extract the sequences from TCDB
  my $tcAccFile = "$seqDir/tcacc.txt";
  getSeqs4hmmtop(\%tcacc, $tcAccFile, "tcdb",  $allSeqFile);

  #Extract the sequences from the genome
  my $gnmAccFile = "$seqDir/gnmacc.txt";
  getSeqs4hmmtop(\%gacc, $gnmAccFile, $blastdb,  $allSeqFile);




  #----------------------------------------------------------------------
  #Run HMMTOP


  my $hmmtopOut = "$blastoutDir/hmmtop.out";
  my $cmd3 = qq(hmmtop -if=$allSeqFile -of=$hmmtopOut);

  system $cmd3 unless (-f $hmmtopOut && !(-z $hmmtopOut));
  die "Problem predicting TMSs: $hmmtopOut" unless (-f $hmmtopOut && !(-z $hmmtopOut));


  #----------------------------------------------------------------------
  #parse HMMTOP output and update hash %filteredBlast, just add the
  #keys: qtms (for TC protein) and stms (for genome protein)


  parseHMMTOP($hmmtopOut, $tms);
}


#==========================================================================
#Given the output of HMMTOP only get the accession and the number of
#TMSs


sub parseHMMTOP {

  my ($hmmtop, $tms) = @_;


  open (my $fh, '<', $hmmtop) || die $!;
  while (<$fh>) {
    chomp;
    next unless (/^\>/);

    my ($leftPart, $center, $rightPart) = ();
    if (/^(.+)(IN|OUT)(.+)$/) {
      ($leftPart, $center, $rightPart) = (/^(.+)(IN|OUT)(.+)$/)? ($1, $2, $3) : undef;
    }

    #Get the accession
    my ($len, $acc) = ($leftPart =~ /^\>HP:\s+(\d+)\s+(\S+)/)? ($1, $2) : undef;
    die "Could not extract the accession: |$_|\n|$leftPart|\n|$rightPart|\n" unless ($acc);

    #Remove version number from accession
    $acc =~ s/\.\d$//;


    #Get the number of TMSs
    my $nTMS = ($rightPart =~ /^\s+(\d+)/)? $1 : undef;
    die "Could not extract TMSs: |$_|\n|$leftPart|\n|$rightPart|" unless (defined $nTMS && $nTMS >= 0);


    $tms->{$acc} = [$nTMS, $len];
  }
  close $fh;

#  print Data::Dumper->Dump([$tms], [qw(*tms )]);
#  exit;
}



#==========================================================================
#Given a set of accessions get their protein sequences in preparation
#to run HMMTOP.


sub getSeqs4hmmtop {

  my ($accHash, $accFile, $db, $seqFile) = @_;

  #Create accessions file
  unless (-f $accFile && !(-z $accFile)) {
    open (my $fh1, ">", $accFile) || die $!;
    print $fh1 join("\n", keys %{ $accHash }), "\n";
    close $fh1;
  }
  die "Could not extract accessions: $accFile" unless (-f $accFile && !(-z $accFile));


  #Extract sequences
  my $cmd = qq($blastBinDir/blastdbcmd -db $db -entry_batch $accFile -target_only >> $seqFile);
  system $cmd;
  die "Failed to extract sequences: $seqFile" unless (-f $seqFile && !(-z $seqFile));


  #Remove accessions file
  my $cmd2 = qq(rm $accFile);
  system $cmd2 if (-f $accFile);
}


#==========================================================================
#Parse the ssearch36 output file, and select top hits


sub filterSSEARCHoutput {

  my ($infile, $filtered, $gblastMCS, $tcdbMCS) = @_;

  my %rawBlast = ();


  my $parser    = new Bio::SearchIO (-format => 'fasta', -file => $infile);
  my $formatTmp = $parser->format();

 RT:while (my $result = $parser->next_result) {

    my $qseqid = $result->query_name;
    my $qlen   = $result->query_length;

#    #For controlled debuging
#    next RT unless ($qseqid =~ /1\.B\.22\.1\.3|2\.A\.6\.1\.8/);

  HIT:while (my $hit = $result->next_hit) {
    HSP:while(my $hsp = $hit->next_hsp) {


	#Alignment parameters
	my $sseqid = $hit->name;
	my $slen   = $hit->length;
	my $evalue = $hsp->evalue;
	my $pident = $hsp->frac_identical('total') * 100;


	#There must be a minimum sequence identity
	next HSP if ($pident < $min_identity);


	#coordinates and sequence
	my $qstart  = $hsp->start('query');
	my $qend    = $hsp->end('query');
	my $sstart  = $hsp->start('subject');
	my $send    = $hsp->end('subject');
	my $qseq    = $hsp->query_string;
	my $sseq    = $hsp->hit_string;
	my $hstr    = $hsp->homology_string;


	#If the alignment has less than $minLength aas, ignore it
        my $qtmp = $qseq; $qtmp =~ s/-//g;
	my $stmp = $sseq; $stmp =~ s/-//g;
	my $length = (sort {$a<=>$b} (length($qtmp), length($stmp)))[0];


	#go to next hit if the alignment is not long enough
	#next HSP if ($length < $min_alignment);


	#Exit unless all data have been parsed properly
	unless ($qseqid && $sseqid && $qlen && $slen && ($evalue >= 0) && $pident && $length && $qseq && $sseq) {
	  print "Error: could not extract all info for hit $qseqid vs ${sseqid}:\n";
	  print Data::Dumper->Dump([$hsp ], [qw(*HSP)]);
	  exit;
	}


	#Calculate coverages
	my $qcov = ($qend - $qstart) / $qlen * 100;
	my $scov = ($send - $sstart) / $slen * 100;


	#Get TCIDs and clean accessions
	my ($tcid, $qacc) = split(/-/,  $qseqid);
	my ($sacc, $ver)  = split(/\./, $sseqid);


#	print "Before:\n", Data::Dumper->Dump([$qseqid, $sseqid, $evalue, $pident, $qcov, $scov], [qw(*qacc *sacc *evalue *pident *qcov *scov )]);
#	<STDIN>;

	#*** Ignore hits with not enough coverage and a minimum alignment length ***
	next HSP if (
		     #Discard if any sequence is covered poorly and has bad E-value
		     (($qcov < $minCovDiscard  || $scov < $minCovDiscard) && $evalue > $badCov_evalue) ||

		     #Discard if both coverages are mediums and the E-value is not acceptable
		     (($scov < $high_coverage && $qcov < $high_coverage) && $evalue > $midCov_evalue)  ||

		     #Discard if both coverages are poor reagrdless of the E-value
		     ($scov <= $min_coverage && $qcov <= $min_coverage)  ||

		     #E-value must be acceptable
		     $evalue > $minEvalDiscard);


#	print "After:\n", Data::Dumper->Dump([$qseqid, $sseqid, $evalue, $pident, $qcov, $scov], [qw(*qacc *sacc *evalue *pident *qcov *scov )]);
#	<STDIN>;



	#Format some varaibles for printing.
	my $qcovStr = sprintf ("%.1f", $qcov);
	my $scovStr = sprintf ("%.1f", $scov);


	#The initial status of this hit
	my $status = undef;
	if (exists $gblastTCIDs{$tcid} && exists $gblastTCIDs{$tcid}{$qacc} && exists $gblastTCIDs{$tcid}{$qacc}{$sacc}) {
	  $status = "GBLAST";
	}
	else {
	  $status = "RawBlast";
	}


	my $hitData = {tcid=>$tcid,
		       qacc=>$qacc,
		       sacc=>$sacc,
		       qlen=>$qlen,
		       slen=>$slen,
		       eval=>$evalue,
		       ident=>$pident,
		       aln=>$length,
		       qcov=>$qcovStr,
		       scov=>$scovStr,
		       qstart=>$qstart,   #query coords for quod
		       qend=>$qend,       #query coords for quod
		       sstart=>$sstart,   #subject coords for quod
		       send=>$send,       #subject coords for quod
		       qseq=>"-", #$qseq,
		       sseq=>"-", #$sseq,
		       status=>$status};


#	print Data::Dumper->Dump([$tcid,$qacc,$sseqid,$qlen,$slen,$evalue,$pident,$length,$qcovStr,$scovStr,$qseq,$sseq],
#				 [qw(*tcid *qacc *sseqid *qlen *slen *evalue *pident *length *qcovs *scov *qseq $sseq)]);
#	<STDIN>;


	#Hits with minimum coverage, and minimum e-value
	#If both proteins have less than mid_coverage, require at most short coverage e_value
	unless ($scov < $mid_coverage && $qcov < $mid_coverage && $evalue > $shortCov_evalue) {
	  push(@{ $rawBlast{$tcid}{$qacc} }, $hitData);
	}



	#------------------------------------------------------------
	#Determine if this is a good hit


	#If one protein is much larger than the other, there must be high coverage of the
	#smallest protein
	my $lenRatio = ($qlen >= $slen)? $qlen/$slen : $slen/$qlen;


	#Very bad coverage in at least one protein (requires better e-value: $badCov_evalue)
	if ($scov < $min_coverage || $qcov < $min_coverage) {

	  #One coverage is bad, at least one protein must have acceptable coverage with
	  #acceptable E-value
	  if ($evalue <= $badCov_evalue && ($scov >= $mid_coverage || $qcov >= $mid_coverage)) {
	    $hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
	    push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
	  }
	}


	#small to mid coverage in at least one protein
	elsif ($scov < $mid_coverage || $qcov < $mid_coverage) {
	  if ($evalue <= $shortCov_evalue) {
	    $hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
	    push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
	  }
	}


	#mid to high coverage in at least one protein
	elsif ($scov < $high_coverage || $qcov < $high_coverage) {
	  if ($evalue <= $midCov_evalue) {
	    $hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
	    push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
	  }
	}


	#At this point both coverages are >= $high_coverage
	elsif ($evalue < $worst_evalue) {
	  $hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
	  push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
	}

#	print Data::Dumper->Dump([$filtered ], [qw(*parsed )]);
#	<STDIN>;

      } #HSP while
    } #HIT while
  }  #RT while



  #Now complement filtered hits with raw blast hits to have as many components as
  #possible identified.
# TCID:foreach my $tcid (keys %{ $gblastMCS }) {
 TCID:foreach my $tcid (keys %{ $tcdbMCS }) {

    #Total tcdb components for system
    my $tcdbComp = $tcdbMCS->{$tcid};

    #Search hits for each component in this system
  ACC:foreach my $acc (@{ $tcdbComp }) {

      if (exists $filtered->{$tcid} && exists $filtered->{$tcid}->{$acc}) {
	next ACC;
      }

      elsif (exists $rawBlast{$tcid} && exists $rawBlast{$tcid}{$acc}) {
	$filtered->{$tcid}->{$acc} = $rawBlast{$tcid}{$acc};
      }
      else {
	$filtered->{$tcid}->{$acc} = [{status=>'NoHit'}];
      }
    }

#    print Data::Dumper->Dump([$tcid, $filtered->{$tcid}], [qw(*tcid *complemented)]);
#    <STDIN>;
  }

}



#==========================================================================
#Parse the blastp output file and select top blastp hits

sub filterBlastOutput {

  my ($infile, $filtered, $gblastMCS, $tcdbMCS) = @_;

  my %rawBlast = ();

  open (my $h2, "<", $infile) || die $!;
  while (<$h2>) {

    chomp;
    next unless ($_);

    my ($qseqid, $sseqid, $qlen, $slen, $evalue, $pident, $length, $qcov, $qstart, $qend, $sstart, $send, $qseq, $sseq) = split (/,/, $_);

    unless ($qseqid && $sseqid && $qlen && $slen && $evalue && $pident && $length && $qcov >= 0 && $qseq && $sseq) {
      die "Error: not a csv blastp output line -> $_";
    }


    #
    #NOTE: Do not use $length to calculate coverages because it includes gaps, this may cause significant
    #      differences in the bars/wedges delimiting alignments within hydropathy plots. Therefore use
    #      alignment coordinates to calculate coverage.
    #

    #The coverage of the query protein
    my ($tcid, $qacc) = split(/-/, $qseqid);
    my $qcovs = ($qend - $qstart) / $qlen * 100;
#    my $qcovs = ($qcov_tmp >= 100.0)? 100.0 : $qcov_tmp;


    #The coverage of the subject protein
    my ($sacc, $ver) = split(/\./, $sseqid);
    my $scov = ($send - $sstart) / $slen * 100;
#    my $scov = ($scov_tmp >= 100.0)? 100.0 : $scov_tmp;


    #*** Ignore hits with not enough coverage and a minimum alignment length ***
    next if ($qcovs < $minCovDiscard  || $scov < $minCovDiscard ||   #At least one sequence is covered very poorly (discard!)
	     $qcovs < $min_coverage   || $scov < $min_coverage  ||   #Both sequences must have the minimum coverage
	     $length < $min_alignment || $evalue > $minEvalDiscard); #Check for the minimum alignment size and max E-evalue


    #Format some varaibles for printing.
    my $qcovStr = sprintf ("%.1f", $qcovs);
    my $scovStr = sprintf ("%.1f", $scov);


    #The initial status of this hit
    my $status = undef;
    if (exists $gblastTCIDs{$tcid} && exists $gblastTCIDs{$tcid}{$qacc} && exists $gblastTCIDs{$tcid}{$qacc}{$sacc}) {
      $status = "GBLAST";
    }
    else {
      $status = "RawBlast";
    }


    my $hitData = {tcid=>$tcid,
		   qacc=>$qacc,
		   sacc=>$sacc,
		   qlen=>$qlen,
		   slen=>$slen,
		   eval=>$evalue,
		   ident=>$pident,
		   aln=>$length,
		   qcov=>$qcovStr,
		   scov=>$scovStr,
		   qstart=>$qstart,   #query coords for quod
		   qend=>$qend,       #query coords for quod
		   sstart=>$sstart,   #subject coords for quod
		   send=>$send,       #subject coords for quod
		   qseq=>$qseq,
		   sseq=>$sseq,
		   status=>$status};


#    print Data::Dumper->Dump([$tcid,$qacc,$sseqid,$qlen,$slen,$evalue,$pident,$length,$qcovStr,$scovStr,$qseq,$sseq],
#			     [qw(*tcid *qacc *sseqid *qlen *slen *evalue *pident *length *qcovs *scov *qseq $sseq)]);
#    <STDIN>;


    #Hits with minimum coverage and any evalue
    push(@{ $rawBlast{$tcid}{$qacc} }, $hitData);



    #------------------------------------------------------------
    #Determine if this is a good hit


    #If one protein is much larger than the other, there must be high coverage of the
    #smallest protein
    my $lenRatio = ($qlen >= $slen)? $qlen/$slen : $slen/$qlen;


    #Very bad coverage in at least one protein.
    #This condition will never be met because now I require that both sequences
    #Should have the minum coverage (see above)
    if ($scov < $min_coverage || $qcovs < $min_coverage) {
      if ($evalue <= $badCov_evalue) {
	$hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
	push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
      }
    }


    #small to mid coverage 30-40% in at least one protein
    elsif ($scov < $mid_coverage || $qcovs < $mid_coverage) {
      if ($evalue <= $shortCov_evalue) {
	$hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
	push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
      }
    }


    #mid to high coverage in at least one protein
    elsif ($scov < $high_coverage || $qcovs < $high_coverage) {
      if ($evalue <= $midCov_evalue) {
	$hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
	push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
      }
    }


    #At this point both coverages are >= $high_coverage
    elsif ($evalue <= $worst_evalue) {
	$hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
	push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
    }



    #based on protein length
    #selecting potential fusions
    # if ($lenRatio >= $fusion_lenRatio) {
    #   if ($scov >= $high_coverage || $qcovs >= $high_coverage) {
    # 	if ($evalue <= $worst_evalue || $pident >= $min_identity) {
    # 	  $hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
    # 	  push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
    # 	}
    #   }
    # }

    # #large proteins
    # elsif ($qlen >= $large_protein || $slen >= $large_protein) {
    #   if ($scov >= $min_coverage || $qcovs >= $min_coverage) {
    # 	if ($evalue <= $max_evalue || $pident >= $min_identity) {
    # 	  $hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
    # 	  push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
    # 	}
    #   }
    # }

    # #mid size proteins
    # elsif ($qlen >= $small_protein || $slen >= $small_protein) {
    #   if ($length >= $min_alignment || $scov >= $mid_coverage || $qcovs >= $mid_coverage) {
    # 	if ($evalue <= $max_evalue || $pident >= $min_identity) {
    # 	  $hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
    # 	  push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
    # 	}
    #   }
    # }

    # #small proteins
    # else {
    #   if ($scov >= $high_coverage || $qcovs >= $high_coverage) {
    # 	if ($evalue <= $max_evalue || $pident > $min_identity + 5.0) {
    # 	  $hitData->{status} = 'Candidate' unless ($hitData->{status} eq 'GBLAST');
    # 	  push (@{ $filtered->{$tcid}->{$qacc} }, $hitData);
    # 	}
    #   }
    # }

  }
  close $h2;



  #Now complement filtered hits with raw blast hits to have as many components as
  #possible identified.
 TCID:foreach my $tcid (keys %{ $gblastMCS }) {

    #Total tcdb components for system
    my $tcdbComp = $tcdbMCS->{$tcid};

    #Search hits for each component in this system
  ACC:foreach my $acc (@{ $tcdbComp }) {

      if (exists $filtered->{$tcid} && exists $filtered->{$tcid}->{$acc}) {
	#$filtered->{$tcid}->{$acc} = $filtered->{$tcid}->{$acc};
	next ACC;
      }

      elsif (exists $rawBlast{$tcid} && exists $rawBlast{$tcid}{$acc}) {
	$filtered->{$tcid}->{$acc} = $rawBlast{$tcid}{$acc};
      }
      else {
	$filtered->{$tcid}->{$acc} = [{status=>'NoHit'}];
      }
    }

#    print Data::Dumper->Dump([$tcid, $filtered->{$tcid}], [qw(*tcid *complemented)]);
#    <STDIN>;
  }
}







#==========================================================================
#Download the sequences of the multicomponent systems and put them all
#in a single file that will be used as a query for blastp searches.

sub downloadSeqsForMCS {

  my ($tcids, $all_tcdb, $outSeqsFile) = @_;

#  print Data::Dumper->Dump([$all_tcdb ], [qw(*all_tcdb )]);
#  exit;

  unless (-f $outSeqsFile && ! (-z $outSeqsFile)) {

    #Save accessions to a file
    my $accFile="$seqDir/mcs_tcid.acc";
    open (my $fh, ">", $accFile) || die $!;
    print $fh join ("\n", keys %{ $tcids }), "\n";
    close $fh;

    #Download systems now
    system "extractTCDB.pl -i $accFile -o $seqDir -f fasta";

    #Verify that each system was subbessfully downloaded
    foreach my $system (keys %{ $tcids }) { #Based on GBLAST output

      my $seqFile = "$seqDir/tcdb-${system}.faa";
      die "Could not download sequences for system:  $system" unless  (-f $seqFile && ! (-z  $seqFile));
    }

    #Put all sequences in a single file
    system "cat $seqDir/tcdb-*.faa > $outSeqsFile";

    unless (-f $outSeqsFile && ! (-z $outSeqsFile)) {
      die "Error: failed to generate sequence file -> $outSeqsFile";
    }
  }
}





#==========================================================================
#Run quod

sub run_quod {

  my ($tcid, $q, $s) = @_;


  my $pdir   = "$plotDir/${q}_vs_$s";
  my $htfile = "$pdir/ssearch_${s}_vs_${q}/report.html";
  my $cmd    = qq(examineGBhit.pl -q $s -s $q -t $tcid -o $pdir -bdb $blastdb);
  system $cmd unless (-f $htfile && !(-z $htfile));
#  print "$cmd\n";
#  exit;


  if (-f $htfile && !(-z $htfile)) {
    return 1;
  }
  else {
    return undef;
  }
}




#===========================================================================
#Read command line and print help


sub read_command_line {

    print_help() unless (@ARGV);

    my $status = GetOptions(
      "gdir|genome-dir=s"       => \&read_genome_dir,
      "gacc|genome-accession=s" => \$gnmAcc,
      "rs|replicon-structure=s" => \$repForm,
      "gb|gblast=s"             => \&read_gblast_file,
      "dbn|blastdb-name=s"      => \&read_blastdb_name,
      "tcs|tcdb-faa=s"          => \$tcdbFaa,
      "bo|tcblast-overwrite!"   => \$blastOverWrite,
      "of|output-format=s"      => \$outFmt,
      "p|prog=s"                => \&read_prog,
      "o|outdir=s"              => \$outdir,
      "d|max-gene-dist=i"       => \$geneDis,
      "f|fusion-len-ratio=f"    => \$fusion_lenRatio,
      "lp|large-protein=i"      => \$large_protein,
      "sp|small-protein=i"      => \$small_protein,
      "hc|high-coverage=f"      => \$high_coverage,
      "c|mid_coverage=f"        => \$mid_coverage,
      "mc|min-coverage=f"       => \$min_coverage,
      "id|identity=f"           => \$min_identity,
      "e|evalue=f"              => \$minEvalDiscard,
      "a|min-aln-length=i"      => \$min_alignment,
      "h|help"                  => sub { print_help(); },
       "<>"                     => sub { die "Error: Unknown argument: $_[0]\n"; });
    exit unless ($status);

    #validate replicon structure
    die "Unknown replicon structure: $repForm -> " unless ($repForm =~ /(circular|linear)/);

    #Make sure there is a gblast file
    die "GBLAST output file (with extension .tsv) is mandatory -> " unless ($gblastFile);

    #Check that the genome directory exists
    die "Genome directory not found: $gnmDir -> " unless (-d $gnmDir);

    #Check that the feature table exists (to check genomic context)
    $gnmFeatTable = "$gnmDir/${gnmAcc}_feature_table.txt.gz";
    die "Genome feature table  file not found: $gnmFeatTable -> " unless (-f $gnmFeatTable);


    #Check that proteome file exists
    $file_proteome = "$gnmDir/${gnmAcc}_protein.faa.gz";
    die "Proteome file not found: $file_proteome -> " unless (-f $file_proteome);

    system "mkdir -p $outdir" unless (-d $outdir);
}



#==========================================================================
#Read option -p 

sub read_prog {

  my ($opt, $value) = @_;

  my $tmp = lc $value;

  unless ($tmp =~ /^(blastp|ssearch36)$/) {
    die "Error: illegal option value ($value) for -p. Use 'blastp' or 'ssearch36' (case insensitive)\n";
  }

  $prog = $tmp;
}


#==========================================================================
#Read the -gblast option

sub read_gblast_file {

  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
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



#==========================================================================
#Read option 

sub read_proteome {
    my ($opt, $value) = @_;

    unless (-f $value) {
      die "Error: file with the proteome not found: -> $value";
    }

    unless ($value =~ /.+\.(faa|fasta|gz|bz2)$/) {
      die "Error: file does not have the right extensions (fasta|faa|.gz|.bz2)";
    }

    $file_proteome = $value;
}



#==========================================================================
#read option

sub read_blastdb_name {
    my ($opt, $value) = @_;

    unless (length($value) >= 3) {
      die "Error: name of blastdb must be at least 3 characters long -> $value";
    }

    $blastdb_name = $value;
}





sub print_help {

    my $help = <<"HELP";

This program performs the analysis of multicomponent transport systems in a genome.
BlastP searches of all components are carried out followed by correlation with
genomic context. Finally the program generated the hydropathy plots for the best 
candidates.

The main inputs are 1) the output file of gblast in tab-delimited format; 2)
the downloaded NCBI folder of the genome of interest containing the Protein 
sequences of the genome; and 3) The sequences

NOTE: This program works currently for genomes with only one replicon!

Options:

-gdir, --genome-dir { path }  (Mandatory)
  Directory with all the assembly files  as downloaded from NCBI for the genome
  under analysis.

-gacc, --genome-accesssion { string } (Mandatory)
  The prefix of all the files in the genome directory. This is normally the name of
  the folder in NCBI. For example:  GCA_000995795.1_ASM99579v1

-rs, --replicon-structure {string}  (Optional; Default: circular)
  The structure of the replicon under analysis. Only the following values
  are accepted: circular or linear

  Note: use 'linear' when dealing with incomplete genomes as the sequence
        will be in multiple scaffolds.

-gb, --gblast { file } (Mandatory)
   The GBLAST output file in tab-separated format.

-dbn, --blastdb-name {string}  (Optional; Default: proteome)
  Name of the blast database that will be generated with the proteome
  of the genome.

-bo, --tcblast-overwrite  (Optional; by default the DB is not overwritten)
  If given, the blast DB with all the protein content in TCDB will be overwritten.
  The location of this database is ~/db/blastdb. By default the database is not
  overwritten, so be careful and see if your current database is not too old.

-of, --output-format { string } (default: tsv)
  The format of the text file that will be generated with the final report.
  By default the output is tab-separated values (tsv), but it can also be
  comma-separated values (csv).

-o, --outdir { path } (Optional; Default: MultiComponentSystemsAnalysis)
  Directory where the results of the program will be stored.

-d, --max-gene-dist { int } (Optional; Default: 20)
  Maximum distance in number of genes to consider two genes as in the
  neighborhood.

-f, --fusion-len-ratio { float } (Optional; Default: 1.8)
  Minimal length ratio beteen the query and subject proteins that specifies
  al least how much larger one protein must be must be relative to the other
  in order to consider the larger protein a fusion protein.
  (This option is experimental and may disappear

-lp, --large-protein { int } (Optional; Default: 200)
  Minimal length to identify large proteins.

-sp, --small-protein { int } (Optional; Default: 100)
  Maximal length to consider a protein small.

-hc, --high-coverage { float } (Optional; Default: 60.0);
  Assumption of good coverage between query and subject to consider
  two proteins homologous (value must be grater than for option -c).

-c, --mid-coverage { float } (Optional; Defaul: 50.0)
  Assumption of OK coverage to consider two proteins homologous
  (value must be greater than for option -mc).

-mc, --min-coverage { float } (Optional;  Default: 40.0)
  Assumption of poor coverage between two proteins but still
  it may be worth taking a look at the alignment
  (Value must be larger than 20.0)

-id, --identity { float }  (Optional;  Default: 18.0)
  Minimal alignment identity to consider a hit in the results.

-e, --evalue { float } (Optional;  Default: 1e-3)
  Worst acceptable evalue for alignments with high coverage.

-a, --min-aln-length { int } (Optional; Default: 55)
  Minimal accepted alignment length.

-h, --help
  Print this help. This option takes precedence over anyother option.

HELP

    print $help;
    exit;
}

