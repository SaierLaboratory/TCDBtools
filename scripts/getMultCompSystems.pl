#!/usr/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use List::Util qw(max min);

use TCDB::Assorted;

###########################################################################
#
# In genome analysis, after running GBLAST there will be multicomponent
# systems with missing components. This is approached by running blastp
# with the the TCDB proteins of the missing components against the genome
# of interest, and by checking whether proteins matching other components
# withing the same system are located nearby in the genome.
#
# This script will run blastp of a series of TCDB proteins expected to be
# part of multicomponent systems. Top hits that will be selected must meet
# the following criteria to be called "Candidates":
#
# 1) id% >= 24%
# 2) coverage of at least 45% of either the query or the subject
# 3) If the smaller protein is less than 100aa long, the coverage muust
#    be at least 70% of the smallest protein
# 4) E-value < 1.0
#
# After the top hits are filtered, and the genome neighborhood is determined
# for each protein, then the program WHAT (quod.py) is run to generate plots
# for the alignment reported by blast and the full sequence of the query
# (TCDB protein) and the subject (protein in the reference genome).
#
#--------------------------------------------------------------------------
#
# Writen by: Arturo Medrano
# Date: 08/2017
#
###########################################################################



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
my $blastOverWrite  = 0;      #Recreate the blastDB of all TCDP proteins
my $outFmt          = "tsv";  #can also be "csv"
my $outdir          = "./MultiComponentSystemsAnalysis";


#Variables that will be used to select good blastp hits (Candidates)
my $geneDis         = 15;    #How many genes away to look for neighbors
my $fusion_lenRatio = 1.8;   #At what point a protein is considered fused to another
my $large_protein   = 200;
my $small_protein   = 100;

my $high_coverage   = 50.0;  #Assumption of good coverage for candidates 55
my $mid_coverage    = 40.0;  #Assumption of ok coverage  for candidates  50
my $min_coverage    = 30.0;  #Assumption of low coverage for candidate   40
my $minCovDiscard   = 20.0;  #don't look at  blast hits if one of the proteins has this or lower coverage

my $worst_evalue    = 1e-3;  #E-value for high coverage >=50%
my $badCov_evalue   = 1e-9;  #E-value for really bad coverage >=20%
my $shortCov_evalue = 1e-7;  #E-value for short coverage alignments >=30%
my $midCov_evalue   = 1e-5;  #E-value for mid coveage proteins >=40%

my $min_identity    = 20.0;  #minimum accepted identity to consider a good hit

my $compStats       = 'F';   #Use composition statistis to filter blast output
my $min_alignment   = 55;    #minimal alignment length





read_command_line();


#print Data::Dumper->Dump([$gnmAcc, $gnmDir, $repForm, $blastdb_name, $file_proteome,
#			  $gnmFeatTable, $gblastFile, $blastOverWrite, $outFmt, $outdir,
#			  $geneDis, $fusion_lenRatio, $large_protein, $small_protein,
#			  $high_coverage, $mid_coverage, $min_coverage, $min_identity,
#			  $max_evalue, $min_alignment],
#			 [qw(*gnmAcc *gnmDir *repForm *blastdb_name *file_proteome
#			  *gnmFeatTable *gblastFile *blastOverWrite *outFmt *outdir
#			  *geneDis *fusion_lenRatio *large_protein *small_protein
#			  *high_coverage *mid_coverage *min_coverage *min_identity
#			  *max_evalue *min_alignment)]);
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

my $blastoutDir = "$outdir/blast_output";
system "mkdir $blastoutDir" unless (-d $blastoutDir);


my $outfile = "$outdir/plot_index.html";



#==========================================================================
#Get the TCIDs of multi-component systems in TCDB. This information will
#be used to extract the TCIDs of all multicomponent systems in TCDB.

print "Downloading TCDB systems.\n";


#Download the TCDB database in fasta format
my $tcdbFaa = "$outdir/tcdb/tcdb.faa";
system "extractFamily.pl -i tcdb -o $outdir/tcdb -f fasta" unless (-f $tcdbFaa);
die "TCDB seqs not found: $tcdbFaa" unless (-f $tcdbFaa);

my $multSystems = getModeSystems($tcdbFaa, 'multi');

#print Data::Dumper->Dump([$multSystems ], [qw(*multSystems )]);
#exit;



#==========================================================================
#Get the multicomponent systems identified by GBLAST

print "Extracting multicomponent systems in TCDB\n";


#First get all TCIDs of multicomponent systems found by GBLAST.
#NOTE: Don't filter by threshold to capture everything that GBLAST found.
my %gblastTCIDs = ();
getGBLASTmultiComponentTCIDs($gblastFile, 10, $multSystems, \%gblastTCIDs);

#print Data::Dumper->Dump([\%gblastTCIDs ], [qw(*gblastTCIDs)]);
#exit;





#==========================================================================
#Download the sequences of the multicomponent systems and put them all
#in a single file that will be used as a query for blastp searches.

print "Downloading sequences for multicomponent systems.\n";


my $file_blastp_query = "$seqDir/systems.faa";

downloadSeqsForMCS(\%gblastTCIDs, $file_blastp_query);





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
  system "extractFamily.pl -i tcdb -o $tcdbBlastDBdir -f blast";
}




#==========================================================================
#Create the blast database with the reference genome. This will
#allow searching for missing components and extracting sequences
#of individual proteins in the genome

print "Creating blast database for genome!\n";



#Prepare command
my $blastdb = "$blastdb_dir/$blastdb_name";
my $bdb_cmd =
  qq(makeblastdb -dbtype prot -input_type fasta -title 'Sequences: $blastdb_name' ) .
  qq(-parse_seqids -hash_index -out $blastdb);


#generate the database
unless (-f "${blastdb}.phr") {

  if ($file_proteome =~ /.+\.gz$/) {
    print qq(zcat $file_proteome | $bdb_cmd -in - \n);
    system qq(zcat $file_proteome | $bdb_cmd -in -);
  }
  elsif ($file_proteome =~ /.+\.bz2$/) {
    system qq(bzcat $file_proteome | $bdb_cmd -in -);
  }
  else {
    system qq($bdb_cmd -in $file_proteome);
  }


  #Verify that the database was successfully generated
  die "Blast databse $blastdb_name was not created" unless (-f "${blastdb}.phr");
}




#==========================================================================
#blastp the TCDB proteins (multicomponent systems) against the genome


print "Blasting multi-component systems against the genome!\n";

my $blast_file = "$blastoutDir/blastOutput.csv";


#Prepare the command
my $outfmt = qq(10 qacc sacc qlen slen evalue pident length qcovs qstart qend sstart send qseq sseq);
my $blast_cmd =
  qq(blastp -db $blastdb -query $file_blastp_query -evalue $worst_evalue -use_sw_tback ) .
  qq( -comp_based_stats $compStats -max_hsps 1 -outfmt '$outfmt' -out $blast_file );


#Run blastp
system $blast_cmd unless (-f $blast_file && !(-z  $blast_file));
die "Error: Blastp output file not found -> $blast_file " unless (-f $blast_file && !(-z  $blast_file));




#==========================================================================
#Parse the blastp output file, select top blastp hits, and get the DNA
#sequences of the aligned regions


print "Extracting best blast hits\n";


my %filteredBlast =();
filterBlastOutput($blast_file, \%filteredBlast, \%gblastTCIDs, $multSystems);


#print Data::Dumper->Dump([\%filteredBlast], [qw(*filterdBlast)]);
#exit;



#==========================================================================
#Now check the genomic context of all multicomponent systems and
#idenfiy those TCDB components that  share genomic neighborhood
#in the genome under analysis.

print "Analyzing genome context\n";

verify_neighborhood(\%filteredBlast, $gnmFeatTable);

#print Data::Dumper->Dump([\%filteredBlast ], [qw(*filteredBlast )]);
#exit;






#==========================================================================
#Create text file with top hits and genome context results


print "Saving results to text and html reports\n";

generate_reports(\%filteredBlast);

print "Done!\n";






###########################################################################
##                                                                       ##
##                        Subroutine definitions                         ##
##                                                                       ##
###########################################################################



#==========================================================================
#Generate the text and html output files with the results of the analysis


sub generate_reports {

  my $systems = shift;

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
  my @hHeader = qw (tcid  query_accession  status  subject_accession  hydropathy
		    query_length  subject_length  evalue  perc_idenity
		    alignment_length  query_coverage  subject_coverage
		    neighbors);
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
  my $noHitColSpan = scalar(@hHeader) - 2;


  #Get the text header
  my @tHeader = qw (tcid  query_accession  status  subject_accession
		    query_length  subject_length  evalue  perc_idenity
		    alignment_length  query_coverage  subject_coverage
		    neighbors);
  my $txtHeader = join ($sep, @tHeader) . "\n";



  open  (my $htmh, ">", $htmlFile) || die $!;
  print $htmh $htmlHeader;




  open (my $outh, ">", $outfile) || die $!;
  print $outh $txtHeader;


  my $rowCnt = 1;

  foreach my $tcid (sort_by_system keys %{ $systems }) {

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
	my $aln       = "";
	my $qcov      = "";
	my $scov      = "";
	my $neighbors = "";


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
	}

	#print the text version of the row
	print $outh "$tcid${sep}$tcAcc${sep}$status${sep}$sacc${sep}$qlen${sep}$slen${sep}$eval${sep}$ident${sep}$aln${sep}$qcov${sep}$scov${sep}$neighbors\n";


	if ($status =~ /(GBLAST|Candidate|RawBlast)/) {

	  #generate hydropathy plots and print the corresponding HTML row.
	  my $plotFile = "${tcAcc}_vs_${sacc}.html";
	  print "  Creating plots comparing:  ${tcAcc} vs ${sacc}\n";


	  my $good = run_quod($tcid, $tcAcc, $sacc, $hit->{qstart}, $hit->{qend}, $hit->{sstart}, $hit->{send}, $hit->{qseq},$hit->{sseq});
	  unless ($good) {
	    print Data::Dumper->Dump([$tcid, $tcAcc, $sacc, $hit->{qstart}, $hit->{qend}, $hit->{sstart}, $hit->{send}, $hit->{qseq},$hit->{sseq} ],
				     [qw(*tcid *tcAcc *sacc *qstart *qend *sstart *send *qseq *sseq )]);
	    die "Could not generate plot for: ${tcid}-$tcAcc vs $sacc -> ";
	  }
	  my @data = ($tcid, $tcAcc, $status, $sacc, $qlen, $slen, $eval, $ident, $aln, $qcov, $scov, $neighbors);
	  my $hitRaw = getMatchRowHTMLstring($sysURL, $accURL, $ncbiURL, $plotFile, $rowColor, \@data);

	  print $htmh $hitRaw;

	}

	#print html row with no hit
	else {
	  my $noHitRaw = getNoHitRowHTMLstring($sysURL, $accURL, $tcid, $tcAcc, $status, $noHitColSpan, $rowColor);
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
#get a string with the HTML code for a row of data

#Print an HTML raw with the info of blastp match
sub getMatchRowHTMLstring {

  my ($sysURL, $accURL, $ncbiURL, $plotFile, $color, $data) = @_;


  my ($tcid, $tcAcc, $status, $sacc, $qlen, $slen, $eval, $ident, $aln, $qcov, $scov, $neighbors) = @$data;

  my $identity = sprintf ("%.1f", $ident);

  $neighbors =~ s/\|/\<br \\>/g;

  my $row = <<"NOHIT";
    <tr>
      <td style='text-align:left;   background-color:$color;'><a href='$sysURL${tcid}' target='_blank'>$tcid</a></td>
      <td style='text-align:center; background-color:$color;'><a href='$accURL${tcAcc}' target='_blank'>$tcAcc</a></td>
      <td style='text-align:center; background-color:$color;'>$status</td>
      <td style='text-align:center; background-color:$color;'><a href='${ncbiURL}$sacc' target='_blank'>$sacc</a></td>
      <td style='text-align:center; background-color:$color;'><a href='html/$plotFile' target='_blank'>plots</a></td>
      <td style='text-align:right;  background-color:$color;'>$qlen</td>
      <td style='text-align:right;  background-color:$color;'>$slen</td>
      <td style='text-align:center; background-color:$color;'>$eval</td>
      <td style='text-align:right;  background-color:$color;'>$identity</td>
      <td style='text-align:right;  background-color:$color;'>$aln</td>
      <td style='text-align:right;  background-color:$color;'>$qcov</td>
      <td style='text-align:right;  background-color:$color;'>$scov</td>
      <td style='text-align:center; background-color:$color;'>$neighbors</td>
    </tr>
NOHIT

  return $row;

}




#==========================================================================
#get the HTML string to print a row that had no hits

sub getNoHitRowHTMLstring {

  my ($tcURL, $accURL, $tcid, $tcAcc, $status, $colSpan, $color) = @_;

  my $row = <<"NOHIT";
    <tr>
      <td style='text-align:left; background-color:$color;'><a href='$tcURL${tcid}' target='_blank'>$tcid</a></td>
      <td style='text-align:center; background-color:$color;'><a href='$accURL${tcAcc}' target='_blank'>$tcAcc</a></td>
      <td style='text-align:left; background-color:$color; padding-left:10px;' colspan='$colSpan'>$status</td>
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

#  print Data::Dumper->Dump([\@cds ], [qw(*cds )]);
#  exit;


  #-----------------------------------------------------------------
  #Get the ranked position of each CDS

  my $pos = 1;
  foreach my $orf (sort { $a->{start} <=> $b->{start} } @cds) {
    $pos2cds{ $pos }         = $orf;
    $cds2pos{ $orf->{acc} } = $pos;
    $pos++;
  }

#  print Data::Dumper->Dump([\%pos2cds], [qw(*pos2cds )]);
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
    #and their position rank in the replicon
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


      #register whether there were neighbors or not
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
#Parse the blastp output file, select top blastp hits, and get the DNA
#sequences of the aligned regions

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
    my $qcov_tmp = ($qend - $qstart) / $qlen * 100;
    my $qcovs = ($qcov_tmp >= 100.0)? 100.0 : $qcov_tmp;


    #The coverage of the subject protein
    my ($sacc, $ver) = split(/\./, $sseqid);
    my $scov_tmp = ($send - $sstart) / $slen * 100;
    my $scov = ($scov_tmp >= 100.0)? 100.0 : $scov_tmp;


    #*** Ignore hits with not enough coverage and a minimum alignment length ***
    next if (($qcovs < $minCovDiscard || $scov < $minCovDiscard) ||  #At least one sequence is covered very poorly (discard!)
	     ($qcovs < $min_coverage  && $scov < $min_coverage)  ||  #At least one sequence must have the minimum coverage for candidates
	     $length < $min_alignment);                              #Check for the minimum alignment size


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


    #Very bad coverage (<30%) in at least one protein
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

  my ($tcids, $outSeqsFile) = @_;

  unless (-f $outSeqsFile && ! (-z $outSeqsFile)) {

    foreach my $system (keys %{ $tcids }) {

      #extract the sequences associated with this system
      my $seqFile = "$seqDir/family-${system}.faa";
      unless (-f $seqFile && ! (-z  $seqFile)) {
	system "extractFamily.pl -i $system -o $seqDir -f fasta";
      }
      die "Could not download sequences for system:  $system" unless  (-f $seqFile && ! (-z  $seqFile));
    }

    #Put all sequences in a single file
    system "cat $seqDir/family-*.faa > $outSeqsFile";

    unless (-f $outSeqsFile && ! (-z $outSeqsFile)) {
      die "Error: failed to generate sequence file -> $outSeqsFile";
    }
  }
}





#==========================================================================
#Run quod

sub run_quod {

  my ($tcid, $q, $s, $qs, $qe, $ss, $se, $qseq, $sseq) = @_;


  #extract sequences for query and subject
  my $tcAcc = "${tcid}-$q";
  extract_full_sequences($tcAcc,$s);



  #-----------------------------------------------------------------
  #Run quod for the alignment

  my $alnFig = "${q}_vs_${s}.png";
  my $cmd1 = qq(quod.py -q -s -l "$q (red) and $s (blue)"  -o $plotDir/$alnFig --xticks 50 --width 15 $qseq $sseq );
  system $cmd1 unless (-f "$plotDir/$alnFig");
  return undef unless (-f "$plotDir/$alnFig");


  #-----------------------------------------------------------------
  #Run quod for the full sequencess of the query and subject proteins

  my $qName = "${q}_vs_${s}_qAln";
  my $cmd2 = qq(quod.py -q -f -l "$tcAcc"  -o $plotDir/${qName}.png --width 15 --color red --xticks 100 -W $qs,2 $qe,-2 --  $seqDir/${tcAcc}.faa);
  system $cmd2 unless (-f "$plotDir/${qName}.png");
  return undef unless (-f "$plotDir/${qName}.png");


  my $sName = "${q}_vs_${s}_sAln";
  my $cmd3 = qq(quod.py -q -f -l "$s"  -o $plotDir/${sName}.png --width 15 --color blue --xticks 100 -W $ss,2 $se,-2 --  $seqDir/${s}.faa);
  system $cmd3 unless (-f "$plotDir/${sName}.png");
  return undef unless (-f "$plotDir/${sName}.png");



  #-----------------------------------------------------------------
  #Generate the html file that will display the three plots

  generate_html_page ("${q}_vs_${s}", $qName, $sName, $tcAcc, $s);

  return 1;
}



#==========================================================================
#Generate the html file  with the three plots for a given blastp hit

sub generate_html_page {

  my ($alnImg, $qImg, $sImg, $q, $s) = @_;


  my $html = <<HTML;
<html>
  <head><title>Hydropaty: $q vs $s</title></head>
  <body>
    <table style="width:100%">
      <tr>
        <td><img src="../plots/${qImg}.png" alt="Query protein"></td>
        <td><img src="../plots/${sImg}.png" alt="Subject protein"></td>
      </tr>
      <tr>
        <td colspan="2" style="text-align: center;"><img src="../plots/${alnImg}.png" alt="Aligned region"></td>
      </tr>
  </body>
</html>
HTML


  my $outPath = "$htmlDir/$alnImg.html";
  open (my $fh, ">", $outPath) || die $!;
  print $fh $html;
  close $fh;
}




#==========================================================================
#Extract the full sequences of the query and subject proteins
#Examples:  1.F.1.1.1-O00161   AKM80767.1



sub extract_full_sequences {

  my ($tcid, $s) = @_;


  my $tcid_seq = "$seqDir/${tcid}.faa";
  my $acc_seq  = "$seqDir/${s}.faa";

  #extract the query secuence from tcdb and the subject from the custom blastdb
  my $cmd1 = qq(blastdbcmd  -db $tcdbBlastDBdir/tcdb  -entry $tcid  -out $tcid_seq);
  system "$cmd1" unless (-f $tcid_seq && !(-z $tcid_seq));
  die "Could not extract sequence for $tcid" unless (-f $tcid_seq && !(-z $tcid_seq));


  my $cmd2 = qq(blastdbcmd  -db $blastdb  -entry ${s}.1  -out $acc_seq);
  system "$cmd2" unless (-f $acc_seq && !(-z $acc_seq));
  die "Could not extract sequence for $s" unless (-f $acc_seq && !(-z $acc_seq));
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
      "bo|tcblast-overwrite!"   => \$blastOverWrite,
      "of|output-format=s"      => \$outFmt,
      "o|outdir=s"              => \$outdir,

      "d|max-gene-dist=i"       => \$geneDis,
      "f|fusion-len-ratio=f"    => \$fusion_lenRatio,
      "lp|large-protein=i"      => \$large_protein,
      "sp|small-protein=i"      => \$small_protein,
      "hc|high-coverage=f"      => \$high_coverage,
      "c|mid_coverage=f"        => \$mid_coverage,
      "mc|min-coverage=f"       => \$min_coverage,
      "id|identity=f"           => \$min_identity,
      "e|evalue=f"              => \$worst_evalue,
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

The main inputs are the output file of gblast in tab-delimited format and
the downloaded NCBI folder of the genome of interest.

NOTE: This program works currently for genomes with only one replicon!

Options:

-gdir, --genome-dir { path }  (Mandatory)
  Directory with all the files as downloaded from NCBI for the genome under analysis.

-gacc, --genome-accesssion { string } (Mandatory)
  The prefix of all the files in the genome directory. This is normally the name of
  the folder in NCBI. For example:  GCA_000995795.1_ASM99579v1

-rs, --replicon-structure {string}  (Optional; Default: circular)
  The structure of the replicon under analysis. Only the following values
  are accepted: circular or linear

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

-d, --max-gene-dist { int } (Optional; Default: 15)
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

-hc, --high-coverage { float } (Optional; Default: 50.0);
  Assumption of good coverage between query and subject to consider
  two proteins homologous (value must be grater than for option -c).

-c, --mid-coverage { float } (Optional; Defaul: 40.0)
  Assumption of OK coverage to consider two proteins homologous
  (value must be greater than for option -mc).

-mc, --min-coverage { float } (Optional;  Default: 30.0)
  Assumption of poor coverage between two proteins but still
  it may be worth taking a look at the alignment
  (Value must be larger than 20.0)

-id, --identity { float }  (Optional;  Default: 20.0)
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

