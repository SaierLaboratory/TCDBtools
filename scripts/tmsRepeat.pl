#!/usr/bin/env perl -w

use warnings;
use strict;
use Data::Dumper;

$Data::Dumper::Deepcopy = 1;
$Data::Dumper::Indent   = 1;
#$Data::Dumper::Purity   = 0;
$Data::Dumper::Sortkeys = 1;

use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;


use TCDB::CheckDependencies;
use TCDB::Assorted;



#==========================================================================
#Check dependencies

my @dependencies = ("water", "ssearch36", "extractFamily.pl", "tmsplit", "quod.py");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;



#==========================================================================
#Read command line options

my $gs_infile      = "";
my $infileFmt      = "hmmtop"; #The other option is 'tms' which is the ID and TMS
my $gs_idFormat    = "";
my $gs_repUnit     = 0;
my $gs_seqDir      = "";
my $gs_tail        = 5;
my $gs_evalue      = 0.1;
my $gs_coverage    = 0.8;
my $gs_identity    = 0.25;
my $gsatShuffles   = 1000;
my $min_gsat_score = 4.0;

my $compStatsFlag  = 1;
my $compStats      = "";
my $outdir         = "repeats";
my $repDir         = "reports";
my $seqDir         = "sequences";
my $alignDir       = "alignments";
my $plotsDir       = "plots";
my $goodHitsOnly   = 1;  #print only significant results, ignore everything else


#all (all sequences in output file)
#each (generate one directory per sequence.. for better organization)
#debug (it will print the contents of the hash table one sequences at a time)
my $mode           = "all";

read_command_line_arguments();

#print Data::Dumper->Dump([$gs_infile, $gs_idFormat, $gs_repUnit, $gs_seqDir,
#			  $gs_tail, $gs_evalue, $gs_coverage, $gs_identity, $gsatShuffles, $compStatsFlag, $compStats],
#			 [qw(*infile *idFormat *repUnit *seqDir *tail *evalue
#			     *coverage *identity $gsatShuffles *compStatFlag *compStats)]);
#exit;

# ssearch36 -p -k 1000 -z 11 -E 1.0 -s BL62 -W 0 4.B.1_4tms_all/sequences/4.B.1.1.2-Q4QLL1_bundle1.faa 4.B.1_4tms_all/sequences/lib_4.B.1.1.2-Q4QLL1_bundle1.faa


#==========================================================================
#Read file with coordinates of TMSs and verify that the sequences are
#available

my %gh_tms = ();

read_tms_coordinates_file($gs_infile, \%gh_tms);

#print Data::Dumper->Dump([ \%gh_tms], [qw(*tms )]);
#exit;


#===========================================================================
#Main Output directory

#Root directory for all results
system "mkdir -p $outdir" unless (-d $outdir);
die "Could not generate output directory: $outdir" unless (-d $outdir);



#==========================================================================
#Search for repeats inside query sequences

my %results = ();
my %origSeqLength = ();  #To calculate x-ticks spacing in hydropathy plots

foreach my $ls_sid (keys %gh_tms) {

  my %gh_bundleSeqs = ();
  my %gh_topHits    = ();


  print "Processing: $ls_sid\n";


  #Clean results if one output directory is generated per input sequence
  %results = () if ($mode eq 'each');


  #Cut sequences in non overlaping regions with as many TMS as the
  #repeat unit we want to find.
  cut_seq_in_tms_regions ($ls_sid, $gs_repUnit, \%gh_tms,  \%gh_bundleSeqs);


#  print Data::Dumper->Dump([\%gh_bundleSeqs ], [qw(*bundleSeqs)]);
#  <STDIN>;


  #run ssearch to find potential repeats.
  align_bundles($ls_sid,\%gh_bundleSeqs, \%gh_topHits);


#  print Data::Dumper->Dump([\%gh_topHits ], [qw(*topHits )]);
#  <STDIN>;


  #Collect results for final table
  $results{$ls_sid} = \%gh_topHits;

  #present results per input sequence to verify everything looks fine.
  if ($mode eq 'debug') {
    print Data::Dumper->Dump([\%gh_topHits], [qw(*topHits)]);
    <STDIN>;
  }

  print_reports(\%results) if ($mode eq 'each');
}





#===========================================================================
#Print final results in summarized or detailed format

#print Data::Dumper->Dump([\%results ], [qw(*results )]);
#<STDIN>;

print_reports(\%results) if ($mode eq 'all');





###########################################################################
#                                                                         #
#                        Subroutine definition                            #
#                                                                         #
###########################################################################


#print final_report

sub print_reports {

  my $res = shift;


  #Get the directory where reports will be saved
  my $reportDir = undef;
  if ($mode eq 'all') {
    $reportDir = getReportsDir();
  }
  else {

    #one id per report
    my @ids = keys %$res;
    my $seqId = $ids[0];

    $reportDir = getReportsDir($seqId);
  }
  die "Error: invalid report dir" unless ($reportDir);


  my $sumFile     = "$reportDir/repeats_summary_report.txt";
  my $detailsFile = "$reportDir/repeats_detailed_report.txt";
  my $htmlFile    = "$reportDir/report.html";


  open (my $htmlfh, ">", $htmlFile) || die $!;

   my $htmlHeader = <<HEADER;
<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />

    <style type="text/css">

.label {
   text-align: right;
   width: 50px;
}

.data {
   text-align: left;
   padding-left: 8px;
   width: 100px;
}

.uline {
   text-decoration: underline;
}

.seq {
   border: 2px solid black;
   height: 70px;
   width:  100%;
   overflow-x: auto;
   overflow-y: hidden;
   margin: 1em 0;
   background: gray;
   color: white;
}

img {
   display: block;
   margin-left: auto;
   margin-right: auto;
   height: 250px;
   width:  auto;
   max-width: 1500px;
   max-height: 300px;
}

    </style>
    <title>Inferring repeats of $gs_repUnit TMS</title>
  </head>
  <br />
  <h1 style='text-align:center'>Inferred Repeats Based On ${gs_repUnit}-TMS Bundles</h1>
  <body>

HEADER

  print $htmlfh $htmlHeader;
  open (my $sumh, ">", $sumFile) || die $!;
  open (my $deth, ">", $detailsFile) || die $!;


  #Header for summary table
  print $sumh "#Accession\tQ_bundle\tS_bundle\tQ_len\tS_len\tE-value\tIdentity\tGSAT\tAln_len\tQ_cov\tS_cov\n";


#  print Data::Dumper->Dump([$res ], [qw(*res )]);
#  <STDIN>;


 P:foreach my $id (sort {$a cmp $b} keys %$res) {

    #Jump to next result if there are NO hits for this protein and
    #ONLY good hits are going to be recorded.
    unless (%{ $res->{$id} }) {
      next P if ($goodHitsOnly);
    }


    print $deth "===========================================================================\n";
    print $htmlfh "   <br /><hr style=\"border-style:solid; border-width:5px; color:black;\"/>\n";

    #There must be results to continue
    unless (%{ $res->{$id} }) {
      print $sumh   "$id\tNo_hits\n";
      print $deth   "$id\tNo_hits\n\n\n";
      print $htmlfh "    <h2 style=\"text-align:center;\">$id</h2>\n    <p><b>No candidate repeats found</b></p>\n";
    }

    print $deth "$id\n\n";
    print $htmlfh "    <h2 style=\"text-align:center;\">$id</h2>\n";




    #get the long bundle names
  BS:foreach my $bundleName (sort {$a cmp $b} keys %{ $res->{$id} }) {

    BN:foreach my $bundleNumber (sort {$a <=> $b} keys %{ $res->{$id}->{$bundleName} }) {

	my $qName = $res->{$id}->{$bundleName}->{$bundleNumber}->{qName};
	my $qLen  = $res->{$id}->{$bundleName}->{$bundleNumber}->{qLen};

	#Each of the hits for this bundle
	my @hits_tmp = @{ $res->{$id}->{$bundleName}->{$bundleNumber}->{hits} };

	#To get rid of a warning when there is only one hit.
	my @hits = (scalar (@hits_tmp) > 1)?
	  sort {$a->{hName} cmp $b->{hName}} @hits_tmp : @hits_tmp;

	foreach my $hit (@hits) {

	  my $hName  = $hit->{hName};
	  my $hLen   = $hit->{hLen};

	  my $evalue = sprintf ("%.1e", $hit->{hEvalue});
	  my $ident  = sprintf ("%.1f", $hit->{hId} * 100);
	  my $sim    = sprintf ("%.1f", $hit->{hSim} * 100);
	  my $gsat   = sprintf ("%.1f", $hit->{gsat});

	  my $alnLen = $hit->{alnLen};
	  my $qCov   = sprintf("%.1f", $hit->{qCov} * 100);
	  my $hCov   = sprintf("%.1f", $hit->{hCov} * 100);


	  #The alignment
	  my $qstart = $hit->{qstart};
	  my $qend   = $hit->{qend};
	  my $sstart = $hit->{sstart};
	  my $send   = $hit->{send};
	  my $qSeq   = $hit->{qSeq};
	  my $homStr = $hit->{homStr};
	  my $sSeq   = $hit->{sSeq};

	  my $plot   = $hit->{plot};

	  #For summary tab-delimitedfile (everything except the alignment)
	  print $sumh "$id\t$qName\t$hName\t$qLen\t$hLen\t$evalue\t$ident\t$gsat\t$alnLen\t$qCov\t$hCov\n";


	  #Detailed report that includes the alignment
	  print $deth "----------\n";
	  print $deth "$qName ($qLen) vs $hName ($hLen)\n\n";
	  print $deth "E-value: $evalue     Identity: ${ident}%     GSAT: $gsat\n";
	  print $deth "Q_cov: ${qCov}%     S_cov: ${hCov}%     Aln_length: $alnLen\n\n";
	  print $deth "Alignment ($qName|${qstart}-$qend vs $hName|${sstart}-$send):\n$qSeq\n$homStr\n$sSeq\n\n\n";


	  #The HTML report (includes alignment and hydropathy image
	  my $repHit = <<HIT;

    <p><b>$qName ($qLen) vs $hName ($hLen)</b></p>

    <table width="600px" border="0" cellspacing="0" cellpadding="2">
      <tr>
         <td class='label'><b>E-value:</b></td>
         <td class='data'>$evalue</td>
         <td class='label'><b>Identity:</b></td>
         <td class='data'>${ident}%</td>
         <td class='label'><b>Similarity:</b></td>
         <td class='data'>${sim}%</td>
         <td class='label'><b>GSAT:</b></td>
         <td class='data'>$gsat</td>
      </tr>
      <tr>
         <td class='label'><b>Aln:</b></td>
         <td class='data'>$alnLen</td>
         <td class='label'><b>Q_cov:</b></td>
         <td class='data'>${qCov}%</td>
         <td class='label'><b>S_cov:</b></td>
         <td class='data'>${hCov}%</td>
         <td class='label'></td>
         <td class='data'></td>
      </tr>
    </table>

    <p><b>Alignment (</b>$qName:<b class="uline">${qstart}-$qend</b> vs $hName:<b class="uline">${sstart}-$send</b><b>):</b></p>
    <div class='seq'>
    <pre>
$qSeq
$homStr
$sSeq
    </pre>
    </div>
    <a href="$plot" target="_blank"><img src="$plot"/></a>
    <br />
    <hr />

HIT

	  print $htmlfh $repHit;

	} #hit
      } #reference bundle number
    } #Reference bundle name
  } #Query protein

  #Close HTML report
  my $closeRep = <<CLOSE;
  </body>
</html>
CLOSE

  print $htmlfh $closeRep;

  close $sumh;
  close $deth;
  close $htmlfh;
}



#==========================================================================
#Run ssearch36 between the different bundles in a sequence

sub align_bundles {

  my ($seqId, $lhr_bundleSeqFiles, $lhr_topHits) = @_;

  %$lhr_topHits = ();

  #Directory where the sequences of TMS bundles are saved
  my $sequencesDir = undef;
  my $alignmentsDir = undef;
  my $hydroPlotsDir = undef;

  if ($mode eq 'all') {
    $sequencesDir  = getSequencesDir();
    $alignmentsDir = getAlignmentsDir();
    $hydroPlotsDir = getPlotsDir();
  }
  else {
    $sequencesDir  = getSequencesDir($seqId);
    $alignmentsDir = getAlignmentsDir($seqId);
    $hydroPlotsDir = getPlotsDir($seqId);
  }
  die "Error: invalid sequences dir" unless ($sequencesDir);
  die "Error: invalid alignments dir" unless ($alignmentsDir);
  die "Error: invalid plots dir" unless ($hydroPlotsDir);


#  print Data::Dumper->Dump([$lhr_bundleSeqFiles ], [qw(*files )]);
#  <STDIN>;


  #The bundle that will be used as reference for the comparison
  REF:foreach my $bundle (sort {$a <=> $b} keys %$lhr_bundleSeqFiles) {

      my $rFile = "$sequencesDir/" . $lhr_bundleSeqFiles->{$bundle}->[0];


      #Id to name ssearch36 output files
      my $id = $lhr_bundleSeqFiles->{$bundle}->[0];
      $id =~ s/\.faa//;


      #For naming GSAT files (ID of system or protein accession)
      my $tcAcc = ($id =~ /(\S+)_bundle.*/)? $1 : undef;
      die "Could not extract accession from $id!" unless ($id);


#      print Data::Dumper->Dump([$id, $tcAcc ], [qw(*id *tcAcc)]);
#      <STDIN>;


      #--------------------------------------------------------------------
      #Get the non-overlapping bundles to compare them against the
      #reference bundle

      my @cmpFiles = ();

      #Initialize the index to the first non-overlapping bundle
      my $next_bundle_idx = $bundle + $gs_repUnit;

    CMP:while (1) {

	#Exit if next bundle is not in bundles hash
	last CMP unless (exists $lhr_bundleSeqFiles->{$next_bundle_idx});

	#Get file name for this non-overlapping bundle
	my $cmpBundle = $sequencesDir . "/" . $lhr_bundleSeqFiles->{$next_bundle_idx}->[0];
	push (@cmpFiles, $cmpBundle);

	#Update the index to the next non-overlapping bundle
	$next_bundle_idx = $next_bundle_idx + $gs_repUnit;
      }

      #go to next reference bundle if there are no non-overlapping bundles.
      next REF unless (@cmpFiles);


#      print Data::Dumper->Dump([\@cmpFiles ], [qw(*cmpFiles )]);
#      <STDIN>;


      #--------------------------------------------------------------------
      #Now run ssearch36 of the reference bundle against all its
      #non-overlapping bundles

      #put all non-overlapping bundles into a file
      my $libFile = "$sequencesDir/lib_$id.faa";
      my $cmd = "cat " . join(" ", @cmpFiles) . " > $libFile";
      system $cmd;


      #run ssearch36 of $rFile vs @cmpFile
      my $ssearchOut = "$alignmentsDir/ssearch_$id.out";
      my $ssearch_params = qq(-p $compStats -E $gs_evalue -s BL62 -W 0 $rFile $libFile > $ssearchOut);
      system "ssearch36 $ssearch_params" unless (-f $ssearchOut);


#      print Data::Dumper->Dump([$ssearchOut ], [qw(*ssearchOut )]);
#      <STDIN>;


      #---------------------------------------------------------------------
      #Estimate here the spacing between x-ticks for hydropathy plots

      my $protLen = $origSeqLength{$seqId};

      my $xticksSpacing = undef;
      if ($protLen <= 500) {
	$xticksSpacing = 25;
      }
      elsif ($protLen <= 1000) {
	$xticksSpacing = 50;
      }
      else {
	$xticksSpacing = 100;
      }



      #--------------------------------------------------------------------
      #parse ssearch36 output. For BioPerl resouces check:
      #http://search.cpan.org/dist/BioPerl/Bio/SearchIO.pm
      #https://classes.soe.ucsc.edu/bme060/Winter07/bptutorial.html

      my $parser = new Bio::SearchIO (-format => 'fasta', -file => $ssearchOut);


      #put hir the top hits
      my %lh_hits = ();


      while (my $result = $parser->next_result) {


	my $qLen = $result->query_length;
	$lh_hits{$bundle}{qName} = $result->query_name;
	$lh_hits{$bundle}{qLen}  = $qLen;
	$lh_hits{$bundle}{hits}  = [];


      HIT:while (my $hit = $result->next_hit) {

	HSP:while(my $hsp = $hit->next_hsp) {


#	    print Data::Dumper->Dump([$hsp ], [qw(*hsp )]);
#	    <STDIN>;


	    my %tmp = ();

	    my $alnLen  = $hsp->hsp_length;
	    my $hLen    = $hit->length;
	    my $hEvalue = $hsp->evalue;
	    my $hId     = $hsp->frac_identical('total');  #identity in the alignment
	    my $hSim    = $hsp->frac_conserved('total');  #similarity in the alignment


	    #coordinates in the alignment to properly calculate coverages
	    my $qstart  = $hsp->start('query');
	    my $qend    = $hsp->end('query');
	    my $sstart  = $hsp->start('subject');
	    my $send    = $hsp->end('subject');


	    #Calculate coverages properly (do not use alignment length as it includes gaps

	    my $qCov_tmp = ($qend - $qstart + 1) / $qLen;
	    my $qCov     = ($qCov_tmp > 1.0)? 1.0 : $qCov_tmp;

	    my $hCov_tmp = ($send - $sstart + 1) / $hLen;
	    my $hCov = ($hCov_tmp > 1.0)? 1.0 : $hCov_tmp;


#	    print Data::Dumper->Dump([$qLen, $qCov, $hLen, $hCov, $gs_coverage, $hEvalue, $gs_evalue, $hId, $gs_identity],
#				     [qw(*qLen *qCov $hLen *hCov *coverageCutoff *evalue *evalCutoff *hId *IDcutoff)]);
#	    <STDIN>;


	    #Before storing hit results check minimum coverage, identity and evalue
	    next HSP unless (($qCov >= $gs_coverage || $hCov >= $gs_coverage) &&
			     ($hEvalue <= $gs_evalue) && ($hId >= $gs_identity));


	    #hit identity
	    $tmp{hName}   = $hit->name;
	    $tmp{hLen}    = $hLen;


	    #hit statistics
	    $tmp{alnLen}  = $alnLen;
	    $tmp{hEvalue} = $hEvalue;
	    $tmp{hId}     = $hId;
	    $tmp{hSim}    = $hSim;
	    $tmp{qCov}    = $qCov;
	    $tmp{hCov}    = $hCov;


	    #The alignment
	    $tmp{qstart}  = $qstart;
	    $tmp{qend}    = $qend;
	    $tmp{sstart}  = $sstart;
	    $tmp{send}    = $send;

	    $tmp{qSeq}    = $hsp->query_string;
	    $tmp{sSeq}    = $hsp->hit_string;
	    $tmp{homStr}  = $hsp->homology_string;


	    #Get the GSAT score
	    my $gsat_outFile = "$alignmentsDir/${tcAcc}_" . $lh_hits{$bundle}{qName} . "_vs_" . $tmp{hName} . ".gsat";


#	    print "gsat.py $tmp{qSeq} $tmp{sSeq} $gsatShuffles > $gsat_outFile\n";
#	    exit;

	    system "gsat.py $tmp{qSeq} $tmp{sSeq} $gsatShuffles > $gsat_outFile" unless (-f $gsat_outFile);

	    my $gsat_score = TCDB::Assorted::get_gsat_score ($gsat_outFile);
	    $tmp{gsat} = $gsat_score;


#	    print Data::Dumper->Dump([\%tmp ], [qw(*matchData )]);
#	    <STDIN>;


	    #GSAT is the last filter
	    next HSP unless ($gsat_score >= $min_gsat_score);

	    #------------------------------------------------------------
	    #Generate quod plot with the repeat

	    my $whole_prot_seq = "$gs_seqDir/${seqId}.faa";
	    die "Protein sequence not found: $whole_prot_seq" unless (-f $whole_prot_seq);


	    my $plotFile  = "$hydroPlotsDir/${seqId}_" . $lh_hits{$bundle}{qName} . "_vs_" . $tmp{hName};
	    my $fileName  = "../$plotsDir/${seqId}_" . $lh_hits{$bundle}{qName} . "_vs_" . $tmp{hName} . ".png";
	    my $plotTitle = $lh_hits{$bundle}{qName} . " vs " . $tmp{hName};

	    #Get hydrophobic peaks coords
	    my $hydroPeaks = $gh_tms{$seqId};
	    die "No hydrophobic peaks found for sequence: $seqId" unless (@{ $hydroPeaks });


	    #format the hydrophobic peaks for quod
	    my @peaks = map { join ("-", @$_) . ":orange" } @$hydroPeaks;
	    my $pstring = join (" ", @peaks);


	    #----------
	    #Calculate the positions of the aligned section of each bundle in the full sequence.

	    my $q_bid = ($lh_hits{$bundle}{qName} =~ /BDL(\d+)/)? $1 : undef;
	    my $s_bid = ( $tmp{hName} =~ /BDL(\d+)/)? $1 : undef;
	    die "Could not extract bundle number for: $lh_hits{$bundle}{qName} or $tmp{hName}" unless ($q_bid && $s_bid);


	    #extract initial positions for both bundles
	    my $qbstart = $lhr_bundleSeqFiles->{$q_bid}->[1];
	    my $qbend   = $lhr_bundleSeqFiles->{$q_bid}->[2]; #$qLen - 1;
	    my $sbstart = $lhr_bundleSeqFiles->{$s_bid}->[1];
	    my $sbend   = $lhr_bundleSeqFiles->{$s_bid}->[2]; #$hLen - 1;
	    die "Could not extract coords for bundle $q_bid" unless ($qbstart && $qbend);
	    die "Could not extract coords for bundle $s_bid" unless ($sbstart && $sbend);


	    #Calculate bundle positions here
	    my $qgp_start = $qbstart + ($qstart - 1);
	    my $qgp_end   = $qbstart + ($qend - 1);

	    my $sgp_start = $sbstart + ($sstart - 1);
	    my $sgp_end   = $sbstart + ($send - 1);


	    #Format the coordinates for the repeats now
	    my $qrep = "${qgp_start}-${qgp_end}:green";
	    my $srep = "${sgp_start}-${sgp_end}:blue";

	    #Format the coordinates for the bar delimiting the bundles
	    my $bars = "-w ${qbstart}-${qbend}::1 ${sbstart}-${sbend}::1";

	    #The quod command line
	    my $cmd = "quod.py $whole_prot_seq  -t png -l '$plotTitle' -o $plotFile -q -r 80 $bars --xticks $xticksSpacing -nt +0 -at ${pstring} ${qrep} ${srep}";

	    my $img = "${plotFile}.png";
	    system $cmd unless (-f $img);
	    die "Could not generate plot: $img" unless (-f $img);

	    $tmp{plot} = $fileName;


	    #load the data into the hits section for this bundle
	    push (@{ $lh_hits{$bundle}{hits} }, \%tmp);


	  } #HSP
	} #HIT
      } #While


      #Add results to the topHits hash
      if (@{ $lh_hits{$bundle}{hits} }) {
	$lhr_topHits->{$id} = \%lh_hits;
      }

    }
}




#==========================================================================
#Given a sequence, its TMS coordinates and a repeat size (rsize), cut the
#sequence in TMS bundles of length rsize.


sub cut_seq_in_tms_regions {

  my ($ls_pid, $ls_repeat, $lhr_tms, $lhr_seqSegs) = @_;


  %$lhr_seqSegs = ();


  #Get the directory where bundle sequences will be saved
  my $sequencesDir = undef;

  if ($mode eq 'all') {
    $sequencesDir = getSequencesDir();
  }
  else {
    $sequencesDir = getSequencesDir($ls_pid);
  }
  die "Error: invalid sequence dir" unless ($sequencesDir);


  #----------------------------------------------------------------------
  #Get the coordinates of the overlapping bundles

  my @la_tms = @{ $lhr_tms->{$ls_pid} };



  #Get the Length of the sequence of the query protein
  my $seqFile = "$gs_seqDir/${ls_pid}.faa";
  my $obj  = Bio::SeqIO->new(-file => $seqFile , -format => "fasta");
  my $seqObj = $obj->next_seq;
  my $qlength = $seqObj->length;
  die "Could not extract protein length." unless ($qlength);

  #Store the length of the original sequence for proper calculation of
  #the x-ticks in the hydropathy plots of the results
  $origSeqLength{$ls_pid} = $qlength;



  #Number of TMS in protein
  my $ls_ntms = scalar (@la_tms);



  for (my $idx=1; $idx <= ($ls_ntms - $ls_repeat + 1); $idx++) {

    #TMS in bundle
    my $left_tms  = $la_tms[$idx - 1];
    my $right_tms = $la_tms[$idx + $ls_repeat - 2];


    #The coordinates of the bundle
    my $left_pos  = (($left_tms->[0]  - $gs_tail) <= 0)? 1 : $left_tms->[0]  -  $gs_tail;
    #my $right_pos = (($right_tms->[1] + $gs_tail) >= $qlength)? $right_tms->[1] : $right_tms->[1] + $gs_tail;
    my $right_pos = (($right_tms->[1] + $gs_tail) >= $qlength)? $qlength - 1 : $right_tms->[1] + $gs_tail;


    #Cut and name the bundles only if bundle file does not exist
    my $outfile = "${ls_pid}_bundle${idx}";
    unless (-f "$sequencesDir/${outfile}.faa") {

      #cutting bundle
      my $args = qq(-if $seqFile -od $sequencesDir -of $outfile -rangeCut -s $left_pos -e $right_pos -t 0);
      system "tmsplit $args > /dev/null";

      #replace protein ID with bundle number to the ID so alignments can be easily identified
      system qq(perl -i -pe 's/>\\S+/>BDL$idx/' $sequencesDir/${outfile}.faa);
    }

    $lhr_seqSegs->{$idx} = ["${outfile}.faa", $left_pos, $right_pos];
  }
}




#==========================================================================
#Read file with the TMS coordinates of the input proteins. The TMS
#must have been validated with WHAT to make sure they are reliable.


sub read_tms_coordinates_file {

  my ($s_coordsFile, $hr_tms) = @_;

  open (my $fileh, "<", $s_coordsFile) || die $!;

  #-----------------------------------------------------------------
  #The format of this file is protein ID followed by pairs of
  #coordinates separated by dash:
  #  2.A.43.1.1-O60931   1-20  25-35  50-68 ....
  if ($infileFmt eq 'tms') {

    while(<$fileh>) {
      chomp;

      #ignore empty lines;
      next unless ($_);

      #extract id and TMSs coordinates
      my ($id, @tms_str) = split(/\s+/, $_);
      my @tms = map { [ split(/-/, $_) ] }  @tms_str;


      #For debugging purposes
#      next unless ($id eq 'WP_100644534');


      $hr_tms->{$id} = \@tms;

      #Verify that the sequence is available for this protein
      unless (-f "$gs_seqDir/${id}.faa" && ! (-z "$gs_seqDir/${id}.faa")) {
	die "Could not find sequence for protein: $id in dir: $gs_seqDir -->";
      }
    } #while
  }

  #Input file is in HMMTOP format
  else {
    while(<$fileh>) {
      chomp;

      #Remove trailing spaces
      s/\s+$//;

      #ignore empty lines
      next unless ($_);


      #parse hmmtop line
      my ($id, $ntms, $tms_str) = (/\S+\s+\d+\s+(\S+).+(IN|OUT)\s+(\d+)\s+([\d\s-]+)/)? ($1, $3, $4) : ();

      #For debugging purposes
#      next unless ($id eq 'WP_100644534');


      if ($id && $ntms && $tms_str) {

	#extract the pairs of coordinates for TMS
	my @coords = split(/\s+/, $tms_str);
	my @tms = ();
	for (my $i=0; $i < $#coords; $i += 2) {
	  push (@tms, [$coords[$i], $coords[$i+1]]);
	}

	$hr_tms->{$id} = \@tms;

      }
      else {
	print "problem parsing HMMTOP line: $_\n";;
	print Data::Dumper->Dump([$id, $ntms, $tms_str ], [qw(*id *ntms *tms_str )]);
	exit;;
      }
    }
  }

  close $fileh;
}



#==========================================================================
#Get the directory where the sequences of bundles will be saved.

sub getSequencesDir {

  my $protId = shift;

  my $dir = undef;

  if ($mode eq 'all') {
    $dir = "$outdir/$seqDir";
  }
  else {
    die "Error: protein accession missing" unless ($protId);
    $dir = "$outdir/$protId/$seqDir";
  }

  system "mkdir -p $dir" unless (-d $dir);
  die "No dir for bundle sequences found: $dir" unless (-d $dir);

  return $dir;
}


#==========================================================================
#Get the directory where the alignments will be saved

sub getAlignmentsDir {

  my $protId = shift;

  my $dir = undef;

  if ($mode eq 'all') {
    $dir = "$outdir/$alignDir";
  }
  else {
    die "Error: protein accession missing" unless ($protId);
    $dir = "$outdir/$protId/$alignDir";
  }

  system "mkdir -p $dir" unless (-d $dir);
  die "No dir for alignments found: $dir" unless (-d $dir);

  return $dir;
}


#==========================================================================
#Get the directory where hydropathy plots will be saved

sub getPlotsDir {

  my $protId = shift;

  my $dir = undef;

  if ($mode eq 'all') {
    $dir = "$outdir/$plotsDir";
  }
  else {
    die "Error: protein accession missing" unless ($protId);
    $dir = "$outdir/$protId/$plotsDir";
  }

  system "mkdir -p $dir" unless (-d $dir);
  die "No dir for plots found: $dir" unless (-d $dir);

  return $dir;
}


#==========================================================================
#Get the directory where the reports will be saved

sub getReportsDir {

  my $protId = shift;

  my $dir = undef;

  if ($mode eq 'all') {
    $dir = "$outdir/$repDir";
  }
  else {
    die "Error: protein accession missing" unless ($protId);
    $dir = "$outdir/$protId/$repDir";
  }

  system "mkdir -p $dir" unless (-d $dir);
  die "No dir for reports found: $dir" unless (-d $dir);

  return $dir;
}





#==========================================================================
#Read command-line arguments

sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $ls_status = GetOptions(
      "i|infile=s"         => \$gs_infile,
      "if|infile-format=s" => \$infileFmt,
      "o|outdir=s"         => \$outdir,
      "f|id-format=s"      => \$gs_idFormat,
      "r|rep-unit=i"       => \$gs_repUnit,
      "t|tail-size=i"      => \$gs_tail,
      "s|seqs=s"           => \$gs_seqDir,
      "e|evalue=f"         => \$gs_evalue,
      "c|coverage=f"       => \$gs_coverage,
      "id|identity=f"      => \$gs_identity,
      "ncs|no-comp-stats!" => \$compStatsFlag,
      "gs|gsat-shuffles=i" => \$gsatShuffles,
      "z|gsat-cutoff=f"    => \$min_gsat_score,
      "m|mode=s"           => \$mode,
      "h|help"        => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($ls_status);

  #----------------------------------------------------------------------
  #Validate command line arguments

  die "Error: argument -i is mandatory.\n" unless ($gs_infile);
  die "Error: argument -r is mandatory and must be greater than 0.0\n" unless ($gs_repUnit > 0);
  die "Error: augument -t must be grater than 0 and less than 16\n" if ($gs_tail > 15 || $gs_tail < 0);
  die "Error: argument -e must be greater than 0\n" unless ($gs_evalue >=0 );
  die "Error: argument -c must be between 0.5 and 1.0\n" unless ($gs_coverage >= 0.0 && $gs_coverage <= 1.0);
  die "Error: argument -id must be between 0.25 and 1.0\n" unless ($gs_identity >= 0.0 && $gs_identity <= 1.0);

  #Option -f
  $gs_idFormat = lc $gs_idFormat;
  unless ($gs_idFormat =~ /^(tc|tca|o)$/) {
    die "Error: There are 3 Valid options for -f (tc, tca, o)\n";
  }


  #option -if
  $infileFmt = lc $infileFmt;
  unless ($infileFmt =~ /^(hmmtop|tms)$/) {
    die "Error: invalid input file format: '$infileFmt' (Valid options: hmmtop, tms).\n";
  }


  #option -m
  $mode = lc $mode;
  unless ($mode =~ /^(all|each|debug)$/) {
    die "Error: invalid mode of operation '$mode'. Valid options are: all, each!\n";
  }


  #Option -s
  unless (-d $gs_seqDir) {
    die "Error: Directory with sequences must exits -> $gs_seqDir\n";
  }


  #Validate GSAT cutoff
  unless ($min_gsat_score >= 0) {
    die "Use GSAT cutoff >= 3.0!\n";
  }


  #option -ncs
  $compStats  = ($compStatsFlag)? "" : "-k 1000 -z 11";
}



sub print_help {

  my $help = <<'HELP';

This script searches for regions of TMSs repeated in a full protein.

-i, --infile {path}
   Input file with id/accession(s) of the protein(s) to analyze and the coordinates
   of the TMSs in that protein(s). Use option -if to specify the format of this
   file.
   (Argument is mandatory).

-if, --infile-format {string} (optional)
  Format of the TMS coordenates. It can be either 'tms' or 'hmmtop'.
  (Default: hmmtop)

-o, --outdir {path}
  Output directory where results will be saved.
  (Default: repeats)

-s, --seqs {path}
  Directory to access the sequences in FASTA format that will be used to 
  search for repeats. One file per sequence, and the name of the file is
  the accession of the protein followed by '.faa'
  (Argument is mandatory)

-f, --id-format {string}
   Acceptable formats for identifiers:
   tc    plain tcdb identifier of a system (e.g., 2.A.1.8.1)
   tca   tcdb id and accession separated by dash (e.g. 2.A.1.8.3-Q9R6U5)
   o     other, it can be refSeq, uniprot or custom, but it is requried
         that is is a single string without spaces.
   (Argument is mandatory)

-r, --repeat-unit {int)
   Size in TMS of the repeat unit to search in the protein.
   (Argument is mandatory)

-t, --tail-size {int}
   Number of residues to add to the beginning and end of TMS regions before
   running comparisons. Value should be less than or equal to 15 residues.
   (Default: 5);

-e, --evalue {float}
   Maximum evalue to consider an alignment between two TMS bundles significant.
   (Default: 0.1);

-ncs, --no-comp-stats {FLAG}
   If present, this flag indicates that  E-values will not be corrected using
   compositional statistics.
   (Default: apply correction).

-c, --coverage {float}
   Minimum alignment coverage of the smallest bundle to consider an alignment
   signifiant.
   (Default: 0.8)

-id, --identity {float}
   Minimum identity, expressed as a float in the 0-1 range, to consider an
   alignment signficant.
   (Defatul: 0.25);

-gs, --gsat-shuffles {int}
  Number of shuffles that will be used to run GSAT on good matches.
  (Default: 1000);

-z, --gsat-cutoff {int}
  Minimum GSAT score cutoff to select good hits. 
  (Default: 4.0)

-h, --help
  Print this help message. It takes precedence to any other option.

HELP

  print $help;
  exit;

}
