#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

$Data::Dumper::Deepcopy = 1;

use Getopt::Long;
use LWP;
use Bio::SeqIO;
use Bio::SearchIO;

#Local libraries
use TCDB::CheckDependencies;
use TCDB::Domain::PfamParser;
use TCDB::Domain::Characterize;
use TCDB::Assorted;


###########################################################################
#
# Comapre two files with fasta sequences and report the alignment parameters
# Along with hydropathy plots and PFAM domains.
#
###########################################################################

#==========================================================================
#Check dependencies

my @dependencies = ("zgrep", "blastp", "ssearch36", "hmmtop", "blastdbcmd",
		    "hmmscan");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;


#This will prevent quod and alnquod from going into interactive mode
$ENV{"MPLBACKEND"} = "agg";


#==========================================================================
#Read command line arguments

my $qfile      = "";
my $qProt      = "";
my $sfile      = "";
my $sProt      = "";
my $qlabel     = "Query";
my $slabel     = "Subject";
my $outdir     = "";
my $prog       = 'ssearch36'; #'blastp';
my $evalue     = 1e-4;
my $identity   = 20.0;
my $coverage   = 40.0;
my $covControl = "X";
my $blastComp  = "F"; #2;
my $segFilter  = 'no';
my $minLength  = 30;  #Min legnth of proteins to analyze (without gaps)
my $subMatrix  = 'BL50';

#this can be used to remove long sequences from results
my $maxProtLength = 100000; #default threshold to allow any length
my $LengthControl  = "N";    #same meaning as $covControl

#internal directories
my $filesDir = "";
my $plotsDir = "";
my $seqDir   = "";
my $blastDir = "";


read_command_line();

#print Data::Dumper->Dump([$qfile, $qProt, $sfile, $sProt, $qlabel, $slabel, $outdir, $prog,
#                          $evalue, $coverage, $covControl, $blastComp, $segFilter, $maxProtLength,
#			  $LengthControl],
#                        [qw(*qfile *qProt *sfile *sProt *qlabel *slabel *outdir *prog
#                         *evalue *coverage *covControl *blastComp *segFilter *maxProtLength
#                         *LengthControl)]);
#exit;



#==========================================================================
#Output files

#The alignment file by blastp or ssearch36
my $alnFile = "$filesDir/${prog}.out";

#The results of running hmmscan
my $pfamFile = "$filesDir/hmmscan.out";

#The results from running hmmtop
my $hmmtopFile = "$filesDir/hmmtop.out";

#The blast database to retrieve sequences for ploting
my $blastdb  = "$blastDir/sequences";



#==========================================================================
#Run the alignment first

print "Running $prog and parsing output....\n";
run_alignment();


my @alnHits = ();
if ($prog eq 'blastp') { parse_blast(\@alnHits); }
elsif ($prog eq 'ssearch36') { parse_ssearch(\@alnHits)}

#print Data::Dumper->Dump([\@alnHits ], [qw($alnHits )]);
#exit;


die "No significant blastHits found!\n" unless (@alnHits);



#==========================================================================
#Run pfam (get clans, hmmtop,  and parse results

my %pfamHits   = ();
my %clans      = ();
my %hmmtopHits = ();
run_pfam_hmmtop(\%pfamHits, \%hmmtopHits,\%clans );

#print Data::Dumper->Dump([\%clans], [qw(*clans)]);
#exit;


#==========================================================================
#Parse the alignment results to make sure there are signficant results,
#get domains for significant hits and plot the corresponding hydropathies.




print "Geneating report...\n";
generate_report();



#==========================================================================
################   Subroutines definition beyond ths point   ##############
#==========================================================================


#==========================================================================
#Generate output for significant hits

sub generate_report {


  #Prepare output files
  my $htmlFile  = "$outdir/report.html";
  my $plotsFile = "$outdir/plots.html";


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

.pfam {
   text-align: center;
   vertical-align: middle;
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


.dom {
   border: 2px solid black;
   height: 100px;
   width:  100%;
   overflow-x: auto;
   overflow-y: auto;
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
    <title>$qlabel vs $slabel</title>
  </head>
  <br />
  <h1 style='text-align:center'>$qlabel vs $slabel</h1>
  <body>

HEADER


  open  (my $outh, ">", $htmlFile) || die $!;
  print $outh $htmlHeader;


  foreach my $hit (sort by_evalue @alnHits) {


    my $qacc   = $hit->{qacc};
    my $qlen   = $hit->{qlen};
    my $qseq   = $hit->{qseq};
    my $qcov   = sprintf("%.1f", $hit->{qcov});
    my $qstart = $hit->{qstart};
    my $qend   = $hit->{qend};

    my $sacc   = $hit->{sacc};
    my $slen   = $hit->{slen};
    my $sseq   = $hit->{sseq};
    my $scov   = sprintf("%.1f", $hit->{scov});
    my $sstart = $hit->{sstart};
    my $send   = $hit->{send};

    my $eval   = sprintf ("%.1e", $hit->{evalue});
    my $id     = sprintf ("%.1f", $hit->{id});
    my $hstr   = $hit->{hstr};

    my $alnHit = <<HIT;

    <br /><hr style=\"border-style:solid; border-width:5px; color:black;\"/>

    <p><b>$qacc ($qlen) vs $sacc ($slen)</b></p>

    <table width="600px" border="0" cellspacing="0" cellpadding="2">
      <tr>
         <td class='label'><b>E-value:</b></td>
         <td class='data'>$eval</td>
         <td class='label'><b>Identity:</b></td>
         <td class='data'>${id}%</td>
         <td class='label'><b>Q_coverage:</b></td>
         <td class='data'>${qcov}%</td>
         <td class='label'><b>S_coverage:</b></td>
         <td class='data'>${scov}%</td>
      </tr>
      <tr>
         <td class='label'><b>Q_aln:</b></td>
         <td class='data'>${qstart}-$qend</td>
         <td class='label'><b>S_aln:</b></td>
         <td class='data'>${sstart}-$send</td>
         <td class='label'></td>
         <td class='data'></td>
         <td class='label'></td>
         <td class='data'></td>
      </tr>
    </table>
    <br />

    <p><b>Alignment:</b></p>
    <div class='seq'>
    <pre>
$qseq
$hstr
$sseq
    </pre>
    </div>

HIT

    print $outh $alnHit;


    #Generate the hydropathy plots
    my $good = run_quod($qacc, $sacc, $qstart, $qend, $sstart, $send, $qseq, $sseq);
    die "Could not generate plots for hit: $qacc vs $sacc" unless ($good);

    my $domData = generate_domain_data($qacc, $sacc);
    my $domHTML = "";
    if ($domData) {
      $domHTML =<<DATA;
    <br />
    <hr />
    <p><b>Pfam info:</b></p>
    <div class='dom'>

$domData

    </div>
DATA
    }

    my $plot_aln = "plots/${qacc}_vs_${sacc}_qs${qstart}_qe${qend}_ss${sstart}_se${send}.png";
    my $qplot    = "plots/${qacc}_vs_${sacc}_qaln_qs${qstart}_qe${qend}.png";
    my $splot    = "plots/${qacc}_vs_${sacc}_saln_ss${sstart}_se${send}.png";

    #now include the plots
    my $prtPlots =<<PLOTS;

    <br />
    <table style="width:100%">
      <tr>
        <td><a href="$qplot" target="_blank"><img src="$qplot" alt="$qacc"></a></td>
        <td><a href="$splot" target="_blank"><img src="$splot" alt="$sacc"></a></td>
      </tr>
      <tr>
        <td colspan="2" style="text-align: center;">
           <a href="$plot_aln" target="_blank"><img src="$plot_aln" alt="$qacc vs $sacc alignment"></a>
        </td>
      </tr>
    </table>

$domHTML

PLOTS

      print $outh $prtPlots;
  }

  #Close HTML report
  my $closeRep = <<CLOSE;
  </body>
</html>
CLOSE

  print $outh $closeRep;

  close $outh;
}


#==========================================================================
#Generate domain data for the html report


sub generate_domain_data {

  my ($q, $s) = @_;

  #Format of PFAM hash:
  #      push (@{ $out->{$qacc}->{$pfamID} },
  #	    {pfamid=> $pfamID, dlen=>$pfamLen,  dstart=>$dstart,  dend=>$dend, evalue=>$eval,
  #	      dname=>$pfamName, def=>$def, qlen=>$qlen, qstart=>$qstart, qend=>$qend });




  my @cols   = qw(Query Domain Clan Dom_length E-value Dom_start Dom_end Q_Start Q_end Dom_Name Dom_Info);
  my $colStr = "      <th>" . join ("</th>\n      <th>", @cols) .
    "</th>\n";

  my $header =<<HEADER;
    <table border='1', style='width:100%'>
    <tr>
$colStr
    </tr>
HEADER


  my $res = "";
  foreach my $prot ($q, $s) {

    if  (exists $pfamHits{$prot}) {
      my @Doms = keys %{ $pfamHits{$prot} };

      foreach my $d (@Doms) {

	my $clan = ($clans{$d})? $clans{$d} : "N/A";

	my @hits = @{ $pfamHits{$prot}{$d} };
	foreach my $hit (@hits) {
	  my $dlen   = $hit->{dlen};
	  my $eval   = $hit->{evalue};
	  my $qstart = $hit->{qstart};
	  my $qend   = $hit->{qend};
	  my $dstart = $hit->{dstart};
	  my $dend   = $hit->{dend};
	  my $name   = $hit->{dname};
	  my $def    = $hit->{def};

	  $res .=<<ROW;
    <tr>
       <td class="pfam">$prot</td>
       <td class="pfam">$d</td>
       <td class="pfam">$clan</td>
       <td class="pfam">$dlen</td>
       <td class="pfam">$eval</td>
       <td class="pfam">$dstart</td>
       <td class="pfam">$dend</td>
       <td class="pfam">$qstart</td>
       <td class="pfam">$qend</td>
       <td class="pfam">$name</td>
       <td class="pfam">$def</td>
    </tr>
ROW
	}
      }
    }
  }

  #Return final result
  if ($res) {
    $header .= $res;
    $header .= "    </table>\n";
    return $header;
  }
  else {
    return $res;
  }
}




#==========================================================================
#Run quod on the query, subject and the alignment.


sub run_quod {

  my ($q, $s, $qs, $qe, $ss, $se, $qseq, $sseq) = @_;


  #extract sequences for query and subject
  extract_full_sequences($q,$s);


  #-----------------------------------------------------------------
  #Run quod for the alignment

  #First save aligned segments to files
  my $qalnFile ="$seqDir/${q}_aln.faa";
  open(my $qfh, '>', $qalnFile) || die $!;
  print $qfh ">$q alignment\n$qseq\n";
  close $qfh;

  my $salnFile ="$seqDir/${s}_aln.faa";
  open(my $sfh, '>', $salnFile) || die $!;
  print $sfh ">$s alignment\n$sseq\n";
  close $sfh;


  #Note alnquod requires to add the extension to the image name
  my $alnFig = "$plotsDir/${q}_vs_${s}_qs${qs}_qe${qe}_ss${ss}_se${se}.png";
  my $cmd1 = qq(alnquod.py --grid -q -l "$q (red) and $s (blue)"  -o $alnFig --xticks 25 --width 15 -- $qalnFile  $seqDir/${q}.faa $salnFile  $seqDir/${s}.faa);
  #print "$cmd1\n\n";
  system $cmd1 unless (-f "${alnFig}");
  return undef unless (-f "${alnFig}");


  #-----------------------------------------------------------------
  #Run quod for the full sequencess of the query and subject proteins


  #Extract TMS coordinates for query
  die "Error: no hmmtop results for: $q" unless (exists $hmmtopHits{$q});
  my $qTMS = "";
  if (scalar @{ $hmmtopHits{$q}{coords} } > 0) {
    $qTMS  = "-at " . join(",", @{ $hmmtopHits{$q}{coords} }) . ":orange";
  }


  #Plot query hydropathy
  my $qPfam = get_pfam_coords_for_quod($q, "red");
  my $qName = "$plotsDir/${q}_vs_${s}_qaln_qs${qs}_qe${qe}";
  my $cmd2  = qq(quod.py --grid -q -l "$q"  -o $qName --width 15 --color red --xticks 25 -w ${qs}-${qe}::1 -t png -nt +0  $qTMS $qPfam --  $seqDir/${q}.faa);
  #print "$cmd2\n\n";
  system $cmd2 unless (-f "${qName}.png");
  return undef unless (-f "${qName}.png");



  #TMS coords for the subject
  die "Error: no hmmtop results for: $s" unless (exists $hmmtopHits{$s});
  my $sTMS = "";
  if (scalar @{ $hmmtopHits{$s}{coords} } > 0) {
    $sTMS  = "-at " . join(",", @{ $hmmtopHits{$s}{coords} }) . ":cyan";
  }

  #Plot Subject hydropaty
  my $sPfam = get_pfam_coords_for_quod($s, "blue");
  my $sName = "$plotsDir/${q}_vs_${s}_saln_ss${ss}_se${se}";
  my $cmd3  = qq(quod.py --grid -q -l "$s"  -o $sName --width 15 --color blue --xticks 25 -w ${ss}-${se}::1 -t png -nt +0 $sTMS $sPfam --  $seqDir/${s}.faa);
  #print "$cmd3\n\n";
  system $cmd3 unless (-f "${sName}.png");
  return undef unless (-f "${sName}.png");


  return 1;
}

#==========================================================================
#Get the string for quod that will plot the PFAM domains

sub get_pfam_coords_for_quod {

  my ($prot, $color) = @_;

  #Format of PFAM hash:
  #      push (@{ $out->{$qacc}->{$pfamID} },
  #	    {pfamid=> $pfamID, dlen=>$pfamLen,  dstart=>$dstart,  dend=>$dend,
  #	     def=>$def, qlen=>$qlen, qstart=>$qstart, qend=>$qend });



  my $str = "";

  if  (exists $pfamHits{$prot}) {
    my @Doms = keys %{ $pfamHits{$prot} };
    my $dcnt = 0;
    $str = "--region-font 12 -ar ";
    foreach my $d (@Doms) {

      my @hits = @{ $pfamHits{$prot}{$d} };
      foreach my $hit (@hits) {
	my $left  = $hit->{qstart};
	my $right = $hit->{qend};

	my $ypos = -2.8 + $dcnt * 0.4;
	$str .= "${left}-${right}:'${d}':${ypos}:$color ";
	$dcnt++;
      }
    }
  }

  return $str;

}


#==========================================================================
#Extract the full sequences of the query and subject proteins
#Examples:  AKM80767.1



sub extract_full_sequences {

  my ($q, $s) = @_;


  my $q_seq = "$seqDir/${q}.faa";
  my $s_seq = "$seqDir/${s}.faa";

  #extract the query secuence from tcdb and the subject from the custom blastdb
  my $cmd1 = qq(blastdbcmd  -db $blastdb  -entry $q -target_only -out $q_seq);
  system "$cmd1" unless (-f $q_seq && !(-z $q_seq));
  die "Could not extract sequence for $q" unless (-f $q_seq && !(-z $q_seq));


  my $cmd2 = qq(blastdbcmd  -db $blastdb  -entry $s -target_only -out $s_seq);
  system "$cmd2" unless (-f $s_seq && !(-z $s_seq));
  die "Could not extract sequence for $s" unless (-f $s_seq && !(-z $s_seq));
}







#==========================================================================
#Sort alignmnet results by E-value

sub by_evalue {
  $a->{evalue} <=> $b->{evalue};
}



#==========================================================================
#Run PFAM, hmmtop and parse results


sub run_pfam_hmmtop {
  my ($pfamOut, $hmmtopOut, $pfamClans) = @_;


  #----------------------------------------------------------------------
  #Generate blast DB for easy sequence retrieval

  print "Generate Blast DB with sequences for fast sequence retrieval...\n";

  #Get the sequences for which hmmscan will run
  my $allSeqsFile = "$seqDir/all_seqs.faa";
  system qq(cat $qfile $sfile > $allSeqsFile) unless (-f $allSeqsFile && !(-z $allSeqsFile));
  die "Could not generate file: $allSeqsFile" unless (-f $allSeqsFile && !(-z $allSeqsFile));


  #Generate blastdb ...assuming there are no duplicate sequences.
  my $cmd1 = qq(makeblastdb -dbtype prot -in $allSeqsFile -title '$qlabel plus $slabel' -parse_seqids -hash_index -out $blastdb);
  print "$cmd1\n";
  system $cmd1 unless (-f "${blastdb}.pin");
  system "rm $allSeqsFile" if (-f $allSeqsFile);


  #----------------------------------------------------------------------
  #Get the accessions of the top hits in the alignments

  #get the accessions with significant hits
  my %accList = ();
  foreach my $hit (@alnHits) {
    $accList{$hit->{qacc}} = 1;
    $accList{$hit->{sacc}} = 1;
  }


  #Save accessions to a file
  my $idFile = "$seqDir/top_hits_accs.txt";
  unless (-f $idFile) {
    open (my $afh, ">", $idFile) || die $!;
    print $afh join("\n", keys %accList), "\n";
    close $afh;
  }

  #----------------------------------------------------------------------
  #Extract full sequences for top hits.

  my $topHitsSeqs = "$seqDir/topHits.faa";
  my $cmdTopHits = qq(blastdbcmd -db $blastdb -entry_batch $idFile -target_only -out $topHitsSeqs);
  system $cmdTopHits unless (-f $topHitsSeqs && !(-z $topHitsSeqs));


  #----------------------------------------------------------------------
  #run hmmscan on all the sequences for both files

  print "\nRunning hmmscan and parsing output....\n";

  my $pfamDB = ($ENV{PFAMDB})? $ENV{PFAMDB} : "$ENV{RESEARCH_DATA}/pfam/pfamdb/Pfam-A.hmm";
  my $cmd2 = qq(hmmscan --cpu 4 --noali --cut_ga -o /dev/null --domtblout $pfamFile $pfamDB  $topHitsSeqs);
  system $cmd2 unless (-f $pfamFile && !(-z $pfamFile));


  #parse Pfam output 
  TCDB::Assorted::parse_pfam($pfamFile, $pfamOut, $pfamClans);
#  print Data::Dumper->Dump([$pfamOut, $pfamClans ], [qw(*pfamOut *pfamClans )]);
#  exit;

  #----------------------------------------------------------------------
  #Extract clans

  TCDB::Assorted::get_clans($pfamClans, $filesDir);
#  print Data::Dumper->Dump([$pfamClans ], [qw(*clans )]);
#  exit;

  #--------------------------------------------------------------------------
  #Run hmmtop on top hits for later hydropathy plots.

  print "Runnign HMMTOP and parsing output...\n";

  my $cmd3 = qq(hmmtop -if=$topHitsSeqs -of=$hmmtopFile -sf=FAS -pi=spred -is=pseudo);
  system $cmd3 unless (-f $hmmtopFile);
  system "rm $topHitsSeqs" if (-f $topHitsSeqs);

  #Parse hmmtop output
  TCDB::Assorted::parse_hmmtop($hmmtopOut, $hmmtopFile);

}


#==========================================================================
#Parse ssearch36 output

sub parse_ssearch {

  my $out = shift;

  my $parser = new Bio::SearchIO (-format => 'fasta', -file => $alnFile);

  my $formatTmp = $parser->format();
#  print Data::Dumper->Dump([$formatTmp ], [qw(*fileFormat )]);
#  exit;

  while (my $result = $parser->next_result) {

    my $qacc = $result->query_name;
    my $qlen = $result->query_length;


  HIT:while (my $hit = $result->next_hit) {
    HSP:while(my $hsp = $hit->next_hsp) {

	#Alignment parameters
	my $sacc = $hit->name;
	my $slen = $hit->length;
	my $eval = $hsp->evalue;
	my $id   = $hsp->frac_identical('total') * 100;

	#coordinates and sequence
	my $qstart  = $hsp->start('query');
	my $qend    = $hsp->end('query');
	my $sstart  = $hsp->start('subject');
	my $send    = $hsp->end('subject');
	my $qseq    = $hsp->query_string;
	my $sseq    = $hsp->hit_string;
	my $hstr    = $hsp->homology_string;


	#Check first that both proteins have the right length
	next HSP if (max_length_violation($qlen, $slen, $maxProtLength, $LengthControl));

	#If the alignment has less than $minLength aas, ignore it
        my $qtmp = $qseq; $qtmp =~ s/-//g;
	my $stmp = $sseq; $stmp =~ s/-//g;
	next HSP if (length($qtmp) < $minLength || length($stmp) < $minLength);

	#Calculate coverages properly (do not use alignment length as it includes gaps
	my $qCov_tmp = ($qend - $qstart + 1) / $qlen * 100;
	my $qcov     = ($qCov_tmp > 100.0)? 100 : $qCov_tmp;

	my $sCov_tmp = ($send - $sstart + 1) / $slen * 100;
	my $scov = ($sCov_tmp > 100.0)? 100 : $sCov_tmp;


	if ($eval <= $evalue && TCDB::Assorted::coverage_ok($qcov, $scov, $coverage, $covControl)) {

	  push(@{ $out }, {qacc=>$qacc, sacc=>$sacc,     qlen=>$qlen, slen=>$slen,     qcov=>$qcov,
			   scov=>$scov, evalue=>$eval,   id=>$id,     qstart=>$qstart, qend=>$qend,
			   sstart=>$sstart, send=>$send, qseq=>$qseq, sseq=>$sseq, hstr=>$hstr});
	}
      } # hsp
    } # hit
  } # result
}



#==========================================================================
#Test whether the lengths of two proteins are withing a predefined
#legnth specified by the user. This are the options for control:
# X:  Either protein is larger than the cutoff
# B:  Both proteins are larger than the cutoff
# Q:  Only the query protein is larger than the cutoff
# S:  Only the subject protein is larger than the cutoff
# N:  No control. Any length is ok.

sub max_length_violation {

  my ($qlen, $slen, $maxLen, $control) = @_;

  if ($control eq "X") {
    (($qlen >= $maxLen) || ($slen >= $maxLen))? return 1 : return 0;
  }

  if ($control eq "B") {
    (($qlen >= $maxLen) && ($slen >= $maxLen))? return 1 : return 0;
  }

  if ($control eq "Q") {
    ($qlen >= $maxLen)? return 1 : return 0;
  }

  if ($control eq "S") {
    ($slen >= $maxLen)? return 1 : return 0;
  }

  if ($control eq "N") {
    return 0;
  }

  die "Unknown control mode: $control";

}




#==========================================================================
#Parse blast output


sub parse_blast {
  my $out = shift;

  open (my $fh, "<", $alnFile) || die $!;
 LINE:while (<$fh>) {
    chomp;
    next unless ($_);
    next if (/^#/);

    #Blast columns: qacc sacc qlen slen evalue pident qstart qend sstart send qseq sseq
    my ($qacc, $sacc, $qlen, $slen, $eval, $id, $qstart, $qend, $sstart, $send, $qseq, $sseq) = split (/\t/, $_);


    if ($eval <= $evalue) {

      my $qcov = ($qend - $qstart + 1) / $qlen * 100;
      my $scov = ($send - $sstart + 1) / $slen * 100;

      if (TCDB::Assorted::coverage_ok($qcov, $scov, $coverage, $covControl)) {

	push(@{ $out }, {qacc=>$qacc, sacc=>$sacc,     qlen=>$qlen, slen=>$slen,     qcov=>$qcov,
			 scov=>$scov, evalue=>$eval,   id=>$id,     qstart=>$qstart, qend=>$qend,
			 sstart=>$sstart, send=>$send, qseq=>$qseq, sseq=>$sseq});
      }
    }
  }
  close $fh;
}





#==========================================================================
#Run the alignemnt between the two files depending on the program
#Selected by the user.

sub run_alignment {

  my $cmd = "";


  if ($prog eq 'blastp') {

    my $compStr = "-comp_based_stats $blastComp";
    my $segStr  = "-seg $segFilter";
    my $outFmt  = qq(-outfmt '7 qacc sacc qlen slen evalue pident qstart qend sstart send qseq sseq');

    #Run blast
    $cmd = qq(blastp -query $qfile -subject $sfile -matrix BLOSUM62 -out $alnFile $outFmt -evalue $evalue -use_sw_tback $compStr $segStr);
    print "$cmd\n";
    system $cmd unless (-f $alnFile && !(-z $alnFile));

    #Append command line to the end of results file
    open (my $fh, ">>", $alnFile) || die $!;
    print $fh "\n# $cmd\n";
    close $fh;
  }
  elsif ($prog eq 'ssearch36') {
    $cmd = qq(ssearch36 -z 11 -k 1000 -s $subMatrix -E $evalue -W 0 -m 0  $qfile $sfile > $alnFile );
    print "$cmd\n";
    system $cmd unless (-f $alnFile && !(-z $alnFile));
  }
  elsif ($prog eq 'glsearch36' || $prog eq 'ggsearch36') {
    $cmd = qq($prog -z 11 -k 1000 -s $subMatrix -E $evalue -m 0  $qfile $sfile > $alnFile );
    print "$cmd\n";
    system $cmd unless (-f $alnFile && !(-z $alnFile));
  }
}




#===========================================================================
#Read command line and print help


sub read_command_line {

  print_help() unless (@ARGV);

  my $status = GetOptions(
	"q|qfile=s"        => \&read_qfile,
	"s|sfile=s"        => \&read_sfile,
	"ql|qlabel=s"      => \&read_qlabel,
        "sl|slabel=s"      => \&read_slabel,
        "p|prog=s"         => \&read_prog,
	"o|outdir=s"       => \$outdir,
	"e|evalue=f"       => \$evalue,
	"c|coverage=f"     => \$coverage,
	"cc|cov-control=s" => \&read_covControl,
	"l|max-len=i"      => \$maxProtLength,
	"lc|len-control=s" => \&read_lenControl,
	"m|sub-matrix=s"   => \&read_subMatrix,
	"scs|seq-comp-stats=s"     => \&read_blastComp,
	"lcf|low-complex-filter=s" => \&read_segFilter,
	"h|help"           => sub { print_help(); },
	"<>"               => sub { die "Error: Unknown argument: $_[0]\n"; });
  exit unless ($status);


  #The output directories
  $outdir = "${prog}_${qlabel}_vs_${slabel}" unless ($outdir);
  system "mkdir -p $outdir" unless (-d $outdir);

  $filesDir = "$outdir/files";
  system "mkdir $filesDir" unless (-d $filesDir);

  $plotsDir = "$outdir/plots";
  system "mkdir $plotsDir" unless (-d $plotsDir);

  $seqDir = "$outdir/seqs";
  system "mkdir $seqDir" unless (-d $seqDir);

  $blastDir = "$outdir/blastdb";
  system "mkdir $blastDir" unless (-d $blastDir);


  #Get sequences if an accession was given by the user.
  if ($qProt) {
    $qfile = "$seqDir/${qProt}.faa";
    get_sequence ('nr', $qProt, $qfile);
  }

  if ($sProt) {
    $sfile = "$seqDir/${sProt}.faa";
    get_sequence ('nr', $sProt, $sfile);
  }

  #Check for incompatibilities and errors
  die "Error: options -q and -s are mandatory!\n" unless ($sfile && $qfile);


}


#==========================================================================
#Extract a sequence from blastdb


sub get_sequence {
  my ($db, $acc, $outfile) = @_;

  my $cmd = qq(blastdbcmd -db $db -entry $acc -target_only > $outfile);
  system $cmd;

  unless (-f $outfile && !(-z $outfile)) {
    die "Could not extract protein ($acc) from blast DB ($db)";
  }
}





#==========================================================================
#Option -q (Query can be a file or a RefSeq accession)

sub read_qfile {
  my ($opt, $value) = @_;

  if (-f $value) {

    unless (!(-z $value)) {
      die "Error in option -$opt: Query file is empty --> $value\n";
    }
    $qfile = $value;
    return;
  }

  #Not a file, assume protein accession
  else {
    $qProt = $value;
  }
}


#==========================================================================
#Option -f2

sub read_sfile {
  my ($opt, $value) = @_;

    if (-f $value) {

    unless (!(-z $value)) {
      die "Error in option -$opt: Subject file is empty --> $value\n";
    }
    $sfile = $value;
    return;
  }

  #Not a file, assume protein accession
  else {
    $sProt = $value;
  }
}


#==========================================================================
#Check if a user provided label chas legal characters

sub check_label {

  my($option, $label) = @_;

  unless ($label =~ /^[\+0-9a-zA-Z_\.\+-]+$/) {
    die "Error in option ${option}: illegal characters in label. Valid characters are 0-9a-zA-Z_.+-: $label\n";
  }
}


#==========================================================================
#Option -ql

sub read_qlabel {
  my ($opt, $value) = @_;

  check_label ("-ql", $value);

  $qlabel = $value;
}


#==========================================================================
#Option -sl

sub read_slabel {
  my ($opt, $value) = @_;

  check_label("-sl", $value);

  $slabel = $value;
}


#==========================================================================
#Option -p

sub read_prog {
  my ($opt, $value) = @_;

  my $tmp = lc $value;
  unless ($tmp =~ /^(ssearch36|blastp)$/) {
    die "Error in option -$opt: illegal program ($value). Valid programs are blastp and ssearch36\n";
  }

  $prog = $tmp;
}



#==========================================================================
#Option -cc

sub read_covControl {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^[XQSB]$/) {
    die "Error in option -$opt: illegal charater ($value). Valid characters are Q,S,B,X\n";
  }

  $covControl = $tmp;
}


#==========================================================================
#Option -lc

sub read_lenControl {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^[XQSBN]$/) {
    die "Error in option -$opt: illegal charater ($value). Valid characters are Q,S,B,X,N\n";
  }

  $LengthControl = $tmp;
}






#==========================================================================
#Option -m (Any matrix supported by ssearch

sub read_subMatrix {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^(BL50|BL62|P250|OPT5|VT200|VT160|P120|VT120|BL80|VT80|MD40|VT40|MD20|VT20|MD10|VT10)$/) {
    die "Error in option -$opt: illegal matrix ($value). Value should be any matrix supported by SSEARCH\n";
  }

  $subMatrix = $tmp;
}



#==========================================================================
#Option -scs

sub read_blastComp {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^[TF]$/) {
    die "Error in option -$opt: illegal charater ($value). Valid characters are T and F\n";
  }

  $blastComp = ($tmp eq 'T')? 2 : 0;
}


#==========================================================================
#Option -lcf

sub read_segFilter {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^[TF]$/) {
    die "Error in option -$opt: illegal charater ($value). Valid characters are T and F\n";
  }


  $segFilter = ($tmp eq 'T')? 'yes' : 'no';
}







#==========================================================================
#option -h


sub print_help {

    my $help = <<'HELP';

Compare two files with sequences in FASTA format. Sequences can be compared
using several programs, and hydropathy plots are generated to show the
alignments and domain content for each sequence.

Options:

-q, --qfile {file} (Mandatory)
   Query sequence file (in FASTA format) or NCBI accession.

-s, --sfile {file} (Mandatory)
   Subject sequence file (in FASTA format) or NCBI accession.

-ql, --qlabel {string} (Optional. Default: Query)
   Short ilabel identifying the query file. Valid characters are numbers, dots,
   dashes and underscores. No spaces allowed.

-sl, --slabel {string} (Optional. Default: Subject)
   Short label identifying the subject file. Valid characters are numbers, dots,
   dashes and underscores. No spaces allowed.

-p, --prog {string} (Optional. Default: ssearch36);
   Program that will be used to align the sequences. 
   Valid options are (case insensitive): blastp|ssearch36

-o, --oudir {path} (Optional, Default: see below)
   Output directory to store the results. Depending on the chosen  alignment
   program (prog), plus the labels for the Query (qlabel) and Subject (slabel)
   files, the default output directory will be:
             ./{prog}_{slabel}_vs_{qlabel}

-e, --evalue {float} (Optional. Default: 1e-4)
   E-value threshold to use in the alignments.

-m, --sub-matrix {string} (Optional. Default: BL50)
   Substitution matrix to use in the alignments. Any matrix supported by
   ssearch36 can be used.

-c, --coverage {float} (Optional. Default: 40.0)
   Threshold for Alignment coverage percentage to use in the alignments.

-cc, --cov-control {char} (Optional. Defaul: X)
   Controls how the coverage threshold will be applied to the alignments
   (case insensitive):
   X:  Coverage applies to at least one protein.
   B:  Coverage applies to both proteins.
   Q:  Coverage applies to the query protein only.
   S:  Coverage applies to the subject protein only.

-l, --max-len {int} (Optional. Default: no restructions)
   Maximum length for either query or subject proteins. This option can
   be used to filter results when long proteins when necessary.

-lc, --len-control {char} (Optional. Default: N)
   Controls how the max-len threshold will be applied to the query/subject
   proteins (case insensitive):
   X:  max-len applies to at least one protein.
   B:  max-len applies to both proteins.
   Q:  max-len applies to the query protein only.
   S:  max-len applies to the subject protein only.
   N:  max-len is ignored. Any length is allowed. This option supersedes
       any value passed through option -l.

-scs, --seq-comp-stats { T/F } (Optional. Default: T)
   Perform sequence composition statistics when calculating E-values with
   BlastP. SSEARCH always corrects for sequence composition bias.

-lcf, --low-complex-filter { T/F } (Optional. Default: F)
    Filter low complexity regions from alignments before calculating
    E-values. Only useful when using blastp for option -p.


HELP

    print $help;
    exit;
}



