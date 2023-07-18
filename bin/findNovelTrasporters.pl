#!/usr/bin/env perl -w

use warnings;
use strict;
use Data::Dumper;

use Getopt::Long;

use TCDB::CheckDependencies;


#==========================================================================
#  This script identifies transporters in a genome that have no
#  significant hits in TCDB.
#
#--------------------------------------------------------------------------
#  Written by: Arturo Medrano
#  Date:       12/05/2017
#==========================================================================



#==========================================================================
#Check dependencies

my @dependencies = ('grep', 'hmmtop', 'blastdbcmd',  'blastp', 'extractTCDB.pl');

my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;


#==========================================================================
#read command line

#Output directory where the data for this analysis will be stored
my $workDir = ".";

my $dataDir = undef;
my $plotsDir = undef;
my $seqPlotsDir = undef;

#Path to the directory with the TCDB blast database
my $blastDir  = undef;
my $tcblastdb  = undef;
my $pblastdb  = undef;

#Path to the proteome that will be analyzed (compressed or not)
my $proteomeFile = undef;


#Output file with the output of BLASTing the genome against TCDB.
my $blastOutputFile = undef;


#Files where the unknown transporters will be stored
my $outfile = undef;
my $htmlReport = undef;

#Minimal number of TMS expected in acceptable distant homologs
my $minTMS = 4;


#For queries that require results higher than this evalue
my $greaterThanEvalue = 1e-3;


#Evalue threshold to consider two remote homologs redundant.
my $redundantEvalue = 1e-5;


#How will the results be sorted (asc, desc, no)
my $sort = 'desc';

read_command_line_arguments();

#print Data::Dumper->Dump([$workDir, $dataDir, $blastDir, $tcblastdb, $pblastdb, $proteomeFile, $plotsDir
#			  $outfile, $minTMS, $greaterThanEvalue, $redundantEvalue, $sort],
#			 [qw(*workDir *dataDir *blastDir *tcblastdb *pblastdb *proteomeFile *plotsDir
#			     *outfile *minTMS *greaterThanEvalue  *redundantEvalue *sort)]);
#exit;



#==========================================================================
#To easily extract sequences from the proteome, generate a blast database.

create_proteome_blastdb();



#==========================================================================
#First run hmmtop in all proteins in the proteome, this way we will only
#compare proteins against TCDB that have at least the minimum number of
#of TMS.

print "Extracting sequence of proteins with at least $minTMS TMSs...\n";

my %all_trans_data  = ();
my $hmmtopOutfile      = "$dataDir/hmmtop.out";
my $seqsCandidatesFile = "$dataDir/seqsCandTransporters.faa";

getProteinsWithMinNumberOfTMS(\%all_trans_data, $hmmtopOutfile, $seqsCandidatesFile);


#print Data::Dumper->Dump([\%all_trans_data], [qw(*all_trans_data)]);
#print "\nCandidatesData: ", scalar(keys %all_trans_data), "\n";
#exit;



#==========================================================================
#Run blast of the candidate transporters against TCDB

print "Blasting candidate transporters against TCDB\n";

die "No sequence file with candidates transporters found: $seqsCandidatesFile" unless (-f $seqsCandidatesFile);

#File with blast output
$blastOutputFile = "$workDir/data/blastOutput.tsv";

my $fmtBlastOut = '7 qseqid sseqid qlen slen evalue pident nident length qcovs qstart qend sstart send';

unless (-f $blastOutputFile) {
  system "blastp -db $tcblastdb -evalue 10 -use_sw_tback -comp_based_stats f -seg no -outfmt '$fmtBlastOut' -query $seqsCandidatesFile -out $blastOutputFile";
}

die "TC Blast output file not found: $blastOutputFile" unless (-f $blastOutputFile);




#==========================================================================
#Parse TCDB blast output and discard all matches with lower e-values


my %candidates = ();

parse_tcblast_output($blastOutputFile, \%candidates, \%all_trans_data);

#verify if candidate transporters were found
unless (keys %candidates) {
  print "\n### No transporters with poor hits to TCDB were found! ####";
  exit;
}

#print Data::Dumper->Dump([\%candidates ], [qw(*candidates )]);
#print "Total candidates: ", scalar( keys %candidates ), "\n";
#exit;




#==========================================================================
#remove redundancy from the the top candidates to unknown transporters

my %bestCandidates = ();

remove_redundant_candidates(\%candidates, $seqsCandidatesFile, \%bestCandidates);

print Data::Dumper->Dump([\%bestCandidates ], [qw(*bestCandidates )]);
exit;



#==========================================================================
#Generate reports with the unknown transporters. Also generate the
#hydropathies for final candidates, this will help remove Candidates
#with bad TMS predictions (e.g. 4 TMS proteins that have fewer TMS).

print "Generating reports\n";



#HTML report
generate_reports();





print "\n### Finished ###\n";




###########################################################################
##                    Subroutine Defintions                              ##
###########################################################################


sub generate_reports {


  #----------------------------------------------------------------------
  #Parse HMMTOP output



  #Text report
  open (my $outh, ">", $outfile) || die $!;

  #HTML report
  open (my $hth, ">", $htmlReport) || die $!;


  my $hHeader =<<HEADER;
<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <title>Novel Transporters</title>
  </head>
  <br />
  <h1 style='text-align:center'>Novel Transporters</h1>
  <body>

HEADER

  print $hth $hHeader;

  print $outh "#Accession\tLength\tTMSs\tDescription\n";
  foreach my $acc (sort by_tmsNlength  keys %bestCandidates) {


    #Ignore accession if the sequence is partial
    next if ($bestCandidates{$acc}{annot} =~ /partial/i);

    #text report
    print $outh "$acc\t", $bestCandidates{$acc}{length}, "\t",
      $bestCandidates{$acc}{tms}, "\t", $bestCandidates{$acc}{annot}, "\n";


    print $hth "  <br /><hr style=\"border-style:solid; border-width:5px; color:black;\"/>\n\n";
    print $hth "  <p><b>$acc  $bestCandidates{$acc}{annot}</b></p>\n\n";

    #HTML report
    my $plotFile = "$plotsDir/${acc}.png";
    print "  Generating plot for: $acc\n";
    generate_plot($acc, $plotFile, $bestCandidates{$acc}{tmspos}) unless (-f $plotFile);

    my $pPath = "plots/${acc}.png";
    my $iPlot =<<PLOT;

  <br />
  <div>
    <a href="$pPath" target="_blank"><img src="$pPath";" alt="$acc hydropathy"></a>
  </div>
  <br />

PLOT

    print $hth $iPlot;

  }

  my $closeRep = <<CLOSE;
  </body>
</html>
CLOSE

  print $hth $closeRep;


  close $outh;
  close $hth;



}


#==========================================================================
#Sort final candidates by number of TMS and by length.

sub by_tmsNlength {

  my $alen = $bestCandidates{$a}{length};
  my $atms = $bestCandidates{$a}{tms};

  my $blen = $bestCandidates{$b}{length};
  my $btms = $bestCandidates{$b}{tms};


  if ($atms == $btms) {
    $blen <=> $alen;
  }
  else {
    $btms <=> $atms;
  }
}



#==========================================================================
#Generate hydropathy plot for a given accession

sub generate_plot {
  my ($id, $pfile, $coords) = @_;

  my $sfile  = "$seqPlotsDir/${id}.faa";

  die "Error: no tms coords found for protein: $id" unless (@{ $coords });

#  print Data::Dumper->Dump([$id, $pfile, $sfile, $coords ],
#			   [qw(*id *pfile *sfile *coords )]);
#  <STDIN>;



  #Extract the sequence of the protein
  my $cmd1 = qq(blastdbcmd -db $pblastdb -target_only -entry $id -out $sfile);
  system $cmd1 unless (-f $sfile && !(-z $sfile));

  my $qTMS = "--add-tms " . join(",", @{ $coords }) . ":cyan";
  my $cmd2 = qq(quod.py -q -l "$id" -o $pfile --width 15 --edgecolor blue --xticks 25 --no-tms +0  $qTMS --  $sfile);
  system $cmd2 unless (-f $pfile && !(-z $pfile));

}


#==========================================================================
#remove redundancy from the resulting sequences



sub remove_redundant_candidates {

  my ($cand, $seqsCand, $best) = @_;


  #-----------------------------------------------------------------
  #Get the sequences of the final candidate transporters

  #Save the IDs of the final candidates to a file
  my $idsCandidatesFile = "$dataDir/idsFinalCandidates.txt";

  unless (-f $idsCandidatesFile && !(-z $idsCandidatesFile)) {
    open (my $fidh, ">", $idsCandidatesFile) || die $!;
    map {print $fidh $_, "\n";} sort {$a cmp $b} keys %{$cand};
    close $fidh;

    die "Ids file not created or empty" unless (-f $idsCandidatesFile && !(-z $idsCandidatesFile));
  }


  #Extract the sequences of the final candidates
  my $seqFile = "$dataDir/finalCandidatesSeqs.faa";
  getSequencesForCandTransporters($cand, $idsCandidatesFile, $seqFile);


  #-----------------------------------------------------------------
  #Blast all vs all candidates to remove redundance

  my $blastOutFile = "$dataDir/finalCandBlast.out";

  unless (-f $blastOutFile && !(-z $blastOutFile)) {
    system "blastp -query $seqFile -subject $seqFile -evalue 10 -comp_based_stats f -seg no -outfmt '$fmtBlastOut' -out $blastOutFile";
  }

  #Verify that the blast output file exists
  die "Blast output file not generated or empty: $blastOutFile" unless (-f $blastOutFile && !(-z $blastOutFile));



  #-----------------------------------------------------------------
  #Identify the redundant candidates. Candidates with very different size
  #are allowed. When 2 redundant candidates have similar size, remove the
  #the smallest one.

  my %redundantCandidates = ();

  open (my $blasth, "<", $blastOutFile) || die $!;
  while (<$blasth>) {

    chomp;  next if (/^#/ || !$_);

    my ($qry, $sbj, $qlen, $slen, $ev, $pIdent, $nIdent, $alen, $qcovs, @kk) = split(/\t/, $_);


    #There must be IDs to continue
    die "Could not extract id for query: $qry" unless ($qry);
    die "Could not extract id for query: $sbj" unless ($sbj);


    #Remove self comparisons
    next if ($qry eq $sbj);


    #When there are redundant hits with very different size keep both.
    my $lenRatio = (sort {$b <=> $a} ($qlen/$slen, $slen/$qlen))[0];


    if ($ev < $redundantEvalue && $lenRatio < 1.8) {

      #When size is comparable, the smallest protein is redundant
      if ($qlen > $slen) {
	$redundantCandidates{$qry} = 1;
      }
      elsif ($qlen < $slen) {
	$redundantCandidates{$sbj} = 1;
      }
      else {

	#Same length: sort by accession and keep the first one
	my $redAcc = (sort {$a cmp $b} ($qry, $sbj))[0];
	$redundantCandidates{$redAcc} = 1;
      }
    }
  }
  close $blasth;



  #-----------------------------------------------------------------
  #Get the available annotations for each protein

  #remove '>' and genome name from fasta headers
  my @annot = map {s/\>//; s/\s+\[.+\]//; $_} split(/\n/, `grep '>' $seqsCand`);

  #Generate a hash with the annotations per protein
  my %func = ();
  map {  $func{$1} = $2 if(/(\S+)\s+(.+)/) } @annot;



  #-----------------------------------------------------------------
  #generate hash with final results

  foreach my $c (keys %{ $cand }) {

    unless (exists $redundantCandidates{$c}) {
      $best->{$c}->{annot}  = $func{$c};
      $best->{$c}->{tms}    = $cand->{$c}->{tms};
      $best->{$c}->{tmspos} = $cand->{$c}->{tmspos};
      $best->{$c}->{length} = $cand->{$c}->{qlen};
    }

  }
}









#==========================================================================
#Now discard all matches that involve lower evalues

sub parse_tcblast_output {

  my ($blastFile, $hits, $data) = @_;


  #Put hear the parsed raw blast output
  my $parsed = {};


  #-----------------------------------------------------------------
  #For candidate transporter get the blast hits against TCDB.

  open (my $blasth, "<", $blastFile) || die $!;
 HIT:while(<$blasth>) {

    chomp;

    #Ignore comment lines or empty lines
    next HIT if (/^#/ || !$_);


    my ($ncbi, $tcid, $qlen, $slen, $ev, $pIdent, $nIdent, $alen, $qcovs, @kk) = split(/\t/, $_);


    #There must be data to process
    die "$blastFile -> Line in does not contain parseable blast output: $_" unless ($tcid);

    #make sure the query protein has the minimum number of TMS... A little bit of
    #redunance to allow parsing of previously existing blast outputs.
    next HIT unless ($data->{$ncbi}->[1] >= $minTMS);


    #Collect all hits for each accession
    push (@{ $parsed->{$ncbi} }, {qid=>$ncbi, tcid=>$tcid, qlen=>$qlen, slen=>$slen, eval=>$ev,
				  pident=>$pIdent, qcovs=>$qcovs, tms=>$data->{$ncbi}->[1],
	                          tmspos=>$data->{$ncbi}->[2]});
  }
  close $blasth;


#  print Data::Dumper->Dump([$parsed ], [qw(*parsed )]);
#  exit;

  #-----------------------------------------------------------------
  #Now sort the blast hits by evalue

  foreach my $acc (keys %{ $parsed }) {

    my @matches = sort by_evalue @{ $parsed->{$acc} };

    if ($matches[0]->{eval} >= $greaterThanEvalue) {
      $hits->{$acc} = $matches[0];
    }
  }



  #-----------------------------------------------------------------
  #Add to the results and Accession that did not have a blast match

  foreach my $acc (keys %{ $data }) {
    unless (exists $parsed->{$acc}) {
      $hits->{$acc} = { qlen=>$data->{$acc}->[0], tms=>$data->{$acc}->[1] };
    }
  }

}




#==========================================================================
#Sort parsed blast matches by E-value

sub by_evalue {
  $a->{eval} <=> $b->{eval}
}







#==========================================================================
#Run hmmtop for all proteins in the proteome, keep proteins with the
#the minimum number of TMS.

sub getProteinsWithMinNumberOfTMS {

  my ($candTransData, $tmsOutfile, $candOutfile) = @_;

  #-----------------------------------------------------------------
  #Run HMMTOP

  print "    ... Running HMMTOP\n";

  my $cmd = "";
  if ($proteomeFile =~ /.+\.gz$/) {
    $cmd = qq(gunzip -c $proteomeFile | hmmtop -if=--);
  }
  elsif ($proteomeFile =~ /.+\.bz2$/) {
     $cmd = qq(bzcat $proteomeFile | hmmtop -if=--);
  }
  else {
    $cmd = qq(hmmtop -if=$proteomeFile );

  }
  $cmd .= " -of=$tmsOutfile -is=pseudo -pi=spred -sf=FAS";

  system $cmd unless (-f $tmsOutfile);


  #-----------------------------------------------------------------
  #Parse HMMTOP and keep proteins with the minimum number of TMS
  #Specified by the user


  print "    ... Parsing HHMTOP output\n";

  parseHMMTOP($candTransData, $tmsOutfile, $minTMS);

#  print Data::Dumper->Dump([$candTransData], [qw(*candiates )]);
#  exit;

  #-----------------------------------------------------------------
  #Generate a fasta file with the sequences of the candidate transporters.
  #These are the sequences that will be blasted against TCDB to search
  #for transporters not yet represented in TCDB.

  print "    ... Generating sequence file with candidate transporters\n";

  my $idFile = "$dataDir/candTransportersIDs.txt";
  getSequencesForCandTransporters($candTransData, $idFile, $candOutfile);

}



#===========================================================================
#Generate a fasta file with the sequences of the candidate transporters.
#These are the sequences that will be blasted against TCDB to search
#for transporters not yet represented in TCDB.

sub getSequencesForCandTransporters {
  my ($candsData, $idFile, $outfile) = @_;


  #Save the accesssions of the putative transporters into a file
  my @ids = sort {$a cmp $b} keys %{ $candsData };


  unless (-f $idFile && !(-z $idFile)) {

    open (my $idh, ">", $idFile) || die $!;
    print $idh join("\n", @ids), "\n";
    close $idh;

    unless (-f $idFile && !(-z $idFile)) {
      die "File with putative transporters ID not found or empty: $idFile";
    }
  }



  #Now retrieve the sequences of the putative transporters
  unless (-f $outfile && !(-z $outfile)) {
    system qq(blastdbcmd -dbtype prot -db $pblastdb -target_only -entry_batch $idFile -out $outfile);
  }


  #Verify that the sequences were properly extacted
  unless (-f $outfile && !(-z $outfile)) {
    die "Sequence file with candidate transporters not found or empty: $outfile";
  }
}





#==========================================================================
#Parse hmmtop output and  get the the proteins with the minimum number
#of TMSs.


sub parseHMMTOP {

  my ($candsData, $tmsFile, $mTMS) = @_;


  open (my $tmsh, "<", $tmsFile) || die $1;
  while(<$tmsh>) {
    chomp;

    my ($len, $acc, $ntms, $location) = (/\s+(\d+)\s+(\S+).+\s(IN|OUT)\s+(\d+)\s+(.+)/)? ($1, $2, $4, $5) : ();

    if ($ntms && ($ntms >= $mTMS)) {

      my @coords = split (/\s+/, $location);

      my @res = ();
      for(my $i=0; $i <= $#coords - 1; $i += 2) {
	push (@res, "$coords[$i]-$coords[$i+1]");
      }

      $candsData->{$acc} = [$len, $ntms, \@res];
    }

  }
  close $tmsh;

#  print Data::Dumper->Dump([ $candsData ], [qw( *candsData )]);
#  exit;

}




#==========================================================================

sub create_proteome_blastdb {


  print "Generating Blast DB of the input proteome\n";

  my $cmd = "";
  if ($proteomeFile =~ /.+\.gz$/) {
    $cmd = qq(gunzip -c $proteomeFile | makeblastdb -in -);
  }
  elsif ($proteomeFile =~ /.+\.bz2$/) {
     $cmd = qq(bzcat $proteomeFile | makeblastdb -in -);
  }
  else {
    $cmd = qq(makeblastdb -in $proteomeFile);

  }
  $cmd .= " -dbtype prot -input_type fasta -title 'Proteome DB' -parse_seqids -hash_index -out $pblastdb";


  my $pblastdbFile = "${pblastdb}.pin";
  system $cmd unless (-f $pblastdbFile);


  die "Could not generate proteome blast database: $pblastdb" unless (-f $pblastdbFile);


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
      "bdb|blastdb=s"          => \$tcblastdb,
      "p|proteome=s"           => \$proteomeFile,
      "tms|min-tms=s"          => \$minTMS,
      "e|min-evalue=f"         => \$greaterThanEvalue,
      "re|redundant-evalue=f"  => \$redundantEvalue,
      "h|help"                 => sub { print_help(); },
      "<>"                     => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);


  #Validate proteome file
  die "Proteome file does not exist ==> $proteomeFile" unless (-f $proteomeFile);


  #make sure working directory exists
  system "mkdir -p $workDir" unless ($workDir && -d $workDir);


  #now create the data directory where intermediate files will be saved (blast, hmmtop, etc.)
  $dataDir = "$workDir/data";
  system "mkdir $dataDir" unless (-d $dataDir);


  #Output file with final Unknown transporters in the input genome
  $outfile    = "$workDir/novelTransporters.txt";
  $htmlReport = "$workDir/novelTransporters.html";

  #Generate directory where hydropathy plots will be saved
  $plotsDir = "$workDir/plots";
  $seqPlotsDir = "$plotsDir/sequences";
  system "mkdir -p $seqPlotsDir" unless (-d $seqPlotsDir);




  #Validate the blast directory
  $blastDir = "$workDir/blastdb";
  system qq(mkdir $blastDir) unless (-d $blastDir);

  #Validate the tcblastdb
  if ( $tcblastdb ) {
    my $tcblastdbFile = "${tcblastdb}.pin";
    die "BlastDB not found in: $tcblastdb" unless (-f $tcblastdbFile);
  }
  else {
    $tcblastdb = "$blastDir/tcdb";
    my $tcblastdbFile = "${tcblastdb}.pin";
    #get the blast database
    print "Extracting TCDB blast database...\n";
    system qq(extractTCDB.pl -i tcdb -f blast -o $blastDir) unless (-f $tcblastdbFile);
  }


  #the name of the proteome blastdb
  $pblastdb = "$blastDir/proteome";
}






sub print_help {

  #
  # $errMsg: The error message to be diplayed
  #
  # $printHelp: boolean value indicating whether the help of the
  #             program will be displayed.
  #

  my $help = <<'HELP';

 This script searches for transporters in a genome that have no
 hits in TCDB.

 Input paramateres

 -wd, --workdir {path}  (Defalut: current directory)
    Path to the directory where the ouput and temporary files will be stored.

 -bdb, --blastdb {path}  (Default: ./blastdb/tcdb)
    Full path to the TCDB blast database that will be used.

 -p, --proteome {path}  (Mandatory)
    File in fasta format with the proteome that will be analyzed
    (the proteome file can be compressed in gz or bz2 format)
    Fasta headers must contain the ID (accession) of the
    protein followed by a space and the functional description
    (e.g. >AKM78775.1 AP4A hydrolase ...).

 -tms, --min-tms {integer}  (Default: 4)
    Minimum number of TMS that candidate transporters should have.

 -e, --min-evalue {float}  (default: 1e-3)
    Proteins with E-value greater than this value will be considered.

 -re,  --redundancy-evalue {float}  (default: 1e-5)
    Evalue threshold to consider two proteins redundant for the
    the purpose of presenting final results.

 -h, --help
    Display this help. This argument takes precedence over any other
    argument.

HELP

  print $help;
  exit;
}
