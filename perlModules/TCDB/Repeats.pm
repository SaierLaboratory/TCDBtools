package TCDB::Repeat;


#use warnings;

#To remove harmeless warnings from overriding accessors
no warnings;

use strict;
use Data::Dumper;
use Ref::Util qw(is_arrayref);
use Scalar::Util qw (looks_like_number);


#Other TCDB useful functions
use TCDB::Assorted;

#Bioperl
use Bio::DB::Fasta;
use Bio::SearchIO;
use Bio::SeqIO;



#For OOP in perl.
use Class::Struct;


#==========================================================================
# The purpose of this class is to provide functions to search for repeats
# under two conditions:
#
# 1. Given a set of sequences with well defined TMS coordinates,
#    Identify repeats within the same sequences.
#
# 2. Given a set of sequences with well defined TMS coordinates,
#    cut sequences at a specified location, and align all sequences
#    to the left of the cutting point to all sequences to the right.
#
# 3. Specify two regions in TMSs or coordinates, and compares all
#    sequences in region 1 with all sequences of region 2.
#    (Note: point 2 is a special case of this functionality).
#
# -------
# Date created:  02/20/2021
# Programmer:    Arturo Medrano
# Institution:   UCSD
#==========================================================================



struct ("TCDB::Repeat" =>
	{

	 ####  user definable variables  ####
	 'tmsFile'             => '$',  #Input file with the TMS locations for all seqs in the analysis
	 'tmsFormat'           => '$',  #Format of the tmsFile: tms or hmmtop (Default: hmmtop)
	 'seqsDir'             => '$',  #Sequences with all proteins to be analyzed (one file per sequence).
	 'seqsFile'            => '$',  #file with all full sequences (When seqsDir is not given)
	 'outDir'              => '$',  #Turn this flag to 1 if dependencies may not be executable.
	 'mode'                => '$',  #all (seqs in input file; default), each (one output dir per seq), debug

	 #These 4 options are incompatible with each other 
	 'repeatUnit'          => '$',  #The expected size of the repeat unit
	 'AAranges2search'     => '@',  #ranges of amino acids to search
	 'TMSranges2search'    => '$',  #ranges of TMSs to search
	 'cutSeqsPoint'        => '$',  #2 TMSs betwen which sequences will be split in two sections.

	 #parametners to control alignments
	 'tails'               => '$',  #when cutting TMSs or regions, this is the tail in aa to leave to the left and right.
	 'alignmentProg'       => '$',  #Alignment program to use (default: ssearch36).
	 'evalueCutoff'        => '$',  #Maximal evalue to accept in alignments (default: 0.01)
	 'coverageCutoff'      => '$',  #Minimal coverage in alignments (default: 0.7)
	 'identityCutoff'      => '$',  #Minimal Identity to accept alignments (default: 0.2)
	 'compStatistics'      => '$',  #Indicates whether E-values will be corrected by compositional biases (Default: 1)
	 'aaSubstMatrix'       => 'S',  #Amino acid substitution matrix (essentially for ssearch alignments)
	 'gsatShuffles'        => '$',  #for GSAT scores (default 1000)
	 'gsatCutoff'          => '$',  #Minimal GSAT score to accept alignments (default 4)
	 'maxHitsPerQuery'     => '$',  #Maximal number of hits to analyze per ssearch query (Default: 20);

	 #This parameters control reports
	 'goodHitsOnly'        => '$',  #Indicates whether or not reports will only contain significant alignments (default: 1)
	 'clobber'             => '$',  #Indicates whehter files will be overwritten (default: 0)


	 ### Internal Variables  ###
	 'compStatisticsStr'   => '$',  #String to use when running alignments depending on the value of compStatistics
	 'blastdbName'         => '$',  #Name of the BlastDB that will be created
	 'blastdbDirName'      => '$',  #The name of the BlastDB dir
	 'blastdbDir'          => '$',  #Full path to the BlastDB dir.
	 'procSeqDirName'      => '$',  #Name of the directory where sequences of cut regions will be stored
	 'procSeqDir'          => '$',  #Full path to the directory where sequences of cut regions will be stored
	 'procAlnDirName'      => '$',  #Name of the directory where alignment files will be stored
	 'procAlnDir'          => '$',  #Full path to the directory where alignment files will be stored
	 'procPlotsDirName'    => '$',  #Name of the directory where plots will be stored
	 'procPlotsDir'        => '$',  #Full path to the directory where plots will be stored
	 'procReportDirName'   => '$',  #Name of the directory where reports will be saved
	 'procReportDir'       => '$',  #Full path to Directory where reports will be saved
	 'coordsTMS'           => '%',  #TMS inferred by HMMTOP or any other method.
	 'topAlignments'       => '%',  #Where alignments with top scores will be stored
	}
       );



###########################################################################
#              Set up constructors for internal varaibles
###########################################################################

#==========================================================================
#Read file with TMS coordinates

sub tmsFile {

  my ($self, $file) = @_;

  if ( $file ) {
    die "Error: File with TMSs coords must exists ($file)." unless (-f $file);
    $self->{'TCDB::Repeat::tmsFile'} = $file;
  }


  return $self->{'TCDB::Repeat::tmsFile'};
}



#==========================================================================
#Read format of file with TMSs coordinates

sub tmsFormat {

  my ($self, $value) = @_;

  my $default_format = 'hmmtop';

  if ($value) {

    my $format = lc $value;

    #validate
    unless  ( $format =~ /hmmtop/ ) {
      die "Error: only hmmtop format is supported for now ($format).";
    }

    $self->{'TCDB::Repeat::tmsFormat'} = $format;
  }

  unless ($self->{'TCDB::Repeat::tmsFormat'}) {
    $self->{'TCDB::Repeat::tmsFormat'} = $default_format;
  }

  return $self->{'TCDB::Repeat::tmsFormat'};
}




#==========================================================================
#Read directory with sequence files, one sequence per protein to be
#analyzed

sub seqsDir {

  my ($self, $dir) = @_;

  if ( $dir ) {
    die "Error: Dir with input seqeunces must exists ($dir)." unless (-d $dir);
    $self->{'TCDB::Repeat::seqsDir'} = $dir;
  }

  return $self->{'TCDB::Repeat::seqsDir'};
}


#==========================================================================
#Read file with input sequences

sub seqsFile {

  my ($self, $file) = @_;

  if ( $file ) {
    die "Error: File with TMSs coords must exists ($file)." unless (-f $file);
    $self->{'TCDB::Repeat::seqsFile'} = $file;
  }


  return $self->{'TCDB::Repeat::seqsFile'};
}




#==========================================================================
#Read Output dir


sub outDir {

  my ($self, $dir) = @_;

  my $default_dir = './Repeats';

  if ( $dir ) {
    $self->{'TCDB::Repeat::outDir'} = $dir;
  }

  unless ($self->{'TCDB::Repeat::outDir'}) {
    $self->{'TCDB::Repeat::outDir'} = $default_dir;
  }

  return $self->{'TCDB::Repeat::outDir'};
}



#==========================================================================
#Read mode of operation for the program: all, each, debug

sub mode {

  my ($self, $value) = @_;

  my $default_mode = 'all';

  if ($value) {

    my $mode = lc $value;

    #validate
    unless  ( $value =~ /all|each|debug/ ) {
      die "Error: unrecognized mdoe $value";
    }

    $self->{'TCDB::Repeat::mode'} = $mode;
  }

  unless ($self->{'TCDB::Repeat::mode'}) {
    $self->{'TCDB::Repeat::mode'} = $default_mode;
  }

  return $self->{'TCDB::Repeat::mode'};
}


#==========================================================================
#Read mode of operation for the program: all, each, debug

sub tails {

  my ($self, $value) = @_;

  my $default = 5;

  if ($value) {

    #validate
    unless  ( $value =~ /^\d+$/ && $value > 0) {
      die "Error: tails should be given an integer > 0";
    }

    $self->{'TCDB::Repeat::tails'} = $value;
  }

  unless ($self->{'TCDB::Repeat::tails'}) {
    $self->{'TCDB::Repeat::tails'} = $default;
  }

  return $self->{'TCDB::Repeat::tails'};
}



#==========================================================================
#Read mode of operation for the program: all, each, debug

sub alignmentProg {

  my ($self, $value) = @_;

  my $default = 'ssearch36';

  if ($value) {

    my $prog = lc $value;

    #validate
    unless  ( $prog =~ /hmmtop/) {
      die "Error: Acceptable values are: ssearch36";
    }

    $self->{'TCDB::Repeat::alignmentProg'} = $prog;
  }

  unless ($self->{'TCDB::Repeat::alignmentProg'}) {
    $self->{'TCDB::Repeat::alignmentProg'} = $default;
  }

  return $self->{'TCDB::Repeat::alignmentProg'};
}


#==========================================================================
#Read E-value cutoff for the analysis


sub evalueCutoff {

  my ($self, $value) = @_;

  my $default = 1e-3;

  if ($value) {

    #validate
    unless  ( Scalar::Util::looks_like_number($value) && $value >= 0) {
      die "Error: Evalue should be a non-negative number.";
    }

    $self->{'TCDB::Repeat::evalueCutoff'} = $value;
  }

  unless ($self->{'TCDB::Repeat::evalueCutoff'}) {
    $self->{'TCDB::Repeat::evalueCutoff'} = $default;
  }

  return $self->{'TCDB::Repeat::evalueCutoff'};
}


#==========================================================================
#Read coverage cutoff for the analysis


sub coverageCutoff {

  my ($self, $value) = @_;

  my $default = 0.85;

  if ($value) {

    #validate
    unless  ( Scalar::Util::looks_like_number($value) && ($value >= 0 && $value <= 1.0)) {
      die "Error: Coverage should be a  number in the range [0, 1].";
    }

    $self->{'TCDB::Repeat::coverageCutoff'} = $value;
  }

  unless ($self->{'TCDB::Repeat::coverageCutoff'}) {
    $self->{'TCDB::Repeat::coverageCutoff'} = $default;
  }

  return $self->{'TCDB::Repeat::coverageCutoff'};
}


#==========================================================================
#Read the identity cutoff for the analysis

sub identityCutoff {

  my ($self, $value) = @_;

  my $default = 0.2;

  if ($value) {

    #validate
    unless  ( Scalar::Util::looks_like_number($value) && ($value >= 0 && $value <= 1.0)) {
      die "Error: Identity should be a number in the range [0, 1].";
    }

    $self->{'TCDB::Repeat::identityCutoff'} = $value;
  }

  unless ($self->{'TCDB::Repeat::identityCutoff'}) {
    $self->{'TCDB::Repeat::identityCutoff'} = $default;
  }

  return $self->{'TCDB::Repeat::identityCutoff'};
}


#==========================================================================
#Read flag indicating whether evalues will be corrected for compositional
#biases

sub compStatistics {

  my ($self, $value) = @_;

  my $default = 1;

  if ($value) {

    #validate
    unless  ( $value =~ /^[01]$/) {
      die "Error: composiiton flag only accepts 0 and 1 values.";
    }

    $self->{'TCDB::Repeat::compStatistics'} = $value;
  }

  unless ($self->{'TCDB::Repeat::compStatistics'}) {
    $self->{'TCDB::Repeat::compStatistics'} = $default;
  }

  return $self->{'TCDB::Repeat::compStatistics'};
}

#==========================================================================
#Amino acid substitution matrix (essentially for ssearch alignments)

sub aaSubstMatrix {

  my ($self, $value) = @_;

  my $default = 'BL50';

  #I'm not validating it because in principle the name can vary between
  #alignment programs and this design will allow for multiple alignment
  #programs.
  $self->{'TCDB::Repeat::aaSubstMatrix'} = $value if ($value);

  unless ($self->{'TCDB::Repeat::aaSubstMatrix'}) {
    $self->{'TCDB::Repeat::aaSubstMatrix'} = $default;
  }

  return $self->{'TCDB::Repeat::aaSubstMatrix'};
}





#==========================================================================
#Read number of shuffles to calculate GSAT scores

sub gsatShuffles {

  my ($self, $value) = @_;

  my $default = 1000;

  if ($value) {

    #validate
    unless  ( $value =~ /^\d+$/ && $value >= 500) {
      die "Error: The number of GSAT shuffles should be an integer >= 500.";
    }

    $self->{'TCDB::Repeat::gsatShuffles'} = $value;
  }

  unless ($self->{'TCDB::Repeat::gsatShuffles'}) {
    $self->{'TCDB::Repeat::gsatShuffles'} = $default;
  }

  return $self->{'TCDB::Repeat::gsatShuffles'};
}



#==========================================================================
#Read the gsat cutoff for alignments

sub gsatCutoff {

  my ($self, $value) = @_;

  my $default = 5;

  if ($value) {

    #validate
    unless  ( Scalar::Util::looks_like_number($value) && $value >= 4) {
      die "Error: The GSAT cutoff should be a number >= 4.";
    }

    $self->{'TCDB::Repeat::gsatCutoff'} = $value;
  }

  unless ($self->{'TCDB::Repeat::gsatCutoff'}) {
    $self->{'TCDB::Repeat::gsatCutoff'} = $default;
  }

  return $self->{'TCDB::Repeat::gsatCutoff'};
}



#
#==========================================================================
#Read max number of hits per ssearch query

sub maxHitsPerQuery {

  my ($self, $value) = @_;

  my $default = 20;

  if ($value) {

    #validate
    unless  ( $value =~/^\d+$/  && $value >= 1) {
      die "Error: Max hits per query should be  >= 1!";
    }

    $self->{'TCDB::Repeat::maxHitsPerQuery'} = $value;
  }

  unless ($self->{'TCDB::Repeat::maxHitsPerQuery'}) {
    $self->{'TCDB::Repeat::maxHitsPerQuery'} = $default;
  }

  return $self->{'TCDB::Repeat::maxHitsPerQuery'};
}


#==========================================================================
#Read flag indicating whether results will include only the hits passing
#thresholds or everything (meaning reporting no_hits for query proteins).

sub goodHitsOnly {

  my ($self, $value) = @_;

  my $default = 1;

  if ($value) {

    #validate
    unless  ( $value =~ /^[01]$/) {
      die "Error: Results format  flag only accepts 0 and 1 values.";
    }

    $self->{'TCDB::Repeat::goodHitsOnly'} = $value;
  }

  unless ($self->{'TCDB::Repeat::goodHitsOnly'}) {
    $self->{'TCDB::Repeat::goodHitsOnly'} = $default;
  }

  return $self->{'TCDB::Repeat::goodHitsOnly'};
}


#==========================================================================
#Read flag indicating whether automatically generated files will be
#overwritten.

sub clobber {

  my ($self, $value) = @_;

  my $default = 0;

  if ($value) {

    #validate
    unless  ( $value =~ /^[01]$/) {
      die "Error: Results format  flag only accepts 0 and 1 values.";
    }

    $self->{'TCDB::Repeat::clobber'} = $value;
  }

  unless ($self->{'TCDB::Repeat::clobber'}) {
    $self->{'TCDB::Repeat::clobber'} = $default;
  }

  return $self->{'TCDB::Repeat::clobber'};
}



#==========================================================================
#Read the name of the directory where blastDBs will be saved

sub blastdbDirName {

  my ($self, $value) = @_;

  my $default = 'blastdb';

  $self->{'TCDB::Repeat::blastdbDirName'} = $value if ($value);


  unless ($self->{'TCDB::Repeat::blastdbDirName'}) {
    $self->{'TCDB::Repeat::blastdbDirName'} = $default;
  }

  return $self->{'TCDB::Repeat::blastdbDirName'};
}


#==========================================================================
#Read format of file with TMSs coordinates

sub blastdbName {

  my ($self, $value) = @_;

  my $default = 'qprot';

  $self->{'TCDB::Repeat::blastdbName'} = $value if ($value);

  unless ($self->{'TCDB::Repeat::blastdbName'}) {
    $self->{'TCDB::Repeat::blastdbName'} = $default;
  }

  return $self->{'TCDB::Repeat::blastdbName'};
}


#==========================================================================
#Read name of directory where sequences of cut regions will be saved.

sub procSeqDirName {

  my ($self, $value) = @_;

  my $default = 'sequences';

  $self->{'TCDB::Repeat::procSeqDirName'} = $value if ($value);


  unless ($self->{'TCDB::Repeat::procSeqDirName'}) {
    $self->{'TCDB::Repeat::procSeqDirName'} = $default;
  }

  return $self->{'TCDB::Repeat::procSeqDirName'};
}


#==========================================================================
#Read name of directory where alignments of cut regions will be saved.

sub procAlnDirName {

  my ($self, $value) = @_;

  my $default = 'alignments';

  $self->{'TCDB::Repeat::procAlnDirName'} = $value if ($value);


  unless ($self->{'TCDB::Repeat::procAlnDirName'}) {
    $self->{'TCDB::Repeat::procAlnDirName'} = $default;
  }

  return $self->{'TCDB::Repeat::procAlnDirName'};
}


#==========================================================================
#Read name of directory where plots will be stored

sub procPlotsDirName {

  my ($self, $value) = @_;

  my $default = 'plots';

  $self->{'TCDB::Repeat::procPlotsDirName'} = $value if ($value);

  unless ($self->{'TCDB::Repeat::procPlotsDirName'}) {
    $self->{'TCDB::Repeat::procPlotsDirName'} = $default;
  }

  return $self->{'TCDB::Repeat::procPlotsDirName'};
}


#==========================================================================
#Read name of directory where reports will be stored

sub procReportDirName {

  my ($self, $value) = @_;

  my $default = 'reports';

  $self->{'TCDB::Repeat::procReportDirName'} = $value if ($value);

  unless ($self->{'TCDB::Repeat::procReportDirName'}) {
    $self->{'TCDB::Repeat::procReportDirName'} = $default;
  }

  return $self->{'TCDB::Repeat::procReportDirName'};
}







###########################################################################
#
#               Define Class functions from here
#
###########################################################################


#==========================================================================
#For this analysis to make sense, all sequences must have the same number
#of TMSs. Ideally, TMSs should be verfified before running this routine
#for highest reliability.

sub findRepeatsTMSranges {

  my $self = shift;


  #The output directory
  my $cmd1 = "mkdir -p  " . $self->outDir;
  system $cmd1 unless (-d $self->outDir);


  #------------------------------------------------------------
  #Verify TMS coordenates (HMMTOP format) were provided and seqsDir

  #Validate TMS file
  unless (-f $self->tmsFile && !(-z $self->tmsFile)) {

    #
    #No tmsFile! If seqsFile was given, run HMMTOP.
    #

    #Check that there was an input file with sequences given.
    die "Error: sequence file is mandatory if no TMS file was givn" unless ($self->seqsFile);

    #Run HMMTOP on input sequences
    my $tmfile = $self->outDir . "/hmmtop.out";
    my $cmd1 = "hmmtop -if=" . $self->seqsFile . " -of=$tmfile -sf=FAS -pi=spred -is=pseudo";
    system $cmd1 unless (-f $tmfile && !(-z $tmfile));
    die "Error: File with TMS coordinates not found -> $tmfile" unless (-f $tmfile);
    $self->tmsFile($tmfile);
  }


  #Validate directory with the input sequences
  unless ($self->seqsDir && -d $self->seqsDir) {

    #
    #No seqsDir! if seqsFile was given, generate sequences directory and fill
    #it with with one sequence per file
    #

    #Check that there was an input file with sequences given.
    die "Error: sequence file is mandatory if no sequence dir (option -d) was given" unless ($self->seqsFile);


    #Create the seqsDir within the output directory
    my $kkDir = $self->outDir . "/sequences";
    system "mkdir -p $kkDir" unless (-d $kkDir);
    die "Error: could not create sequence directory -> $kkDir" unless (-d $kkDir);
    $self->seqsDir($kkDir);


    #Split sequences and put them in the seqsDir
    $self->splitSeqs();

  }

  #The ranges of TMS to cut
  $self->checkTMSrangeFormat();



  #------------------------------------------------------------
  #parse coordinates of TMSs

  print "Reading TMS coordinates...\n";

  if ($self->tmsFormat eq 'hmmtop') {
    $self->parseHMMTOP();
  }
  else {
    die "Unsupported format for TMS coordinates.";
  }

#  print Data::Dumper->Dump([$self->coordsTMS], [qw(*coordsTMS)]);
#  exit;


  #------------------------------------------------------------
  #Generate a blast DB with the relevant sequences. This will
  #facilitate the extraction of protein regions with the
  #command: blastdbcmd

#  print "Making BlastDB...\n";
#  $self->makeBlastDB();


  #------------------------------------------------------------
  #Cut TMS ranges and add their coordinates to $self->coordsTMS

  print "Cutting sequence regions...\n";


  #The directory where cut sequences will be stored
  my $defDir = $self->outDir . "/" . $self->procSeqDirName;
  $self->procSeqDir($defDir) unless ($self->procSeqDir);
  system "mkdir -p " . $self->procSeqDir unless (-d $self->procSeqDir);
  unless (-d $self->procSeqDir) {
    die "Error: no sequences directory found: " . $self->procSeqDir;
  }

  $self->cutTMSranges();
#  print Data::Dumper->Dump([$self->coordsTMS ], [qw(*coordsTMS)]);
#  exit;


  #------------------------------------------------------------
  #Align sequences in different regions

  print "Aligning sequence regions and extracting results...\n";
  $self->alignSequences();


  #------------------------------------------------------------
  #Generate reports

  print "Generating reports...\n";
  $self->generateReports();

  print "\nFinished!\n";

}


#==========================================================================
#Take a file with sequences in fasta format, and generate 1 file per
#sequence

sub splitSeqs {

  my $self = shift;

#  print Data::Dumper->Dump([$self->seqsFile, $self->seqsDir ],
#			   [qw(*seqsFile *tmsFile *seqsDir )]);
#  exit;


  my $obj  = Bio::SeqIO->new(-file => $self->seqsFile , -format => "fasta");

  while(my $seqObj = $obj->next_seq) {

    my $id  = $seqObj->primary_id;
    my $seq = $seqObj->seq;

    die "Error: no sequence data for -> $id" unless ($id && $seq);

    my $outfile = $self->seqsDir. "/${id}.faa";
    open (my $outh, ">", $outfile) || die $!;
    print $outh ">$id\n$seq\n";
    close $outh;

    die "Error: problem generating file -> $outfile" unless (-f $outfile && !(-z $outfile));
  }
}


#==========================================================================
#Generate reports with the top alignments (one report per pair of
#compared regions).
#
#Two reports will be generated, in HTML and tsv format, per pair of
#compared regions.


sub generateReports {
  my $self = shift;

  #
  #Independently of what alignment program was used, the format of
  #$self->topAlignments is the same.
  #

  #Generate the directory where reports will be stored
  my $oDir = $self->outDir . "/" . $self->procReportDirName;
  $self->procReportDir($oDir) unless ($self->procReportDir);
  system "mkdir -p " . $self->procReportDir unless (-d $self->procReportDir);
  unless (-d $self->procReportDir) {
    die "Error: directory to save reports was not found: " . $self->procReportDir;
  }



 CMP:foreach my $rCmpID (keys %{ $self->topAlignments }) {

    #Report files for this comparison
    my $sumFile     = $self->procReportDir . "/summary_${rCmpID}.txt";
    my $htmlFile    = $self->procReportDir . "/report_${rCmpID}.html";

    open (my $htmlh, ">", $htmlFile)    || die $!;
    open (my $sumh,  ">", $sumFile)     || die $!;

    #header for HTML report
    print $htmlh $self->getHTMLheader($rCmpID);

    #Header for summary tsv report
    print $sumh "#Query\tSubject\tQ_len\tS_len\tE-value\tBitscore\tGSAT\tIdentity\tAln_len\tQ_cov\tS_cov\n";


    my $topResults = $self->topAlignments->{$rCmpID};

#    print Data::Dumper->Dump([$topResults], [qw(*topResults)]);
#    exit;

  QRY:foreach my $resh (@$topResults) {


#      print Data::Dumper->Dump([$resh ], [qw(*resh )]);
#      <STDIN>;


      #extract the query protein (key of hash)
      my $qName = (keys %$resh)[0];
      die "Error: could not extract query name from: \$resh" unless ($qName);

      my $qAcc  = (exists $resh->{$qName}->{qAcc})? $resh->{$qName}->{qAcc} :
	die "Error: could not extract query accession from: \$resh";

      #Unexpected error: there are NO hits for this protein
      die "Error: no hits for $qName... check what's going on!" unless (@{ $resh->{$qName}->{hits} });


      my ($qRegStart, $qRegEnd) = @{ $resh->{$qName}->{qRegion} };

      my $qFullLength = ($self->coordsTMS->{ $resh->{$qName}->{qAcc} }->{length})?
	$self->coordsTMS->{ $resh->{$qName}->{qAcc} }->{length} :
	die "Could not extract full legnth for query: $qAcc";

      my $qRegLength  =  $qRegEnd - $qRegStart + 1;


    HIT:foreach my $hit (@{ $resh->{$qName}->{hits} }) {


	my $hName = $hit->{hName};
	my $hAcc  = $hit->{hAcc};
	my $hFullLength = ($self->coordsTMS->{ $hit->{hAcc} }->{length})?
	  $self->coordsTMS->{ $hit->{hAcc} }->{length} :
	  die "Could not extract full legnth for hit: $hAcc";
	my ($hRegStart, $hRegEnd) = @{ $hit->{hRegion} };
	my $hRegLength  = $hRegEnd - $hRegStart + 1;


#	print Data::Dumper->Dump([$qName, $qAcc, $qFullLength, $qRegStart, $qRegEnd,
#				  $hName, $hAcc, $hFullLength, $hRegStart, $hRegEnd],
#				 [qw(*qName *qAcc *qFullLength *qRegStart *qRegEnd
#				     *hName *hAcc *hFullLength *hRegStart *hRegEnd)]);
#	<STDIN>;

	#Only print the first HSP;
	my $hspCnt = 1;
      HSP:foreach my $hsp (@{ $hit->{HSPs} }) {

	  last HSP if ($hspCnt > 1);


	  #print summary data
	  my $idStr   = sprintf ("%.2f", $hsp->{hId} * 100);
	  my $gsatStr = sprintf ("%.2f", $hsp->{gsat});
	  my $line = "$qName\t$hName\t$qRegLength\t$hRegLength\t" . $hsp->{hEvalue} . "\t";
	  $line .= $hsp->{hBits}  . "\t" . $idStr       . "\t"    .  $gsatStr . "\t";
	  $line .= $hsp->{alnLen} . "\t" . $hsp->{qCov} . "\t"    . $hsp->{hCov} . "\n";
	  print $sumh $line;


	  #Print HTML alignment parameters
	  my $alignmentParams = $self->htmlAlignmentParams($qName, $qFullLength, $qRegStart, $qRegEnd, $hName, $hFullLength, $hRegStart, $hRegEnd, $hsp);
	  print $htmlh $alignmentParams;

	  #print HTML for hydropathy plots
	  my $quodPlots = $self->htmlQuodPlots($qName, $hName, $hsp);
	  print $htmlh $quodPlots;


#	  print Data::Dumper->Dump([$qName, $hsp ], [qw(*qName *hsp )]);
#	  exit;

	  $hspCnt++;
	} #HSP
      } #HIT
    } #QRY


    #Close report files
    my $closeHTMLfile = $self->closeHTMLreport();
    print $htmlh $closeHTMLfile;

    close $htmlh;
    close $sumh;
  } #CMP

}


#==========================================================================
#close html report

sub closeHTMLreport {

  my $self = shift;

  my $closeRep = <<CLOSE;
  </body>
</html>
CLOSE

  return $closeRep;

}



#==========================================================================
#Generate HTML code to display the quod plots for a given hit

sub htmlQuodPlots {

  my ($self, $qname, $hname, $hsp) = @_;

  my $qplot = $hsp->{qPlot};
  my $hplot = $hsp->{hPlot};
  my $aPlot = $hsp->{alnPlot};

  my $quodPlots =<<PLOTS;

    <br />
    <table style="width:100%">
      <tr>
        <td><a href="$qplot" target="_blank"><img src="$qplot" alt="$qname"></a></td>
        <td><a href="$hplot" target="_blank"><img src="$hplot" alt="$hname"></a></td>
      </tr>
      <tr>
        <td colspan="2" style="text-align: center;">
           <a href="$aPlot" target="_blank"><img src="$aPlot" alt="$qname vs $hname alignment"></a>
        </td>
      </tr>
    </table>

PLOTS

  return $quodPlots;

}



#==========================================================================
#get HTML code for printing the parameters of the alignment

sub htmlAlignmentParams {

  my ($self, $qname, $qfLength, $qregStart, $qregEnd, $hname, $hfLength, $hregStart, $hregEnd, $hsp) = @_;


  #Format float numbers to be printed
  my $idStr   = sprintf ("%.2f", $hsp->{hId} * 100);
  my $gsatStr = sprintf ("%.2f", $hsp->{gsat});
  my $qcov    = sprintf ("%.2f", $hsp->{qCov} * 100);
  my $hcov    = sprintf ("%.2f", $hsp->{hCov} * 100);

  my $qrlength = $qregEnd - $qregStart + 1;
  my $hrlength = $hregEnd - $hregStart + 1;

  my $eval   = $hsp->{hEvalue};
  my $bits   = $hsp->{hBits};
  my $qstart = $hsp->{qAlnStart};
  my $qend   = $hsp->{qAlnEnd};
  my $hstart = $hsp->{hAlnStart};
  my $hend   = $hsp->{hAlnEnd};
  my $alnLen = $hsp->{alnLen};
  my $qseq   = $hsp->{qSeq};
  my $hseq   = $hsp->{hSeq};
  my $homStr = $hsp->{homStr};


  my $alnHit = <<HIT;

    <br /><hr style=\"border-style:solid; border-width:5px; color:black;\"/>

    <p><b>$qname ($qfLength) vs $hname ($hfLength)</b></p>

    <table width="900px" border="0" cellspacing="0" cellpadding="2">
      <tr>
         <td class='label'><b>E-value:</b></td>
         <td class='data'>$eval</td>
         <td class='label'><b>Bits:</b></td>
         <td class='data'>$bits</td>
         <td class='label'><b>GSAT:</b></td>
         <td class='data'>$gsatStr</td>
         <td class='label'><b>Identity:</b></td>
         <td class='data'>${idStr}%</td>
      </tr>
      <tr>
         <td class='label'><b>Q_cov:</b></td>
         <td class='data'>${qcov}%</td>
         <td class='label'><b>H_cov:</b></td>
         <td class='data'>${hcov}%</td>
         <td class='label'><b>Q_regLen:</b></td>
         <td class='data'>$qrlength</td>
         <td class='label'><b>H_regLen:</b></td>
         <td class='data'>$qrlength</td>
      </tr>
      <tr>
         <td class='label'><b>AlnLen:</b></td>
         <td class='data'>$alnLen</td>
         <td class='label'><b>Q_aln:</b></td>
         <td class='data'>${qstart}-$qend</td>
         <td class='label'><b>H_aln:</b></td>
         <td class='data'>${hstart}-$hend</td>
         <td class='label'></td>
         <td class='data'></td>
      </tr>
    </table>
    <br />

    <p><b>Alignment:</b></p>
    <div class='seq'>
    <pre>
$qseq
$homStr
$hseq
    </pre>
    </div>

HIT

  return $alnHit

}



#==========================================================================
#get the HTML code for the  header of the report

sub getHTMLheader {

  my ($self, $title) = @_;

  $title =~ s/_vs_/ vs /;


   my $header = <<HEADER;
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
    <title>$title</title>
  </head>
  <br />
  <h1 style='text-align:center'>Inferred Repeats for $title</h1>
  <body>

HEADER

  return $header;

}


#==========================================================================
#Align cut regions with user selected program. So far, only ssearch is
#supported.

sub alignSequences {

  my $self = shift;

  #Directory where alignments will be saved must exist
  my $alnDir = $self->outDir . "/" . $self->procAlnDirName;
  $self->procAlnDir($alnDir) unless ($self->procAlnDir);
  system "mkdir -p " . $self->procAlnDir unless (-d $self->procAlnDir);
  unless (-d $self->procAlnDir) {
    die "Error: no sequences directory found: " . $self->procAlnDir;
  }


  #Directory where plots will be saved must exist
  my $pltDir = $self->outDir . "/" . $self->procPlotsDirName;
  $self->procPlotsDir($pltDir) unless ($self->procPlotsDir);
  system "mkdir -p " . $self->procPlotsDir unless (-d $self->procPlotsDir);
  unless (-d $self->procPlotsDir) {
    die "Error: no plots directory found: " . $self->procPlotsDir;
  }


 IDX1:foreach my $idx1 (1 .. scalar @{ $self->TMSranges2search }) {

    my $scfile1 = "Region_${idx1}";

  IDX2:foreach my $idx2 ($idx1 + 1 .. scalar @{ $self->TMSranges2search }) {

      my $scfile2 = "Region_${idx2}";

      print "  Workign with:    $scfile1 vs $scfile2\n";
      if ($self->alignmentProg eq 'ssearch36') {
	$self->ssearchAlign($scfile1, $scfile2);
      }
      else {
	die "Unexpected error: " . $self->alignmentProg;
      }
    }
  }
}



#==========================================================================
#Perform sequence alignments using ssearch36

sub ssearchAlign {
  my ($self, $fname1, $fname2) = @_;

  #Verify that sequence files exist
  my $file1 = $self->procSeqDir . "/${fname1}.faa";
  my $file2 = $self->procSeqDir . "/${fname2}.faa";
  die "File with sequence regions not found: $file1" unless (-f $file1);
  die "File with sequence regions not found: $file2" unless (-f $file2);

  die "Aligment dir must exist: " . $self->procAlnDir unless (-d $self->procAlnDir);

  #File where alignment will be saved
  my $ofile = $self->procAlnDir . "/ssearch_${fname1}_vs_${fname2}.out";

  if ($self->compStatistics) {
    $self->compStatisticsStr('-z 11 -k 1000');
  }

  #Build ssearch comand.
  #-C 25 generates buggy behavior inserting spaces in the sequence
  my $cmd = "-p " . $self->compStatisticsStr . " -E " . $self->evalueCutoff . " -W 0 ";
  $cmd .= "-s " . $self->aaSubstMatrix . " $file1 $file2 > $ofile";

  #run ssearch
  system "ssearch36 $cmd" unless (-f $ofile && !(-z $ofile));

  #Parse refults file for good hits
  $self->parseSSEARCHoutput($fname1, $fname2, $ofile);
}



#==========================================================================
#Parse the output of ssearch for the comparisons of two sequence regions

sub parseSSEARCHoutput {

  my ($self, $qregion, $sregion, $alnfile) = @_;

  #To name results for the two files being aligned
  my $alnID = "${qregion}_vs_$sregion";

  my $parser = new Bio::SearchIO (-format => 'fasta', -file => $alnfile);

  my %topHits = ();

 RES:while (my $result = $parser->next_result) {

    my %hits = ();

    my $qRegLen  = $result->query_length;
    my $qName    = $result->query_name;
    my $qDesc    = $result->query_description;

    #Get key to extract the coordinates of the cut region in the query
    my ($qAcc, $qReg) = ($qName =~ /(\S+)_region(\d+)/)? ($1, $2) :
      die "Could not extract accesion and region from: $qName";

    my $qKeyRegion = "Region_$qReg";
    my $qRegion = $self->coordsTMS->{$qAcc}->{$qKeyRegion};

    my $qFullLength = ($self->coordsTMS->{$qAcc}->{length})?
      $self->coordsTMS->{$qAcc}->{length} :
      die "Could not extract full legnth for query: $qAcc";


    my $nhits = 1;


  HIT:while (my $hit = $result->next_hit) {

      #Only analyze the hits indicated by the user
      last HIT if ($nhits > $self->maxHitsPerQuery);

      my %hit = ();

      my $hName  = $hit->name();
      my $hDesc  = $hit->description();
      my $hspCnt = 0;
      my @HSPs   = ();

      #Get the key to extract the coordinates of the cut region in hit
      my ($hAcc, $hReg) = ($hName =~ /(\S+)_region(\d+)/)? ($1, $2) :
	die "Could not extract accesion and region from: $hName";
      my $hKeyRegion = "Region_$hReg";
      my $hRegion = $self->coordsTMS->{$hAcc}->{$hKeyRegion};
      my $hFullLength = ($self->coordsTMS->{$hAcc}->{length})?
	$self->coordsTMS->{$hAcc}->{length} :
	die "Could not extract full legnth for hit: $hAcc";


    HSP:while(my $hsp = $hit->next_hsp) {

#	print Data::Dumper->Dump([$hsp], [qw(*hsp )]);
#	<STDIN>;

	$hspCnt++;
	my %tmp = ();

	my $alnLen  = $hsp->hsp_length;
	my $hRegLen = $hit->length;
	my $hEvalue = $hsp->evalue;
	my $hBits   = $hsp->bits;
	my $hId     = $hsp->frac_identical('total');  #identity in the alignment
	my $hSim    = $hsp->frac_conserved('total');  #similarity in the alignment

	#coordinates in the alignment to properly calculate coverages
	my $qstart  = $hsp->start('query');
	my $qend    = $hsp->end('query');
	my $sstart  = $hsp->start('subject');
	my $send    = $hsp->end('subject');


	#Calculate coverages properly (do not use alignment length as it includes gaps
	my $qCov     = ($qend - $qstart + 1) / $qRegLen;
	my $hCov     = ($send - $sstart + 1) / $hRegLen;

	#Before storing hit results check minimum coverage, identity and evalue
	next HSP unless (($qCov >= $self->coverageCutoff || $hCov >= $self->coverageCutoff) &&
			 ($hEvalue <= $self->evalueCutoff) && ($hId >= $self->identityCutoff));

#	print Data::Dumper->Dump([$qName, $qDesc, $Nname, $hDesc, $qLen, $qCov, $hLen, $hCov, $self->coverageCutoff,
#				  $hEvalue, $self->evalueCutoff, $hBits, $hId, $self->identityCutoff],
#				 [qw(*qName *qDesc *hName *hDesc *qLen *qCov $hLen *hCov *coverageCutoff *evalue *evalCutoff *bitscore *hId *IDcutoff)]);
#	<STDIN>;


#	print "   Hit: $nhits    Processing: $qName vs $hName    HSP: $hspCnt\n";

	#hit identity
	$tmp{hName}   = $hName;
	$tmp{hDesc}   = $hDesc;
	$tmp{hRegLen} = $hRegLen;

	#hit statistics
	$tmp{hsp}     = $hspCnt;
	$tmp{alnLen}  = $alnLen;
	$tmp{hEvalue} = $hEvalue;
	$tmp{hBits}   = $hBits;
	$tmp{hId}     = $hId;
	$tmp{hSim}    = $hSim;
	$tmp{qCov}    = $qCov;
	$tmp{hCov}    = $hCov;

	#The alignment
	$tmp{qstart}  = $qstart;
	$tmp{qend}    = $qend;
	$tmp{hstart}  = $sstart;
	$tmp{hend}    = $send;

	$tmp{qSeq}    = $hsp->query_string;
	$tmp{hSeq}    = $hsp->hit_string;
	$tmp{homStr}  = $hsp->homology_string;


	#------------------------------------------------------------
	#Get the GSAT score

	my $gsat_outFile = $self->procAlnDir . "/${qName}_vs_${hName}_hsp${hspCnt}.gsat";

	#remove dashes and spaces from sequences in order to run GSAT (this is becaue labels lengeth
	my $gsatCmd = "gsat.py $tmp{qSeq} $tmp{hSeq} " . $self->gsatShuffles . " > $gsat_outFile";
#	print "$gsatCmd\n";

	unless (-f $gsat_outFile && !(-z $gsat_outFile)) {
	  system "gsat.py $tmp{qSeq} $tmp{hSeq} " . $self->gsatShuffles . " > $gsat_outFile";
	}

	my $gsat_score = TCDB::Assorted::get_gsat_score ($gsat_outFile);
	$tmp{gsat} = $gsat_score;

	#delete gsat file to save disk space
#	system "rm $gsat_outfile" if (-f $gsat_outFile);

	#------------------------------------------------------------
	#Generate 3 hydropathy plots:
	#   a) full length query and hit proteins with aligned regions marked with bars and wedges
	#      plusq shaded regions denoting the segment covered by the alignment.
	#   b) alignment of the 2 sequence regions.


	#Calculate coordinates for the aligned regions by ssearch
	my $qAln_start  = $qRegion->[0] + ($qstart - 1);
	my $qAln_end    = $qRegion->[0] + ($qend   - 1);
	$tmp{qAlnStart} = $qAln_start;
	$tmp{qAlnEnd}   = $qAln_end;

	my $hAln_start = $hRegion->[0] + ($sstart - 1);
	my $hAln_end   = $hRegion->[0] + ($send - 1);
	$tmp{hAlnStart} = $hAln_start;
	$tmp{hAlnEnd}   = $hAln_end;


	#Format the coordinates for the bar/wedges delimiting of cut regions
	my $qBars = "-w " . $qRegion->[0] . "-" . $qRegion->[1] . "::1 ";
	my $hBars = "-w " . $hRegion->[0] . "-" . $hRegion->[1] . "::1 ";

	#Format the coordinates for the repeats now
	my $qAln = "${qAln_start}-${qAln_end}:red";
	my $hAln = "${hAln_start}-${hAln_end}:blue";


	#Format TMS bars for quod
	my @qPeaks = map { join ("-", @$_) . ":orange" } @{ $self->coordsTMS->{$qAcc}->{tms} };
	my $qTMS   = join (" ", @qPeaks);

	my @hPeaks = map { join ("-", @$_) . ":cyan" } @{ $self->coordsTMS->{$hAcc}->{tms} };
	my $hTMS   = join (" ", @hPeaks);


	#format the xtics for the plots
	my $qXtics   = $self->xtics4plot($qFullLength);
	my $hXtics   = $self->xtics4plot($hFullLength);
	my $alnXtics = $self->xtics4plot($alnLen);


	#Sequence and plot files for full protein plots.
	#
	#NOTE: Need to include the alignment info in the name because full protein plots include the
	#      region aligned with the other protein.
	my $qFullSeq  = $self->seqsDir . "/${qAcc}.faa";
	my $qPlotName = "${qName}_vs_${hName}_fullProt_${qName}_hsp${hspCnt}";
	my $qPlotFile = $self->procPlotsDir . "/$qPlotName";

	my $hFullSeq  = $self->seqsDir . "/${hAcc}.faa";
	my $hPlotName = "${qName}_vs_${hName}_fullProt_${hName}_hsp${hspCnt}";
	my $hPlotFile = $self->procPlotsDir . "/$hPlotName";

	#The quod command line for full query protein
	my $qCmd = "quod.py -c red -t png -l '$qName' -o $qPlotFile -q -r 80 $qBars --xticks $qXtics -nt +0 -at $qTMS ${qAln} -- $qFullSeq";
	my $qImg = "${qPlotFile}.png";
	system $qCmd unless (-f $qImg);
	die "Could not generate plot: $qImg" unless (-f $qImg);
	$tmp{qPlot} = "../" . $self->procPlotsDirName . "/${qPlotName}.png";


	#The quod command line for full hit protein
	my $hCmd = "quod.py -c blue -t png -l '$hName' -o $hPlotFile -q -r 80 $hBars --xticks $hXtics -nt +0 -at $hTMS ${hAln} -- $hFullSeq";
	my $hImg = "${hPlotFile}.png";
	system $hCmd unless (-f $hImg);
	die "Could not generate plot: $hImg" unless (-f $hImg);
	$tmp{hPlot} = "../" . $self->procPlotsDirName . "/${hPlotName}.png";


	#Sequence files for the alignments
	my $qAlnFile = $self->procSeqDir . "/${qName}_vs_${hName}_alnSeq_${qName}_hsp${hspCnt}.faa";
	unless (-f $qAlnFile && !(-z $qAlnFile)) {
	  open(my $qh, ">", $qAlnFile) || die $!;
	  print $qh ">$qName alignment\n", $hsp->query_string, "\n";
	  close $qh;
	  die "Error: Problem generating alignment sequence file: $qAlnFile" unless (-f $qAlnFile && !(-z $qAlnFile));
	}

	my $hAlnFile = $self->procSeqDir . "/${qName}_vs_${hName}_alnSeq_${hName}_hsp${hspCnt}.faa";
	unless (-f $hAlnFile && !(-z $hAlnFile)) {
	  open(my $hh, ">", $hAlnFile) || die $!;
	  print $hh ">$hName alignment\n", $hsp->hit_string, "\n";
	  close $hh;
	  die "ErrorK Problem generating alignment sequence file: $hAlnFile" unless (-f $hAlnFile && !(-z $hAlnFile));
	}

	#alignment plot file
	my $alnPlotName = "${qName}_vs_${hName}_hsp${hspCnt}.png";
	my $alnPlotFile =  $self->procPlotsDir . "/$alnPlotName";


	#Run quod for the alignment
	my $alnCmd = qq(alnquod.py --grid -q -l "$qName (red) and $hName (blue)" -o $alnPlotFile --xticks $alnXtics --width 15 -- $qAlnFile $qFullSeq $hAlnFile $hFullSeq);
	system $alnCmd unless (-f $alnPlotFile);
	return "Could not generate alignment plot: $alnPlotFile" unless (-f $alnPlotFile);
	$tmp{alnPlot} = "../" . $self->procPlotsDirName . "/$alnPlotName";

#	print Data::Dumper->Dump([\%tmp ], [qw(*hsp )]);
#	print "$alnCmd\n\n$alnPlotFile\n\n";
#	<STDIN>;


	push (@HSPs, \%tmp);
      } #HSP


      #If there were significant alignments save the hit
      if (@HSPs) {

	$hits{$qName}{qName}   = $qName;
	$hits{$qName}{qLength} = $qFullLength;
	$hits{$qName}{qRegLen} = $qRegLen;
	$hits{$qName}{qDesc}   = $qDesc;
	$hits{$qName}{qAcc}    = $qAcc;
	$hits{$qName}{qRegion} = $qRegion;

	$hit{hName}   = $hName;
	$hit{hLength} = $hFullLength;
	$hit{hDesc}   = $hDesc;
	$hit{hAcc}    = $hAcc;
	$hit{hRegion} = $hRegion;
	$hit{HSPs}    = \@HSPs;

	push (@{ $hits{$qName}{hits} }, \%hit);

#	print Data::Dumper->Dump([\%hits ], [qw(*hits)]);
#	<STDIN>;

	$nhits++;
      }


    } #HIT


    #Add results to 
    if ( exists $hits{$qName} && @{ $hits{$qName}{hits} }) {
      push(@{ $self->topAlignments->{$alnID} },  \%hits);
    }

  } #RES

}


#==========================================================================
#Calculate the xtics specing for a quod plot given the length of the
#input sequence

sub xtics4plot {

  my ($self, $protLen) = @_;

  my $xticksSpacing = undef;

  if ($protLen <= 100) {
    $xticksSpacing = 10;
  }
  elsif ($protLen <= 300) {
    $xticksSpacing = 20;
  }
  elsif ($protLen <= 500) {
    $xticksSpacing = 25;
  }
  elsif ($protLen <= 500) {
    $xticksSpacing = 25;
  }
  elsif ($protLen <= 1000) {
    $xticksSpacing = 50;
  }
  else {
    $xticksSpacing = 100;
  }

  return $xticksSpacing;
}





#==========================================================================
#For each input sequences extract the regions of TMS defined by the user.

sub cutTMSranges {

  my $self = shift;

  #Verify whether files with cut regions already exist
  my $fileExist = 0;
  foreach my $rangeN (1 .. scalar @{ $self->TMSranges2search }) {
    my $scfile = $self->procSeqDir . "/Region_${rangeN}.faa";
    $fileExist++ if (-f $scfile && !(-z $scfile));
  }


  #Don't cut sequences if all expected files already exist and are
  #not empty
  my $cutRegions = ($fileExist == scalar @{ $self->TMSranges2search })? 0 : 1;

  #If regions will be cut, remove any sequence files that might exist.
  #This is because the following code will append sequences to files and
  #the process needs to start clean.
  my $cmd =  "rm " . $self->procSeqDir . "/Region_*.faa &>/dev/null";
  system $cmd if ($cutRegions);


  #index db to access and cut sequences
  my $db = Bio::DB::Fasta->new($self->seqsDir);
#  print Data::Dumper->Dump([ $db ], [qw(*fastaDB )]);


  #Open files where cut sequences will be added.
 RANGE:foreach my $rangeN (1 .. scalar @{ $self->TMSranges2search }) {

    my $outfile = $self->procSeqDir . "/Region_${rangeN}.faa";

    #Open file only if regions will be cut
    my $fileh = undef;
    if ($cutRegions) {
      open($fileh, ">>", $outfile) || die "Could not open: $outfile -> $!";
    }

    #Cut and save sequence regions
  PROT:foreach my $acc (keys %{ $self->coordsTMS }) {

      #Get the ranges of TMSs to compare
      my @range = @{ $self->TMSranges2search->[$rangeN - 1] };


      #extract the coordinates of the TMS for this protein
      my @tms = ($self->coordsTMS->{$acc} && @{ $self->coordsTMS->{$acc}->{tms} })?
	@{ $self->coordsTMS->{$acc}->{tms} } : undef;
      next PROT unless (@tms);


      #Get the coordinates that will be used to cut sequences for this range
      my $left  = $tms[$range[0] - 1]->[0] - $self->tails;
      $left = 1 if ($left <= 0);

      my $right = $tms[$range[1] - 1]->[1] + $self->tails;
      if ($right > $self->coordsTMS->{$acc}->{length}) {
	$right = $self->coordsTMS->{$acc}->{length};
      }

      #Add coordinates of this region to $self->coordsTMS, this will allow
      #adding this information to top alignment results.
      my $regionLabel = "Region_${rangeN}";
      $self->coordsTMS->{$acc}->{$regionLabel} = [$left, $right];


      #Cut regions here if sequence files were not found.
      my $subseq   = "";
      my $fastaSeq = "";
      if ($cutRegions) {

	#Extracting sequence region
	$subseq = $db->seq($acc, $left, $right);

	#subsequence in fasta format
	$fastaSeq = ">${acc}_region$rangeN aa:${left}-$right\n$subseq\n";

	print $fileh $fastaSeq;
      }

#      print Data::Dumper->Dump([$acc, $rangeN, \@range, \@tms, $left, $right, $subseq ],
#			       [qw(*acc *rangeN *range *tms *left *right *subseq)]);
#      <STDIN>;

    } #PROT

    close $fileh if ($cutRegions);
  } #RANGE

}



#==========================================================================
#Verify that the TMS ranges were given in the right format

sub checkTMSrangeFormat {

  my $self = shift;

  unless (Ref::Util::is_arrayref($self->TMSranges2search)) {
    die "Error: ranges should be given as an array reference!";
  }

  unless (@{ $self->TMSranges2search }) {
    die "Error: Ranges of TMSs to compare must be given!";
  }
  unless (scalar @{ $self->TMSranges2search } > 1) {
    die "Error: At least 2 TMS ranges shoul be provied!";
  }


  #First make sure that the TMS ranges have the right format
  foreach my $range (@{ $self->TMSranges2search }) {

    unless (Ref::Util::is_arrayref($range)) {
      die "Error: each range should be an array of two integers";
    }

    unless (scalar (@{ $range }) == 2) {
      die "Error: TMS ranges should be given as pairs of TMS numbers per range!";
    }

    foreach my $TMSn (@{ $range }) {
      unless ($TMSn =~ /^\d+$/ && $TMSn > 0) {
	die "Error: Numbers in ranges should be integers > 0 -> ($TMSn).";
      }
    }
  }
}








#==========================================================================
#Generate Blast DB with relevant sequences. This is to facilitate the
#extraction of protein regions with the command: blastdbcmd

sub makeBlastDB {

  my $self = shift;

  #generate directory where blastDB will be created.
  my $tmpDir = $self->outDir . "/" . $self->blastdbDirName;
  $self->blastdbDir($tmpDir) unless ($self->blastdbDir);

  #Create BlastDB directory
  my $cmd1 = "mkdir -p " . $self->blastdbDir;
  system $cmd1 unless (-d $self->blastdbDir);


  #command to generate the BlastDB
  my $seqsF =  $self->blastdbDir . "/" . $self->blastdbName . "_all.faa";
  my $cmd2 = "find " . $self->seqsDir . " -type f -name \"*.faa\" | xargs cat >> $seqsF";
#  print "$cmd2\n";
  system $cmd2 unless (-f $seqsF);

  my $path2db = $self->blastdbDir . "/" . $self->blastdbName;
  my $cmd3 = "makeblastdb -in $seqsF -dbtype prot -blastdb_version 5 ";
  $cmd3 .= "-input_type fasta -parse_seqids  -hash_index ";
  $cmd3 .= "-title " . $self->blastdbName . " -out $path2db";
#  print "$cmd3\n";
  system $cmd3 unless (-f "${path2db}.pin");

  die "Error: not able to generate BlastDB: $path2db" unless (-f "${path2db}.pin");
}







#==========================================================================
#Parse the TMS coordinates in HMMTOP format for each protein.

sub parseHMMTOP {

  my $self = shift;


  open(my $fh, "<", $self->tmsFile) || die $!;

 PROT:while(<$fh>) {

    chomp;

    #Remove trailing spaces
    s/^\s+//;
    s/\s+$//;

    #ignore empty lines
    next PROT unless ($_);


    #parse hmmtop line
    my ($len, $id, $ntms, $tms_str) = (/\S+\s+(\d+)\s+(\S+).+(IN|OUT)\s+(\d+)\s+([\d\s-]+)/)? ($1, $2, $4, $5) : ();


    #For debugging purposes
#    next unless ($id =~ /WP_100644534/);


    #Jump to next protein is sequence file does not exist
    my $seqF = $self->seqsDir . "/${id}.faa";
    next PROT unless (-f $seqF);


    if ($id && $ntms && $tms_str) {

      #extract the pairs of coordinates for TMS
      my @coords = split(/\s+/, $tms_str);
      my @tms = ();
      for (my $i=0; $i < $#coords; $i += 2) {
	push (@tms, [$coords[$i], $coords[$i+1]]);
      }

      $self->coordsTMS->{$id} = {length=>$len, ntms=>$ntms, tms=>\@tms};

    }
    else {
      print "problem parsing HMMTOP line: $_\n";;
      print Data::Dumper->Dump([$id, $ntms, $tms_str ], [qw(*id *ntms *tms_str )]);
      exit;
    }
  }
}



# #==========================================================================
# # This function verifies that the list of dependencies are installed
# # in the systems $path variable and/or in a set of user supplied
# # directories.
# #


# sub checkDependencies {

#     my ($self, $ar_progList) =  @_;

#     my $output = {};

#     #Check if a list of programs were given as argument to this function.
#     #If the list of programs is not detected it is assumed that it was given
#     #to the object when it was initialized.
#     if ($ar_progList) {
# 	$self->dependencies_list($ar_progList);
#     }


#     #Check if the list of dependencies to search is not available
#     unless (@{ $self->dependencies_list } || $self->allow_no_dependencies) {
#       die "Error: you must specify a list of dependencies for this script.\n\n",
# 	"Source -->";
#     }


#     #Verify that every program is installed here
#   PROG:foreach my $prog (@{ $self->dependencies_list }) {


#       #Check if program is in variable $path
#       my $whichOutput = $self->iWhich($prog);

#       #Dependency in $path and is executable
#       if (! $whichOutput ) {
# 	  next PROG;
#       }


#       #Dependency in $path but is NOT executable
#       elsif ($whichOutput == 102) {
# 	  push (@{ $output->{102} }, $prog);
# 	  next PROG;
#       }


#       #To this point the program is not in variable $path, check if the user
#       #supplied additional dirs to search
#       foreach my $upath (@{ $self->user_dirs }) {

# 	  my $uProg = "$upath/$prog";

# 	  #Dependency is in user-supplied directory and is executable
# 	  if (-x $uProg) {
# 	      next PROG;
# 	  }

# 	  #Dependency is in user-supplied directory but is NOT executable.
# 	  elsif (-e $uProg) {
# 	      push (@{ $output->{102} }, $uProg);
# 	      next PROG;
# 	  }
#       }


#       #Dependency was not detected at all
#       push (@{ $output->{101} }, $prog);
#     }


#     #If dependency errors where found, report them to user and exit;
#     if ( %{ $output } ) {

#       my $errorMsg = "";

#       foreach my $error (keys %{ $output }) {

# 	if ($error == 101) {
# 	  $errorMsg .= "Error: Make sure the following program(s) is(are) installed in the system: \n  " .
# 	    join ("\n  ", @{ $output->{$error} }) . "\n\n";
# 	}	
# 	elsif ($error == 102) {
# 	  unless ($self->allow_nonexecutables) {
# 	    $errorMsg .= "Error: the following program(s) is(are) not executable: \n  " .
# 	      join ("\n  ", @{ $output->{$error} }) . "\n\n";
# 	  }
# 	}
# 	else {
# 	  die "Unknown error code: $error -->";
# 	}
#       }

#       die $errorMsg, "Source --> " if ($errorMsg);
#     }
# }



# #==========================================================================
# #This function verifies that a dependency is in the environment variable $path.
# #Fuction returns an integer, which is a code indication the success or failure in
# #finding the dependency. These are the  return codes:
# #
# # 100: Empty dependency
# # 101: dependency not detected
# # 102: dependency was found but NOT  executable


# sub iWhich {

#     my ($self, $dependency) = @_;


#     #Dependency can't be empty
#     return 100 unless ($dependency);


#     my @dirs = split (/:/, $ENV{'PATH'});


#   PATH:foreach my $dir (@dirs) {

#       my $path = "$dir/$dependency";

#       #Found and executable
#       if (-x $path) {
# 	  return 0;
#       }

#       #Found but not executable
#       elsif (-e $path) {
# 	  return 102;
#       }
#   }

#     #Dependency not found
#     return 101;
# }




#==========================================================================
#The help for this class

sub help {

  my $self = shift;

  my $help = << 'HELP';

==========================================================================
 The purpose of this class is to verify that all third party programs
 required by a script in order to run properly are installed and
 accessible.  Programs will be searched in:

 a) the invironment table $path via the 'which' command
 b) a list of user supplied directories.

 If a dependency is not properly validated, the dependencies triggering
 the error will be reported and the program will abort.

 The class has the following variables:

 dependencies_list
    Reference to an array of dependencies to be verified.

 user_dirs
    Reference to an array of additional directories to search for
    dependencies. This is necessary when some programs are
    installed in directories not registered in the environment
    variable $PATH.

 allow_no_dependencies
    Flag that if true (1) won't generate an error if the list
    of dependencies is empty. Useful when programs have no
    other dependencies. Default is undef;

 allow_nonexecutables
    Flag that if true (1) won't generate an error if a dependency
    is found but the program is not executable. Default is undef.



 The function that verifies dependencies is: checkDependencies()


 -------


 GENERAL GUIDE FOR USAGE:

 Load class:
 use TCDB ::CheckDependencies;


 #Object Constructor:
 my $myObj = new TCDB::CheckDependencies();



 #Get help about this class:
 my $myObj = new TCDB::CheckDependencies()->help;


 $myObj->help;



 #Initiallizing variables within the constructor:
 my $myObj = new TCDB::CheckDependencies(
    dependencies_list => ['dep1', 'dep2', 'dep3'...],
    user_dirs => ['dir1', 'dir2', 'dir3', ... ],
    allow_no_dependencies => 1,
    allow_nonexecutables => 1);



 #Initializing variables through accessors:
 my $myObj = new TCDB::CheckDependencies();
 $myObj->dependencies_list(['dep1', 'dep2', 'dep3'...]);
 $myObj->user_dirs(['dir1', 'dir2', 'dir3', ... ]);
 $myObj->allow_no_dependencies(1);
 $myObj->allow_nonexecutables(1);



 #Once dependencies are defined, verify that all dependencies exist in
 #in the environment:
 $myObj->checkDependencies;



 -------
 Date created:  11/13/2015
 Programmer:    Arturo Medrano
 Institution:   UCSD, WLU
==========================================================================

HELP

  print $help;
  exit;
}




1;

