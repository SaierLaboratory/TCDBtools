#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;

$Data::Dumper::Deepcopy = 1;

#To read command line options
use Getopt::Long;

#to check dependencies
use TCDB::CheckDependencies;




#==========================================================================
#Check dependencies

my @dependencies = ("sort", "grep", "hvordan.py", "quod.py");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;




#==========================================================================
#Read command line

my $fxpand_dir     = "";
my $proto2_root    = "";
my $proto2_dir     = "";
my $tscores_dir    = "";  #Directory with the results of running list_top_gsat_hits
my $outdir         = "";  #Where the output of the program will be saved
my $nTopHits       = 15;
my $p2_low_cutoff  = 15;
my $p2_high_cutoff = 10000;


read_command_line();

#print Data::Dumper->Dump([$fxpand_dir, $proto2_root, $proto2_dir, $outdir, $nTopHits, $p2_low_cutoff, $p2_high_cutoff],
#                         [qw(*fxpand_dir *proto2_root *proto2_dir *outidr *nTopHits *p2_low_cutoff *p2_high_cutoff)]);
#exit;




#==========================================================================
#Get all the output files from script list_top_gsat_hits

opendir(my $dirh, $tscores_dir);
my @scoreFiles = grep { /^top_gsat/ && -f "$tscores_dir/$_" } readdir($dirh);
closedir $dirh;

#print Data::Dumper->Dump([\@scoreFiles ], [qw(*scoreFiles )]);
#exit;



#==========================================================================
#For top scores file get the genes B and C of the top hits


foreach my $sFile (@scoreFiles) {

  #Read hits file
  open (my $tmph, "<", "$tscores_dir/$sFile") || die $!;
  chomp (my @content = <$tmph>);
  close $tmph;


  #Extract the BC hits
  my @bc = grep { /^HIT/ } @content;
  next unless (@bc);


  #Get the families being compared from the file name
  my $fam1 = "";
  my $fam2 = "";
  if ($sFile =~ /top_gsat_([A-Z0-9\.]+)_vs_([A-Z0-9\.]+)\.scores/) {
    $fam1 = $1;
    $fam2 = $2;
  }
  die "Could not get the family names from file: $sFile" unless ($fam1 && $fam2);

  my $famDir = "$outdir/${fam1}_vs_${fam2}";
  system "mkdir $famDir" unless (-d $famDir);


  #For the top n hits run hvordan.py
  my $cnt = 1;
 HIT:foreach my $hit (@bc) {

    #Exit loop if the right number of hits was already processed.
    last HIT if ($cnt > $nTopHits);

    #Get accessions and run hvordan.py
    if ($hit =~ /HIT B-C: (\S+) vs (\S+)\s+\(([0-9\.-]+)\)/) {
      my $pid1  = $1;
      my $pid2  = $2;
      my $score = $3;

      my $pairDir    = "${pid1}_vs_${pid2}_h${cnt}_s${score}";
      my $pairProts  = "${pid1}_vs_${pid2}";
      my $htmlFile   = "$famDir/$pairDir/html/${pairProts}.html";

      #If this pair of proteins was already analyzed, go to next hit.
      if (-f $htmlFile) {
	print "Output file exists: $htmlFile\n";
	$cnt++;
	next HIT;
      }


      #Sleep for a couple of seconds just to avoid hammering the NCBI servers
      sleep (2);


      #Run hvordan.py
      my $args = "--p1d $fxpand_dir --p2d $proto2_root --outdir $famDir/$pairDir --fams $fam1 $fam2 -p $pid1 $pid2 -z $p2_low_cutoff";
      print "\n\nhit: $cnt\nhvordan.py $args\n\n";
      system "hvordan.py $args";
      $cnt++;

    }
  }
}




###########################################################################
################   Subroutines definition beyond ths point   ##############
###########################################################################


#===========================================================================
#Read command line and print help


sub read_command_line {

  print_help() unless (@ARGV);

  my $status = GetOptions(
	"fxd|fxpand-dir=s"           => \&read_fxpand_dir,
	"p2r|proto2-root-dir=s"      => \&read_proto2root_dir,
	"p2d|proto2-dir=s"           => \&read_proto2_dir,
	"tsd|top-scores-dir=s"       => \&read_top_scores_dir,
	"o|outdir=s"                 => \&read_output_dir,
	"n|top-hits=i"               => \$nTopHits,
	"p2lt|proto2-low-cutoff=f"   => \$p2_low_cutoff,
	"p2ht|proto2-high-cutoff=f"  => \$p2_high_cutoff,
	"h|help"                 => sub { print_help(); },
	"<>"                     => sub { die "Error: Unknown argument: $_[0]\n"; });
  exit unless ($status);


  if ($proto2_root && $proto2_dir) {
    die "Error: options -p2r and -p2d are incompatible.\n";
  }

  unless ($nTopHits || $p2_low_cutoff || $p2_high_cutoff) {
    die "Error: one option -n -p2lt -p2ht must give it";
  }

  unless ($outdir) {
    die "Error: output directory is mandatory";
  }

  unless ($tscores_dir) {
    if ($proto2_root) {
      $tscores_dir = $proto2_root;
    }
  }
}



#==========================================================================
#Read directory where the results of this programs will be stored

sub read_output_dir {
  my ($opt, $value) = @_;

  if ($value) {
    $outdir = $value;
    system "mkdir -p $outdir" unless (-d $outdir);
  }
  else {
    die "Error: directory passed to option --$opt must exist.\n";
  }
}





#==========================================================================
# Read directory where the results of running list_top_gsat_hits

sub read_top_scores_dir {
  my ($opt, $value) = @_;

  if ($value && -d $value) {
    $tscores_dir = $value;
  }
  else {
    die "Error: directory passed to option --$opt must exist.\n";
  }
}



#==========================================================================
# Read root diretory for all protocol2 directories in the analysis

sub read_proto2_dir {
  my ($opt, $value) = @_;

  if ($value && -d $value) {
    $proto2_root = $value;
  }
  else {
    die "Error: directory passed to option --$opt must exist.\n";
  }
}



#==========================================================================
# Read root diretory for all protocol2 directories in the analysis

sub read_proto2root_dir {
  my ($opt, $value) = @_;

  if ($value && -d $value) {
    $proto2_root = $value;
  }
  else {
    die "Error: directory passed to option --$opt must exist.\n";
  }
}



#==========================================================================
# Read famXpander directory

sub read_fxpand_dir {
  my ($opt, $value) = @_;

  if ($value && -d $value) {
    $fxpand_dir = $value;
  }
  else {
    die "Error: directory passed to option --$opt must exist.\n";
  }
}



#==========================================================================
#The help of the program

sub print_help {

    my $help = <<'HELP';

After running scripts areFamiliesHomologous and list_top_gsat_hits, use this script 
to generate plots that will help you make the biological interpretation of the results.


Input parameters:

-fxd, --fxpand-dir  { path }  (Mandatory)
  Directory with the results of runing famXpander on the families
  under analysis (mandatory).

-p2r, --proto2-root-dir { path }  (Optional but see -p2d)
  This is the main directory where the the results of running protocol2
  for multiple pairs of families are found. This option can't be used
  in combination with option -p2d.

-p2d, --proto2-dir {path}  (Optional but see -p2r)
  The directory with the results of running  protocol2 for a couple
  of families. This is the directory where file report.html is located.
  This option cannot  be used in combination with option -p2r

-tsd, --top-scores-dir  { path }  (Optional)
  Directory where the results of running list_top_gsat_hits are stored.
  By default it takes the value passed to option -p2r.

-n, --num-top-hits { int }  (Optional; Default: 15)
  Number of significant top hits for which hvordan.py will be run.

-p2lt, --proto2-low-cutoff {float}
  Minimum threshold to identify significant hits from Protocol_2 and GSAT
  (default is 14.0)

-p2ht, --proto2-high-cutoff {float}
  Maximum threshold to identify significant hits from Protocol_2 and GSAT.
  Use this option along with -p2lt to focus the analysis on a range of scores.
  (default is unlimited)

-o, --outdir { path }  (Mandatory)
  Directory with the hydropathy and blast results will be saved.

-h, --help
   Display this help. Also displayed if script is run without arguments.

HELP

    print $help;
    exit;
}
