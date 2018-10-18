#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

$Data::Dumper::Deepcopy = 1;

use Getopt::Long;
use LWP;
use Bio::SeqIO;

#Local libraries
use TCDB::CheckDependencies;
use TCDB::Assorted;


###########################################################################
#
# Parse GBLAST output file and return the sequence of the proteins that
# meet the user's criteria.
#
###########################################################################

#==========================================================================
#Check dependencies

my @dependencies = ("grep", "gblast3.py", "hmmtop", "blastdbcmd");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;




#==========================================================================
#Read command line arguments

my $infile   = undef;
my $outdir   = undef;
my $blastdb  = undef;
my $proteins = undef;
my $nTMS     = undef;
my $evalue   = 1e-10;
my $coverage = 70;
my $maxSeqs  = 400;
my $covControl = "X";

read_command_line();

#print Data::Dumper->Dump([$infile, $outdir, $blastdb, $nTMS,
#                          $evalue, $coverage],
#                         [qw(*infile *outfile *blastdb *nTMS
#                          *evalue *coverage)]);
#exit;



#==========================================================================
#Get the accessions of the proteins that meet the user's parameters

my %proteins = ();
parse_gblast(\%proteins);

#print Data::Dumper->Dump([\%proteins ], [qw(*proteins )]);
#exit;


#==========================================================================
#Get the sequences of the proteins


getSequences(\%proteins);



#==========================================================================
################   Subroutines definition beyond ths point   ##############
#==========================================================================


sub getSequences {

  my $prots = shift;

  my $outfile = "";
  if ($nTMS) {
    $outfile = "$outdir/gblast_seqs_with_${nTMS}tms.faa";
  }
  else {
    $outfile = "$outdir/gblast_seqs.faa";
  }


  #Save accessions to file
  my $tmpFile = "$outdir/accs.txt";
  open (my $fh, ">", $tmpFile) || die $!;
  print $fh join ("\n", keys %$prots), "\n";
  close $fh;


  #Blast DB provided by user
  my $cmd = "";
  if ($blastdb) {

    #Extract sequences
    $cmd = qq(blastdbcmd -db $blastdb -entry_batch $tmpFile -target_only -outfmt '\%f' -out $outfile);
    system $cmd unless (-f $outfile);
  }


  #Extract sequences from NCBI
  else {
    $cmd = qq(blastdbcmd -db nr -entry_batch $tmpFile -target_only -outfmt '\%f' -out $outfile);
    system $cmd unless (-f $outfile);
  }


  #Remove the version and everything else from the sequence file
  my $cmd2 = qq(perl -i.bkp -pe 's/^\\>(\\w+)\.\*\$/\\>\$1/;' $outfile);
#  print "$cmd2\n";
#  exit;
  system ($cmd2) unless (-f "${outfile}.pkp");
}





#==========================================================================
#Extract the accessions of the proteins that meet user's parameters

sub parse_gblast {

  my $out = shift;

  #format query proteins for search... if available
  my $protStr = ($proteins)? join ("|", split(/,/, $proteins)) : undef;


  #To count the number of hits being added to the final results.
  my $hitCnt = 0;

  open (my $fh, "<", $infile) || die $!;
 HIT:while (<$fh>) {

    chomp;
    next HIT if (/^#/);

    #Identify query proteins if given
    next HIT if ($proteins && !( /$protStr/ ));


    #Columns in GBLAST file
    #0.  Query Accession
    #1.  UniProt Id of TCDB hit
    #2.  tcid of hit
    #3.  Description of TC hit
    #4.  Alignment length
    #5.  E-value
    #6.  Identity %
    #7.  Length of query
    #8.  Length of hit.
    #9.  Coverage of query
    #10. Coverage of hit
    #11. TMS in query
    #12. TMS in subject
    #13. TMS overlap.
    #14. family abbreviation of hit
    #15. Substrate
    #16. Row number.
    my @data = split(/\t/);

    my $tmsOk = ($nTMS)? (($data[11] == $nTMS)? 1 : 0) : 1;

#    print Data::Dumper->Dump([$hitCnt, $data[5], $data[9], $data[10], $coverage, $covControl, $tmsOk],
#			     [qw(*row *evalue *qcov *scov *covCutoff *covControl *tmsOk)]);
#   <STDIN>;


    if ($tmsOk && ($data[5] <= $evalue) && TCDB::Assorted::coverage_ok($data[9], $data[10], $coverage, $covControl)) {

      last HIT if ($hitCnt >= $maxSeqs);

      my ($pid, $v) = split (/\./, $data[0]);

      $out->{$pid} = {hit=>[$data[1], $data[2]], eval=>$data[5], id=>$data[6],
		      qcov=>$data[9], scov=>$data[10]};
      $hitCnt++;
    }

  }
  close $fh;

}






#===========================================================================
#Read command line and print help


sub read_command_line {

  print_help() unless (@ARGV);

  my $status = GetOptions(
	"i|infile=s"       => \&read_infile,
	"o|outdir=s"       => \$outdir,
	"bdb|blastdb=s"    => \&read_blastdb,
        "p|proteins=s"     => \$proteins,
	"e|evalue=f"       => \$evalue,
	"c|coverage=f"     => \$coverage,
        "cc|cov-control=s" => \&read_covControl,
	"n|tms=i"          => \$nTMS,
	"m|max-seqs=i"     => \$maxSeqs,
	"h|help"           => sub { print_help(); },
	"<>"               => sub { die "Error: Unknown argument: $_[0]\n"; });
  exit unless ($status);


  #Check for incompatibilities and errors
  unless (-f $infile && !(-z $infile)) {
    die "Error: gblast output file must exist and not be empty! -> $infile";
  }

  #Default value for output directory
  $outdir = "./gblast_filtered_seqs" unless ($outdir);
  system "mkdir -p $outdir" unless (-d $outdir);
}


#==========================================================================
#Option -i

sub read_infile {
  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
    die "Error in option -i: File with sequences does not exist or is empty!\n";
  }

  $infile = $value;
}


#==========================================================================
#Option -bdb

sub read_blastdb {
  my ($opt, $value) = @_;

  my $tmpFile = "${value}.pin";

  unless (-f $tmpFile && !(-z $tmpFile)) {
    die "Error in option -bdb: Blast DB does not exist! -> $value";
  }

  $blastdb = $value;

}

#==========================================================================
#Option -cc

sub read_covControl {
  my ($opt, $value) = @_;

  my $tmp = uc $value;
  unless ($tmp =~ /^[XQSB]$/) {
    die "Error in option -cc: illegal charater ($value). Valid characters are Q,S,B,X\n";
  }

  $covControl = $tmp;
}




#==========================================================================
#option -h


sub print_help {

    my $help = <<'HELP';

Tyee the help here.

HELP

    print $help;
    exit;
}
