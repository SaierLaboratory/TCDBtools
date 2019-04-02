#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

$Data::Dumper::Deepcopy = 1;

use Getopt::Long;
use Bio::SearchIO;
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

my @dependencies = ("glsearch36", "blastdbcmd");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;




#==========================================================================
#Read command line arguments

my $fragment    = undef;
my $accession   = undef;
my $accFile     = undef
my $outdir      = undef;
my $blastdb     = undef;
my $evalue      = 1e-2;
my $subMatrix   = 'BL50';
my $quiet       = 0;
my $interactive = 0;

read_command_line();

#print Data::Dumper->Dump([$fragment, $accession, $accFile, $outdir, $blastdb, $evalue, $quiet],
#                         [qw(*fragment *accession *accFile *outdir *blastdb *evalue *quiet)]);
#exit;




#==========================================================================
#Get the sequences for the fragment and full proteins in files


my ($fragFile, $protFile) = @{ getSequences($fragment, $accession) };
#print Data::Dumper->Dump([$fragFile, $protFile ], [qw(*fragFile *protFile)]);

die "Both files must exist:\n  $fragFile\n  $protFile" unless (-f $fragFile && -f  $protFile);

#==========================================================================
#align fragment to full protein

my $alignFile = "$outdir/${accession}_glsearch.out";
my $cmd = qq(glsearch36 -s BL62 -z 21 -k 10000 -E $evalue -s $subMatrix $fragFile $protFile > $alignFile);
system $cmd unless (-f  $alignFile);


#==========================================================================
#Parse glsearch output and run quod


run_quod ($alignFile, $protFile);





#==========================================================================
################   Subroutines definition beyond ths point   ##############
#==========================================================================


sub run_quod {

  my ($alignment, $sequence) = @_;


  my $parser = new Bio::SearchIO (-format => 'fasta', -file => $alignment);

  my @res = ();

  #----------------------------------------------------------------------
  #Parse glsearch output

  my $hsp_cnt = 1;
  while (my $result = $parser->next_result) {

  HIT:while (my $hit = $result->next_hit) {
    HSP:while(my $hsp = $hit->next_hsp) {

	my %data = ();
	my $key  = "hsp_$hsp_cnt";
	my $hId  = $hsp->frac_identical('total');  #identity in the alignment
	my $hSim = $hsp->frac_conserved('total');  #similarity in the alignment

	#Alignment parameters
	$data{hsp}    = $key;
	$data{evalue} = $hsp->evalue;
	$data{id}     = $hId;
	$data{sim}    = $hSim;

	#coordinates in the alignment to plot bars
	$data{qstart} = $hsp->start('query');
	$data{qend}   = $hsp->end('query');
	$data{sstart} = $hsp->start('subject');
	$data{send}   = $hsp->end('subject');

	push (@res, \%data);
	$hsp_cnt++;
      }
    }
  }

  die "Error: no match found between fragment and sequence" unless (@res);
#  print Data::Dumper->Dump([\@res ], [qw(*res )]);


  #----------------------------------------------------------------------
  #Generate quod plot

  #Format string for the regions
  my $regions = "-at ";
  my $coords  = "";
  foreach my $hit (@res) {

    $coords = $hit->{sstart} . "-" . $hit->{send};
    $regions .= "${coords}:green";

    #only the best HSP is required to be plotted
    last;
  }

  my $qstring = ($quiet)? "-q" : "";
  my $iString = ($interactive)? "--show" : "";

  my $cmd = qq(quod.py $qstring $iString -l "$accession ($coords)" --xticks 25 --grid  $regions -- $sequence);
#  print "$cmd\n";
  system $cmd;
}



sub getSequences {

  my ($frag, $acc) = @_;

  my $accSeq = "$outdir/${acc}.faa";


  #Save fragment to file
  my $tmpFile = "$outdir/${acc}_frag.faa";
  unless (-f $tmpFile) {
    open (my $fh, ">", $tmpFile) || die $!;
    print $fh ">${acc} fragment\n$frag\n";
    close $fh;
  }

  unless (-f $accSeq) {
    #Blast DB to be used
    my $db = ($blastdb)? $blastdb : 'nr';

    my $cmd = qq(blastdbcmd -db $db -entry $acc -target_only -outfmt '\%f' -out $accSeq);
    system $cmd;

    #Remove the version and annotations from the sequence file
    my $cmd2 = qq(perl -i.bkp -pe 's/^\\>(\\w+)\.\*/\\>\$1/;' $accSeq);
    system $cmd2 unless (-f "${accSeq}.pkp");
  }

  #Return files to be aligned
  return [$tmpFile, $accSeq];
}




#===========================================================================
#Read command line and print help


sub read_command_line {

  print_help() unless (@ARGV);

  my $status = GetOptions(
	"i|acc-file=s"     => \&read_accFile,
	"a|accession=s"    => \&read_accession,
	"f|fragment=s"     => \&read_fragment,
	"o|outdir=s"       => \$outdir,
	"bdb|blastdb=s"    => \&read_blastdb,
	"e|evalue=f"       => \$evalue,
	"m|sub-matrix=s"   => \&read_subMatrix,
	"t|interactive!"   => \$interactive,
	"q|quiet!"         => \$quiet,
	"h|help"           => sub { print_help(); },
	"<>"               => sub { die "Error: Unknown argument: $_[0]\n"; });
  exit unless ($status);


  die "Error: option -f is mandatory." unless ($fragment);
  die "Error: options -i or -a are mandatory." unless ($accession || $accFile);
  die "Error: flags -t and -q cannot be set at the same time!" if ($quiet && $interactive);

  #Default value for output directory
  $outdir = "." unless ($outdir);
  system "mkdir -p $outdir" unless (-d $outdir);
}


#==========================================================================
#Option -i

sub read_accFile {
  my ($opt, $value) = @_;

  unless (-f $value && !(-z $value)) {
    die "Error in option -i: File with sequences does not exist or is empty!\n";
  }

  $accFile = $value;
}


#==========================================================================
#Option -a

sub read_accession {
  my ($opt, $value) = @_;

  #Remove version number if any
  $value =~ s/\.\d+$//;

  $accession = $value;
}


#==========================================================================
#Option -f

sub read_fragment {
  my ($opt, $value) = @_;

  #Remove dashes, dots and any other symbol that may represent gaps
  $value =~ s/[\.\s-]+//g;

  $fragment = $value;
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
#option -h


sub print_help {

    my $help = <<'HELP';

Analyze any given GBLAST match. Plot hydropathy and inferred domains.

Options

-a, --accession {String} (Mandatory)
   NCBI accession of query protein

-f, --fragment {String} (Mandatory)
   Sequence fragment to locate within the full protein

-o, --outdir {PATH} (Optional. Default: ./)
   Path to output directory.

-bdb {PATH} (Optional. Default: nr)
   Path to the BLAST database where accessions will be extracted from.

-e, --evalue {FLOAT} (Optional. Default: 0.01)
   E-value cut off when comparing full proteins

-m, --sub-matrix {STRIN} (Optional. Default: BL50)
   Amino acid substitution matrix to use in comparisons.

-t, --interactive [flag] (Optional)
   Generate the plot but open it in interactive mode,
   The user will be able to identify specific residue
   positions. This option is incompatible with -q.

-q, --quiet [flag] (Optional)
    Generate plots but do not attempt to open them with
    an image viewer. Incompatible with option -i.

-h, --help
   Print this help.

HELP

    print $help;
    exit;
}
