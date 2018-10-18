#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

use HTML::Table;
use Getopt::Long;

use TCDB::DrawSequenceObjects;
use TCDB::CheckDependencies;

#==========================================================================
#Input variables

my $hmmtopOutput = undef;
my $mastOutput   = undef;
my $outdir       = "./plots";
my %iMotifs      = ();

read_command_line_arguments();

#print Data::Dumper->Dump([$hmmtopOutput, $outdir, \%iMotifs], [qw(*hmmtopOutput *outdir *iMotifs)]);
#exit;



#==========================================================================
#Parse HMMTOP output in order to get the TMS coordinates of each sequence

my %hmmtop = ();
my %length = ();
parse_hmmtop_output($hmmtopOutput, \%hmmtop, \%length);

#print Data::Dumper->Dump([ \%hmmtop], [qw( *hmmtop )]);
#exit;




#==========================================================================
#Parse all MAST output files, extract the coordinates of the motifs,
#and plot them with the respective TMS.



#Parse MAST and get the positions of motifs
my %mast = ();
my @protIds = ();
parse_mast_output(\%mast, \@protIds);

#print Data::Dumper->Dump([\%mast ], [qw(*mast )]);
#exit;



#Create HTML table with the plots of TMS vs MAST motifs
my $tbl_name = "tms_motifs_plots.html";
generate_html_table($outdir, $tbl_name, \%hmmtop, \%mast, \%length, \@protIds);




###########################################################################
#
#             Functions are declared below this point
#
###########################################################################



sub generate_html_table {

  my ($dir, $tblName, $htms, $hmotifs, $hlen, $gabyIDs) = @_;

  #make directory where the images will be stored
  my $imgDir = "$dir/images";
  mkdir ($imgDir) unless (-d $imgDir);


#  print Data::Dumper->Dump([$hmotifs], [qw(*hmotif )]);
#  exit;


  #Generate the header of the HTML document
  my $table = new HTML::Table;


  foreach my $seqId (@{ $gabyIDs }) {

    my $aLines = $htms->{$seqId};
    my $aRects = $hmotifs->{$seqId};
    my $length = $hlen->{$seqId};

    die "Check tms, motivs and length for: $seqId in $dir" unless ($aLines && $aRects && $length);


    #Generate the image for this sequence
    my $imgName = "${seqId}.png";


    my $plotObj = new R2::DrawSequenceObjects ();
    $plotObj->seqLength($length);
    $plotObj->outImgFile("$imgDir/$imgName");
    $plotObj->outFormat('png');
    $plotObj->lines($aLines);
    $plotObj->rectangles($aRects);

    $plotObj->plot_lines_vs_rectangles();


    #Generate the link to the image
    my $imgLink = "<img src='./images/$imgName'>";

    $table->addRow($seqId, $imgLink);
  }


  #Generate HTML output
  open (my $outh, '>', "$dir/$tblName") || die $!;
  print $outh "<html>\n<body>\n";
  print $outh $table->getTable;
  print $outh "\n</body>\n</html>\n";
  close $outh;

}







#==========================================================================
#Parse mast ouput


sub parse_mast_output {

  my ($hout, $aIDs) = @_;


  #Determine whether just a few or all motifs will be plotted
  my $allMotifs = (%iMotifs)?  0 : 1;


  open(my $masth, "<", $mastOutput) || die $!;


  my %motifs = ();
  while (<$masth>) {

    chomp;

    #Trim spaces at the beginning and end of lines
    s/^\s+//;
    s/\s+$//;


    #Ignore empty lines;
    next unless ($_);


    #Read the length of each motif found
    if (/^#/) {
      my @motifLength = split (/\s+/, $_);
      $motifs{$motifLength[1]} = $motifLength[2];
      next;
    }


    #Read the MAST motifs and their locations in the sequence
    my ($id, $evalue, $mPositions) = split(/\s+/, $_);
    die "No ID, evalue or positions: $_" unless ($id && $evalue && $mPositions);


    push(@{ $aIDs }, $id);

    #Translate motif positions to sequence coordinates
    my @posData = split(/\-/, $mPositions);
    my %tmp = ();
    my $left = 0;
    for (my $i = 0; $i < $#posData; $i++) {

      #Determine wether datum is an offset or a site id.
      my $datum = $posData[$i];
      if ($datum =~ /^[0-9]+$/) {
	$left = $posData[$i] + $left;
      }

      #Datum is a site motif id
      elsif ($posData[$i] =~ /^\[(\d+)\]$/) {

	my $motivo = $1;


	#Get the length of the motif
	my $length = $motifs{$motivo};
	die "No length for motif: $motivo" unless ($length);


	my $start = $left + 1;
	my $end   = $start + $length - 1;
	$left = $end;


	#Only add this motif if it will be plotted
	if (exists $iMotifs{ $motivo } ||  $allMotifs) {
	  push (@{ $hout->{$id} }, [$motivo, $start, $end]);
	}

	#	print "$mPositions\n";
#	print Data::Dumper->Dump([ $out ], [qw(*motifs )]);
#	<STDIN>;
      }
    }
  }
  close $masth;
}



#==========================================================================
#Parse hmmtop and extract the sequence ID, length, and TMS positions

sub parse_hmmtop_output {
  my ($hmmtopPath, $hmmtopHash, $lengthHash) = @_;

  open (my $tmsh, "<", $hmmtopPath) || die $!;
  while (<$tmsh>) {
    chomp;
    my ($header, $orientation, $tmsData) = split (/\s+(IN|OUT)\s+/, $_);
    die "No header or TMS info found: $_" unless ($header && $tmsData);


    #Extract the length and the ID of the sequence
    my ($kk1, $length, $id, @kk) = split(/\s+/, $header);
    die "No seq ID or length: $header" unless ($length && $id);


    #Store sequence length
    $lengthHash->{$id} = $length;


    #Extract the TMS data
    my @positions = split (/\s+/, $tmsData);
    my @coords = ();
    for (my $i = 1; $i < $#positions; $i += 2) {
      push (@coords, [$positions[$i], $positions[$i + 1]]);
    }


    #Store the TMS positions
    $hmmtopHash->{$id} = \@coords;


#    print "$id\t($length)\t$tmsData\n";
#    print Data::Dumper->Dump([\@coords ], [qw( *tmsCoords)]);
#    <STDIN>;

  }
  close $tmsh;
}





#==========================================================================
#Read command line arguments

sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  my $tmpMotifs = [];

  my $status = GetOptions(
      "t|hmmtop-output=s"     => \$hmmtopOutput,
      "m|mast-output=s"       => \$mastOutput,
      "o|outdir=s"            => \$outdir,
      "p|plot-motifs=s@"      => \$tmpMotifs,
      "h|help"                => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);


  #HMMTOP output file must exist and not be empty
  die "File with HMMTOP output must exist and not be empty: $hmmtopOutput" unless (-f $hmmtopOutput && !( -z  $hmmtopOutput));


  #MAST output files must exist and not be empty
  die "File with MAST output must exist and not be empty: $mastOutput" unless (-f $mastOutput && !( -z $mastOutput));


  #Validate motif lists and remove potentially duplicated motif IDs. This are the motifs that will
  #be included in the plot
  if (@{ $tmpMotifs }) {

    @{ $tmpMotifs } = split (/,+/, join (",", @{ $tmpMotifs }));

    foreach my $mID (@{ $tmpMotifs }) {
      die "-m option requires comma-separated integer values (MEME motif identifiers): $mID" unless ($mID =~ /[0-9]+/);
      $iMotifs{$mID} = 1;
    }
  }

  #Output directory must exist or should be created
  system "mkdir -p $outdir" unless (-d $outdir);
}





#==========================================================================
#Print the help to run this script

sub print_help {

  my $help = <<'HELP';

This script reads the preprocessed output of MAST and plots for each
sequence the TMS and the motifs found by MEME. An HTML file will
be generated showing the plots.

-t, --hmmtop-output={string}
   File with the results of running hmmtop on the set of sequences
   analyzed by MEME.

-m, --mast-output={string}
   File with the results of running mast on the set of sequences
   analyzed by MEME. This files is postprocessed and only has the
   lines that indicate the length of each motif, and the lines
   with the coordintates where each motif is located in each
   sequence.


-o, --outdir={string}
    Directory where the plots will be located.

-p, --plot-motifs {comma-separated integers}
   List of comma separated motif identifiers to plot. This allows to select which of
   the MEME motifs to include in the plot. By default all motifs are
   included.

-h, --help
    Print this help.

HELP

  print $help;
  exit;
}
