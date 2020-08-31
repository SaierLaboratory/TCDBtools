#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

use Getopt::Long;
use Bio::SeqIO;

use TCDB::CheckDependencies;



#==========================================================================
#Variables globales

#Archivo con sequencies de la familia
my $fam_seq_file = undef;


#Archivo con hmmtop results
my $hmmtop_file  = undef;


#Output directory where results should be placed
my $outdir = "./OrganizedSequences";



#==========================================================================
#Check the dependencies for this program

my @dependencies = ("mkdir", "hmmtop");
my $CheckDep_obj = new TCDB::CheckDependencies();
$CheckDep_obj -> dependencies_list(\@dependencies);
$CheckDep_obj -> checkDependencies;





#==========================================================================
#Read Command line arguments

read_command_line_arguments();

#print Data::Dumper->Dump([ $hmmtop_file, $fam_seq_file, $outdir], [qw(*hmmtop_file *fam_seq_file *outdir)]);
#exit;



#==========================================================================
#run hmmtop on the input sequences

system "hmmtop -if=$fam_seq_file -of=$hmmtop_file -is=pseudo -pi=spred -sf=FAS";
die "Problem found when running hmmtop" unless (-f $hmmtop_file && !( -z $hmmtop_file));



#==========================================================================
#Extraer resultados de hmmtop

my %hmmtop = ();
parse_hmmtop($hmmtop_file, \%hmmtop);
#print Data::Dumper->Dump([\%hmmtop], ["*hmmtop"]);
#exit;


#==========================================================================
#Extraer la sequencias fasta

my %seqs = ();
read_fasta_sequences($fam_seq_file, \%seqs);
#print Data::Dumper->Dump([\%seqs], ["*seqs"]);



#Averiguar cuantos diferentes TMS hay en los resultados
#my %nTMS = map { $_ => 1 } values %hmmtop;
#print Data::Dumper->Dump([\%nTMS], ["*nTMS"]);




#==========================================================================
#Clasificar las secuencias de acuerdo al numero de TMS

foreach my $id (keys %hmmtop) {

  my $tms  = $hmmtop{$id};

  print "Classifying sequence: $id has $tms TMS\n";

  my $file = "$outdir/tms_${tms}.faa";


  open (my $out, ">>", $file) || die $!;
  print $out $seqs{$id}, "\n";
  close $out;
}






#==========================================================================
#Read all the fasta sequences from input fasta file

sub read_fasta_sequences {

  my ($infile, $seq_hash) = @_;

  my $seqin  = Bio::SeqIO->new(-file => $infile , -format => "fasta");


  while(my $obj = $seqin->next_seq) {

    my $seq =  ">" . $obj->primary_id . " " . $obj->desc . "\n" . $obj->seq . "\n";

    $seq_hash->{ $obj->primary_id} = $seq;
  }
}





#==========================================================================
#Parse the output of HMMTOP

sub parse_hmmtop {

  my ($infile, $tms_hash) = @_;


  open(my $file, $infile) || die $!;
  while (<$file>) {
    chomp;

    if (/^\>HP:\s+\d+\s+(\S+).+(IN|OUT)\s+(\d+)/) {
      $tms_hash->{$1} = $3;
    }
  }
  close $file;

}




#==========================================================================
#Read command line arguments


sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "s|seqs-file=s"         => \$fam_seq_file,
      "o|outdir=s"            => \$outdir,
      "h|help"                => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);


  #sequence file should exist and not be empty
  die "File with input sequences must exist and not be empty: $fam_seq_file" unless (-f $fam_seq_file || !( -z  $fam_seq_file));


  #Output directory must exist or should be created
  system "mkdir -p $outdir" unless (-d $outdir);


  #The name of the hmmtop output file
  $hmmtop_file = "$outdir/hmmtop.out";
}




#==========================================================================
#Print the help to run this script

sub print_help {

  my $help = <<'HELP';

This script takes a file with multiple protein sequences in FASTA fromat
and generates multiple files were each files has proteins with a given
number of HMMTOP-predictes TMSs.

-s, --seqs-file={string}
    Fasta file with the proteins sequences that will be organized.
    Argument is mandatory.

-o, --outdir={string}
    Directory where the output files will be located. If the directory
    does not exist, it will be created.

-h, --help
    Print this help.

HELP

  print $help;
  exit;
}

