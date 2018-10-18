#!/usr/bin/env perl -w

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;



my $blast_output_file     = "";
my $blast_output_format   = "csv";
my $system_sequences_file = "";
my $outfile               = "";


#==========================================================================
#read comand line arguments

read_command_line();


#print Data::Dumper->Dump([$blast_output_file, $blast_output_format, $system_sequences_file, $outfile],
#			 [qw(*blast_output_file *blast_output_format *system_sequences_file *outfile)]);
#exit;


#==========================================================================
#Extract the IDs of all the components in this system


#Get the family for which sequences are being extracted
my $family = ( $system_sequences_file =~ /family-([a_zA-Z0-9\.]+).faa/ )? $1 : undef;
die "Could not extract family from: $system_sequences_file" unless ($family);



my @components = ();
open (my $fh1, "<", $system_sequences_file) || die $!;

while (<$fh1>) {
  chomp;
  if (/^\>(\S+)/) {
    push (@components, $1);
  }
}
close $fh1;

#print Data::Dumper->Dump([\@components], [qw(*components )]);
#exit;




#==========================================================================
#Search each component in the blast results and select all components
#without blast matches


#Read all lines in blasst output
open (my $fh2, "<", $blast_output_file) || die $!;
chomp (my @hits = <$fh2>);
close $fh2;


my @notFound = ();

foreach my $cmp (@components) {

  my @found = grep {/$cmp/} @hits;

  unless (@found) {
    push (@notFound, $cmp);
  }
}

#print Data::Dumper->Dump([\@notFound ], [qw(*notFound )]);
#exit;


#==========================================================================
#Save output to file

open (my $fh3, ">", $outfile) || die $!;

my $ratio = scalar(@notFound) / scalar (@components) * 100;

print $fh3 "#missing components from $family: ", scalar(@notFound), "/", scalar (@components), " (${ratio}%)\n";
print $fh3 join ("\n", @notFound), "\n";

close $fh3;








###########################################################################
####          Functions are defined below this point                   ####
###########################################################################


#===========================================================================
#Read command line and print help


sub read_command_line {

    print_help() unless (@ARGV);

    my $status = GetOptions(
	"b|blast-file=s"         => \&read_blast_file,
	"f|blast-format=s"       => \&read_blast_format,
	"s|sys-seqs-file=s"      => \&read_seqs_file,
	"o|outfile=s"            => \$outfile,
	"h|help"                 => sub { print_help(); },
	"<>"                     => sub { die "Error: Unknown argument: $_[0]\n"; });
    exit unless ($status);
}



sub read_blast_file {
    my ($opt, $value) = @_;

    unless (-f $value) {
      die "Error: blast file does not exist -> $value";
    }

    $blast_output_file = $value;
}




sub read_blast_format {
    my ($opt, $value) = @_;

    unless ($value =~/(csv|tsv)/) {
      die "Error: unknown format for blast file --> $value"
    }

    $blast_output_format = $value;
}



sub read_seqs_file {
    my ($opt, $value) = @_;

    unless (-f $value) {
      die "Error: file with sequences not found --> $value"
    }

    $system_sequences_file = $value;
}




sub print_help {

    my $help = <<'HELP';

When running blastp of a multicomponent system in TCDB against a query genome, this script extracts
all components that were not found by blastp.

-b, --blast-file {path}
   File with the ouput of blastp

-f, --blast-format {string}
   Format of the blast output file, the accepted arguments are:
     csv  Comma-delimited file (default)
     tsv  Tab-delimited file.

-s, --sys-seqs-file {path}
   File with the sequences of all components in a specific transport system
   in fasta format.

-o, --outifle {path}
   File where the results will be saved.

-h, --help
   Display this help. Also displayed if script is run without arguments.

HELP

    print $help;
    exit;
}
