#!/usr/bin/perl

## This program should be able to:
##
## 1. Blastp the genome sequences
## 2. Find reciprocal best hits
## 3. Find ortholog-higher-than-paralog
## 4. Find unidirectional best hits
## 5. Find paralogs (remaining hits after orthologs)

use strict;
use Getopt::Long;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

### checking options
my $cpuNumber
    = qx(sysctl -a | grep 'cpu.thread_count')
        =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
    : qx(sysctl -a 2>/dev/null | grep 'max-threads')
        =~ m{\.max-threads\s+=\s+(\d+)} ? $1
    : 2;

my $inputSeqDir = "";
my (
    $outputFolder, $maxSubject, $Evalue,
    $minCoverage, $cpus
) = ( 'RBH', 10000, "1e-7",
      0.6, $cpuNumber);
my $options = GetOptions(
    "i=s" => \$inputSeqDir,
    "o=s" => \$outputFolder,
    "n=f" => \$maxSubject,
    "e=s" => \$Evalue,
    "c=f" => \$minCoverage,
    "a=s" => \$cpus,
);

my $ownName = $0;
$ownName =~ s{.*/}{};
if ( $inputSeqDir eq "" ) {
    print "usage: " . $ownName . " [options]\n";
    print "\noptions:\n";
    print "   -i sequence files in fasta format, required\n";
    print "   -o output folder, default RBH\n";
    print "   -n max number of aligned sequences to keep, default 10000\n";
    print "   -e evalue threshold, default 1e-7\n";
    print "   -c minimum alignment coverage of aligned sequence,\n"
        . "       default 0.6\n";
    print "   -a number of cpus to use, default in this machine: $cpuNumber\n";
    print "\n";
    exit;
}

my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
### check that the files exist and have sequences (?)
my $checked = checkFastas($inputSeqDir);

### program here



### sub routines

sub checkFastas {
    my $inputSeqDir = @_;
}

sub signalHandler {
    if( -d "$tempFolder" ) {
        print "\n\tcleaning up ...\n";
        system "rm -r $tempFolder";
        die  "    done!\n\n";
    }
    else {
        print "\n\ttemp files cleared out\n";
        die  "    done!\n\n";
    }
}
