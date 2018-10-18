#!/usr/bin/perl
#########################################################################
#									#
#	Author: Gabo Moreno-Hagelsieb           			#
#	Date first raw version: Aug 4, 2017    				#
#									#
#########################################################################

use strict;
use Getopt::Long;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

### ensure that external programs exist
my @missingsoft = ();
for my $xsoftware ( qw( ssearch36 ggsearch36 segmasker shuffleseq ) ) {
    if( my $sfwpath = qx(which $xsoftware) ) {
        chomp($sfwpath);
        #print "$sfwpath<-here\n";
    }
    else {
        print "\tprogram $xsoftware not found\n";
        push(@missingsoft,$xsoftware);
    }
}
my $cntMissing = @missingsoft;
if( $cntMissing > 0 ) {
    die "\tcan't proceed because of missing software\n\t"
        . join("\n\t",@missingsoft) . "\n\n";
}

### check if there's more than one processor or assume there's 2.
my $cpuNumber
    = qx(sysctl -a | grep 'cpu.thread_count')
        =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
    : qx(sysctl -a 2>/dev/null | grep 'max-threads')
        =~ m{\.max-threads\s+=\s+(\d+)} ? $1
    : 2;

my (
    $queryFile, $subjectFile, $outputFolder, $maxSubject, $Evalue, $cpus
) = ( '', '', 'sandwichD', "1e-3", $cpuNumber );
my $options = GetOptions(
    "i=s" => \$queryFile,
    "s=s" => \$subjectFile,
    "o=s" => \$outputFolder,
    "e=s" => \$Evalue,
    "a=s" => \$cpus,
);

my $ownName = $0;
$ownName =~ s{.*/}{};
if ( !$queryFile || !$subjectFile) {
    print "usage: " . $ownName . " [options]\n";
    print "\noptions:\n";
    print "   -i query filename in fasta format, required\n";
    print "   -s subject filename in fasta format, required\n";
    print "   -o output folder, default: sandwichD\n";
    print "   -e evalue threshold, default 1e-3\n";
    print "   -a number of cpus to use, default in this machine: $cpuNumber\n";
    print "\n";
    exit;
}

####### avoid running this program because it's still under development
print "   Program still under development\n\n";
exit;

print "   Using $cpus cpu threads\n";

my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
print "   working in temp folder:\n\t$tempFolder\n";
unless( -d $outputFolder ) {
    system("mkdir $outputFolder");
}

my $querySeqRef = checkFastaFile($queryFile);
my $subjectSeqRef = checkFastaFile($subjectFile);
my @qids = sort keys %$querySeqRef;
my $toRun = @qids;
print "   preparing $toRun sequences\n";
my $tempSeqFile = "$tempFolder/query.seqTemp";
if( open( my $OS,">","$tempSeqFile" ) ) {
    for my $seqID ( @qids ) {
        my $sequence = $querySeqRef->{"$seqID"};
        print {$OS} $sequence;
    }
    close ($OS);
}
else {
    system "rm -r $tempFolder";
    die "\n\tcould not open $tempSeqFile\n\n";
}

print "   will be comparing $toRun sequences\n";
#my $blastOutputHashRef = runBlast($toRun);

print  "\n\tcleaning up ...";
system "rm -r $tempFolder";
print  "\tdone!\n\n";

sub checkFastaFile {
    my $inputSeqFile = $_[0];
    my $seqHashRef     = {};
    my $seqCount       = {};
    my $currentName    = '';
    my $totalSeqs      = 0;
    open( my $IS,"<",$inputSeqFile );
    while (<$IS>) {
        if ( $_ =~ m/^\s+\n/ ) { next; }
        if ( $_ =~ m/\>(\S+)\s/ ) {
            $currentName = $1;
            $currentName =~ s{^(gnl|lcl)\|}{};
            $currentName =~ s{\|$}{}g;
            $currentName =~ s{\|}{-}g;
            $seqCount->{$currentName}++;
            $totalSeqs++;
            s{^(gnl|lcl)\|}{};
        }
        $seqHashRef->{$currentName} .= $_;
    }
    close($IS);
    my $redundantNum = 0;
    while ( my ( $key, $value ) = each %$seqCount ) {
        $redundantNum = ( $value > $redundantNum ) ? $value : $redundantNum;
        if ( $value > 1 ) {
            print "sequence name $key has $value copies\n";
        }
    }
    if ( $redundantNum > 1 ) {
        system "rm -r $tempFolder";
        my $errmessage = "\nPlease remove redundant sequences\n\n";
        die $errmessage;
    }
    if( $totalSeqs < 1 ) {
        system "rm -r $tempFolder";
        die " no sequences or no fasta format in $inputSeqFile\n\n";
    }
    return $seqHashRef;
}

sub signalHandler {
    print "\n\tcleaning up ...";
    system "rm -r $tempFolder";
    die  " done!\n\n";
}

sub progressLine {
    my($done,$toGo,$decim) = @_;
    my $columns = qx(tput cols);
    chomp($columns);
    if( $columns > 50
            && $toGo >= 10
            && -t STDOUT
            && -t STDIN ) {
        if( $decim > 2 ) {
            ### better to count than show percent
            my $buffer     = " " x ( length($toGo) - length($done) + 1 );
            my $counter    = "[$buffer" . $done . " ]";
            my $countSpace = length($counter);
            my $pbwidth    = $columns - ( $countSpace + 3 );
            my $nhashes    = int($done / $toGo * $pbwidth);
            printf("\r% -${pbwidth}s% ${countSpace}s",
                   '#' x $nhashes,$counter);
        }
        else {
            my $percent    = sprintf("%.${decim}f",(100 * $done / $toGo));
            my $maxPC      = sprintf("%.${decim}f",100);
            my $buffer     = " " x ( length($maxPC) - length($percent) + 1 );
            my $counter    = "[$buffer" . $percent . "% ]";
            my $countSpace = length($counter);
            my $pbwidth    = $columns - ( $countSpace + 3 );
            my $nhashes    = int($done / $toGo * $pbwidth);
            printf("\r% -${pbwidth}s% ${countSpace}s",
                   '#' x $nhashes,$counter);
        }
        if( $done == $toGo ) {
            print "\n";
        }
    }
}
