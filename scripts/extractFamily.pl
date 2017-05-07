#!/usr/bin/env perl
#########################################################
#						       	#
#	Author : Gabo Moreno-Hagelsieb         		#
#	Date (first version): Aug 05, 2016     	      	#
#	       						#
#########################################################


## a program to retrieve particular protein sequence sets
## from TCDB using a TCDB ID

use Getopt::Long;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

my ( $tcdbFamID, $outputFolder, $format )
    = ( '', 'Families', 'fasta' );
my $options = GetOptions(
    "i=s" => \$tcdbFamID,
    "o=s" => \$outputFolder,
    "f=s" => \$format,
);

my $ownName = $0;
$ownName =~ s{.*/}{};
if ( !$tcdbFamID ) {
    print "usage: " . $ownName . " [options]\n";
    print "\noptions:\n";
    print "   -i TCDB family ID, required. Example: 1.C.39\n"
        . qq(      "-i tcdb", "-i all", or "-i full" will bring the\n)
        . qq(      complete TCDB database\n);
    print "   -o output folder, default Families\n";
    print "   -f output format (fasta|column|blast), default fasta\n";
    print "\n";
    exit;
}

my $tcdbFamIDtest
    = $tcdbFamID =~ m{^\d+\.[A-Z]+\.\d+} ? $tcdbFamID
    : $tcdbFamID  eq  "tcdb"             ? "tcdb"
    : $tcdbFamID  eq  "full"             ? "tcdb"
    : $tcdbFamID  eq  "all"              ? "tcdb"
    : "none";

if( $tcdbFamIDtest eq "none" ) {
    die "    $tcdbFamID is not a valid TCDB family ID\n\n";
}
else {
    $tcdbFamID = $tcdbFamIDtest;
}

my $format = $format =~ m{fasta|column|blast}i ? $format : "fasta";
my $filename = $tcdbFamID;
my $matcher  = $tcdbFamID;
print "   looking for $tcdbFamID in TCDB\n";
print "   I will print the family sequences into a $format format\n";
if ( $matcher =~ m{\d+$} ) {
    $matcher =~ s{\.}{\\.}g;
    my $cntdots = $matcher =~ tr{.}{.};
    if( $cntdots < 4) {
        $matcher .= "\\.";
    }
    else {
        $matcher .= "\\b";
    }
}

my $ownName = $0;
$ownName =~ s{.*/}{};
my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
print "   working in temp folder:\n\t$tempFolder\n";
mkdir("$outputFolder") unless( -d "$outputFolder" );
my $tcdbref = bringTCDB();
my @keys
    = ( $tcdbFamID eq "tcdb" ) ? keys %$tcdbref
    : grep { m{$matcher} }  keys %$tcdbref;

my $rootName = ( $tcdbFamID eq "tcdb" ) ? "tcdb" : "family-$filename";
my $ending
    = $format eq "column" ? "clm"
    : "faa";

my $tmpFile = "$tempFolder/$rootName.$ending";
my $outFile = "$outputFolder/$rootName.$ending";
my $countPrinted = 0;
open( my $FAMFL,">","$tmpFile" );
for my $id ( sort @keys ) {
    my $print_id = $id;
    $print_id =~ s{lcl\|}{};
    if( $format eq "column" ) {
        print {$FAMFL} join("\t",$print_id,$tcdbref->{"$id"}),"\n";
    }
    elsif( $format eq "blast" ) {
        print {$FAMFL} ">lcl|",$id,"\n",$tcdbref->{"$id"},"\n";
    }
    else {
        print {$FAMFL} ">",$id,"\n",$tcdbref->{"$id"},"\n";
    }
    $countPrinted++;
}
close($FAML);
if( $countPrinted > 0 ) {
    if( $format eq "blast" ) {
        my $date = qx(date);
        chomp($date);
        my $makeDbCmd
            = qq(makeblastdb -in $tmpFile -dbtype prot )
            . qq( -parse_seqids -hash_index -out $tempFolder/$filename )
            . qq( -title "$filename $date");
        system("$makeDbCmd >&/dev/null");
        system("mv $tempFolder/$filename.p* $outputFolder/");
        print  "    cleaning up ...\n";
        system "rm -r $tempFolder";
        print "    $tcdbFamID is saved as blast database at:\n"
            ,"\t$outputFolder/$filename\n\n";
    }
    else {
        system("mv $tmpFile $outFile 2>/dev/null");
        print  "    cleaning up ...\n";
        system "rm -r $tempFolder";
        print "    $tcdbFamID is saved as:\n\t$outFile\n\n";
    }
}
else {
    print "    no member of family $tcdbFamID were found at TCDB\n\n";
    print  "\n\tcleaning up ...\n";
    system "rm -r $tempFolder";
}

sub bringTCDB {
    system("wget -N http://www.tcdb.org/public/tcdb -O $tempFolder/tcdb >&/dev/null");
    my $id = "";
    my %seq = ();
    open( my $TCDBI,"<","$tempFolder/tcdb" );
    while(<$TCDBI>) {
        if( m{^>gnl\|TC-DB\|(\S+?)\s*\|(\S+)} ) {
            $id = join("-",$2,$1);
        }
        else {
            chomp;
            $seq{"$id"} .= $_;
        }
    }
    close($TCDBI);
    my $count = keys %seq;
    if( $count > 0 ) {
        return(\%seq);
    }
    else {
        return();
    }
}

sub signalHandler {
    print "\n\tcleaning up ...";
    system "rm -r $tempFolder";
    die  " done!\n\n";
}
