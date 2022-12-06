#!/usr/bin/perl
#########################################################
#						       	#
#	Author : Gabo Moreno-Hagelsieb         		#
#	Date (first version): Aug 05, 2016     	      	#
#	       						#
#########################################################


## a program to retrieve particular protein sequence sets
## from TCDB using a TCDB ID

use Getopt::Long;
use strict;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

my $tcdbID   = '';
my $outDir   = 'TCDB';
my $format   = 'fasta';
my $database = 'tcdb';

my @formats = qw(
                    fasta
                    column
                    blast
                    diamond
            );
my $formats = join("|",@formats);

my $ownName = $0;
$ownName =~ s{.*/}{};
my $help
    = qq(usage:\n)
    . qq(    $ownName [options]\n\n)
    . qq(\examples:\n)
    . qq(    $ownName -i 1.C -o TCDB1C\n)
    . qq(    $ownName -i 1.C.39 -o Subset -s TCDB1C/tcdb-1.C.faa\n)
    . qq(    $ownName -i 2 -o Class2 -f blast\n\n)
    . qq(options:\n)
    . qq(   -i TCDB ID, required. Example: 1.C.39\n)
    . qq(      "-i tcdb", "-i all", or "-i full" will bring the\n)
    . qq(      complete TCDB database\n)
    . qq(   -o output folder, default $outDir\n)
    . qq(   -f output format ($formats), default $format\n)
    . qq(   -s source database (a fasta file with TCDB entries),\n)
    . qq(      default: $database (online database)\n)
    . qq(\n);

my $options = GetOptions(
    "i=s" => \$tcdbID,
    "o=s" => \$outDir,
    "f=s" => \$format,
    "d=s" => \$database,
) or die $help;

if ( !$tcdbID ) {
    print $help;
    exit;
}

my $tcdbIDtest
    = $tcdbID =~ m{^\d+}               ? $tcdbID
    : $tcdbID =~ m{^(tcdb|full|all)$}i ? "tcdb"
    : "none";

if( $tcdbIDtest eq "none" ) {
    die "    $tcdbID is not a valid TCDB ID\n\n";
}
else {
    $tcdbID = $tcdbIDtest;
}

my $format   = $format =~ m{^($formats)$}i ? lc($format) : "fasta";
my $filename = $tcdbID;
my $matcher  = $tcdbID;
print "   looking for $tcdbID in TCDB\n";
print "   I will format the sequences into $format\n";
$matcher =~ s{\.}{\\.}g;
my $cntdots = $matcher =~ tr{.}{.};
if( $cntdots < 4) {
    $matcher .= "\\.";
}
else {
    $matcher .= "\\b";
}

my $ownName = $0;
$ownName =~ s{.*/}{};
my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
print "   working in temp folder:\n\t$tempFolder\n";
mkdir("$outDir") unless( -d "$outDir" );

#### decide if we use TCDB online or a file with "frozen"
#### tcdb sequences
my $tcdbref = bringTCDB("$database")
    or signalHandler("something is wrong with your database choice:\n $database");
my @keys
    = ( $tcdbID eq "tcdb" ) ? keys %{$tcdbref}
    : grep { m{^$matcher} } keys %{$tcdbref};

my $rootName = ( $tcdbID eq "tcdb" ) ? "tcdb" : "tcdb-$filename";
my $ending
    = $format eq "column" ? "clm"
    : "faa";

my $tmpFile = "$tempFolder/$rootName.$ending";
my $outFile = "$outDir/$rootName.$ending";
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
close($FAMFL);

if( $countPrinted > 0 ) {
    if( $format eq "blast" ) {
        my $date = qx(date);
        chomp($date);
        my $makeDbCmd
            = qq(makeblastdb -in $tmpFile -dbtype prot )
            . qq( -parse_seqids -hash_index -out $tempFolder/$rootName )
            . qq( -title "$rootName $date");
        system("$makeDbCmd >&/dev/null");
        system("mv $tempFolder/$rootName.* $outDir/ 2>/dev/null");
        print  "    cleaning up ...\n";
        system "rm -r $tempFolder";
        print "    $tcdbID is saved as a blast database at:\n"
            ,"\t$outDir/$rootName\n\n";
    }
    elsif( $format eq "diamond" ) {
        my $date = qx(date);
        chomp($date);
        my $makeDbCmd
            = qq(diamond makedb --in $tmpFile )
            . qq( --db $tempFolder/$rootName );
        system("$makeDbCmd >&/dev/null");
        system("mv $tempFolder/$rootName.* $outDir/ 2>/dev/null");
        print  "    cleaning up ...\n";
        system "rm -r $tempFolder";
        print "    $tcdbID is saved as a diamond database at:\n"
            ,"\t$outDir/$rootName\n\n";
    }
    else {
        system("mv $tmpFile $outFile 2>/dev/null");
        print  "    cleaning up ...\n";
        system "rm -r $tempFolder";
        print "    $tcdbID is saved as:\n\t$outFile\n\n";
    }
}
else {
    print "    no member with identity $tcdbID were found at TCDB\n\n";
    print  "\n\tcleaning up ...\n";
    system "rm -r $tempFolder";
}

sub bringTCDB {
    my $inputDB = $_[0];
    my $tcdbURL = "http://www.tcdb.org/";
    if( $inputDB eq "tcdb" ) {
        print "   working with online TCDB database\n";
        system("wget -N $tcdbURL/public/tcdb -O $tempFolder/tcdb >&/dev/null");
    }
    elsif( -f "$inputDB" ) {
        print "   working with local TCDB fasta file\n";
        system("cp $inputDB $tempFolder/tcdb >&/dev/null");
    }
    else {
        return();
    }
    if( -f "$tempFolder/tcdb" ) {
        my $id = "";
        my %seq = ();
        open( my $TCDBI,"<","$tempFolder/tcdb" );
        while(<$TCDBI>) {
            if( m{^>gnl\|TC-DB\|(\S+?)\s*\|(\S+)} ) {
                $id = join("-",$2,$1);
                #### remove version from identifier
                $id =~ s{\.\d+$}{};
            }
            elsif( m{>(\S+)} ) {
                $id = $1;
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
    else {
        return();
    }
}

sub signalHandler {
    my $message = $_[0];
    if( length($message) ) {
        print $message,"\n";
    }
    if( -d "$tempFolder" ) {
        print "\n\tcleaning up ...";
        system "rm -r $tempFolder";
    }
    die  " done!\n\n";
}
