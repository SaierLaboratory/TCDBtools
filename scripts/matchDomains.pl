#!/usr/bin/perl

#########################################################
#                                                       #
#       Author : Gabo Moreno-Hagelsieb                  #
#                                                       #
#########################################################

use strict;
use Getopt::Long;
use Pod::Text;
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);
my $columns = qx(tput cols);
chomp($columns);
my $width = $columns >= 100 ? 80 : $columns - 3;
my $parser = Pod::Text->new (sentence => 0, width => $width, margin => 1);

my $ownName = $0;
$ownName =~ s{.*/}{};

my @domainSets = qw(
                       cog
                       cdd
                       pfam
                       tigrfam
               );
my $matchDom = join('|',@domainSets);

my @programs = qw(
                     hmmscan
                     rpsblast
                     mmseqs
             );

my %defaultProg = (
    'cog'     => 'rpsblast',
    'cdd'     => 'rpsblast',
    'pfam'    => 'hmmscan',
    'tigrfam' => 'hmmscan',
);

my $matchProg   = join('|',@programs);
my $defEvalue   = 1e-3;
my $ncbiEvalue  = 1e-2;
my @queries     = ();
my $domainDB    = '';
my $scanProgram = '';
my $resultsDir  = 'scanDomains';
my $cpus        = 1;
my $Evalue      = '';
my $cluster     = 'F';
my $baseDBdir
    = -d $ENV{"DOMAINDB"}          ? $ENV{"DOMAINDB"}
    : -d "/usr/local/domainDBs"    ? "/usr/local/domainDBs"
    : -d "/ResearchData/domainDBs" ? "/ResearchData/domainDBs"
    : -d $ENV{"GENOMEDB"}          ? $ENV{"GENOMEDB"} . "/domainDBs"
    : ".";
### choosing database
my $defaultNA = "none";
my $cddDB
    = -d "$baseDBdir/cddDB" ? "$baseDBdir/cddDB"
    : $defaultNA;
my $xfamDB
    = -d "$baseDBdir/xfamDB" ? "$baseDBdir/xfamDB"
    : $defaultNA;
my $mmseqsDB
    = -d "$baseDBdir/mmseqsDB" ? "$baseDBdir/mmseqsDB"
    : $defaultNA;

my %progDBdir = (
    'rpsblast' => "$cddDB",
    'hmmscan'  => "$xfamDB",
    'mmseqs'   => "$mmseqsDB",
);

### check if there's more than one processor or assume there's 2.
my $cpu_count
    = qx(sysctl -a | grep 'cpu.thread_count')
    =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
    : qx(sysctl -a 2>/dev/null | grep 'max-threads')
    =~ m{\.max-threads\s+=\s+(\d+)} ? $1
    : 2;

my $podUsage
    = qq(=pod\n\n)
    . qq(=head1 NAME\n\n)
    . qq($ownName - Scan protein sequences against domain profiles\n\n)
    . qq(=head1 SYNOPSIS\n\n)
    . qq($ownName -q [fastaFile] -f [$matchDom] [options]\n\n)
    . qq(=head1 EXAMPLE\n\n)
    . qq($ownName -q GCF_000005845.faa.gz -f pfam -p mmseqs -o PFAM\n\n)
    . qq($ownName -q fastaFiles/*.faa.gz -f cdd -o CDD\n\n)
    . qq(=head1 OPTIONS\n\n)
    . qq(=over\n\n)
    . qq(=item B<-q>\n\n)
    . qq(query fasta file(s), required\n\n)
    . qq(=item B<-f>\n\n)
    . qq(domain database: file or [$matchDom], required\n\n)
    . qq(=item B<-p>\n\n)
    . qq(program [$matchProg], except for mmseqs, can be guessed from\n)
    . qq(database:\n\n)
    . qq(=over\n\n)
    . qq(=item -\n\n)
    . qq(rpsblast for cog and cdd\n\n)
    . qq(=item -\n\n)
    . qq(hmmscan for Pfam and TIGRFAM\n\n)
    . qq(=item -\n\n)
    . qq(mmseqs has to be specified in command line\n\n)
    . qq(=back\n\n)
    . qq(=item B<-o>\n\n)
    . qq(output folder, default: scanDomains\n\n)
    . qq(=item B<-e>\n\n)
    . qq(e-value threshold, default $defEvalue (NCBI uses $ncbiEvalue),\n)
    . qq(scientific notation acceptable (e.g. 1e-3)\n\n)
    . qq(=item B<-x>\n\n)
    . qq(number of CPUs to use, default: 1 (max: $cpu_count)\n\n)
    . qq(=item B<-c>\n\n)
    . qq(running in computer cluster [T|F], default 'F'\n\n)
    . qq(=back\n\n)
    . qq(=head1 DESCRIPTION\n\n)
    . qq(This program scans protein sequences against [$matchDom] domains\n)
    . qq( using either of hmmscan, rpsblast, or mmseqs\n\n)
    . qq(=cut\n\n)
    ;

GetOptions(
    "q=s{,}" => \@queries,
    "f=s"    => \$domainDB,
    "p=s"    => \$scanProgram,
    "o=s"    => \$resultsDir,
    "e=s"    => \$Evalue,
    "x=i"    => \$cpus,
    "c=s"    => \$cluster,
) or podhelp();

if ( !$queries[0] || !$domainDB ) {
    podhelp("I need a protein fasta file and a domain family:");
}
$domainDB =~ s{\.(psq|hmm)}{};
$domainDB =~ s{^\S+/}{};
$domainDB =~ s{\.\S+}{};

my $cpus = $cpus > 0 && $cpus <= $cpu_count ? $cpus : 1;
print "  Using $cpus cpu threads\n";

### working in cluster?
my $cluster = $cluster =~ m{(T|F)}i ? uc($1) : 'F';
print "  working in cluster: $cluster\n";

### choosing database
### (will need better ways fo figure out the available databases
###  and avoit genomeTools)

my $scanProgram
    = $scanProgram =~ m{^($matchProg)$}i ? $1
    : exists $defaultProg{"$domainDB"}   ? $defaultProg{"$domainDB"}
    : "none";

if( $scanProgram eq "none" ) {
    podhelp("I need a program to make these comparisons [$matchProg]")
}

my $Evalue
    = $Evalue > 0 && $Evalue < 1 ? $Evalue
    : $scanProgram eq 'hmmscan'  ? '--cut_ga'
    : $defEvalue;
print "  Using an E-value of $Evalue\n";

my @queries = do { my %seen; grep { !$seen{$_}++ } @queries };
my $countQueries = @queries;
my $foundQueries = 0;
my @missing      = ();
my %faaFile      = ();
for my $query ( @queries ) {
    if( -f "$query" ) {
        my $faaFile = $query;
        $query =~ s{\S+/}{};
        $query =~ s{\.(faa|fasta)\S*}{};
        $faaFile{"$query"} = $faaFile;
        $foundQueries++;
    }
    else {
        push(@missing,$query);
    }
}
if( $foundQueries == $countQueries ) {
    print "  found $countQueries query fasta files\n";
}
else {
    my $missing = @missing;
    my $error = "missing $missing query file(s):\n".join("\n",@missing)."\n";
    die $error;
}

my $dbPath = $progDBdir{"$scanProgram"};
####### test for domain database at the appropriate path:
opendir( my $DBDIR,"$dbPath" );
my @dbFiles
    = sort { -s $b <=> -s $a } grep { m{($domainDB)\.}i } readdir($DBDIR);
closedir($DBDIR);
my $dbName = $dbFiles[0] =~ m{($domainDB)}i ? $1 : 'none';
my $fullDB = $dbPath . "/" . $dbName;
if( $dbName eq 'none' ) {
    die "there is no $domainDB database for $scanProgram at:\n$dbPath\n\n";
}
else {
    if( $scanProgram eq 'hmmscan' ) {
        $fullDB .= ".hmm";
    }
    print "working with:\n $fullDB\n";
}

### fields for rpsblast
my @tableFields = qw(
                        qseqid
                        sseqid
                        evalue
                        bitscore
                        qstart
                        qend
                        qlen
                        sstart
                        send
                        slen
                 );

my $tableFields = qq(') . join(" ","7",@tableFields) . qq(');

### fields for mmseqs
my @mmseqsfields = qw(
                         query
                         target
                         evalue
                         bits
                         qstart
                         qend
                         qlen
                         tstart
                         tend
                         tlen
                 );

my $mmseqsfields = join(",",@mmseqsfields);

system "mkdir -p $resultsDir" unless( -d "$resultsDir" );
my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");

for my $query ( sort keys %faaFile ) {
    my $faaFile = $faaFile{"$query"};
    my ($cat,$queryFile) = figureCompression("$faaFile");
    my $file_name  = "$query." . lc($domainDB) . ".$scanProgram";
    my $tmpFile    = "$tempFolder/$file_name";
    my $outFile    = "$resultsDir/$file_name";
    my $tmpQuery   = $queryFile;
    $tmpQuery      =~ s{\S+/}{};
    $tmpQuery      = "$tempFolder/$tmpQuery";
    
    print "   running $scanProgram $query vs $domainDB\n";
    if( $cluster eq 'T' ) {
        print "    preparing database\n";
        system("cp $fullDB* $tempFolder/");
        print "    preparing query\n";
        system("cp $queryFile $tmpQuery");
    }
    my $dbFile
        = $cluster eq 'T' ? "$tempFolder/$domainDB"
        : "$fullDB";
    if( $scanProgram eq "rpsblast" ) {
        my $input
            = $cluster eq 'T' ? qq($cat $tmpQuery)
            : qq($cat $queryFile);
        my $rpsblast_command
            = qq( rpsblast -query - )
            . qq( -db $dbFile )
            . qq( -seg yes -soft_masking true )
            . qq( -num_threads $cpus )
            . qq( -evalue $Evalue )
            . qq( -parse_deflines )
            . qq( -comp_based_stats 0 )
            . qq( -outfmt $tableFields );
        my $logStuff
            = qx($input | $rpsblast_command | bzip2 -9 > $tmpFile.bz2 2>&1);
        if( verifyResults("$tmpFile.bz2") ) {
            system("mv $tmpFile.bz2 $outFile.bz2 2>/dev/null");
        }
        else {
            print "   no $domainDB matches for $query\n";
        }
    }
    elsif( $scanProgram eq "mmseqs" ) {
        unless( -f "$tmpQuery" ) {
            print "    preparing query\n";
            system("cp $queryFile $tmpQuery");
        }
        my $input = $tmpQuery;
        my $mmseqs_command
            = qq( mmseqs easy-search $input $dbFile $tmpFile $tempFolder)
            . qq( -e $Evalue )
            . qq( --alt-ali 5 )
            . qq( --threads $cpus )
            . qq( --comp-bias-corr 0 )
            . qq( --format-output "$mmseqsfields" );
        print "    running mmseqs\n";
        my $logStuff = qx($mmseqs_command 2>&1);
        if( verifyResults("$tmpFile") ) {
            system("bzip2 -f -9 $tmpFile");
            system("mv $tmpFile.bz2 $outFile.bz2 2>/dev/null");
        }
        else {
            print "   no $domainDB matches for $query\n";
            my $logFile = "$resultsDir/$scanProgram.log";
            open( my $LOG,">","$logFile" );
            print {$LOG} $logStuff;
            close($LOG);
            print "   output from $scanProgram in $logFile\n";
        }
    }
    else { #### default is hmmscan
        my $input
            = $cluster eq 'T' ? qq($cat $tmpQuery)
            : qq($cat $queryFile);
        my $threshold = $Evalue =~ m{cut_ga} ? $Evalue : "-E $Evalue";
        my $hmmscan_command
            = qq( hmmscan --cpu $cpus --noali $threshold -o /dev/null)
            . qq( --domtblout $tmpFile $dbFile - );
        my $logStuff = qx($input | $hmmscan_command 2>&1);
        if( verifyResults("$tmpFile") ) {
            system("bzip2 -f -9 $tmpFile");
            system("mv $tmpFile.bz2 $outFile.bz2 2>/dev/null");
        }
        else {
            print "   no $domainDB matches for $query\n";
            my $logFile = "$resultsDir/$scanProgram.log";
            open( my $LOG,">","$logFile" );
            print {$LOG} $logStuff;
            close($LOG);
            print "   output from $scanProgram in $logFile\n";
        }
    }
}

if( -d "$tempFolder" ) {
    print "\tcleaning up ...";
    system "rm -rf $tempFolder";
}
print "\n    done with $ownName\n\n";

sub signalHandler {
    if( -d "$tempFolder" ) {
        print "\n\tcleaning up ...";
        system "rm -rf $tempFolder";
        die  "    done!\n\n";
    }
    else {
        print "\n\ttemp files cleared out\n\n";
        die  "    done!\n\n";
    }
}

sub verifyResults {
    my $tmpFile = $_[0];
    my $openTest
        = $tmpFile =~ m{\.bz2$} ? "bzip2 -qdc $tmpFile" : "cat $tmpFile";
    my $verify = 0;
    open( my $TESTF,"-|","$openTest" );
    while(<$TESTF>) {
        if( $_ !~ m{^#} ) {
            $verify++;
        }
    }
    close($TESTF);
    return($verify);
}

sub podhelp {
    my $extraMessage = $_[0];
    open( my $PIPE,"|-","less -wiseMR" );
    if( length "$extraMessage" > 2 ) {
        print {$PIPE} "    ",$extraMessage,"\n\n";
    }
    $parser->output_fh($PIPE);
    $parser->parse_string_document($podUsage);
    exit;
}

sub figureCompression {
    my $rootName = $_[0];
    $rootName =~ s{\.(gz|bz2|Z)$}{};
    my $fullName
        = ( -f "$rootName.gz" )  ? "$rootName.gz"
        : ( -f "$rootName.Z" )   ? "$rootName.Z"
        : ( -f "$rootName.bz2" ) ? "$rootName.bz2"
        : ( -f "$rootName" )     ? "$rootName"
        : "none";
    my $catProg
        = $fullName =~ m{\.(gz|Z)$} ? "gzip -qdc"
        : $fullName =~ m{\.bz2$}    ? "bzip2 -qdc"
        : "cat";
    if( $fullName eq "none" ) {
        return();
    }
    else {
        return("$catProg","$fullName");
    }
}
