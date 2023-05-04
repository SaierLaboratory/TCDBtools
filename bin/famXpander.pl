#!/usr/bin/perl
#########################################################################
#									#
#	Author : Gabo Moreno-Hagelsieb          			#
#	Date of first draft: Dec 18, 2015                               #
#                                                                       #
#       Modified for parallelization by Arturo Medrano-Soto             #
#									#
#########################################################################

use strict;
use Getopt::Long;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

my $commandLine = join(" ",$0,@ARGV);

# Get the IP and computer name (for MacOS). For debugging purposes, this
# will identifiy the computer where a particular job ran.
chomp (my $compIP   = qx(ipconfig getifaddr en0));
chomp (my $compName = qx(hostname));

my @blasters = qw(
                     psiblast
                     blastp
                     diamond
                     mmseqs
             );
my $matchBlasters = join("|",@blasters);
my $matchNCBI     = join("|",'psiblast','blastp');

my %dbdir = (
    'psiblast' => $ENV{"BLASTDB"},
    'blastp'   => $ENV{"BLASTDB"},
    'diamond'  => $ENV{"DMDB"},
    'mmseqs'   => $ENV{"MMDB"},
);

my @dbs = qw(
                nr
                trembl
                uniref90
        );
my $matchDBs = join("|",@dbs);

### ensure that external programs (psiblastp, makeblastdb, cd-hit) exist
my @missingsoft = ();

chomp(my $host = qx(hostname));
print "Running in: $host\n";

my $pathBin = "/usr/local/bin";
my %wpath = ();
for my $xsoftware ( @blasters,qw(blastdbcmd cd-hit) ) {
    if( my $sfwpath = qx(which $xsoftware) ) {
        chomp($sfwpath);
        $wpath{"$xsoftware"} = $sfwpath;
    }
    elsif( -x "$pathBin/$xsoftware" ) {
        $wpath{"$xsoftware"} = "$pathBin/$xsoftware";
    }
    else {
        print "$host:\tprogram $xsoftware not found\n";
        push(@missingsoft,$xsoftware);
    }
}
my $cntMissing = @missingsoft;
if( $cntMissing > 0 ) {
    die "$host:\tcan't proceed because of missing software\n\t"
        . join("\n\t",@missingsoft) . "\n\n";
}

####### default values for options
my $defcpu     = 2;
my $defcov     = 80;
my $defblaster = 'psiblast';
my $defEvalue  = 1e-7;
my $defDB      = 'nr';
my $defRW      = 'F';
my $defHSP     = 1;
####### preset values for options
my $inputSeqFile = '';
my $nrDB         = $defDB;
my $outputFolder = 'famXOut';
my $maxSubject   = 10000;
my $maxHSP       = $defHSP;
my $Evalue       = $defEvalue;
my $iEvalue      = $defEvalue * 0.1;
my $iters        = 1;
my $cutRange     = 'T';
my $minCoverage  = $defcov;
my $eitherCov    = 'F';
my $minSize      = 0.8;
my $maxSize      = 1.25;
my $filter       = 0.8;
my $cpus         = $defcpu;
my $rewrite      = $defRW;
my $remote       = 'F';
my $blaster      = $defblaster;

### check if there's more than one processor or assume there's 2.
my $cpuNumber
    = qx(sysctl -a | grep 'cpu.thread_count')
    =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
    : qx(sysctl -a 2>/dev/null | grep 'max-threads')
    =~ m{\.max-threads\s+=\s+(\d+)} ? $1
    : $defcpu;


my $ownName = $0;
$ownName =~ s{.*/}{};
my $helpMsg
    = qq(usage:\n)
    . qq(    $ownName [options]\n)
    . qq(warning:\n)
    . qq(    Comparisons based on mmseqs aren't working yet.\n)
    . qq(\noptions:\n)
    . qq(   -i input filename in fasta format, required\n)
    . qq(   -d non-redundant database [$matchDBs], default $defDB\n)
    . qq(   -b problam for sequence comparisons [$matchBlasters],\n)
    . qq(      default $defblaster\n)
    . qq(   -o output folder, default $outputFolder\n)
    . qq(   -n max number of aligned sequences to keep, default $maxSubject\n)
    . qq(   -e evalue threshold, default $defEvalue\n)
    . qq(   -t psiblast|mmseqs iterations, default $iters\n)
    . qq(   -f psiblast|mmseqs iterations evalue threshold, default $iEvalue\n)
    . qq(   -h keep only aligned region [T/F], default $cutRange\n)
    . qq(   -c minimum alignment coverage of query sequence,\n)
    . qq(      default $defcov (it's a percent)\n)
    . qq(   -x coverage applies to either sequence, default $eitherCov\n)
    . qq(   -s minimal subject seq length relative to query seq length,\n)
    . qq(      default $minSize (ignored if -h T)\n)
    . qq(   -l maximal subject seq length relative to query seq length,\n)
    . qq(      default $maxSize (ignored if -h T)\n)
    . qq(   -r identity redundancy threshold (for cd-hit), default $filter\n)
    . qq(   -a number of cpus to use, default: $defcpu (max $cpuNumber)\n)
    . qq(   -w overwrite previous psiblast.tbl (if it exists) [T/F],\n)
    . qq(      default $defRW\n)
    . qq(   -p run psiblast remotely (at ncbi) [T/F], default $remote\n)
    . qq(\n)
    . qq(warning:\n)
    . qq(    Comparisons based on mmseqs aren't working yet.\n)
    . qq(\n);


my $options
    = GetOptions(
        "i=s" => \$inputSeqFile,
        "d=s" => \$nrDB,
        "b=s" => \$blaster,
        "o=s" => \$outputFolder,
        "n=i" => \$maxSubject,
        "e=f" => \$Evalue,
        "f=f" => \$iEvalue,
        "t=i" => \$iters,
        "h=s" => \$cutRange,
        "c=f" => \$minCoverage,
        "x=s" => \$eitherCov,
        "s=f" => \$minSize,
        "l=f" => \$maxSize,
        "r=f" => \$filter,
        "a=i" => \$cpus,
        "w=s" => \$rewrite,
        "p=s" => \$remote,
        "k=i" => \$maxHSP,
    ) or die $helpMsg;

if ( !$inputSeqFile ) {
    print $helpMsg;
    exit;
}

### make sure that some options are well declared
$remote    = $remote    =~ m{^(T|F)$}i ? uc($1) : "F";
$cutRange  = $cutRange  =~ m{^(T|F)$}i ? uc($1) : "T";
$rewrite   = $rewrite   =~ m{^(T|F)$}i ? uc($1) : $defRW;
$eitherCov = $eitherCov =~ m{^(T|F)$}i ? uc($1) : "F";
$blaster   = $blaster   =~ m{($matchBlasters)}i ? lc($1) : $defblaster;
$maxHSP    = $maxHSP > 0 ? $maxHSP : 1;
####### check remote shit #########
if( $remote eq "T" ) {
    if( $blaster eq "psiblast" ) {
        if( $nrDB ne "nr" ) {
            die "$nrDB is not available at NCBI for remote psiblast runs\n\n";
        }
        elsif( $iters > 2 ) {
            print qq(    can run only 2 iterations at NCBI, resetted to "2"\n);
        }
        print "   Running psiblast at NCBI\n";
    }
    else {
        die "cannot run $blaster remotely at NCBI (don't use -p T)\n\n";
    }
}
else {
    print "   Using $cpus cpu threads\n";
}

my @heading = qw(
                    qseqid
                    sseqid
                    bitscore
                    evalue
                    pident
                    qstart
                    qend
                    qlen
                    qcov
                    sstart
                    send
                    slen
                    scov
                    qseq
                    sseq
            );

my $fields = getFields($blaster);

#### make sure coverage is a percent
my $pc_min_cover
    = ( $minCoverage >= 1 && $minCoverage <= 100 ) ? $minCoverage
    : $minCoverage > 100                           ? $defcov
    : 100 * $minCoverage;
print "   Coverage: $pc_min_cover%\n";

my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
print "   working in temp folder:\n\t$tempFolder\n";
unless( -d "$outputFolder" ) {
    system("mkdir -p $outputFolder");
}

#Save original command line and computer where program is being
#executed
open( my $COMMANDLINE,">>","$outputFolder/command.line" );
print {$COMMANDLINE} "$compName ($compIP)\n";
print {$COMMANDLINE} $commandLine,"\n";
close($COMMANDLINE);


my $querySeqHashRef = checkFastaFile($inputSeqFile);
my @qids = sort keys %{ $querySeqHashRef };
my $toRun = @qids;
print "   preparing $toRun sequences\n";
my $tempSeqFile = "$tempFolder/query.seqTemp";
if( open( my $OS,">","$tempSeqFile" ) ) {
    for my $seqID ( @qids ) {
        my $sequence = $querySeqHashRef->{"$seqID"};
        print {$OS} $sequence;
    }
    close ($OS);
}
else {
    signalHandler("\tcould not open $tempSeqFile");
}

print "   will be ". $blaster . "ing $toRun sequences\n";
my $blastOutputHashRef = runBlast($toRun);
saveSeqs2File( $blastOutputHashRef,$cutRange );

print  "\n\tcleaning up ...";
if( -d $tempFolder ) {
    system "rm -r $tempFolder";
}
print  "\tdone!\n\n";

###########################################################################
###########################################################################
################# subroutines #############################################
###########################################################################
###########################################################################

sub saveSeqs2File {
    my ( $blastOutputHashRef, $cutRange ) = @_;
    my @ids = sort keys %{ $blastOutputHashRef };
    my $tmp_file = "$tempFolder/results.faa";
    my $out_file = "$outputFolder/results.faa";
    open( my $OF, ">","$tmp_file" );
    my $printed_seqs = 0;
    my $seqNum       = @ids;
    if( $seqNum > 0 ) {
        #### using entry_batch to retrieve full seqs:
        if ( $cutRange eq 'F' ) {
            my $entryList = "$tempFolder/entry.list";
            open( my $ENTRYLS,">","$entryList" );
            print {$ENTRYLS} join("\n",sort @ids ),"\n";
            close($ENTRYLS);
            my $get_seqs
                = qq($wpath{"blastdbcmd"} -db $nrDB )
                . qq(-entry_batch $entryList -target_only )
                . qq(-outfmt "%a %i %t %s");
            print "   extracting $seqNum full sequences from $nrDB database:\n";
            print "   $get_seqs\n";
            for my $seqLine ( qx($get_seqs) ) {
                $seqLine =~ s{^(\w+)\.\d+}{$1};
                $seqLine =~ s{(\S+)\n}{\n$1\n};
                print {$OF} ">",$seqLine;
                $printed_seqs++;
                ### feedback to terminal
                progressLine($printed_seqs,$seqNum,0);
                ### end feedback
            }
        }
        else { # $cutRange eq 'T'
            print "   saving $seqNum sequence segments:\n";
            for my $refID ( @ids ) {
                my $hr = $blastOutputHashRef->{$refID};
                #my $fullname = $hr->{'sn'};
                my $range = join("-",$hr->{'ss'},$hr->{'se'});
                #    . $hr->{'full'} . ":$range [" . $fullname . "]\n"
                my $sequence
                    = ">" . $refID . " "
                    . $hr->{'full'} . ":$range\n"
                    . $hr->{'seq'}  . "\n";
                print {$OF} $sequence;
                $printed_seqs++;
                ### feedback to terminal
                progressLine($printed_seqs,$seqNum,0);
                ### end feedback
            }
        }
    }
    close($OF);
    ##### filter redundancy out
    if( $printed_seqs > 1 && $filter < 1 ) {
        my $cdhit_command
            = qq($wpath{"cd-hit"} -i $tmp_file -o $tmp_file.cdhit -c $filter);
        system("$cdhit_command >& /dev/null");
        print "       cleaning redundancy out ($filter)\n";
        system("mv $tmp_file.cdhit $out_file 2>/dev/null");
        # Include cd-hit cluster file in results directory
        system("mv $tmp_file.cdhit.clstr  $outputFolder 2>/dev/null");
    }
    else {
        system("mv $tmp_file $out_file 2>/dev/null");
    }
}

sub runBlast {
    my $totalQueries = $_[0];
    my $outPsiFile   = "$outputFolder/psiblast.tbl";
    my $tmpPsiFile   = "$tempFolder/psiblast.tbl";
    if( -f "$outPsiFile" && $rewrite eq "F" ) {
        print "    using pre-run psiblast results in $outPsiFile\n";
        my( $psiCount,$refIDBlastResultRef ) = parseBlast("$outPsiFile");
        if( $psiCount > 0 ) {
            return($refIDBlastResultRef);
        }
        else {
            my $error
                = "    $outPsiFile has no results\n"
                . "    (rm $outPsiFile and run again,\n"
                . "    or use the '-w T' option)\n";
            signalHandler("$error");
        }
    }
    else {
        ##### a bit redundant, but might be good when
        ##### performing psiblast iterations
        if( $remote eq "T" ) {
            my $rootCmd = prepareCommand("$nrDB");
            print "   ". $blaster . "ing at NCBI (might take a while):\n";
            #### each remote run should contain only one query
            #### for network and wait reasons
            my $pCnt = 0;
            for my $query ( separateQueries() ) {
                $pCnt++;
                my $tmpFile = "$tmpPsiFile.$query";
                print "      "
                    . join(";  ",
                           "Query: $pCnt",
                           "Iteration: 1",
                           "ID: $query"
                       ),"\n";
                my $blastCmdR = $rootCmd . qq( -remote);
                $blastCmdR =~ s{query\s+$tempSeqFile}{query $tempFolder/$query};
                open( my $PSIFL,">","$tmpFile" );
                print {$PSIFL} "# Iteration: 1\n";
                print {$PSIFL} "# ",join("\t",@heading),"\n";
                for my $remoteLine ( qx($blastCmdR 2>/dev/null) ) {
                    if( $remoteLine =~ m{^#} ) {
                        print {$PSIFL} $remoteLine;
                    }
                    elsif( my $cleanLine = cleanLine("$remoteLine") ) {
                        print {$PSIFL} $cleanLine;
                    }
                }
                close($PSIFL);
                #### now second iteration:
                if( $iters > 1 ) {
                    print "        preparing for second iteration\n";
                    my $sendCmd = $rootCmd;
                    $sendCmd
                        =~ s{query\s+$tempSeqFile}{query $tempFolder/$query};
                    if( my $pssm = prepareRemoteIter("$tmpFile","$sendCmd") ) {
                        print "        running second iteration ($query)\n";
                        my $blastCmd2 = $rootCmd . qq( -remote);
                        my $runPssm = "$tempFolder/$pssm";
                        $blastCmd2
                            =~ s{query\s+$tempSeqFile}{in_pssm $runPssm};
                        open( my $PSIFL2,">>","$tmpFile" );
                        print {$PSIFL2} "# Iteration: 2\n";
                        print {$PSIFL2} "# ",join("\t",@heading),"\n";
                        for my $remoteLine2 ( qx($blastCmd2 2>/dev/null) ) {
                            if( $remoteLine2 =~ m{^#} ) {
                                print {$PSIFL2} $remoteLine2;
                            }
                            elsif( my $cleanLine = cleanLine("$remoteLine2") ) {
                                print {$PSIFL2} $cleanLine;
                            }
                        }
                        close($PSIFL2);
                    }
                    else {
                        print "        no need to run second iteration\n";
                    }
                }
            }
            system("cat $tmpPsiFile.* > $tmpPsiFile");
            system("rm $tmpPsiFile.*");
        }
        else {
            #### if running locally, check BLASTDB and database (ex nr)
            my $dbok    = checkDBs();
            my $rootCmd = prepareCommand("$dbok");
            my $addCmd
                = $blaster eq 'psiblast'
                ? qq( -num_threads $cpus -num_iterations $iters )
                : $blaster eq 'blastp'
                ? qq( -num_threads $cpus )
                : qq( -p $cpus );
            my $fullCmd = join( " ",$rootCmd,$addCmd );
            $fullCmd =~ s{\s+}{ }g;
            print "   ". $blaster . "ing now (might take a while):\n";
            runPlusSave("$tmpPsiFile","$fullCmd","$totalQueries");
        }
        my( $psiCount,$refIDBlastResultRef ) = parseBlast("$tmpPsiFile");
        if( $psiCount > 0 ) {
            system("mv $tmpPsiFile $outPsiFile 2>/dev/null");
            return($refIDBlastResultRef);
        }
        else {
            signalHandler("   no $blaster results to report\n");
        }
        unlink($tempSeqFile);
    }
}

sub parseBlast {
    my $psiBlastFile = $_[0];
    my $refIDBlastResultRef = {};
    my $psiCount = 0;
    open( my $PSIBL,"<","$psiBlastFile" );
  BLASTRESULT:
    while( my $blastResult = <$PSIBL> ) {
        chomp $blastResult;
        if( $blastResult =~ m{^\s*$} || $blastResult =~ m{^#} ) {
            next BLASTRESULT;
        }
        my @items = split(/\t/,$blastResult);
        my $nitems = @items;
        if( $nitems < 15 ) {
            next BLASTRESULT;
        }
        my( $qID,$sid,
            $bit_score,$e_value,$identity,
            $q_start,$q_end,$qlen,$qcov,
            $s_start,$s_end,$slen,$scov,
            $qseq,$sseq
        ) = @items;
        if( $eitherCov eq 'T' ) {
            if( $qcov < $pc_min_cover && $scov < $pc_min_cover ) {
                next BLASTRESULT;
            }
        }
        else {
            if( $qcov < $pc_min_cover ) {
                next BLASTRESULT;
            }
        }
        # make sure that the line has results
        if( $qlen < 5 ) {
            next BLASTRESULT;
        }
        my $rel_size = $slen/$qlen;
        if( $rel_size < $minSize ) {
            next BLASTRESULT;
        }
        if( $rel_size > $maxSize ) {
            if( $cutRange eq "F" ) {
                next BLASTRESULT;
            }
        }
        # only counting those that will get all the way to the fasta file
        $psiCount++;
        if( exists $refIDBlastResultRef->{$sid} ) {
            my $oldstart = $refIDBlastResultRef->{$sid}->{'ss'};
            my $oldend   = $refIDBlastResultRef->{$sid}->{'se'};
            my $oldLn = $oldend - $oldstart;
            my $newLn = $s_end  - $s_start;
            if( $oldLn > $newLn ) {
                #print "jumping new $sid result\n";
                next BLASTRESULT;
            }
        }
        if( $cutRange eq "T" ) {
            $sseq =~ s{\-+}{}g;
            my $tempRef = {
                'ss'   => $s_start,
                'se'   => $s_end,
                #'sn'   => $sciname,
                'full' => $sid,
                'seq'  => $sseq,
            };
            $refIDBlastResultRef->{$sid} = $tempRef;
        }
        else {
            my $tempRef = {
                'ss' => $s_start,
                'se' => $s_end,
                #'sn' => $sciname,
            };
            $refIDBlastResultRef->{$sid} = $tempRef;
        }
    }
    close($PSIBL);
    return($psiCount,$refIDBlastResultRef);
}

sub extractName {
    my $fullname = $_[0];
    $fullname =~ s{\[|\(|\)|\]}{}g;
    if( $fullname =~ m{^([A-Z])\S*\s([a-z])([a-z])} ) {
        my $fix = $1 . $2 . $3;
        return("$fix");
    }
    elsif( $fullname =~ m{^([A-Z][a-z][a-z])} ) {
        my $fix = $1;
        return("$fix");
    }
    else {
        return("Unk");
    }
}

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
        signalHandler("Please remove redundant sequences\n");
    }
    if( $totalSeqs < 1 ) {
        signalHandler("no sequences or no fasta format in $inputSeqFile\n");
    }
    return $seqHashRef;
}

sub signalHandler {
    my $msg = $_[0];
    print "\n\tcleaning up ...\n";
    if( -f "$tempFolder" ) {
        system "rm -rf $tempFolder";
    }
    if( length("$msg") > 1 ) {
        die  "$msg\n done!\n\n";
    }
    else {
        die  " done!\n\n";
    }
}

sub prepareRemoteIter {
    my ($rawFile,$blastPssm) = @_;
    my $pssmFile = "$tempFolder/pssm";
    system("rm ${pssmFile}* &>/dev/null");
    my $accList  = "$tempFolder/accList";
    my $dbFile   = "$tempFolder/dbFile";
    ### first extract the identifiers from the raw psiblast file
    #my %listFor = ();
    my %cAcc    = ();
    open( my $RAWPSI,"<","$rawFile" );
    while(<$RAWPSI>) {
        next if( m{^#} );
        chomp;
        my @items = split;
        my $cItems = @items;
        if( $cItems > 14 ) {
            #$listFor{"$items[0]"}++;
            $cAcc{"$items[1]"}++;
        }
    }
    close($RAWPSI);
    my @accs = keys %cAcc;
    my $cAcc = @accs;
    if( $cAcc < 10 ) {
        print "        no results to produce pssms\n";
        return();
    }
    else {
        open( my $LIST,">","$accList" );
        print {$LIST} join("\n",@accs),"\n";
        close($LIST);
        ### now build a blast database
        my $buildDBCmd
            = qq($wpath{"blastdbcmd"} -entry_batch $accList -target_only 2>/dev/null)
            . qq( | makeblastdb -dbtype prot -out $dbFile -title "dbFile");
        system("$buildDBCmd &>/dev/null");
        ### now run psiblast to build pssm files as necessary
        ### here it cannot use a specific nr because it runs with
        ### databases at NCBI
        $blastPssm =~ s{db\s+nr\s+}{db $dbFile };
        $blastPssm
            .= qq( -save_pssm_after_last_round )
            .  qq( -out_pssm $pssmFile -save_each_pssm);
        system("$blastPssm >& /dev/null");
        opendir( my $TMPDIR,"$tempFolder");
        my @pssms = grep { m{^pssm} } readdir($TMPDIR);
        my $cPssm = @pssms;
        if( $cPssm > 0 ) {
            my @returned
                = sort { length($a) <=> length($b) || $a cmp $b } @pssms;
            return($returned[0]);
        }
        else {
            print "        did not produce pssms\n";
            return();
        }
    }
}

sub runPlusSave {
    my($tmpPsiFile,$compCmd,$qtoRun) = @_;
    my $numLn    = length($qtoRun);
    my $iter     = 0;
    my $query    = "NA";
    my $queryCnt = 0;
    open( my $PSITMP,">","$tmpPsiFile");
    print {$PSITMP} "# ",join("\t",@heading),"\n";
    open( my $GETPSIRES,"-|","$compCmd 2>/dev/null" );
  GETTINGBLASTED:
    while(<$GETPSIRES>) {
        if( m{^#} ) {
            print {$PSITMP} $_;
        }
        elsif( my $cleanLine = cleanLine("$_") ) {
            print {$PSITMP} $cleanLine;
        }
        if( m{Iteration:\s+(\d+)} ) {
            $iter = $1;
        }
        if( m{Query:\s+(\S+)} ) {
            my $newquery = $1;
            if( $newquery ne $query ) {
                $queryCnt++;
            }
            $query = $newquery;
            my $pCnt = " " x ($numLn - length($queryCnt)) . $queryCnt;
            print "      "
                .join(";  ",
                      "Query: $pCnt",
                      "Iteration: $iter",
                      "ID: $query"
                  ),"\n";
        }
    }
    close($GETPSIRES);
    close($PSITMP);
    ### return something to verify that this ran all right
    if( $queryCnt > 0 ) {
        return($queryCnt);
    }
    else {
        return();
    }
}

sub cleanLine {
    ##### clean result line
    my $line = $_[0];
    chomp($line);
    my( $qID,$sID,
        $bitscore,$evalue,$identity,
        $qstart,$qend,$qlen,
        $sstart,$send,$slen,
        $qseq,$sseq
    ) = split(/\t/,$line);
    if( $qlen eq '' ) {
        return();
    }
    elsif( $sID =~ m{pdb\|} ) {
        ## eliminate pdb
        return();
        ## pdb hereby eliminated
    }
    else {
        my $qcov = calcCoverage($qstart,$qend,$qlen);
        my $scov = calcCoverage($sstart,$send,$slen);
        my $fixed =
            join("\t",cleanID($qID),cleanID($sID),
                 $bitscore,$evalue,$identity,
                 $qstart,$qend,$qlen,$qcov,
                 $sstart,$send,$slen,$scov,
                 $qseq,$sseq) . "\n";
        return($fixed);
    }
}

sub separateQueries {
    my @queries = ();
    for my $seqID ( @qids ) {
        my $sequence = $querySeqHashRef->{"$seqID"};
        #print "separating $seqID into $tempFolder/$seqID\n";
        open( my $SINGLE,">","$tempFolder/$seqID");
        print {$SINGLE} $sequence;
        close($SINGLE);
        push(@queries,"$seqID");
    }
    return(@queries);
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

sub calcCoverage {
    my( $start,$end,$ln ) = @_;
    my $coverage = sprintf("%.1f",100 * ( ( $end - $start + 1 ) / $ln ));
    return($coverage);
}

sub cleanID {
    my $ugly = $_[0];
    my $clean = $ugly;
    $clean =~ s{^\S*?\|}{};
    $clean =~ s{\|\S*}{};
    my $cleanID = length($clean) > 1 ? $clean : $ugly;
    return($cleanID);
}

sub getFields {
    my $blaster = $_[0];
    my @columns
        = $blaster =~ m{$matchNCBI|diamond}
        ? qw(
                qseqid
                sseqid
                bitscore
                evalue
                pident
                qstart
                qend
                qlen
                sstart
                send
                slen
                qseq
                sseq
        )
        : qw(
                query
                target
                bits
                evalue
                pident
                qstart
                qend
                qlen
                tstart
                tend
                tlen
                qaln
                taln
        );

    my $fields
        = $blaster =~ m{$matchNCBI} ? q(").join(" ","6",@columns).q(")
        : $blaster eq 'diamond'     ? join(" ","6",@columns)
        : join(",",@columns);
    return($fields);
}

sub checkDBs {
    my $dbok  = 'no';
    if( $blaster =~ m{$matchNCBI}
        || $cutRange eq 'F' ) {
        print "   checking NCBI formatted databases\n";
        if( $blaster !~ m{$matchNCBI} ) {
            print "   needed to extract full sequences\n";
        }
        my $blastdbok = checkSpecificDBs("psiblast");
        if( $blaster =~ m{$matchNCBI} ) {
            $dbok = $blastdbok;
        }
    }
    else {
        $dbok = checkSpecificDBs("$blaster");
    }
    return($dbok);
}

sub checkSpecificDBs {
    my $testblaster = $_[0];
    my $otherdbok = 'no';
    my $dbdir     = $dbdir{"$blaster"};
    my $dbfile
        = $blaster =~ m{$matchNCBI} ? "$nrDB.pal"
        : $blaster eq 'diamond'     ? "$nrDB.dmnd"
        : "$nrDB";
    #### check current dir:
    if( -f "$dbfile" ) {
        $otherdbok = $nrDB;
        print "   working with $otherdbok\n";
    }
    elsif( length($dbdir) > 0 ) {
        my @otherdbs  = split(/:/,$dbdir);
        my @trueDBs   = ();
        my $verifyDBs = 0;
        my $verifyNR  = 0;
        for my $testDir ( @otherdbs ) {
            if( -d $testDir ) {
                push(@trueDBs,$testDir);
                my $sureNR = $testDir . "/" . $dbfile;
                if( -f "$sureNR" ) {
                    $verifyNR++;
                    if( $verifyNR == 1 ) {
                        $otherdbok = $testDir . "/" . $nrDB;
                        print "   working with $sureNR\n";
                    }
                }
            }
        }
        my $foundDBs = @trueDBs;
        if( $foundDBs == 0 ) {
            my $error
                = qq(\tDatabase directory not found:\n$dbdir\n);
            signalHandler("$error");
        }
        elsif( $verifyNR == 0 ) {
            my $error
                = qq(\tno $nrDB database in:\n$dbdir\n);
            signalHandler("$error");
        }
    }
    else {
        my $dbsetup
            = $blaster =~ m{$matchNCBI} ? q($ENV{"BLASTP"})
            : $blaster eq 'diamond'     ? q($ENV{"DMDB"})
            : q($ENV{"MMDB"});
        signalHandler(qq(\tDatabase directory is not set up: $dbsetup));
    }
    if( $blaster =~ m{$matchNCBI} ) {
        return($nrDB);
    }
    else {
        return($otherdbok);
    }
}

sub prepareCommand {
    my $trueDB = $_[0];
    my $init
        = $blaster =~ m{$matchNCBI} ? $wpath{"$blaster"}
        : $blaster eq 'diamond'     ? $wpath{"$blaster"} . " blastp"
        : $wpath{"$blaster"} . " easy-search";
    my $moreCmd
        = $blaster =~ m{$matchNCBI|diamond}
        ? qq( -query $tempSeqFile )
        . qq( -db TRUEDB )
        . qq( -max_target_seqs $maxSubject )
        . qq( -max_hsps $maxHSP )
        . qq( -evalue $Evalue )
        . qq( -outfmt $fields)
        : qq(mmseqsthingie);
    if( $blaster eq 'psiblast' ) {
        $moreCmd .= qq( -inclusion_ethresh $iEvalue )
    }
    elsif( $blaster eq 'diamond' ) {
        $moreCmd =~ s{\s+(-\w)}{ -$1}g;
        $moreCmd =~ s{(\w)_(\w)}{$1-$2}g;
    }
    $moreCmd =~ s{\s+TRUEDB\s+}{ $trueDB };
    return( join(" ",$init,$moreCmd) );
}
