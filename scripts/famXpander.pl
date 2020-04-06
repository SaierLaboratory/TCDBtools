#!/usr/bin/perl
#########################################################################
#									#
#	Author : Gabo Moreno-Hagelsieb          			#
#	Date of first draft: Dec 18, 2015              			#
#									#
#########################################################################

use strict;
use Getopt::Long;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

my $commandLine = join(" ",$0,@ARGV);

### ensure that external programs (psiblastp, makeblastdb, cd-hit) exist
my @missingsoft = ();
for my $xsoftware ( qw( psiblast blastdbcmd cd-hit ) ) {
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

####### default values for options
my $inputSeqFile = '';
my $nrDB         = 'nr';
my $outputFolder = 'famXOut';
my $maxSubject   = 10000;
my $Evalue       = "1e-7";
my $iEvalue      = "1e-5";
my $iters        = 1;
my $cutRange     = 'T';
my $minCoverage  = 0.8;
my $eitherCov    = 'F';
my $minSize      = 0.8;
my $maxSize      = 1.25;
my $filter       = 0.8;
my $cpus         = $cpuNumber;
my $rewrite      = 'T';
my $remote       = 'F';

my $ownName = $0;
$ownName =~ s{.*/}{};
my $helpMsg
    = qq(usage: " . $ownName . " [options]\n)
    . qq(\noptions:\n)
    . qq(   -i input filename in fasta format, required\n)
    . qq(   -d non-redundant database, default $nrDB\n)
    . qq(   -o output folder, default $outputFolder\n)
    . qq(   -n max number of aligned sequences to keep, default $maxSubject\n)
    . qq(   -e evalue threshold, default $Evalue\n)
    . qq(   -f psiblast evalue threshold, default $iEvalue\n)
    . qq(   -t psiblast iterations, default $iters\n)
    . qq(   -h keep only aligned region [T/F], default $cutRange\n)
    . qq(   -c minimum alignment coverage of query sequence,\n)
    . qq(       default $minCoverage\n)
    . qq(   -x coverage applies to either sequence, default $eitherCov\n)
    . qq(   -s minimal subject seq length relative to query seq length,\n)
    . qq(       default $minSize (ignored if -h T)\n)
    . qq(   -l maximal subject seq length relative to query seq length,\n)
    . qq(       default $maxSize (ignored if -h T)\n)
    . qq(   -r identity redundancy threshold (for cd-hit), default $filter\n)
    . qq(   -a number of cpus to use, default in this machine: $cpuNumber\n)
    . qq(   -w overwrite previous psiblast.tbl (if it exists) [T/F],\n)
    . qq(       default $rewrite\n)
    . qq(   -p run remotely (at ncbi) [T/F], default $remote\n)
    . qq(\n);


my $options
    = GetOptions(
        "i=s" => \$inputSeqFile,
        "d=s" => \$nrDB,
        "o=s" => \$outputFolder,
        "n=f" => \$maxSubject,
        "e=s" => \$Evalue,
        "f=s" => \$iEvalue,
        "t=s" => \$iters,
        "h=s" => \$cutRange,
        "c=f" => \$minCoverage,
        "x=s" => \$eitherCov,
        "s=f" => \$minSize,
        "l=f" => \$maxSize,
        "r=f" => \$filter,
        "a=s" => \$cpus,
        "w=s" => \$rewrite,
        "p=s" => \$remote,
    );

if ( !$inputSeqFile ) {
    print $helpMsg;
    exit;
}

### make sure that some options are well declared
$remote    = $remote    =~ m{^(T|F)$}i ? uc($1) : "F";
$cutRange  = $cutRange  =~ m{^(T|F)$}i ? uc($1) : "T";
$rewrite   = $rewrite   =~ m{^(T|F)$}i ? uc($1) : "T";
$eitherCov = $eitherCov =~ m{^(T|F)$}i ? uc($1) : "F";
if( $remote eq "T" && $iters > 2) {
    if( $nrDB ne "nr" ) {
        die "$nrDB is not available at NCBI for remote psiblast runs\n\n";
    }
    print qq(    can run only two iterations at NCBI, -t reset to "2"\n);
}

my @columns = qw(
                    qseqid
                    sacc
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
                    ssciname
                    qseq
                    sseq
);

if( $remote eq "T" ) {
    print "   Running psiblast at NCBI\n";
}
else {
    print "   Using $cpus cpu threads\n";
}

#### make sure coverage is a percent
my $pc_min_cover
    = ( $minCoverage >= 1 && $minCoverage <= 100 ) ? $minCoverage
    : $minCoverage > 100 ? 80
    : 100 * $minCoverage;
print "   Coverage: $pc_min_cover%\n";

my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
print "   working in temp folder:\n\t$tempFolder\n";
unless( -d $outputFolder ) {
    system("mkdir $outputFolder");
}

open( my $COMMANDLINE,">>","$outputFolder/command.line" );
print {$COMMANDLINE} $commandLine,"\n";
close($COMMANDLINE);

my $querySeqHashRef = checkFastaFile($inputSeqFile);
my @qids = sort keys %$querySeqHashRef;
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
    system "rm -r $tempFolder";
    die "\n\tcould not open $tempSeqFile\n\n";
}

print "   will be psiblasting $toRun sequences\n";
my $blastOutputHashRef = runBlast($toRun);
PasteSeqToFiles( $blastOutputHashRef,$cutRange );

print  "\n\tcleaning up ...";
system "rm -r $tempFolder";
print  "\tdone!\n\n";

sub PasteSeqToFiles {
    my ( $blastOutputHashRef, $cutRange ) = @_;
    my $refIDArrayRef = [];
    push @{ $refIDArrayRef }, keys %$blastOutputHashRef;
    my $tmp_file = "$tempFolder/results.faa";
    my $out_file = "$outputFolder/results.faa";
    open( my $OF, ">","$tmp_file" );
    my $printed_seqs = 0;
    my $seqNum       = (@$refIDArrayRef);
    if( $seqNum > 0 ) {
        #### using entry_batch to retrieve full seqs:
        if ( $cutRange eq 'F' ) {
            my $entryList = "$tempFolder/entry.list";
            open( my $ENTRYLS,">","$entryList" );
            print {$ENTRYLS} join("\n",sort @{ $refIDArrayRef }),"\n";
            close($ENTRYLS);
            my $get_seqs
                = qq(blastdbcmd -entry_batch $entryList -target_only )
                . qq(-outfmt "%a %i %t %s");
            print "   extracting $seqNum full sequences from $nrDB database:\n";
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
            for my $refID ( sort @{ $refIDArrayRef } ) {
                my $hr = $blastOutputHashRef->{$refID};
                my $fullname = $hr->{'sn'};
                my ( $s_start, $s_end ) = ( $hr->{'ss'}, $hr->{'se'} );
                my $range .= qq($s_start-$s_end);
                my $sequence
                    = ">" . $refID . " "
                    . $hr->{'full'} . ":$range [" . $fullname . "]\n"
                    . $hr->{'seq'}  . "\n";
                #print "      saving seq: $refID ($s_start-$s_end)\n";
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
            = qq(cd-hit -i $tmp_file -o $tmp_file.cdhit -c $filter);
        system("$cdhit_command >& /dev/null");
        print "       cleaning redundancy out ($filter)\n";
        system("mv $tmp_file.cdhit $out_file 2>/dev/null");

        #Include cd-hit cluster file in results directory
        system("mv $tmp_file.cdhit.clstr  $outputFolder 2>/dev/null"); 
    }
    else {
        system("mv $tmp_file $out_file 2>/dev/null");
    }
}

sub runBlast {
    my $totalQueries = $_[0];
    my $tempSeqFile = "$tempFolder/query.seqTemp";
    my $outPsiFile = "$outputFolder/psiblast.tbl";
    my $tmpPsiFile = "$tempFolder/psiblast.tbl";
    if( -f "$outPsiFile" && $rewrite eq "F" ) {
        print "    using pre-run psiblast results in $outPsiFile\n";
        my( $psiCount,$refIDBlastResultRef ) = parseBlast("$outPsiFile");
        if( $psiCount > 0 ) {
            return($refIDBlastResultRef);
        }
        else {
            print "    $outPsiFile has no results\n"
                . "     (rm $outPsiFile and run again)\n";
            signalHandler();
        }
    }
    else {
        #### if running locally, check for BLASTDB and nr within it
        if( $remote eq "F" ) {
            if( -f "$nrDB.pal" ) {
                print "   working with $nrDB\n";
            }
            elsif( length($ENV{"BLASTDB"}) > 0 ) {
                my @blastdbs = split(/:/,$ENV{"BLASTDB"});
                my $trueDBs  = @blastdbs;
                my $verifDir = 0;
                my $vefirNR  = 0;
                for my $testDir ( @blastdbs ) {
                    if( -d $testDir ) {
                        $verifDir++;
                        my $nr_sure = $testDir . "/" . $nrDB . ".pal";
                        if( -f "$nr_sure" ) {
                            $vefirNR++;
                        }
                    }
                }
                unless( $verifDir == $trueDBs ) {
                    system "rm -r $tempFolder";
                    my $diefor
                        = qq(\n\tproblems with BLASTDB directories:\n)
                        . qq(\t$ENV{"BLASTDB"}\n\n);
                    die "$diefor";
                }
                if( $vefirNR < 1 ) {
                    system "rm -r $tempFolder";
                    die qq(\n\tno $nrDB database in:\n\t$ENV{"BLASTDB"}\n\n);
                }
            }
            else {
                system "rm -r $tempFolder";
                die qq(\n\tBLASTDB is not set up\n\n);
            }
        }
        my $blastRootCmd
            = qq(psiblast -query $tempSeqFile -db $nrDB )
            . qq( -max_target_seqs $maxSubject )
            . qq( -max_hsps 1 )
            . qq( -evalue $Evalue -inclusion_ethresh $iEvalue )
            . qq( -outfmt ') . join(" ","7",@columns) . qq(');
        ##### a bit redundant, but might be good when
        ##### performing psiblast iterations
        if( $eitherCov eq 'F' ) {
            $blastRootCmd .= qq( -qcov_hsp_perc $pc_min_cover );
        }
        if( $remote eq "T" ) {
            print "   psiblasting at NCBI (might take a while):\n";
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
                my $blastCmdR = $blastRootCmd . qq( -remote);
                $blastCmdR =~ s{query\s+$tempSeqFile}{query $tempFolder/$query};
                open( my $PSIFL,">","$tmpFile" );
                print {$PSIFL} "# Iteration: 1\n";
                for my $remoteLine ( qx($blastCmdR 2>/dev/null) ) {
                    my @items = split(/\s+/,$remoteLine);
                    if( scalar(@items) > 3 || $remoteLine =~ m{^#} ) {
                        print {$PSIFL} $remoteLine;
                    }
                }
                close($PSIFL);
                #### now second iteration:
                if( $iters > 1 ) {
                    print "        preparing for second iteration\n";
                    my $sendCmd = $blastRootCmd;
                    $sendCmd
                        =~ s{query\s+$tempSeqFile}{query $tempFolder/$query};
                    if( my $pssm = prepareRemoteIter("$tmpFile","$sendCmd") ) {
                        print "        running second iteration ($query)\n";
                        my $blastCmd2 = $blastRootCmd . qq( -remote);
                        my $runPssm = "$tempFolder/$pssm";
                        $blastCmd2
                            =~ s{query\s+$tempSeqFile}{in_pssm $runPssm};
                        open( my $PSIFL2,">>","$tmpFile" );
                        print {$PSIFL2} "# Iteration: 2\n";
                        for my $remoteLine2 ( qx($blastCmd2 2>/dev/null) ) {
                            my @items2 = split(/\s+/,$remoteLine2);
                            if( scalar(@items2) > 3 || $remoteLine2 =~ m{^#} ) {
                                print {$PSIFL2} $remoteLine2;
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
            my $blastCmd = $blastRootCmd
                . qq( -num_iterations $iters -num_threads $cpus);
            print "   psiblasting now (might take a while):\n";
            runPlusSave("$tmpPsiFile","$blastCmd","$totalQueries");
        }
        my( $psiCount,$refIDBlastResultRef ) = parseBlast("$tmpPsiFile");
        if( $psiCount > 0 ) {
            system("mv $tmpPsiFile $outPsiFile 2>/dev/null");
            return($refIDBlastResultRef);
        }
        else {
            print "   no psiblast results to report\n\n";
            signalHandler();
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
        my( $qID,$acc,$full_s_id,
            $bit_score,$e_value,$identity,
            $q_start,$q_end,$qlen,
            $s_start,$s_end,$slen,
            $sciname,
            $qseq,$sseq
        ) = @items;
        my $qcov = 100 * ( ( $q_end - $q_start + 1 ) / $qlen );
        my $scov = 100 * ( ( $s_end - $s_start + 1 ) / $slen );
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
        ## eliminate pdb
        if( $full_s_id =~ m{pdb\|} ) {
            next BLASTRESULT;
        }
        ## pdb hereby eliminated
        # only counting those that will get all the way to the fasta file
        $psiCount++;
        if( exists $refIDBlastResultRef->{$acc} ) {
            my $oldstart = $refIDBlastResultRef->{$acc}->{'ss'};
            my $oldend   = $refIDBlastResultRef->{$acc}->{'se'};
            my $oldLn = $oldend - $oldstart;
            my $newLn = $s_end  - $s_start;
            if( $oldLn > $newLn ) {
                #print "jumping new $acc result\n";
                next BLASTRESULT;
            }
        }
        if( $cutRange eq "T" ) {
            $sseq =~ s{\-+}{}g;
            my $tempRef = {
                'ss'   => $s_start,
                'se'   => $s_end,
                'sn'   => $sciname,
                'full' => $full_s_id,
                'seq'  => $sseq,
            };
            $refIDBlastResultRef->{$acc} = $tempRef;
        }
        else {
            my $tempRef = {
                'ss' => $s_start,
                'se' => $s_end,
                'sn' => $sciname,
            };
            $refIDBlastResultRef->{$acc} = $tempRef;
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
            = qq(blastdbcmd -entry_batch $accList -target_only 2>/dev/null)
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
    my($tmpPsiFile,$blastCmd,$qtoRun) = @_;
    my $numLn    = length($qtoRun);
    my $iter     = 0;
    my $query    = "NA";
    my $queryCnt = 0;
    open( my $PSITMP,">","$tmpPsiFile");
    open( my $GETPSIRES,"-|","$blastCmd 2>/dev/null" );
    while(<$GETPSIRES>) {
        print {$PSITMP} $_;
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
