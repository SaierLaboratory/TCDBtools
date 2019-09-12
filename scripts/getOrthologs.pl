#!/usr/bin/perl

## a little program to extract orthologs from a homologs file
## using the reciprocal best hits definition
## an example line is printed when no arguments are given
use strict;
use Getopt::Long;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

my @pwProgs = qw(
                      blastp
                      diamond
            );
my $matchProg = join("|",@pwProgs);

my $cov1       = 60;
my $cov2       = 99;
my $defCov     = 60;
my $defAln     = "F";
my $defPW      = 'blastp';
###### assignable:
my $faaQuery   = '';
my $faaSubject = '';
my $pwDir      = "compRuns";
my $rbhDir     = "RBHs";
my $minCov     = $defCov;
my $maxOverlap = 0.1; # maximum overlap in fused proteins/genes
my $maxEvalue  = 1e-6;
my $alnSeqs    = $defAln;
my $pwProg     = $defPW;

my $options = GetOptions(
    "i=s" => \$faaQuery,
    "s=s" => \$faaSubject,
    "m=s" => \$pwProg,
    "p=s" => \$pwDir,
    "o=s" => \$rbhDir,
    "c=s" => \$minCov,
    "f=s" => \$maxOverlap,
    "e=s" => \$maxEvalue,
    "a=s" => \$alnSeqs,
);

my $ownName = $0;
$ownName =~ s{.*/}{};
my $helpBrief
    = qq(about:\n)
    . qq(  This program produces files with reciprocal best hits\n\n)
    . qq(usage: $ownName [options]\n)
    . qq(example:\n)
    . qq(       $ownName -i [filename1.faa] -s [filename2.faa]\n)
    . qq(\noptions:\n)
    . qq(   -i query file in fasta format, required\n)
    . qq(       (files can be compressed with gzip or bzip2)\n)
    . qq(   -s subject file in fasta format, required\n)
    . qq(       (files can be compressed with gzip or bzip2)\n)
    . qq(   -m program for pairwise comparisons [$matchProg],\n)
    . qq(       default: $defPW\n)
    . qq(   -p comparison results directory, default compRuns\n)
    . qq(   -o RBH directory, default RBHs\n)
    . qq(   -c minimum coverage of shortest sequence [60 - 99],\n)
    . qq(       default $defCov\n)
    . qq(   -f maximum overlap between fused sequences [0 - 0.2],\n)
    . qq(       default 0.1\n)
    . qq(   -e maximum e-value, default 1e-6\n)
    . qq(   -a include aligned sequences in comparison results [T/F],\n)
    . qq(       default $defAln\n)
    . qq(\n)
    . qq(requirements:\n)
    . qq(  This program requires either blastp or diamond\n)
    . qq(  to cpmpare protein sequences, or appropriately\n)
    . qq(  formatted comparison results for the query/subject\n)
    . qq(  genomes\n\n)
    ;


if( !$faaQuery || !$faaSubject ) {
    print $helpBrief;
    exit;
}
unless( -f "$faaQuery" ) {
    print "Error: There's no $faaQuery file\n\n",$helpBrief;
}
unless( -f "$faaQuery" ) {
    print "Error: There's no $faaSubject file\n\n",$helpBrief;
}

#### directories where to find files:
### in case we need to run blastp
### temporary working directory:
my $queryGnm   = nakedName($faaQuery);
my $subjectGnm = nakedName($faaSubject);
my $cwd         = qx(pwd);
my $tempFolder  = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
my $pwDB        = $tempFolder . "/$subjectGnm";
my $alnFile
    = $pwDir      . "/$queryGnm.$subjectGnm." . $pwProg . ".bz2";
### directory for RBHs:
unless( -d "$rbhDir" ) {
    mkdir("$rbhDir");
}
my $pwProg = $pwProg =~ m{^($matchProg)$} ? lc($1) : $defPW;
print "comparing sequences with $pwProg\n";
my $minCov = $minCov >= $cov1 && $minCov <= $cov2 ? $minCov : $defCov;
print "minimum coverage of shortest sequence: $minCov\n";
my $maxOverlap
    = $maxOverlap >= 0 && $maxOverlap <= 0.2
    ? $maxOverlap : 0.1;
print "maximum overlap for fusions: $maxOverlap\n";
my $alnSeqs = $alnSeqs =~ m{^(T|F)$}i ? uc($1) : $defAln;
print "include aligned seqs in blast results: $alnSeqs\n";

##### the table format below is more informative than the default
my @pwTbl = qw(
                     qseqid
                     sseqid
                     evalue
                     bitscore
                     score
                     pident
                     ppos
                     qstart
                     qend
                     qlen
                     sstart
                     send
                     slen
             );
##### add alignments to blast results table?
if( $alnSeqs eq "T" ) {
    push(@pwTbl,"qseq","sseq");
}

my $pwTbl
    = $pwProg eq "blastp" ? join(" ","7",@pwTbl)
    : join(" ","6",@pwTbl);
my $blastOptions
    = qq( -query - -db $pwDB -evalue $maxEvalue -max_hsps 1 )
    . qq( -seg yes -soft_masking true )
    . qq( -outfmt '$pwTbl' );
my $diamondOptions
    = qq( --db $pwDB --evalue $maxEvalue --masking 0 --sensitive )
    . qq( --quiet --outfmt $pwTbl );
unless( -f "$alnFile" ) {
    if( $pwProg eq "blastp" ) {
        runBlastp("$faaQuery","$faaSubject");
    }
    else {
        runDiamond("$faaQuery","$faaSubject");
    }
}

#########################################################################
######### reciprocal best hits
#########################################################################
my $orthQuery   = "$queryGnm.$subjectGnm.rbh.bz2";
my $orthSubject = "$subjectGnm.$queryGnm.rbh.bz2";
my $tmpQuery    = $tempFolder . "/" . $orthQuery;
my $tmpSubject  = $tempFolder . "/" . $orthSubject;

## learn top homologs
print "   learning homologies\n";
##### first learn lines and their scores from query genome:
my %pair_query   = (); # obvious
my %pair_subject = (); # obvious
my %pair_bitsc   = (); # capture bit score
my %pair_E_value = (); # capture E-value
my %pair_line    = (); # capture line
my %ident_bitsc  = (); # capture bit score for identical proteins
my %ident_line   = (); # capture line for identical proteins
##### open original homologs
open( my $IN,"-|","bzip2 -qdc $alnFile" );
HOMOLOGLINE:
while(<$IN>) {
    next HOMOLOGLINE if( m{^#} );
    chomp;
    my(
        $query,$subject,
        $evalue,$bitscore,$score,$pident,$ppos,
        $qstart,$qend,$qlen,
        $sstart,$send,$slen,
        $qalnseq,$salnseq
    ) = split;
    ## percent coverages
    my $qcov = calcCoverage($qstart,$qend,$qlen);
    my $scov = calcCoverage($sstart,$send,$slen);
    if( ( $qcov >= $minCov ) || ( $scov >= $minCov ) ) {
        my $stats = join("\t",
                         $evalue,$bitscore,
                         $qstart,$qend,$qcov,
                         $sstart,$send,$scov
                     );
        if( ( $query eq $subject )
                && ( $bitscore > $ident_bitsc{"$query"} ) ) {
            $ident_bitsc{"$query"} = $bitscore;
            $ident_line{"$query"}
                = join("\t",$query,$subject,$stats);
        }
        my $pair = join("::",$query,$subject);
        if( $bitscore > $pair_bitsc{"$pair"} ) {
            $pair_line{"$pair"}    = join("\t",$query,$subject,$stats);
            $pair_query{"$pair"}   = $query;
            $pair_subject{"$pair"} = $subject;
            $pair_bitsc{"$pair"}   = $bitscore;
            $pair_E_value{"$pair"} = $evalue;
        }
    }
}
close($IN);

#### now sort, first bit-score (highest to lowest),
#### then E value (lowest to highest), then query, then subject
my @query_lines = ();
my @query_pairs
    = sort {
        $pair_bitsc{"$b"} <=> $pair_bitsc{"$a"}
            || $pair_E_value{"$a"} <=> $pair_E_value{"$b"}
            || $pair_query{"$a"} <=> $pair_query{"$b"}
            || $pair_query{"$a"} cmp $pair_query{"$b"}
            || $pair_subject{"$a"} <=> $pair_subject{"$b"}
            || $pair_subject{"$a"} cmp $pair_subject{"$b"}
        } keys %pair_query;

for my $pair ( @query_pairs ) {
    push(@query_lines,$pair_line{"$pair"});
}

#### now sort reciprocals
my @reciprocal_lines = ();
my @pretend_reciprocal
    = sort {
        $pair_bitsc{"$b"} <=> $pair_bitsc{"$a"}
            || $pair_E_value{"$a"} <=> $pair_E_value{"$b"}
            || $pair_subject{"$a"} <=> $pair_subject{"$b"}
            || $pair_subject{"$a"} cmp $pair_subject{"$b"}
            || $pair_query{"$a"} <=> $pair_query{"$b"}
            || $pair_query{"$a"} cmp $pair_query{"$b"}
        } keys %pair_query;

for my $pair ( @pretend_reciprocal ) {
    my (
        $query,$subject,
        $evalue,$bit_score,
        $q_start,$q_end,$qcov,
        $s_start,$s_end,$scov
    ) = split(/\t/,$pair_line{"$pair"});
    my $oppLine = join("\t",
                       $subject,$query,
                       $evalue,$bit_score,
                       $s_start,$s_end,$scov,
                       $q_start,$q_end,$qcov
                   );
    push(@reciprocal_lines,$oppLine);
}

###### now orthology
print "   extracting orthologs\n";
###### making this into a subroutine to allow working both ways:
produceRBH($tmpQuery,\@query_lines,\@reciprocal_lines);
if( -s "$tmpQuery" ) {
    system( "mv $tmpQuery $rbhDir/$orthQuery 2>/dev/null" );
}

produceRBH($tmpSubject,\@reciprocal_lines,\@query_lines);
if( -s "$tmpSubject" ) {
    system( "mv $tmpSubject $rbhDir/$orthSubject 2>/dev/null" );
}

#########################################################################
######### finish
#########################################################################
if( -d "$tempFolder" ) {
    print "   cleaning up\n";
    system("rm -r $tempFolder");
}
print "      Done with $ownName!\n\n";

sub produceRBH {
    my ($outFile,$rqLines,$rsLines) = @_;
    my $tmpFile = $outFile . ".tmp";
    my %query_best_hits  = ();
    my %query_bitsc      = ();
    my %query_E_value    = ();
    my %init_subject     = ();
    my %end_subject      = ();
    my %print_line       = ();
  QLINE:
    for my $line ( @{$rqLines} ) {
        my ( $query,$subject,$evalue,$bit_score,
             $q_start,$q_end,$qcov,
             $s_start,$s_end,$scov
         ) = split(/\s+/,$line);
        next QLINE if( $ident_bitsc{"$query"}   > 0 ); ### jump identicals
        my $print_line = join("\t",
                              $query,$subject,
                              $evalue,$bit_score,
                              $q_start,$q_end,$qcov,
                              $s_start,$s_end,$scov
                          );
        if( length($query_best_hits{"$query"}) > 0 ) {
            if( $query_bitsc{"$query"} == $bit_score ) {
                $query_best_hits{"$query"}        .= "," . $subject;
                $print_line{"$query.$subject"}     = $print_line;
                $init_subject{"$query.$subject"}   = $s_start;
                $end_subject{"$query.$subject"}    = $s_end;
            }
        }
        else {
            $query_best_hits{"$query"}         = $subject;
            $print_line{"$query.$subject"}     = $print_line;
            $init_subject{"$query.$subject"}   = $s_start;
            $end_subject{"$query.$subject"}    = $s_end;
            $query_bitsc{"$query"}             = $bit_score;
            $query_E_value{"$query"}           = $evalue;
        }
    }
    #####
    my %subject_best_hits = ();
    my %subject_bitsc     = ();
    my %subject_E_value   = ();
  SLINE:
    for my $line ( @{$rsLines} ) {
        my ( $query,$subject,$evalue,$bit_score,
             $q_start,$q_end,$qcov,
             $s_start,$s_end,$scov
         ) = split(/\s+/,$line);
        next SLINE if( $ident_bitsc{"$query"}   > 0 ); ### jump identicals
        if( length $subject_best_hits{"$query"} > 0 ) {
            if( $subject_bitsc{"$query"} == $bit_score ) {
                #    && $subject_E_value{"$query"} == $evalue ) {
                $subject_best_hits{"$query"} .= "," . $subject;
            }
        }
        else {
            $subject_best_hits{"$query"} = $subject;
            $subject_bitsc{"$query"}     = $bit_score;
            $subject_E_value{"$query"}   = $evalue;
        }
    }
    #####
    my %orth           = ();
    my %share_this_hit = ();
    ##### sort the queries by bit score to automatically sort the
    ##### best hits by their bit scores (I guess)
    my @queries = sort {
        $query_bitsc{"$b"} <=> $query_bitsc{"$a"}
            || $query_E_value{"$a"} <=> $query_E_value{"$b"}
            || $a <=> $b
            || $a cmp $b
        } keys(%query_best_hits);
    ##### now do the finding of orthologs
    for my $query ( @queries ) {
        for my $query_best_hit ( split(/\,/,$query_best_hits{"$query"}) ) {
            for my $subject_best_hit (
                split(/\,/,$subject_best_hits{"$query_best_hit"})
            ) {
                ## Reciprocal best hits:
                if( $subject_best_hit eq $query ) {
                    $orth{"$query"} .= ",$query_best_hit";
                    $orth{"$query"} =~ s/^\,+//;
                }
                else {
                    $share_this_hit{"$query_best_hit"} .= ",$query";
                    $share_this_hit{"$query_best_hit"} =~ s/^\,+//;
                }
            }
        }
    }
    ##### count printed lines (to ensure content)
    my $printedOrths = 0;
    open( my $ORTHS,"|-","bzip2 -9 >$tmpFile" );
    print {$ORTHS} "#",join("\t",
                            "Query","Subject",
                            "Evalue","bitScore",
                            "queryStart","queryEnd","queryCoverage",
                            "subjectStart","subjectEnd","subjectCoverage",
                            "Qualifier"
                        ),"\n";
    ###### first print identicals:
    for my $query ( sort { $a <=> $b || $a cmp $b } keys %ident_bitsc ) {
        print {$ORTHS} join("\t",$ident_line{"$query"},"Identical"),"\n";
        $printedOrths++;
    }
    ######
    ###### now reciprocal best hits, and fusions:
    for my $query ( sort { $a <=> $b || $a cmp $b } keys %orth ) {
      ORTH:
        for my $orth ( split(/,+/,$orth{"$query"}) ) {
            #### print reciprocal best hit ortholog
            print {$ORTHS} $print_line{"$query.$orth"},"\tRBH\n";
            $printedOrths++;
            ### find fusions from query to subject by testing
            ### other queries that find the same subject as their top
            ### hit, but not reciprocal
            my @share_this_hit = split(/\,/,$share_this_hit{"$orth"});
            my $shared = @share_this_hit;
            next ORTH unless( $shared > 0 );
            my ($i_query,$f_query)
                = ($init_subject{"$query.$orth"},$end_subject{"$query.$orth"});
            my $rbh_coords = join(":",$i_query,$f_query);
            my @coords = ("$rbh_coords");
          PARTNER:
            for my $partner ( @share_this_hit ) {
                #print $print_line{"$partner.$orth"},"\n";
                my $q_start = $init_subject{"$partner.$orth"};
                my $q_end   = $end_subject{"$partner.$orth"};
                ##### allowing a bit of overlap for fusions:
                my $q_max_o
                    = sprintf("%.0f",$maxOverlap * ($q_end - $q_start));
              COORDS:
                for my $coords ( @coords ) {
                    my($cover_start,$cover_end) = split(/:/,$coords);
                    ##### allowing a bit of overlap:
                    my $cover_max_o
                        = sprintf("%.0f",
                                  $maxOverlap*($cover_end - $cover_start));
                    my $min_overlap
                        = $cover_max_o > $q_max_o ? (-1 * $q_max_o)
                        : (-1 * $cover_max_o);
                    my $overlap1 = $q_start - ( $cover_end + 1);
                    my $overlap2 = $cover_start - ( $q_end + 1);
                    #print "OVERLAPS: $overlap1, $overlap2, $min_overlap\n";
                    if( $overlap1 < $min_overlap
                            && $overlap2 < $min_overlap ) {
                        next PARTNER;
                    }
                    #print "PASSED\n";
                }
                my $new_coords = join(":",$q_start,$q_end);
                push(@coords,$new_coords);
                if( $q_start < $i_query ) {
                    print {$ORTHS}
                        $print_line{"$partner.$orth"}
                        ,"\tLEFT of " . $query . "\n";
                    $printedOrths++;
                }
                else {
                    print {$ORTHS}
                        $print_line{"$partner.$orth"}
                        ,"\tRIGHT of " . $query . "\n";
                    $printedOrths++;
                }
            } # for partner
        } # for orth
    } # for query
    close($ORTHS);
    ### rename resulting file, then check if it has content
    ### erase otherwise
    rename($tmpFile,$outFile);
    if( -z "$outFile" || $printedOrths == 0 ) {
        print "   found no orthologs\n";
        unlink("$outFile");
    }
}

#################################################################
######################## subroutines ############################
#################################################################
sub nakedName {
    my $inName  = $_[0];
    my $outName = $inName;
    $outName =~ s{^\S+/}{};
    $outName =~ s{\.(Z|gz|bz2)$}{};
    $outName =~ s{\.faa\S*}{};
    $outName =~ s{\.fasta\S*}{};
    if( length($outName) > 0 ) {
        return("$outName");
    }
    else {
        print "there's a problem with $inName -> $outName\n";
        return();
    }
}

sub nameDB {
    my $inName  = $_[0];
    my $outName = nakedName("$inName");
    if( length($outName) > 0 ) {
        return("$tempFolder/$outName");
    }
    else {
        print "there's a problem with $inName -> $outName\n";
        return();
    }
}

sub signalHandler {
    if( length($tempFolder) > 1 && -d "$tempFolder" ) {
        print "\n\tcleaning up ...\n";
        system "rm -r $tempFolder";
    }
    else {
        print "\n\tquitting $ownName\n";
    }
    die  "\tdone!\n\n";
}

sub how2open {
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

sub inflateFile {
    my $origfile = $_[0];
    my( $opener ) = how2open($origfile);
    my $rootName = nakedName("$origfile");
    my $inflated = "$tempFolder/$rootName.faa";
    unless( -f "$inflated" ) {
        system("$opener $origfile > $inflated");
    }
    return("$inflated");
}

sub formatDB {
    my($file,$dbfile,$dbType) = @_;
    my($opener) = how2open("$file");
    if( $dbType eq "blast" ) {
        if( -f "$dbfile.nsq" ) {
            print "  the blast DB file is already there\n";
        }
        else {
            print "producing blastDB: $dbfile\n";
            my $mkblastdb
                = qq($opener $file |)
                . qq( makeblastdb -parse_seqids )
                . qq( -dbtype prot )
                . qq( -title $dbfile )
                . qq( -out $dbfile );
            system("$mkblastdb 1>/dev/null");
        }
    }
    elsif( $dbType eq "diamond" ) {
        if( -f "$dbfile.dmnd" ) {
            print "  the diamond DB file is already there\n";
        }
        else {
            print "producing diamondDB: $dbfile\n";
            my $mkdiamondDB
                = qq($opener $file |)
                . qq( diamond makedb -d $dbfile --quiet);
            system("$mkdiamondDB 1>/dev/null");
        }
    }
    else { ### lastal for now
        if( -f "$dbfile.bck" ) {
            print "  the lastal DB is already there\n";
        }
        else {
            print "producing lastDB: $dbfile\n";
            my $mklastdb
                = qq( $opener $file | )
                . qq(lastdb $dbfile );
            system("$mklastdb 1>/dev/null");
        }
    }
}

sub calcCoverage {
    my($start,$end,$ln) = @_;
    if( $ln < 1 ) {
        return();
    }
    else {
        my $coverage = 100 * ( ( $end - $start + 1 ) / $ln );
        my $rounded  = sprintf( "%.1f", $coverage );
        return($rounded);
    }
}

sub runBlastp {
    my($queryFile,$subjectFile) = @_;
    my $dbfile = nameDB("$subjectFile");
    formatDB("$subjectFile","$dbfile","blast");
    print "running blastp:\n   ",nakedName("$queryFile"),
        " vs ",nakedName("$subjectFile"),"\n";
    unless( -d "$pwDir" ) {
        mkdir("$pwDir");
    }
    my $catFile
        = $queryFile =~ m{\.bz2$} ? "bzip2 -qdc $queryFile"
        : $queryFile =~ m{\.gz$}  ? "gzip  -qdc $queryFile"
        : "cat $queryFile";
    my $blcommand = qq($catFile | blastp $blastOptions);
    my $lineCount = 0;
    open( my $PWBLAST,"|-","bzip2 -9 > $tempFolder/blast.bz2" );
    for my $blastLine ( qx($blcommand) ) {
        $lineCount++;
        print {$PWBLAST} $blastLine;
    }
    close($PWBLAST);
    if( $lineCount > 0 ) {
        system qq(mv $tempFolder/blast.bz2 $alnFile &>/dev/null);
        print "   ",nakedName("$queryFile"),
            " vs ",nakedName("$subjectFile")," done\n";
    }
    else {
        print "   blastp failed (empty file)\n";
        signalHandler();
    }
}

sub runDiamond {
    my($queryFile,$subjectFile) = @_;
    my $dbfile = nameDB("$subjectFile");
    formatDB("$subjectFile","$dbfile","diamond");
    print "running diamond:\n   ",nakedName("$queryFile"),
        " vs ",nakedName("$subjectFile"),"\n";
    unless( -d "$pwDir" ) {
        mkdir("$pwDir");
    }
    my $catFile
        = $queryFile =~ m{\.bz2$} ? "bzip2 -qdc $queryFile"
        : $queryFile =~ m{\.gz$}  ? "gzip  -qdc $queryFile"
        : "cat $queryFile";
    my $dmcommand = qq($catFile | diamond blastp $diamondOptions);
    my $lineCount = 0;
    open( my $PWDMD,"|-","bzip2 -9 > $tempFolder/diamond.bz2" );
    for my $dmLine ( qx($dmcommand) ) {
        $lineCount++;
        print {$PWDMD} $dmLine;
    }
    close($PWDMD);
    if( $lineCount > 0 ) {
        system qq(mv $tempFolder/diamond.bz2 $alnFile &>/dev/null);
        print "   ",nakedName("$queryFile"),
            " vs ",nakedName("$subjectFile")," done\n";
    }
    else {
        print "   diamond failed (empty file)\n";
        signalHandler();
    }
}
