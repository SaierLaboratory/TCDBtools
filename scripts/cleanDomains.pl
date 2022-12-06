#!/usr/bin/env perl

#########################################################
#						       	#
#	Author : Gabo Moreno-Hagelsieb         		#
#	       						#
#########################################################

use strict;
use Getopt::Long;
### use Pod::Usage qw(pod2usage);
use Pod::Text;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
#use sigtrap qw(handler signalHandler normal-signals);
my $columns = qx(tput cols);
chomp($columns);
my $width = $columns >= 100 ? 80 : $columns - 3;
my $parser = Pod::Text->new (sentence => 0, width => $width, margin => 1);

my $ownName = $0;
$ownName =~ s{.*/}{};

my @domainSets = qw(
                       cog
                       cd
                       cdd
                       pfam
                       tigrfam
                       superfamily
                       VFDB
                       Toxins
               );

my $matchDom = join('|',@domainSets);

my @progs = qw(
                  hmmscan
                  rpsblast
                  mmseqs
          );
my $mProgs = join("|",@progs);

####### defaults:
my $domFile      = '';
my $domFamily    = '';
my $outputFolder = 'cleanDomains';
my $minCover     = 0.60;
my $maxOverlap   = 0.15;
my $appendAnn    = 'F';
my $reference    = '';
my $equivDir     = 'equiv-Digest';

my $podUsage
    = qq(=pod\n\n)
    . qq(=head1 NAME\n\n)
    . qq($ownName - cleaning up results from matchDomains.pl\n\n)
    . qq(=head1 SYNOPSIS\n\n)
    . qq($ownName -q [matchDomainFile] [options]\n\n)
    . qq(=head1 OPTIONS\n\n)
    . qq(=over\n\n)
    . qq(=item B<-q>\n\n)
    . qq(file with domain matches from matchDomains.pl\n)
    . qq([e.g. GCF_000005845.cdd.mmseqs.bz2], required\n\n)
    . qq(=item B<-f>\n\n)
    . qq(domain family, normally [$matchDom], make it explicit if not part\n)
    . qq(of the name of the file with match results\n)
    . qq([e.g. GCF_000005845.cdd.mmseqs.bz2]\n\n)
    . qq(=item B<-o>\n\n)
    . qq(output folder, default: cleanDomains\n\n)
    . qq(=item B<-c>\n\n)
    . qq(minimum coverage of domain model, default 0.60\n\n)
    . qq(=item B<-v>\n\n)
    . qq(maximum overlap between domains, default 0.15\n\n)
    . qq(=item B<-r>\n\n)
    . qq(reference file with all scanned sequences [e.g. GCF_000005845.faa.gz]\n\n)
    . qq(=item B<-a>\n\n)
    . qq(append annotations (T|F). If 'T' the program will produce a file\n)
    . qq(with annotations appended (only works with cdd|cog|cd),\n)
    . qq(default $appendAnn\n\n)
    . qq(=back\n\n)
    . qq(=head1 DESCRIPTION\n\n)
    . qq(B<This program> extracts results from hmmscan, rpsblast, or mmseqs\n)
    . qq(normally obtained using matchDomains.pl\n\n)
    . qq(=cut\n\n)
    ;

GetOptions(
    "q=s" => \$domFile,
    "f=s" => \$domFamily,
    "o=s" => \$outputFolder,
    "c=f" => \$minCover,
    "v=f" => \$maxOverlap,
    "r=s" => \$reference,
    "a=s" => \$appendAnn,
) or podhelp();
#pod2usage(1);

if ( !$domFile ) {
    podhelp("I need a file with resuls from matchDomains.pl")
}

### the default is not to append annotations
my $appendAnn = $appendAnn =~ m{^(t|f)$}i ? uc($1) : "F";

unless( -f "$domFile" ) {
    podhelp("there's no file: $domFile");
    #pod2usage(
    #    -exitval => 0,
    #    -verbose => 2,
    #    -message => "there's no file: $domFile\n",
    #);
}
if( length("$reference") > 1 ) {
    unless( -f "$reference" ) {
        podhelp("$reference file not found");
    }
}

my $extractDir = $domFile =~ m{(\S+)/} ? $1 : "./";
my $mainName
    = $domFile =~ m{$extractDir/(\S+)\.($matchDom)}i ? $1
    : $domFile =~ m{$extractDir/(\S+?)\.} ? $1
    : $domFile =~ m{$extractDir/(\S+)}   ? $1
    : $domFile =~ m{(\S+)} ? $1
    : "none";
my $prog = $mainName =~ s{\.($mProgs)}{} ? $1 : "none";
my $infile = $domFile;
my $xfam
    = $domFamily =~ m{^($matchDom)$}i         ? uc($1)
    : $domFile   =~ m{\.($matchDom)\.}i       ? uc($1)
    : $domFile   =~ m{\.(\S+?)\.($mProgs)\.}i ? lc($1)
    : "none";
if( $xfam eq "none" ) {
    podhelp("I can only work with [$matchDom]");
    #pod2usage(
    #    -exitval => 0,
    #    -verbose => 2,
    #    -message => "I can only work with [$matchDom]\n",
    #);
}
print "  working with $xfam\n";
my $exten = lc($xfam);

my $genomeDB
    = exists $ENV{"GENOMEDB"} ? $ENV{"GENOMEDB"} : '.';
my $equivDir
    = -d "$equivDir"    ? "$equivDir"
    : -d "equiv-Digest" ? "equiv-Digest"
    : $genomeDB . "/equiv-Digest";


mkdir("$outputFolder") unless( -d "$outputFolder" );

#### COG specific
my $ncbiDir
    = -d "$genomeDB//DownLoad/ncbi" ? "$genomeDB//DownLoad/ncbi" : ".";
my $cddDir
    = -d "$ncbiDir/cdd" ? "$ncbiDir/cdd" : "none";
###### learn mini-functions for COGs and other CDDs
my %function  = ();
my %translate = ();
my %desc      = ();
my %priority  = ();
if( $xfam =~ m(^(COG|CD|CDD)$) ) {
    my $cogsDir  = $ncbiDir . "/COG/COG2014/";
    print "  using COGs directory:\n  $cogsDir\n";
    print "  learning COG functions\n";
    open( my $COGF,"<","$cogsDir/data/cognames2003-2014.tab" )
        or die "\tno COG functions file\n\n";
    while(<$COGF>) {
        next if( m{^#} );
        my($cog,$func,@descrip) = split;
        if( $cog =~ m{^COG} ) {
            $function{"$cog"} = $func;
        }
        else {
            print "\ttrouble with $_\n";
        }
    }
    close($COGF);
    ###### learn CDD identifiers
    print "  learning $xfam identifiers from the CDD database\n";
    open( my $CDDF,"-|","gzip -qdc $cddDir/cddid_all.tbl.gz" )
        or die "\tno CDD functions file\n\n";
    while(<$CDDF>) {
        chomp;
        my($cdd,$collection,$secondary,$desc,$length) = split(/\t/,$_);
        $translate{"$cdd"} = $collection;
        $translate{"$collection"} = $collection;
        if( $appendAnn eq "T" ) {
            $desc{"$collection"} = $desc;
        }
        if( $collection =~ m{^$xfam\d+}i ) {
            $priority{"$collection"}++;
        }
    }
    close($CDDF);
}

### open equivalents file to be able to match swissprot homologs
### to the organism under consideration
my($ref2ids,$ref2ln) = learnIDequivs("$mainName");

###### open domain file
###### If working with COGs substitute identifiers for COGs with
###### functions
if( $appendAnn eq "T" and $xfam =~ m{^(cog|cd|cdd)$}i ) {
    open( TRANSLATE,">","$outputFolder/$mainName.$exten.annot" );
}
my $translated = 0;
my $original_format = 0;
my %count
    = ( -f "$reference" ) ? learnReference("$reference")
    : ();
my %checked = ();
print "  working with $infile\n";
my ($cat,$trueFile) = figureCompression("$infile");
open( my $DOMF,"-|","$cat $trueFile" );
DOMLINES:
while(<$DOMF>) {
    if( m{^#} ) {
        $original_format++;
        next DOMLINES;
    }
    my( $query,$dom_id,$score,
        $qStart,$qEnd,$q_ln,
        $domStart,$domEnd,$dom_len
    ) = findFields($original_format,$_);
    my $cover = ( $domEnd - $domStart + 1 ) / $dom_len;
    next DOMLINES if( $cover < $minCover );
    if( $xfam eq "COG" ) {
        my $cog
            = length($translate{"$dom_id"}) > 1 ? $translate{"$dom_id"}
            : "NA";
        if( length($function{"$cog"}) > 0 ) {
            my $remember = join("\t",$query,$cog,$qStart,$qEnd,$score);
            push( @{$checked{"$query"}},$remember );
            $count{"$query"}++;
            $translated++;
            if( $appendAnn eq "T" ) {
                my $newLine = $_;
                my $append
                    = "\t" . $translate{"$dom_id"} . $function{"$cog"}
                    . "\t" . $desc{"$cog"};
                $newLine =~ s{\n}{$append};
                print TRANSLATE $newLine,"\n";
            }
        }
        else {
            # commented out because we know that some cogs
            # are missing in the new version for good reasons
            #print "\tproblems with $cog:\n$_";
        }
    }
    elsif( $xfam =~ m{^(CD|CDD)$} ) {
        if( length($translate{"$dom_id"}) > 2 ) {
            my $remember = join("\t",$query,$translate{"$dom_id"},
                                $qStart,$qEnd,$score);
            push( @{$checked{"$query"}},$remember );
            $count{"$query"}++;
            $translated++;
            if( $appendAnn eq "T" ) {
                my $newLine = $_;
                my $id = $translate{"$dom_id"};
                my $append
                    = "\t" . $id
                    . "\t" . $desc{"$id"};
                $newLine =~ s{\n}{$append};
                print TRANSLATE $newLine,"\n";
            }
        }
    }
    else {
        my $remember = join("\t",$query,$dom_id,$qStart,$qEnd,$score);
        push( @{$checked{"$query"}},$remember );
        $count{"$query"}++;
        $translated++;
    }
}
close($DOMF);
if( $appendAnn eq "T" ) {
    close(TRANSLATE);
}

if( $translated > 0 ) {
    if( $appendAnn eq "T" ) {
        system("bzip2 -9 $outputFolder/$mainName.$exten.annot");
    }
    open( my $TMPF,">","$outputFolder/$mainName.$exten" );
    print {$TMPF} join("\t","#Accs","${xfam}s"),"\n";
  CLEANQUERY:
    for my $query (
        sort { $ref2ids->{"$a"} cmp $ref2ids->{"$b"} || $a cmp $b }
            keys %count
        ) {
        if( $count{"$query"} == 0 ) {
            if( length( $ref2ids->{"$query"} ) > 0 ) {
                print {$TMPF} join("\t",$ref2ids->{"$query"},"NA"),"\n";
            }
            else {
                print {$TMPF} join("\t",$query,"NA"),"\n";
            }
            next CLEANQUERY;
        }
        my @lines = @{$checked{"$query"}};
        my %score = ();
        my %start = ();
        my %end   = ();
        for my $line ( @lines ) {
            my ($query,$domain,$qStart,$qEnd,$score) = split(/\s+/,$line);
            my $coords = join(":",$qStart,$qEnd);
            my $full_domain
                = ( length($function{"$domain"}) > 0 )
                ? $domain . $function{"$domain"} . "($coords)"
                : $domain . "($coords)";
            $score{"$full_domain"} = $score;
            $start{"$full_domain"} = $qStart;
            $end{"$full_domain"}   = $qEnd;
            $priority{"$full_domain"} = $priority{"$domain"};
        }
        my @acceptedCoords  = ();
        my @acceptedDomains = ();
      FULLDOMAIN:
        for my $domain (
            sort {
                $priority{"$b"} <=> $priority{"$a"}
                    || $score{"$b"} <=> $score{"$a"}
                    || $start{"$a"} <=> $start{"$b"}
                    || $end{"$a"}   <=> $end{"$b"}
                    || $a <=> $b
                    || $a cmp $b
                }
                keys %score ) {
            my $testStart = $start{"$domain"};
            my $testEnd   = $end{"$domain"};
            ##### allowing a bit of overlap:
            my $test_max_o = sprintf("%.0f",
                                  $maxOverlap * ($testEnd - $testStart + 1));
            for my $coords ( @acceptedCoords ) {
                my($clearedStart,$clearedEnd) = split(/:/,$coords);
                ##### allowing a bit of overlap:
                my $cleared_max_o
                    = sprintf("%.0f",
                              $maxOverlap * ($clearedEnd - $clearedStart + 1));
                my $acceptOverlap
                    = $cleared_max_o > $test_max_o ? $test_max_o
                    : $cleared_max_o;
                #### start of overlap will be maximum start
                my $coverStart
                    = $testStart >= $clearedStart ? $testStart
                    : $clearedStart;
                #### end of overlap will be minimal end
                my $coverEnd
                    = $testEnd <= $clearedEnd ? $testEnd : $clearedEnd;
                my $overlap
                    = $testEnd < $clearedStart ? 0
                    : $clearedEnd < $testStart ? 0
                    : $coverEnd - $coverStart + 1;
                next FULLDOMAIN if( $overlap > $acceptOverlap );
            }
            ### add coords to accepted coords
            my $addCoords = join(":",$testStart,$testEnd);
            push(@acceptedCoords,$addCoords);
            ### add domain to good, non-overlapping domains to list
            ### of accepted domains
            push(@acceptedDomains,$domain);
        }
        my @sorted_domains
            = sort { $start{"$a"} <=> $start{"$b"} } @acceptedDomains;
        my $print_domains = join(";",@sorted_domains);
        if( length( $ref2ids->{"$query"} ) > 0 ) {
            print {$TMPF} join("\t",$ref2ids->{"$query"},$print_domains),"\n";
        }
        else {
            print {$TMPF} join("\t",$query,$print_domains),"\n";
            #print "problems with $query and $domains\n";
        }
    }
    close($TMPF);
    system("bzip2 -f -9 $outputFolder/$mainName.$exten");
}
else {
    print "   $domFile has no domain results\n";
    unlink("$outputFolder/$mainName.$exten.annot");
}

print "\tdone with $0\n\n";

sub findFields {
    my ($format,$line) = @_;
    if( $xfam =~ m{^(CDD|COG|CD)$}
            || $prog =~ m{^(rpsblast|mmseqs)$} ) {
        my( $digest,$dom_id,$eval,$bitscore,
            $qStart,$qEnd,$q_ln,
            $domStart,$domEnd,$dom_len
        ) = split(/\s+/,$line);
        $digest =~ s{lcl\|}{};
        $dom_id =~ s{gnl\|CDD\|}{};
        return( $digest,$dom_id,$bitscore,
                $qStart,$qEnd,$q_ln,
                $domStart,$domEnd,$dom_len);
    }
    else {
        if( $format > 0 ) { ### it's an hmmscan file
            my( $dom_name,$dom_id,$dom_len,$digest,$caca1,$q_ln,
                $e1,$tscore,$tbias,$n1,$n2,
                $c_eval,$i_eval,$score,$bias,
                $domStart,$domEnd,
                $alnStart,$alnEnd,
                $qStart,$qEnd,$acc
            ) = split(/\s+/,$line);
            $digest =~ s{lcl\|}{};
            my $dom
                = $dom_id eq "-" ? $dom_name
                : $dom_id;
            return( $digest,$dom,$score,
                    $qStart,$qEnd,$q_ln,
                    $domStart,$domEnd,$dom_len);
        }
        else { ### it's an extracted / cleaned file
            my( $digest,$dom_id,$eval,$score,
                $qStart,$qEnd,$q_ln,
                $domStart,$domEnd,$dom_ln
            ) = split(/\s+/,$line);
            $digest =~ s{lcl\|}{};
            return( $digest,$dom_id,$score,
                    $qStart,$qEnd,$q_ln,
                    $domStart,$domEnd,$dom_ln);
        }
    }
}

#sub signalHandler {
#    print "\n\tcleaning up ...";
#    system "rm -r $tempFolder";
#    die  " done!\n\n";
#}

sub podhelp {
    my $extraMessage = $_[0];
    open( my $PIPE,"|-","less" );
    if( length "$extraMessage" > 2 ) {
        print {$PIPE} "    ",$extraMessage,"\n\n";
    }
    $parser->output_fh($PIPE);
    $parser->parse_string_document($podUsage);
    exit;
}

sub learnReference {
    my $inref = $_[0];
    my %ids   = ();
    if( -f "$inref" ) {
        my $openref
            = $inref =~ m{\.bz2}    ? "bzip2 -qdc $inref"
            : $inref =~ m{\.(gz|Z)} ? "gzip -qdc $inref"
            : "cat $inref";
        open( my $INREF,"-|","$openref" );
        while(<$INREF>) {
            if( m{^>(\S+)} ) {
                $ids{"$1"} = 0;
            }
        }
        close($INREF);
        my $countids = keys %ids;
        if( $countids > 0 ) {
            return(%ids);
        }
        else {
            die "reference file $inref had no sequence identifiers\n\n";
        }
    }
    else {
        die "I did not find the reference sequence file: $inref\n\n";
    }
}

sub learnIDequivs {
    my $genomeID = $_[0];
    my $equiv_file = $equivDir . "/" . $genomeID . ".equiv.bz2";
    if( -f "$equiv_file" ) {
        my %seq_ln = ();
        my %equivs = ();
        open( my $EQ,"-|","bzip2 -qdc $equiv_file" );
        while(<$EQ>) {
            next if( m{^#} );
            my( $digest,$seq_ln,$count,@represented ) = split;
            $seq_ln{"$digest"} = $seq_ln;
            for my $full_id ( @represented ) {
                my $represented = $full_id =~ m{gi\|(\d+)} ? $1 : $full_id;
                if( length($equivs{"$digest"}) > 0 ) {
                    $equivs{"$digest"} .= ";" . $represented;
                }
                else {
                    $equivs{"$digest"} = $represented;
                }
            }
        }
        close($EQ);
        my @keys = keys %equivs;
        my $n_keys = @keys;
        if( $n_keys > 0 ) {
            return(\%equivs,\%seq_ln);
        }
        else {
            return();
        }
    }
    else {
        return();
    }
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
