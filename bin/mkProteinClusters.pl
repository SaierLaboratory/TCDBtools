#!/usr/bin/perl

#########################################################################
#                                                                       #
#       Author: Gabo Moreno-Hagelsieb                                   #
#       Date first raw version: Aug 11, 2017                            #
#                                                                       #
#########################################################################
use Getopt::Long;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);
use strict;

### this is to be able to add an option for the agglomerative
### clustering method. Right now only ward because it's the very best
### (as far as our experience indicates).
my @aggMethods = qw(
                        average
                        complete
                        single
                        ward
                        weighted
                );
my $aggMethodsMatch = join("|",@aggMethods);

### This is to add an option for choosing the program to compare the
### sequences, thus obtaining the scores that will later be
### transformed into distances for clustering
my @pairWise = qw(
                     blastp
                     fasta36
                     ssearch36
             );
my $pairWiseMatch = join("|",@pairWise);


### Amino acid substitution matrix to be used with ssearch36
my $ssearchMatrix    = undef; #ssearch36 option -s
my $ssearchAlgorithm = undef; #ssearch36 option -z
my $sserachShuffles  = undef; #ssearch36 option -k


### check if there's more than one processor or assume there's 2.
my $cpuNumber = 2;
#my $cpuNumber
#     = qx(sysctl -a | grep 'cpu.thread_count')
#         =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
#     : qx(sysctl -a 2>/dev/null | grep 'max-threads')
#         =~ m{\.max-threads\s+=\s+(\d+)} ? $1
#     : 2;

my $queryFile    = '';
my $outputFolder = 'Clusters';
my $pwprogram    = 'ssearch36';
my $aggMethod    = 'ward';

my $options = GetOptions(
    "i=s" => \$queryFile,
    "o=s" => \$outputFolder,
    "p=s" => \$pwprogram,
    "c=s" => \$aggMethod,
    "s=s" => \$ssearchMatrix,
    "z=s" => \$ssearchAlgorithm,
    "k=n" => \$sserachShuffles,
);

my $ownName = $0;
$ownName =~ s{.*/}{};
if ( !$queryFile ) {
    print "usage: " . $ownName . " [options]\n";
    print "\noptions:\n";
    print "   -i query filename in fasta format, required\n";
    print "   -o output folder. Default: Clusters\n";
    print "   -p program for pairwise comparisons\n"
      . "        [$pairWiseMatch]. Default: ssearch36\n";
    print "   -s Amino acid substitution matrix that wil be used\n";
    print "      by ssearch36. Default: ssearch36 default for option -s\n";
    print "   -z Algorithm to be used by ssearch36 to calculate E-values,\n";
    print "      default: ssearch36 default for option -z \n";
    print "   -k Number of shuffles to be used by ssearch36 in the\n";
    print "      calculation of E-values. Default: ssearch36 default for option -k\n";
    print "   -c agglomerative clustering method\n"
        . "        [$aggMethodsMatch].\n"
        . "         Default: ward\n";
    print "\n";
    print "requirements:\n"
        . qq(   ) . qq(blastp from NCBI's blast suite\n)
        . qq(   ) . qq(fasta36, ssearch36 from Pearson's fasta36 suite\n)
        . qq(   ) . qq(R from r-project.org\n)
        . qq(   ) . qq(   ) . qq(R packages: cluster, MCMCpack and ape\n\n);
    exit;
}


### Prepare the arguments for ssearch36. If no -z, -k, -s options are given,
### alignments will use the pre-established defaults of ssearch36.
my $ssearchOptions = "";
if ($ssearchMatrix) {
    $ssearchOptions .= "-s $ssearchMatrix ";
}
if ($ssearchAlgorithm) {
    $ssearchOptions .= "-z $ssearchAlgorithm ";
}
if ($sserachShuffles) {
    $ssearchOptions .= "-k $sserachShuffles ";
}


### test that the method is spelled correctly:
if( $aggMethod =~ m{^($aggMethodsMatch)$}i ) {
    $aggMethod = lc($1);
    print "   using $aggMethod agglomerative method\n";
}
else {
    die "   the agglomerative method [$aggMethod] does not exist\n"
        . "   try any of [$aggMethodsMatch] (default: ward)\n\n";
}


### test that the method is spelled correctly:
if( $pwprogram =~ m{^($pairWiseMatch)$}i ) {
    $pwprogram = lc($1);
    print "   using $pwprogram for pairwise alignments\n";
}
else {
    die "   the program [$pwprogram] does not exist\n"
        . "   try any of [$pairWiseMatch] (default: blastp)\n\n";
}

my %cmd = ();
$cmd{"blastp"}
    = qq(blastp -outfmt 7 -max_hsps 1 -parse_deflines -use_sw_tback -out );
$cmd{"fasta36"}   = qq(fasta36   -m 8C );
$cmd{"ssearch36"} = qq(ssearch36 -m 8C $ssearchOptions );

my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
print "   working in temp folder:\n\t$tempFolder\n";
unless( -d "$outputFolder" ) {
    system("mkdir -p $outputFolder");
}

my $family = $queryFile;
$family =~ s{.+\/}{};
$family =~ s{\.\w+$}{};
print "   root name is $family\n";

##### learn all the identifiers
print "      learning IDs\n";
my @ids = ();
open( my $FAA,"<","$queryFile" )
    or die "\tno $queryFile file\n\n";
while(<$FAA>) {
    if( m{^>(\S+)} ) {
        my $id = $1;
        $id =~ s{lcl\|}{};
        push(@ids,$id);
    }
}
close($FAA);

##### establish default score at 10 (arbitrary,
##### taken from SuperFamilyTree)
my %score = ();
print "      establishing default distance\n";
for my $id1 ( @ids ) {
    for my $id2 ( @ids ) {
        my $pair = join("\t",sort($id1,$id2));
        $score{"$pair"} = 10;
    }
}

### now learn actual scores
### to check that we have pairwise alignment results
my $pwFile = "$outputFolder/$family.$pwprogram.bz2";
my $countResults = 0;
unless( -f "$pwFile" ) {
    runPairWise("$queryFile","$pwprogram","$family");
}
open( my $BLASTR,"-|","bzip2 -qdc $pwFile" );
BLASTR:
while(<$BLASTR>) {
    next if( m{^#} );
    chomp;
    my($id1,$id2,@stats) = split;
    $id1 =~ s{lcl\|}{};
    $id2 =~ s{lcl\|}{};
    my $pair = join("\t",sort($id1,$id2));
    if( exists($score{"$pair"}) ) {
        $countResults++;
        if( $stats[-1] > $score{"$pair"} ) {
            $score{"$pair"} = $stats[-1];
        }
    }
}
close($BLASTR);
if( $countResults > 0 ) {
    print "   read $countResults pairwise scores\n";
}
else {
    print "Deleting $pwFile because it seems to be empty\n";
    unlink("$pwFile");
    print "please check your fasta file and try again\n";
    signalHandler();
}

### now save distances into R file
print "   calculating distances\n";
open( my $RFILE,"|-","bzip2 -9 > $outputFolder/$family.$pwprogram.distances.bz2" );
print {$RFILE} join("\t","ID1","ID2","Dist"),"\n";
for my $pair ( sort keys %score ) {
    my $dist = sprintf("%.6f",100/$score{"$pair"});
    my( $id1,$id2 ) = split(/\t/,$pair);
    print {$RFILE} join("\t",$id1,$id2,$dist),"\n";
}
close($RFILE);

## @aggMethods
print "   clustering\n";
for my $clustMethod ( $aggMethod ) {
    my $clustFile = "$outputFolder/$family.$pwprogram.$clustMethod.nw.tree";
    build_clust_R("$clustMethod","$clustFile");
    my $Rcmd = qq(R -f )
        . qq($outputFolder/$family.$pwprogram.$clustMethod.R )
        . qq(>& $outputFolder/$family.$pwprogram.$clustMethod.log);
    print "      running R\n";
    system("$Rcmd");
    print "   The cluster has been saved in Newick format:\n"
        . "      $clustFile\n";
}

sub build_clust_R {
    my($method,$outFile) = @_;
    print "      working with $method"
        . " ($outputFolder/$family.$pwprogram.$method.R)\n";
    open( my $CLUSTR,">","$outputFolder/$family.$pwprogram.$method.R" );
    print {$CLUSTR} << "CLUSTR";
### this R script will read the distance matrix and produce
### a cluster with the $method clustering method
library("cluster")
library("MCMCpack")
library("ape")

fampairw <- read.table("$outputFolder/$family.$pwprogram.distances.bz2",
        header=T,sep="\\t")

#### make into a matrix
number.fam<-length(levels(as.factor(fampairw[,2])))
fam.names<-fampairw[1:number.fam,2]
fammat<-xpnd(fampairw[,3])
rownames(fammat)<-fam.names
colnames(fammat)<-fam.names

#### save matrix
write.table(fammat,
        bzfile("$outputFolder/$family.$pwprogram.matrix.bz2",compression=9),
        row.names=T,col.names=T,quote=F)

## create a dendrogram using the agglomerative (bottom up) approach
## and the distance matrix
aggTree<-agnes(as.dist(fammat),diss=T,method="$method")
aggDend<-as.hclust(aggTree)

## Output the tree in the Newick format
phy <- as.phylo(aggDend)
write.tree(phy, file="$outFile")

## calculate the agglomerative coefficient
coeff.agg<-coef(aggTree)
coeff.agg
coeff.print<-formatC(coeff.agg,digits=3,format="f")
write(paste("$method", coeff.print),
      "$outputFolder/$family.$pwprogram.coeffs",sep = "\\t", append=T)

CLUSTR
    close($CLUSTR);
}

print "\tcleaning up ...";
system "rm -r $tempFolder";
print  "\n        $ownName done!\n\n";

sub runPairWise {
    my ($fastaFile,$program,$familyName) = @_;
    my $tmpFile = "$tempFolder/$familyName.$program";
    my $outFile = "$outputFolder/$familyName.$program.bz2";
    my $rootCmd = $cmd{"$pwprogram"};
    my $pwCmd
        = $pwprogram eq "blastp"
        ? $rootCmd . "$tmpFile -query $fastaFile -subject $fastaFile"
        : $pwprogram eq "fasta36"
        ? $rootCmd . "$fastaFile $fastaFile > $tmpFile"
        : $pwprogram eq "ssearch36"
        ? $rootCmd . "$fastaFile $fastaFile > $tmpFile"
        : "none";
    if( $pwCmd eq "none" ) {
        print "there's something wrong with the $pwprogram\n";
        signalHandler();
    }
    else {
        print "      running $program\n";
        print "      $pwCmd\n";
        system("$pwCmd");
        system("bzip2 -f -9 $tmpFile");
        system("mv $tmpFile.bz2 $outFile 2>/dev/null");
    }
}

sub signalHandler {
    print "\n\tcleaning up ...";
    system "rm -r $tempFolder";
    die  " done!\n\n";
}
