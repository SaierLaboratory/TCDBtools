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

### check if there's more than one processor or assume there's 2.
my $cpuNumber = 2;
# my $cpuNumber
#     = qx(sysctl -a | grep 'cpu.thread_count')
#         =~ m{\.cpu\.thread_count:\s+(\d+)} ? $1
#     : qx(sysctl -a 2>/dev/null | grep 'max-threads')
#         =~ m{\.max-threads\s+=\s+(\d+)} ? $1
#     : 2;

my (
    $queryFile, $outputFolder, $pwprogram, $aggMethod
) = ( '', 'Clusters', 'blastp', 'ward' );

my $options = GetOptions(
    "i=s" => \$queryFile,
    "o=s" => \$outputFolder,
    #"p=s" => \$pwprogram,
    "c=s" => \$aggMethod,
);

my $ownName = $0;
$ownName =~ s{.*/}{};
if ( !$queryFile ) {
    print "usage: " . $ownName . " [options]\n";
    print "\noptions:\n";
    print "   -i query filename in fasta format, required\n";
    print "   -o output folder, default: Clusters\n";
    #print "   -p program for pairwise comparisons\n"
    #    . "        [$pairWiseMatch], default: blastp\n";
    print "   -c clustering method\n"
        . "        [$aggMethodsMatch],\n"
        . "         default: ward\n";
    print "\n";
    print "requirements:\n"
        . qq(   ) . qq(blastp from NCBI's blast suite\n)
        . qq(   ) . qq(R from r-project.org\n)
        . qq(   ) . qq(   ) . qq(R packages: cluster, MCMCpack and ape\n\n);
    exit;
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

my $tempFolder = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
print "   working in temp folder:\n\t$tempFolder\n";
unless( -d "$outputFolder" ) {
    mkdir("$outputFolder");
}

my $family = $queryFile;
$family =~ s{.+\/}{};
$family =~ s{\.\w+$}{};
print "    root name is $family\n";

##### learn all the identifiers
my @ids = ();
open( my $FAA,"<","$queryFile" )
    or die "\tno $queryFile file\n\n";
while(<$FAA>) {
    if( m{^>(\S+)} ) {
        my $id = $1;
        $id =~ s{lcl\|}{};
        if( $id =~ m{\.\d+\.} ) {
            push(@ids,$id);
        }
    }
}
close($FAA);

### establish default score at 10 (arbitrary, coming from SuperFamilyTree)
for my $id1 ( @ids ) {
    for my $id2 ( @ids ) {
        my $pair = join("\t",sort($id1,$id2));
        $score{"$pair"} = 10;
    }
}

### now learn actual scores
unless( -f "$outputFolder/$family.$pwprogram.bz2" ) {
    runPairWise("$queryFile","$pwprogram","$family");
}
open( my $BLASTR,"-|","bzip2 -qdc $outputFolder/$family.$pwprogram.bz2" );
BLASTR:
while(<$BLASTR>) {
    next if( m{^#} );
    my($id1,$id2,@stats) = split;
    $id1 =~ s{lcl\|}{};
    $id2 =~ s{lcl\|}{};
    my $pair = join("\t",sort($id1,$id2));
    if( exists($score{"$pair"}) && ( $stats[-1] > $score{"$pair"} ) ) {
        $score{"$pair"} = $stats[-1];
    }
}
close($BLASTR);

### now save distances into R file
print "    calculating distances\n";
open( my $RFILE,"|-","bzip2 -9 > $outputFolder/$family.$pwprogram.pairwise.bz2" );
print {$RFILE} join("\t","ID1","ID2","Dist"),"\n";
for my $pair ( sort keys %score ) {
    my $dist = sprintf("%.6f",100/$score{"$pair"});
    my( $id1,$id2 ) = split(/\t/,$pair);
    print {$RFILE} join("\t",$id1,$id2,$dist),"\n";
}
close($RFILE);

## @aggMethods
print "    clustering\n";
for my $clust_method ( $aggMethod ) {
    build_clust_R("$clust_method");
    my $Rcmd = qq(R --vanilla < )
        . qq($outputFolder/$family.$pwprogram.$clust_method.R )
        . qq(>& $outputFolder/$family.$pwprogram.$clust_method.log);
    print "      running R\n";
    system("$Rcmd");
    print "    The cluster has been saved in Newick format:\n"
        . "       $outputFolder/$family.$pwprogram.$method.tree\n";
}

sub build_clust_R {
    my($method) = @_;
    print "      working with $method"
        . " ($outputFolder/$family.$pwprogram.$method.R)\n";
    open( my $CLUSTR,">","$outputFolder/$family.$pwprogram.$method.R" );
    print {$CLUSTR} << "CLUSTR";
### this R script will read the distance matrix and produce
### a cluster with the $method clustering method
library(cluster)
library(MCMCpack)
library(ape)

fampairw <- read.table("$outputFolder/$family.$pwprogram.pairwise.bz2",
        header=T)

#### make into a matrix
number.fam<-length(levels(fampairw[,2]))
fam.names<-fampairw[1:number.fam,2]
fammat<-xpnd(fampairw[,3])
rownames(fammat)<-fam.names
colnames(fammat)<-fam.names

#### save matrix
write.table(fammat,
        bzfile("Clusters/$family.$pwprogram.matrix.bz2",compression=9),
        row.names=T,col.names=T,quote=F)

## create a dendrogram using the agglomerative (bottom up) approach
## and the distance matrix
aggTree<-agnes(as.dist(fammat),diss=T,method="$method")
aggDend<-as.hclust(aggTree)

## Output the tree in the Newick format
phy <- as.phylo(aggDend)
write.tree(phy, file="$outputFolder/$family.$pwprogram.$method.tree")

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
    my $tmpFile = "$tempFolder/$familyName.$program.bz2";
    my $outFile = "$outputFolder/$familyName.$program.bz2";
    my $runCmd
        = qq(blastp -query $fastaFile -subject $fastaFile)
        . qq( -outfmt 7 -max_hsps 1 -parse_deflines )
        . qq( | bzip2 -9 > $tempFolder/$familyName.$program.bz2);
    print "    running $program\n";
    system("$runCmd");
    system("mv $tmpFile $outFile 2>/dev/null");
}

sub signalHandler {
    print "\n\tcleaning up ...";
    system "rm -r $tempFolder";
    die  " done!\n\n";
}
