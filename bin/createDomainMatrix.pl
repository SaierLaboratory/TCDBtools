#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

use Set::Scalar;
use TCDB::Assorted;

use Getopt::Long;
use XML::Writer;

#==========================================================================
#                           Approach 1
#
# Based on comparisons of domain content between pairs of families,
# build a matrix to determine which families might actually form
# a super family.
#
# The idea is to create an approach independent from the protocol1/protocol2
# approach, to strengthen the hypothesis of homology between pairs of
# families
#
#==========================================================================




#==========================================================================
#Global variables

#Query families for which the matrix will be built
my $refFams    = [];
my $cmpFams    = [];
my $refFile    = "";
my $cmpFile    = "";
my $inDir      = undef;
my $outPrefix  = "matrix";
my $format     = 'matrix';
#my $clanFile   = "/Users/amedrano/Projects/tcdb/MFS/frozenDatabases/Pfam-A.clans.tsv.gz";
my $clanFile   = "/Users/amedrano/Projects/tcdb/TOG/Domains/Pfam-A.clans.tsv.gz";


my %evidenceRank = ('noHit' => 0, 'Rescued2' => 1, 'Rescued1' => 2, 'DirectHit' => 3);


#Read command line topology
read_command_line_arguments();


#print Data::Dumper->Dump([$refFams, $cmpFams, $refFile, $cmpFile, $inDir, $outPrefix, $clanFile],
#			 [qw(*refFams *cmpFmas *refFile *cmpFile *inDir *outPrefix *clanFile)]);
#exit;



#==========================================================================
#Read the clan to PFAM translations

my %pfam2clan = ();
parse_clans_file (\%pfam2clan);

#print Data::Dumper->Dump([ \%pfam2clan ], [qw(*pfam2clan)]);
#exit;




#==========================================================================
#Parse the Pfam Domains identified/rescued between each pair of families,
#and extract the sets of shared domains and their clans


my %matrixDoms  = ();
my %matrixClans = ();

generateMatrix();

#print Data::Dumper->Dump([\%matrixDoms, \%matrixClans], [qw(*matrixDoms *matrixClans)]);
#exit;


#==========================================================================
#Save matrices of Domains and clans to a files


my $refSet = Set::Scalar->new(@$refFams);
my $cmpSet = Set::Scalar->new(@$cmpFams);

if ($refSet == $cmpSet) {
  #save_symetrical_matrices_to_disk();
  save_ref_vs_cmp_matrces_to_disk();
}
else {
  save_ref_vs_cmp_matrces_to_disk();
}





###########################################################################
##                         Functions                                     ##
###########################################################################



#==========================================================================
#Save the results of the comparison matrix when symmetrical

sub save_ref_vs_cmp_matrces_to_disk {

  my $clanFile = "./${outPrefix}_clans_asym.tsv";

  if ($format eq 'matrix') {
    print_ref_vs_cmp_matrix();
  }
  elsif ($format eq 'gephi') {
    print_ref_vs_cmp_gephi();
  }

}


#==========================================================================
#Print ref vs cmp in gephi format


sub print_ref_vs_cmp_gephi {

  my $clanFile = "./${outPrefix}_clans_gephi.xml";

#  print Data::Dumper->Dump([\%matrixClans], [qw(*matrixClans )]);
#  exit;

#  print Data::Dumper->Dump([ $matrixClans{'2.A.20'}{'2.A.58'},  $matrixClans{'2.A.58'}{'2.A.20'}],
#			   [qw(*20_58 *58_20)]);
#  exit;


  my $edgeCnt = 0;

  #----------------------------------------------------------------------
  #Get the inital order of the nodes (sorted by family)


  my @all_fams = (@$refFams,  @$cmpFams);
  my %nodesTmp1  = map {$_ => 1} @all_fams;
  my @nodesTmp2 = sort_by_family (keys %nodesTmp1);

  my %nodeOrder = ();
  for (my $i = 0; $i <= $#nodesTmp2; $i++) {
    my $f = $nodesTmp2[$i];
    $nodeOrder{$f} = $i;
  }

#  print Data::Dumper->Dump([\%nodeOrder ], [qw( *nodeOrder )]);
#  exit;





  #--------------------------------------------------------------------------
  #Generate data for gephi

  my %nodes = ();
  my %edgesp = ();
  my @edges  = ();

  #Get the nodes and edges for the graph
 REF:foreach my $idx1 (0 .. $#all_fams - 1) {

    my $ref = $all_fams[$idx1];

  CMP:foreach my $idx2 ($idx1 + 1 .. $#all_fams) {

      my $cmp = $all_fams[$idx2];

      next CMP if ($ref eq $cmp);

      #get the shared clans
      my $clanStr = undef;

      #Determine if there will be an edge between this pair of families
      if (exists $matrixClans{$ref} && exists $matrixClans{$ref}{$cmp}) {
	$clanStr = $matrixClans{$ref}{$cmp};
      }
      elsif (exists $matrixClans{$cmp} && exists $matrixClans{$cmp}{$ref}) {
	$clanStr = $matrixClans{$cmp}{$ref};
      }
      else {
	next CMP;
      }

      next CMP if ($clanStr eq 'None');


      #Store node information
      unless (exists $nodes{$ref}) {
	$nodes{$ref} = $nodeOrder{$ref};
      }

      unless (exists $nodes{$cmp}) {
	$nodes{$cmp} = $nodeOrder{$cmp};
      }


      #store edge inforamtion
      unless ((exists $edgesp{$ref} && exists $edgesp{$ref}{$cmp}) ||
	      (exists $edgesp{$cmp} && exists $edgesp{$cmp}{$ref})) {
	$edgesp{$ref}{$cmp} = 1;

	push (@edges, {id => $edgeCnt, source=>$nodes{$ref}, target=>$nodes{$cmp}, label=>$clanStr});
	$edgeCnt++;
      }
    }
  }

#  exit;
#  print Data::Dumper->Dump([\%nodes, \@edges], [qw(*nodes *edges )]);
#  exit;





  #Print matrix in gexf form
  my $writer = XML::Writer->new( OUTPUT => 'self', DATA_MODE => 1);
  $writer->xmlDecl("UTF-8");

  #The start of the document
  $writer->startTag('gexf', 'xmlns' => 'http://www.gexf.net/1.2draft',
		    'xmlns:xsi' => "http://www.w3.org/2001/XMLSchema-instance",
		    'xsi:schemaLocation' => "http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd",
		    'version' => '1.2');

  #Author info
  $writer->startTag('meta', 'lastmodifieddate' => '2018-03-19');
  $writer->dataElement( creator => 'Arturo Medrano' );
  $writer->dataElement( description => 'TOG network data for Gephi' );
  $writer->endTag('meta');


  #The graph starts here
  $writer->startTag('graph', mode => 'static', defaultedgetype => 'undirected');


  #The attributes for nodes and edges
#  $writer->startTag('attributes', 'class' => 'edge');
#  $writer->emptyTag('attribute', id => '0', title => 'domains', type => 'string');
#  $writer->endTag('attributes');


  #----------------------------------------------------------------------
  #The nodes in the graph

  $writer->startTag('nodes');

  #foreach my $node (sort { $nodes{$a} <=> $nodes{$b} } keys %nodes) {
  foreach my $node (sort { $nodeOrder{$a} <=> $nodeOrder{$b} } keys %nodeOrder) {
    #$writer->emptyTag('node', id=>$nodes{$node}, label => $node);
    $writer->emptyTag('node', id=>$nodeOrder{$node}, label => $node);
  }

  $writer->endTag('nodes');



  #----------------------------------------------------------------------
  #The edges in the graph

  $writer->startTag('edges');

  foreach my $edge (sort { $a->{id} <=> $b->{id} } @edges) {
    $writer->emptyTag('edge', id => $edge->{id}, source => $edge->{source}, target => $edge->{target}, label => $edge->{label});
  }

  $writer->endTag('edges');


  #Colosing tags for graph and file
  $writer->endTag('graph');
  $writer->endTag('gexf');

  print $writer->to_string;

}





#==========================================================================
#print matrix format

sub print_ref_vs_cmp_matrix {

    my $clanFile = "./${outPrefix}_clans_asym.tsv";


  #Now generate the comparison matrix (just one half);
  open (my $outh, ">", $clanFile) || die $!;


  #print header
  print $outh "-\t", join("\t", @$cmpFams), "\n";

 REF:foreach my $ref (@$refFams) {

    #reference family
    print $outh "$ref";

  CMP:foreach my $cmp (@$cmpFams) {

      #Identifiy main diagonal
      if ($ref eq $cmp) {
	print $outh "\t-";
	next CMP;
      }


      #Right of diagonal
      else {

	#get the shared clans
	my $clanStr = $matrixClans{$ref}{$cmp};


	if ($clanStr eq "None") {
	  print $outh "\t";
	}
	else {
	  print $outh "\t$clanStr";
	}
      }
    }

    print $outh "\n";
  }

  close $outh;
}





#==========================================================================
#Save the results of the comparison matrix when symmetrical


sub save_symetrical_matrices_to_disk {

  my $clanFile = "./${outPrefix}_clans_sym.tsv";


  #Now generate the comparison matrix (just one half);
  open (my $outh, ">", $clanFile) || die $!;


  #print header
  print $outh "-\t", join("\t", @$cmpFams), "\n";

 REF:foreach my $ref (@$refFams) {

    #reference family
    print $outh "$ref";

    my $left_of_diagonal = 1;

  CMP:foreach my $cmp (@$cmpFams) {

      #Identifiy main diagonal
      if ($ref eq $cmp) {
	print $outh "\t-";
	$left_of_diagonal = 0;
	next CMP;
      }


      #If left of diagonal just print a tab
      if ($left_of_diagonal) {
	print $outh "\t";
      }


      #Right of diagonal
      else {

	#get the shared clans
	my $clanStr = $matrixClans{$ref}{$cmp};


	if ($clanStr eq "None") {
	  print $outh "\t";
	}
	else {
	  print $outh "\t$clanStr";
	}
      }
    }

    print $outh "\n";
  }

  close $outh;
}


sub print_symetrical_matrix_gephi {


}

#==========================================================================
#Generate de Matrix with the domain comparisons between pairs of proteins


sub generateMatrix {


 FAM1:foreach my $f1 (@$refFams) {
  FAM2:foreach my $f2 (@$cmpFams) {

      next FAM2 if ($f1 eq $f2);

      my %sharedDomains = ();
      my %sharedClans   = ();
      identifySharedDomains(\%sharedDomains, \%sharedClans, $f1, $f2);


      #Add shared domains to matrix
      if (%sharedDomains) {
	$matrixDoms{$f1}{$f2} = join("|", keys %sharedDomains);
	$matrixDoms{$f2}{$f1} = join("|", keys %sharedDomains);
      }
      else {
	$matrixDoms{$f1}{$f2} = "None";
	$matrixDoms{$f2}{$f1} = "None";
      }


      #Add clans involved
      if (%sharedClans) {

	my @clans = ();
	#append domains to respective clans
	foreach my $c (keys %sharedClans) {
	  my $d = join ("|", @{ $sharedClans{$c} });
	  push (@clans, "$c|$d");
	}

	$matrixClans{$f1}{$f2} = join (":", @clans);
	$matrixClans{$f2}{$f1} = join (":", @clans);
      }
      else {
	$matrixClans{$f1}{$f2} = "None";
	$matrixClans{$f2}{$f1} = "None";
      }
    }
  }
}



#==========================================================================
#Given a pair of families identify what domains are shared, do it at
#the level of Pfam IDs and Clans


sub identifySharedDomains {

  my ($commonDoms, $clans, $f1, $f2) = @_;

  #Directories for this pair of families
  my $workDir1 = "$inDir/${f1}_vs_$f2";
  my $workDir2 = "$inDir/${f2}_vs_$f1";
  die "Error: directory not found -> $workDir1" unless (-d $workDir1);
  die "Error: directory not found -> $workDir2" unless (-d $workDir2);


  #Domains shared between $f1 and $f2
  my @dom1 = findCommonDomains($f1, $f2, $workDir1);
  my @dom2 = findCommonDomains($f2, $f1, $workDir2);


  #Compute the union of domains
  %{ $commonDoms } = map { $_ => 1 } @dom1, @dom2;


  if (%{ $commonDoms }) {

    #Extract the clans for shared domains
    foreach my $dom (keys %{ $commonDoms }) {
      if (exists $pfam2clan{$dom}) {
	push (@{ $clans->{ $pfam2clan{$dom} }}, $dom);
      }
      else {
	push (@{ $clans->{ NA }}, $dom);
      }
    }
  }

#  print Data::Dumper->Dump([\@dom1, \@dom2, $commonDoms, $clans ], [qw(*dom1 *dom2 *commonDoms *clans)]);
#  exit;
}




#==========================================================================
#Given a pair of families in a specific order, identify shared domains
#and clans.

sub findCommonDomains {

  my ($rfam, $cfam, $dir) = @_;

  my $rescueFile = "$dir/$rfam/reports/${rfam}_rescuedDomains.tsv";
  die "Error: no rescue file --> $rescueFile" unless (-f  $rescueFile && !(-z $rescueFile));


  #-----------------------------------------------------------------
  #parse Domain inferences for this file

  my %res = ();

  #Parse rescue file
  open (my $fh, "<", $rescueFile) || die $!;
  while (<$fh>) {

    chomp;
    next if (/^#/);

    my ($hit, $prot, @domainData) = split (/\s+/);

#   print Data::Dumper->Dump([$hit, $prot, \@domainData], [qw(*hit *prot *domainData)]);
#   <STDIN>;


    #identify family with hit
    my $family = ($prot =~ /$rfam/)? $rfam : $cfam;



    foreach my $pfam (@domainData) {

      my @components = split(/\|/, $pfam);

      #Ignore domains without hits
      next if ($components[-1] eq 'Nohit');

      #determine the highest reliability for this domain
      my $r1    = $evidenceRank{ $components[-1] };


      #keep the highest rank possible for this domain in this family
      if (exists $res{ $family } && exists  $res{ $family }{ $components[0] }) {

	my $r2 = $res{ $family }{ $components[0] };
	$res{ $family }{ $components[0] } = $r1 if ($r1 > $r2);
      }
      else {
	$res{$family}{$components[0]} = $r1;
      }
    }


  }
  close $fh;


#  print Data::Dumper->Dump([\%res ], [qw(*res )]);
#  <STDIN>;



  #-----------------------------------------------------------------
  #identify shared domains

  my @shared = grep { $res{$rfam}{$_} } keys %{ $res{$cfam} };

  return @shared;
}






#==========================================================================
#Parse Pfam clans file


sub parse_clans_file {

  my $c2p = shift;

  open (my $fh, "zcat $clanFile |") || die $!;
 LINE:while (<$fh>) {
    chomp;
    my ($pfam, $clan, @names) = split(/\t/);

    #There must be a clan, otherwise ignore line
    next LINE unless ($clan);

    $c2p->{$pfam} = $clan;
  }
  close $fh;
}






#==========================================================================
#Read command line arguments


sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "rf|ref-fams=s"     => \&readRefFams,      #TCIDs that will be used as reference
      "cf|cmp-fams=s"     => \&readCmdFams,      #TCIDs that will that will compared
      "rfile|ref-file=s"  => \&readRefFile,      #File with reference TCIDs
      "cfile|cmp-file=s"  => \&readCmpFile,      #File with comparison TCIDs
      "c|pfam-clans=s"    => \&read_clans_file,  #PFAM to CLAN mappings in TCDB
      "d|indir=s"         => \&read_indir,       #Directory with all pairwise domain comparisons
      "p|prefix=s"        => \&read_outPrefix,   #Ouput root directory
      "f|format=s"        => \&read_format,      #output format (matrix,gephi,csv)
      "h|help"  => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);

  #----------------------------------------------------------------------
  #Validate command line arguments


  die "Error: input directory is mandatory!\n" unless ($inDir && -d $inDir);

  #Download clans if no imput file was given
  unless (-f $clanFile) {
    system "wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz";
    $clanFile = "./Pfam-A.clans.tsv.gz";
  }

  die "Error: Pfam clans file not found: $clanFile!\n" unless (-f $clanFile);


  #----------------------------------------------------------------------
  #Validate command line arguments


  if (@$refFams && $refFile) {
    die "Error: options -rf and -rfile are mutually exclusive. please, provide only one\n";
  }

  if (@$cmpFams && $cmpFile) {
    die "Error: options -cf and -cfile are mutually exclusive. please, provide only one\n";
  }

  #At this point is guaranteed that only of ref and/or cmp families were given
  $refFams = TCDB::Assorted::readFamilyFile($refFile) if ($refFile);
  $cmpFams = TCDB::Assorted::readFamilyFile($cmpFile) if ($cmpFile);


  die "Error: at least -rf, -cf, -rFile, or -cFile must be given !\n" unless (@$refFams || @$cmpFams);


  #when only one type of families are given
  if (@$cmpFams && !(@$refFams)) {
    $refFams = $cmpFams;
  }
  elsif (@$refFams && !(@$cmpFams)) {
    $cmpFams = $refFams;
  }

}


#==========================================================================
#Read the -rf option

sub readRefFams {

  my ($opt, $value) = @_;

  $refFams = TCDB::Assorted::readStringFams($value);
}



#==========================================================================
#Read the -cf option

sub readCmpFams {

  my ($opt, $value) = @_;

  $cmpFams = TCDB::Assorted::readStringFams($value);
}



#==========================================================================
#Read the -rfile option

sub readRefFile {

  my ($opt, $value) = @_;

  $refFile = $value;
}



#==========================================================================
#Read the -cfile option

sub readCmpFile {

  my ($opt, $value) = @_;

  $cmpFile = $value;
}




#==========================================================================
#Read the -c option

sub read_clans_file {

  my ($opt, $value) = @_;

  die "Error: file with PFAM clans not found!\n" unless (-f $value && !(-z $value));

  $clanFile = $value;
}





#==========================================================================
#Read the -d option

sub read_indir {

  my ($opt, $value) = @_;

  die "Error: input directory must exist" unless (-d $value);

  $inDir = $value;
}



#==========================================================================
#Read the -o option

sub read_outPrefix {

  my ($opt, $value) = @_;

  $outPrefix = $value;
}


#==========================================================================
#Read output format option (-f)


sub read_format {

  my ($opt, $value) = @_;

  unless ($value =~ /^(matrix|gephi|csv)$/) {
    die "Error: not a valied ouput format => $value\n";
  }


  $format = $value;

}




#==========================================================================
#Print help

sub print_help {

  my $help = <<'HELP';

Based on the domains shared between pairs of proteins, build a matrix
indicating pairs of families that shared domains or the domains are in
the same clan.


-rf, --ref-fams {string} (optional)
  List of comma-separated tcdb families that will be used as reference.
  This options is imcompatible with -rfile.

-cf, --cmp-fams {string} (optional)
  List of comma-separated tcdb familes that will be used to compare against
  the set of families given through option -rf. This option is incompatible
  with -cfile.

-rfile, --ref-file {file} (optional)
  File with families that will be used a reference (i.e. --ref-fams). There
  must be only one family per line.

-cfile, --cmp-file {file} (optional)
  File with families that will be compared against the reference set
  (i.e. --cmp-familes). There must be only one family per line.

 -c, --pfam-clans {file} (Optional)
  File with all the clans to PFAM IDs mappings. If this argument is not given,
  and attempt will be made to download the file from the pfam website:
  ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz

 -d, --indir {directory} (Mandatory)
  Directory with all the domain comparisons between pairs of families. These
  directories must contains the results of running script:
            getDomainTopology.pl

 -p, --prefix {string} (optional)
  Prefix that will be included in the output file names.
  Default:  matrix

HELP


  print $help;
  exit;
}
