#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;
#use List::Util qw(sum);

use TCDB::Assorted;
use TCDB::Domain::PfamParser;
use TCDB::Domain::Characterize;

use Getopt::Long;

#==========================================================================
#Global variables

#Query family or families
my @fams = ();

#This is an option for TCDB::Assorted::getSystemAccessions()
my $treatAsSuperfamily   = 0;

#Options for TCDB::Domain::PfamParser
my $domain_cov = 0.7;
my $prot_cov   = 0.1;
my $evalue     = 1e-5;
my $prop_prots_w_domain = 0.05;

#Options for TCDB::Domain::Characterize
my $rootDir              = ".";

#To extract the TCIDs of refernece families
my $tcdbSeqsFile         = "$ENV{RESEARCH_DATA}/pfam/download/tcdb.faa";
my $pfamFile             = "$ENV{RESEARCH_DATA}/pfam/tcdb.pfam-a.hmmscan.bz2";
my $blastdb              = "$ENV{HOME}/db/blastdb/tcdb";
my $prog                 = "ssearch36";
my @candProjProts        = ();
my $analysisLevel        = 'system';



#Read command line topology
read_command_line_arguments();


die "TCDB sequences file not found or empty --> $tcdbSeqsFile\n" unless (-f $tcdbSeqsFile && !(-z $tcdbSeqsFile));
die "TCDB hmmscan output file not found --> $pfamFile\n" unless (-f $pfamFile && !(-z $pfamFile));


#print Data::Dumper->Dump([\@fams, $treatAsSuperfamily,  $rootDir, $tcdbSeqsFile, $pfamFile, $blastdb, $prog, $domain_cov,
#			  $prot_cov, $evalue, $prop_prots_w_domain, \@candProjProts],
#			 [qw(*fams *treatAsSuperfamily *rootDir *tcdbSeqsFile *pfamFile *blastdb *prog *domain_cov
#			     *prot_cov *evalue *prop_prots_w_domain *candProjProts)]);
#exit;


#==========================================================================
#Split tcdb systems into single-component multi-component


if ($treatAsSuperfamily) {

    my $tcids = getSystemAccessions($tcdbSeqsFile, 'both', $analysisLevel, \@fams, $treatAsSuperfamily);

#    print Data::Dumper->Dump([$tcids ], [qw(*tcids )]);
#    exit;



    #==========================================================================
    #Setup the thresholds for parsing the PFAM output


    my $obj = new TCDB::Domain::PfamParser();
    $obj->pfamFile($pfamFile);
    $obj->analysisLevel($analysisLevel);
    $obj->domCovCutoff($domain_cov);
    $obj->tcCovCutoff($prot_cov);
    $obj->evalueCutoff($evalue);
    $obj->minProtsDom($prop_prots_w_domain);
    $obj->treatAsSuperfamily($treatAsSuperfamily);


    my %domFreq   = ();
    my %domCoords = ();
    $obj->getDomainStatsForUserFamilies(\@fams, $tcids, \%domFreq, \%domCoords);

#    print Data::Dumper->Dump([ \%domFreq, \%domCoords ], [qw( *domFreq *domCoords )]);
#    exit;





    #==========================================================================
    #Attempt to rescue the domains that were not recognized by PFAM in some
    #Family members


    my $rescueObj = new TCDB::Domain::Characterize();
    $rescueObj->rootDir($rootDir);
    $rescueObj->tcdbFaa($tcdbSeqsFile);
    $rescueObj->domCoords(\%domCoords);
    $rescueObj->domFreq(\%domFreq);
    $rescueObj->tcids($tcids);
    $rescueObj->searchWith($prog);
    $rescueObj->blastdb($blastdb);
    $rescueObj->evalue($evalue);
    $rescueObj->treatAsSuperfamily($treatAsSuperfamily);


    $rescueObj->rescueDomains(\@fams);

}
else {

  foreach my $fam (@fams) {

    my $tcids = getSystemAccessions($tcdbSeqsFile, 'both', $analysisLevel, [$fam], $treatAsSuperfamily);

#    print Data::Dumper->Dump([$tcids ], [qw( *tcids )]);
#    exit;



    #==========================================================================
    #Setup the thresholds for parsing the PFAM output


    my $obj = new TCDB::Domain::PfamParser();
    $obj->pfamFile($pfamFile);
    $obj->analysisLevel($analysisLevel);
    $obj->domCovCutoff($domain_cov);
    $obj->tcCovCutoff($prot_cov);
    $obj->evalueCutoff($evalue);
    $obj->minProtsDom($prop_prots_w_domain);


    my %domFreq   = ();
    my %domCoords = ();
    $obj->getDomainStatsForUserFamilies([], $tcids, \%domFreq, \%domCoords);

    #print Data::Dumper->Dump([ \%domFreq, \%domCoords ], [qw( *domFreq *domCoords )]);
    #exit;





    #==========================================================================
    #Attempt to rescue the domains that were not recognized by PFAM in some
    #Family members


    my $rescueObj = new TCDB::Domain::Characterize();
    $rescueObj->rootDir($rootDir);
    $rescueObj->tcdbFaa($tcdbSeqsFile);
    $rescueObj->domCoords(\%domCoords);
    $rescueObj->domFreq(\%domFreq);
    $rescueObj->tcids($tcids);
    $rescueObj->searchWith($prog);
    $rescueObj->blastdb($blastdb);
    $rescueObj->evalue($evalue);
    $rescueObj->domCovCutoff($domain_cov);
    $rescueObj->treatAsSuperfamily($treatAsSuperfamily);

    $rescueObj->rescueDomains();
  }
}




###########################################################################
##                         Functions                                     ##
###########################################################################



sub read_command_line_arguments {

  #if no arguments are given print the help
  if (! @ARGV) {
    print_help();
  }

  #----------------------------------------------------------------------
  #Parse command line arguments

  my $status = GetOptions(
      "f|family=s"           => \&read_fams,           #TCIDs of families to analyze (comma separated)

      #Options for TCDB::Domain::PfamParser
      "dc|domain-cov=f"      => \$domain_cov,
      "pc|protein-cov=f"     => \$prot_cov,
      "e|evalue=f"           => \$evalue,
      "m|prots-w-domain=f"   => \$prop_prots_w_domain,

      #Options for TCDB::Domain::Characterize
      "pt|proj-targets=s"    => \&read_proj_targets,    #Target Proteins, NOT in TCDB, to project domains onto
      "o|outdir=s"           => \&read_root_dir,        #Ouput root directory
      "s|tcdb-seqs=s"        => \&read_tcdb_seqs,       #File with all sequences in TCDB
      "sf|superfamily!"      => \$treatAsSuperfamily,	#File with the sequences of the reference family
      "pfam=s"               => \&read_pfam,            #hmmscan output file for whole TCDB
      "b|blastdb=s"          => \&read_blastdb,         #Full path of blastdb to extract sequences
      "p|rescue-prog=s"      => \&read_prog,            #Read the program that will be used to rescue domains (blastp|ssearch36)
      "h|help"               => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);

  #----------------------------------------------------------------------
  #Validate command line arguments


  die "Error: Options -f and -pt are incompatible" if (@fams && @candProjProts);
  die "Error: either -f or -pt must be given" unless (@fams || @candProjProts);

  if (@candProjProts) {
    prepare_seqs_for_projection();
  }

}


#==========================================================================
#Setup the environment for projection of domains onto sequences that are
#not present in TCDB.

sub prepare_seqs_for_projection {

  my $tcdbDir = "$rootDir/tcdb";
  system "mkdir -p $tcdbDir" unless (-d $tcdbDir);

  #to prevent modifying the original files, here I'll save the input
  #sequences with the artificial TCIDs
  my $origInfilesDir = "$rootDir/inputFiles";
  system "mkdir -p $origInfilesDir" unless (-d $origInfilesDir);


  #generate an empty "TCDB sequence file" that will contains proteins not in TCDB
  my $new_tcdbSeqsFile = "$tcdbDir/tcdb.faa";
  system "cat /dev/null > $new_tcdbSeqsFile";


  #----------------------------------------------------------------------
  #generate the new TCDB database relevant for the projection

  foreach my $pair (@candProjProts) {

    my $tcid = $pair->[0];
    my $tgtF = $pair->[1];

    my @comp = split(/\//, $tgtF);
    my $tgtFileName = $comp[-1];

    #Add family to the main array (as if provided by the -f commandline option)
    push (@fams, $tcid);


    #extracts the tcids of the systems under reference family
#    my $tcdbSeqs = $tcdbSeqsFile; #"$ENV{HOME}/db/blastdb/tcdb.faa";
#    die "TCDB sequences not found: $tcdbSeqs" unless (-f $tcdbSeqs);

    my $sysHash = TCDB::Assorted::getSystemAccessions($tcdbSeqsFile, 'both', 'system', [$tcid], 0);

#    print Data::Dumper->Dump([$sysHash, $tcdbSeqs], [qw(*sysHash *tcdbSeqs)]);
#    <STDIN>;

    #determine the TCID that will be used as reference for the target sequences
    my @systems = @{ $sysHash->{$tcid} };
    die "Could not find TCIDs for $tcid in $tcdbSeqsFile" unless (@systems);

    my $tgtTC = $systems[-1]->[0];
    $tgtTC =~ s/\.\d+$/\.10000/;


    #Replace the TCID in the file corresponding to the target proteins
    my $cmd1 = qq(perl -pe 's/\\>([a-zA-Z0-9_-]+).*/\\>${tgtTC}-\$1/;' $tgtF > $origInfilesDir/$tgtFileName);
    system $cmd1 unless (-f "$origInfilesDir/$tgtFileName");


    #Extract sequences for reference family
    my $outFile = "$tcdbDir/tcdb-${tcid}.faa";
    my $cmd2 = qq(extractTCDB.pl -i $tcid -o $tcdbDir -d $tcdbSeqsFile);
    system $cmd2 unless (-f $outFile);
    die "Could not generate sequence file: $outFile" unless (-f $outFile);


    #Add family and target sequences to the new TCDB family
    my $cmd3 = qq(cat $outFile $origInfilesDir/$tgtFileName >> $new_tcdbSeqsFile);
    system $cmd3;
  }

  #----------------------------------------------------------------------
  #Generate the PFam database

  my $pfamD   = "$rootDir/pfam";
  system "mkdir -p $pfamD" unless (-d $pfamD);

  my $pfamTMPfile = "$pfamD/tcdb_pfam.out";
  $pfamFile = "${pfamTMPfile}.bz2";

  #run Pfam
  my $cmd4 = qq (hmmscan --cpu 4 --noali --cut_ga -o /dev/null --domtblout $pfamTMPfile $ENV{RESEARCH_DATA}/pfam/pfamdb/Pfam-A.hmm $new_tcdbSeqsFile);
  system $cmd4 unless (-f $pfamTMPfile || -f $pfamFile);

  #compress pfam file
  my $cmd5 = qq(bzip2 $pfamTMPfile);
  system $cmd5 unless (-f $pfamFile);


  #----------------------------------------------------------------------
  #now generate the blast DB

  my $blastD = "$rootDir/blastdb";
  system "mkdir -p $blastD" unless (-d $blastD);

  $blastdb = "$blastD/tcdb";

  my $cmd6 = qq(extractTCDB.pl -i tcdb -o $blastD -f blast -d $new_tcdbSeqsFile);
  system $cmd6 unless (-f "${blastdb}.pin");


  #For all purposes update the TCDB sequence file to point to the new "customized" file.
  $tcdbSeqsFile = $new_tcdbSeqsFile;
}





#==========================================================================
#Read the -pt option. It is expected that the user provides the family to which
#the target proteins are expected to belong. Example format should is:
# -pt {tcid_1},{file with target sequences 1}:{tcid_2},{file with target sequences 2}.
#
#NOTE: This option is incompatible with -f

sub read_proj_targets {

  my ($opt, $value) = @_;

  my @pairs = split (/:/, $value);
  die "No significant argument passed to option -pt" unless (@pairs);

  foreach my $pair (@pairs) {
    my ($tc, $file) = split (/,/, $pair);
    die "Error: not a valid {tcid},{file} pair: $pair" unless ($tc && $file);

    TCDB::Assorted::validate_tcdb_id([$tc]);

    unless (-f $file && !(-z $file)) {
      die "File with projection targets for $tc was not found or empty: $file";
    }

    push (@candProjProts, [$tc, $file]);

  }
}



#==========================================================================
#Read the -f option

sub read_fams {

  my ($opt, $value) = @_;

  @fams = split (/,/, $value);

  TCDB::Assorted::validate_tcdb_id(\@fams);
}


#==========================================================================
#Read the -d option

sub read_root_dir {

  my ($opt, $value) = @_;

  system "mkdir -p $value" unless (-d $value);

  $rootDir = $value;
}


#==========================================================================
#Read the -s option

sub read_tcdb_seqs {

  my ($opt, $value) = @_;

  die "File with TCDB sequences must exist and not be empty: $value" unless (-f $value && !(-z $value));

  $tcdbSeqsFile = $value;
}



#==========================================================================
#Read the -pfam option

sub read_pfam {

  my ($opt, $value) = @_;

  die "File with TCDB sequences must exist and not be empty: $value" unless (-f $value && !(-z $value));

  $pfamFile = $value;
}


#==========================================================================
#Read the -b option

sub read_blastdb {

  my ($opt, $value) = @_;

  my $tmpFile = "${value}.phd";
  die "Input is not a path to blast database: $value" unless (-f $tmpFile && !(-z $tmpFile));

  $blastdb = $value;
}



#==========================================================================
#Read the - option

sub read_prog {

  my ($opt, $value) = @_;

  die "Unrecognized program: $value" unless ($value =~ /(ssearch36|blastp)/);

  $prog = $value;
}




#==========================================================================
#Print help

sub print_help {

  my $help = <<'HELP';

Identify the main protein domains in a family.

 -f, --family {string} 
  TCID of the family for which Pfam domain analysis will be carried out.
  If multiple TCIDs are given, they whould be comma-separated. The analysis
  will be performed individually for each family, unless the flag -sf is
  given, in which case all TCIDs will be treated as a superfamily.
  This option is incompatible with option -pt. But either -f or -pt
  must be given.

-pt, --proj-targets {string}
  Project the characteristic domains of a reference family onto
  candidate protein(s) not in TCDB that might belong to the family.
  The format indicates pairs tcid,targets separated by ':'. That is:

  {tcid_1},{seq_file_1}:{tcid_2},{seq_file_2}:... {tcid_n},{seq_file_n}

  where tcid_n is the reference family that will project its domains onto
  the sequences in seq_file_n
  This option is incompatible with option -f. But either -f or -pt
  must be given.

 -dc, --domain-cov {float} (Optional; Default: 0.7)
  Minimum coverage of the Pfam domain to consider it a match. If coverage
  is less than the specified threshold, the coverage must apply to the
  the query protein to consider the domain hit significant.

 -pc, --protein-cov {float} (Optional; Default: 0.1)
  Minimum coverage of the query protein sequence per domain.

 -e, --evalue {float} (Optional; Default: 1e-5)
  Maximum evalue threshold to consider a Pfam domain hit.

 -m, --prots-w-domain {float) (Optional; Default: 0.05)
  Minimum proportion of the proteins in the input family that should
  contain a domain, in order to consider the domain as part of the
  family for the purpose of this analyis.

 -o, --outdir {path} (Default: .)
  Directory where results and intermediary files will be saved.

 -s, --tcdb-seqs {file} (Optional; Default: $RESEARCH_DATA/pfam/download/tcdb.faa)
  FASTA file with all sequences in TCDB (as generated by program
  extractTCDB.pl). This file will be used to extract sequences
  and TCIDs for all members of the input family. This allows to
  freeze TCDB contents at a specific date, for the purpose defining
  a project.

 -pfam {file} (Optional; Default: $RESEARCH_DATA/pfam/tcdb.pfam-a.hmmscan.bz2)
  The output of running hmmscan against all TCDB, or at least the input
  family(-ies).

 -b, --blastdb {path} (Default: $HOME/db/blastdb/tcdb)
  Full path to the blast database containing all the sequences in
  TCDB, this is to easily extract segments of proteins that match
  Pfam domains.

 -sf, --superfamily {flag} (default: negated as --no-superfamily)
  This indicates that all families passed to option -f or -pt will be
  treated as a superfamily (i.e. as a single family) for the purpose of
  the Pfam domain  analysis. This is useful when a superfamily is
  composed of two or more different family TCIDs.

 -p, --rescue-prog (Default: ssearch36)
  Program that will be used to project or rescue domains in
  proteins that did have direct hits with Pfam domains. So far,
  only ssearch36 is supported.

 -h, --help
  Display this help. If present, this option takes precedence over any
  other option.

HELP


  print $help;
  exit;
}
