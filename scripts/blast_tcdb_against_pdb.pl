#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

use Bio::SeqIO;
use Getopt::Long;
use TCDB::Assorted;


#==========================================================================
#
# This program compares all proteins in a TCDB family against all protein
# chains in PDB. It will extract the PDB IDs of all significant hits and
# keep those that contain at least one TMS.
#
# The structures reported can be part of TCDB or focus only on those
# structures not yet reported in TCDB.
#
#==========================================================================


#==========================================================================
# Global varaibles

my $inFamily    = "";    #family that will be analyzed
my $overwrite   = 0;     #overwrite downloded files if present
my $evalue      = 1e-5;  #Default E-value for blast searches
my $id_cutoff   = 30;    #minimum ID % for blast alignments
my $cov_cutoff  = 50;    #minimum coverage % of query or subject seqs
my $workdir     = "";    #Ouput dir where relevant files are saved
my $minTMS      = 1;     #Minimum number of TMS in the PDB structures
my $i_pdb_TCDB  = 0;     #Ignore PDBs already available in TCDB


#Read command line arguments
read_command_line_arguments();


#print Data::Dumper->Dump([$inFamily, $overwrite, $evalue, $id_cutoff, $cov_cutoff, $workdir, $minTMS, $ignoreTCDB ],
#			 [qw(*inFamily *overwrite *evalue *identity *coverage *workdir *minTMS *ignoreTCDB)]);
#exit;




#==========================================================================
#Download all the PDB IDs currently available in TCDB. Use this ID
#to discard blast matches in PDB with structures already in TCDB.

my $tcdb_pdb = "$workdir/pdb.tsv";
system "wget -O $tcdb_pdb  http://www.tcdb.org/public/pdb.tsv" if ($overwrite || !(-f $tcdb_pdb));

my %tcdb_structures = ();
open (my $tcdbh, "<", $tcdb_pdb) || die $!;

while (<$tcdbh>) {

  chomp;
  my ($pdb, $fam, $description) = split (/\t/, $_);

  next unless ($pdb && $fam);

  $tcdb_structures{$pdb} = $fam;

}
close $tcdbh;


print "... Downloaded current structures in TCDB (", scalar keys %tcdb_structures, ")\n";
#exit;





#==========================================================================
#Download family sequences from TCDB

my $tcdb_seqs_file = "";

if (lc $inFamily eq "tcdb" || lc $inFamily eq "all" || lc $inFamily eq "full") {
  $tcdb_seqs_file = "$workdir/tcdb.faa";
}
else {
  $tcdb_seqs_file = "$workdir/family-${inFamily}.faa";
}

system "extractFamily.pl -i $inFamily -o  $workdir -f fasta" if ($overwrite || !(-f $tcdb_seqs_file));

print "... TCDB sequences were extracted.\n";
#exit;





#==========================================================================
#Blast all TCDB sequences against PDB

my $blast_hits_clm_file = "$workdir/blast_tcdb_vs_pdb.clm";

my $fmtBlastOut = '6 qseqid sseqid qlen slen evalue pident nident length qcovs qstart qend sstart send sallseqid';
if ($overwrite || !(-f  $blast_hits_clm_file)) {
  system "blastp -db pdbaa -evalue $evalue -outfmt '$fmtBlastOut' -query $tcdb_seqs_file -out $blast_hits_clm_file";
}

die "No blast output file found: $blast_hits_clm_file" unless (-f $blast_hits_clm_file);

print "... Blasted TCDB sequences against PDB.\n";
#exit;




 
#==========================================================================
#Extract the GIs for all hits and make sure they have at least one TMS

my $blast_hits_gi_file = "$workdir/top_hits_gis.txt";

my %redPDB = ();
my %redGI  = ();


open(my $gish, ">", $blast_hits_gi_file) || die $!;

open (my $blasth, "<", $blast_hits_clm_file) || die;
while (<$blasth>) {

  chomp;
  my ($qid, $sid, $qlen, $slen, $ev, $pident, $nident, $alen, $qcov, $qstart, $qend, $sstart, $send, $sredids) = split(/\t/, $_);


  #Alignment coverage for the subject
  my $scov = ($send - $sstart)/$slen * 100;


  #alignment must satisfy the requirements for coverage and identity
  next unless ($pident >= $id_cutoff && ($qcov >= $cov_cutoff || $scov >= $cov_cutoff));


  #Extract TCDB ID of query and the gi and PDB ID of the subject
  my $tcdb_id = "";
  if (lc $inFamily eq "tcdb" || lc $inFamily eq "all" || lc $inFamily eq "full") {
    $tcdb_id = (split (/\|/, $qid))[3];
  }
  else {
    $tcdb_id = (split (/\-/, $qid))[0];
  }
  my ($sgi, $spdb)  = ($sid =~ /gi\|(\d+)\|pdb\|(\w+)/)? ($1, $2) : ();


  #Make sure this line was parsed correctly
  die "could not parse query id: $qid" unless ($tcdb_id && $sgi && $spdb);



#  print "$_\n";
#  print Data::Dumper->Dump([$qid, $sid, $qlen, $slen, $ev, $pident, $nident, $alen, $qcov, $scov, $qstart, $qend, $sstart, $send, $sredids, $tcdb_id, $spdb ],
#			   [qw(*qid, *sid, *qlen, *slen, *ev, *pident, *nident, *alen, *qcov, $scov, *qstart, *qend, *sstart, *send, *sredids *tcdb_id *spdb)]);
#  <STDIN>;



  #Now extract the PDB ids and GIs of the redundant structures
 RED:while($sredids =~ /gi\|(\d+)\|pdb\|(\w+)/g) {

    my ($gi, $pdb) = ($1, $2);

    die "Could not parse blast line: $sredids" unless ($gi && $pdb);

    #If specified, ignore PDB entry if the structure is already in TCDB
    if ($i_pdb_TCDB) {
      next RED if (exists $tcdb_structures{$pdb});
    }

    print $gish "$gi\n" unless (exists $redPDB{$pdb});

    $redGI{$gi}  = [$tcdb_id, $pdb, $ev, $pident, $qcov, $scov, $qlen, $slen];
    $redPDB{$pdb}{$gi} = [$tcdb_id, $ev, $pident, $qcov, $scov, $qlen, $slen];
  }

#  print Data::Dumper->Dump([$spdb, \%redPDB, \%redGI, $sredids ], [qw(*spdb *redPDBs *redGI *sredids)]);
#  <STDIN>;

}

close $blasth;
close $gish;


die "File with GIs of top blast hits does not exist: $blast_hits_gi_file" unless (-f $blast_hits_gi_file);

print "... Extracted the GIs of the top structures returned by blast.\n";
#exit;






#==========================================================================
#Not extract the sequences for top blast matches so HMMTOP can be run
#in the next step

my $blast_hits_seqs_file = "$workdir/top_hits_seqs.faa";
my $get_seq  = qq(blastdbcmd -db pdbaa -entry_batch $blast_hits_gi_file -target_only -out $blast_hits_seqs_file);
system $get_seq if ($overwrite || !( -f  $blast_hits_seqs_file ));

die "No sequences were retrieved for blast hits: $blast_hits_seqs_file" unless ( -f  $blast_hits_seqs_file  );


print "... Extracted the protein sequences of top blast hits against PDB\n";
#exit;





#==========================================================================
#Run hmmtop in the sequences of matching chains, the purpose is
#to keep only structures with at least 1 TMS.

my $hmmtop_top_hits = "$workdir/top_hits_tms.txt";

system "hmmtop -if=$blast_hits_seqs_file -of=$hmmtop_top_hits -is=pseudo -pi=spred -sf=FAS" if ( $overwrite || !( -f $hmmtop_top_hits ));

die "No hmmtop output file found:  $hmmtop_top_hits" unless (-f $hmmtop_top_hits);

print "... Ran HMMTOP on the top blast hits against PDB\n";
#exit;





#==========================================================================
#Parse hmmtop to get the PDBs of structures that have TMSs

my @results = ();
open (my $tmsh, "<", $hmmtop_top_hits) || die $1;
while(<$tmsh>) {
  chomp;

  my ($len, $gi, $pdb, $ntms) = (/\s+(\d+)\s+gi\|(\d+)\|pdb\|(\w+).+\s(IN|OUT)\s+(\d+)/)? ($1, $2, $3, $5) : ();
  next unless ($ntms >= $minTMS);


  #Get the TCDB ID of the transporter matching this protein.
  my ($tcdb_id, $nPDB, $ev, $pident, $qcov, $scov, $qlen, $slen) = @{ $redGI{$gi} };
  die "No TCDB info for: $gi, $pdb" unless ($nPDB && $tcdb_id && $pdb eq $nPDB);

  #print STDERR "Length of hmmtop is different from blast [$tcdb_id, $pdb, $gi]: hmmtop($len), Blast($slen)\n" unless ($len == $slen);


  push (@results, [$pdb, $tcdb_id, $gi, $ntms, $ev, $pident, $qcov, $scov, $qlen, $slen]);
}
close $tmsh;


print "... Filtered out structures that had no TMS\n";
#exit;







#==========================================================================
#Print structures that should be added to TCDB

my $outfile = "";

if (lc $inFamily eq "tcdb" || lc $inFamily eq "all" || lc $inFamily eq "full") {
  $outfile = "$workdir/new_pdb_for_tcdb.txt";
}
else {
 $outfile = "$workdir/new_pdb_for_${inFamily}.txt";
}


open (my $outh, ">", $outfile) || die $!;

map { print $outh join ("\t", @{ $_ }), "\n"; } sort by_tms_len @results;

close $outh;


print "DONE!! Generated the list of structures to add to TCDB for: $inFamily  (Found ", scalar @results, " structures)\n";





#==========================================================================
#Sort results by number of TMS and then by length of subject


sub by_tms_len {

  if ($a->[3] == $b->[3]) {
    $b->[9] <=> $a->[9];
  }
  else {
    $b->[3] <=> $a->[3];
  }
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
      "f|family=s"       => \&read_fam,
      "o|overwrite!"     => \$overwrite,
      "t|ignore-tcdb!"	 => \$i_pdb_TCDB,
      "d|workdir=s"      => \$workdir,
      "e|evalue=f"       => \$evalue,
      "c|coverage=f"     => \&read_coverage,
      "i|identity=f"     => \&read_identity,
      "m|min-tms=i"      => \&read_min_tms,
      "h|help"           => sub { print_help(); },

      #For arguments that do not look like valid options
      "<>"              => sub { die "Error: Unknown argument: $_[0]\n"; }
  );
  die "\n" unless ($status);



  #----------------------------------------------------------------------
  #Validate command line arguments


  #Input TCDB family
  die "Error: argument -f is mandatory.\n" unless ($inFamily);


  #Set the working directory
  unless ($workdir) {
    $workdir = $inFamily;
  }
  system "mkdir -p $workdir";
}



#==========================================================================
#read the minimum number of TMS in target structures

sub read_min_tms {

  my ($opt, $value) = @_;

  die "Artument passed to option -m must be an integer >= 0." if ($value < 0);

  $minTMS = $value;

}



#==========================================================================
# Read the IDentity threshold

sub read_identity {

  my ($opt, $value) = @_;

  if ($value < 25) {
    die "Identity must be at list 25%";
  }

  $id_cutoff = $value;
}



#==========================================================================
# Read the coverage threshold

sub read_coverage {

  my ($opt, $value) = @_;

  if ($value <= 25.0) {
    die "Coverage must be at list 30%";
  }

  $cov_cutoff = $value;
}



#==========================================================================
#Read the -f option

sub read_fam {

  my ($opt, $value) = @_;

  unless (lc $value eq "tcdb" || lc $value eq "all" || lc $value eq "full") {
    TCDB::Assorted::validate_tcdb_id([$value]);
  }

  $inFamily = $value;
}



#==========================================================================
#Print help

sub print_help {

  my $help = <<'HELP';

 -f, --family {string}
   TCDB ID of the family for which the PDB IDs will be extracted.
   If family is: all, full or tcdb, the whole TCDB database
   will be processed. Argument is mandatory.

 -d, -workdir {string}
  Directory where results will be saved. If not given, by default
  the family ID passed with the -f option will be used instead.

 -o, --overwrite (negate with --no-overwrite) 
  Flag that inticates whether or not to overwrite files that were 
  downloaded, including blast comparisons and TMS predictions
  (Default --no-overwrite).

 -t, --ignore-tcdb  (Negate with --no-ignore-tcdb)
  Flag that indicates whether pdb structures already in TCDB should
  be excluded from the results. By defaults all structures are
  reported: --no-ignore-tcdb.

 -e, --evalue {float}
  E-value threshold that will be used in blast against PDB.

 -c, --coverage {float}
  Minimum alignment coverage of either the query and the subject
  proteins in blasts against TCDB. Should be a number between
  30 and 100 (Default 50.0).

 -i, --identity {float}
  Minimum identity requried in the aligned region.
  Should be a number between 25 and 100 (Default 30.0)

 -m, --min-tms {int}
  Minimum number of TMS in structures that will be considered.
  (Default 1).

 -h, --help
   Display the help of this program. Help will also be
   printed if no arguments are passed to the program.

HELP


  print "$help\n";
  exit;
}



