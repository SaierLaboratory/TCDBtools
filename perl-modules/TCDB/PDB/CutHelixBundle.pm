package TCDB:PDB::CutHelixBundle;


use warnings;
no warnings;  #To remove harmeless warnings from overriding accessors

use strict;
use Data::Dumper;


use R2::PDB::Parser;

use Class::Struct;


#==========================================================================
# The purpose of this class is cut PDB files in subsets of structures 
# containing a predefined number of contiguous alpha-helices. For example,
# if there are 5 helices in the structre and we want sets of 3-helix bundles
# per pdb files, the script will generate 3 files with the helices: 
# (1,2,3), (2,3,4) and (3,4,5).
#
# NOTE:
#  Mapping original ATOM positions to a position relative to SEQRES
#  (to easily map atoms to TMS positions) will not always work
#  because in some structures the ATOM sequence is not identical
#  to the SEQRES sequence (e.g. 4PGW). For this reason I decided
#  to use the ATOM sequence to run hmmtop.
#
# -------
#
#
# GENERAL GUIDE FOR USAGE:
#
# Load class:
# use R2::PDB::CutHelixBundle;
#
#
# Constructor:
#
# my $myObj = R2::PDB::CutHelixBundle->new(
#                         pdb_to_cut => $myPDBfile,
#                         resolution_cutoff => $myThreshold,
#                         helix_bundle_size => $myHelixBundle,
#                         pdb => {infile => $myPDBfile, min_helix_size => $myHelixSize}
#                       );
#
#
#
# Get help about this class:
#
# $myObj->help;
#
#
#
# Initiallizing object through constructor:
#
# my $myObj = new R2::PDBparser(infile => '/path/to/pdb/file');
#
#
# Initializing object through accessors:
#
# my $myObj = new R2::PDBparser();
# $myObj->infile('/path/to/pdb/file');
#
#
#
# Once object is initialized parse the PDB file as:
#
# $myObj->parse_pdb();
#
#
#
# -------
# Date created:  11/13/2015
# Programmer:    Arturo Medrano
# Institution:   UCSD, WLU
#==========================================================================



struct ("R2::PDB::CutHelixBundle" =>
	{
	 'pdb_to_cut'          => '$',  #Input PDB file to cut
	 'outdir'              => '$',  #Directory where cut structures will be placed.
	 'helix_bundle_size'   => '$',  #Cut the helices in budles of this number
	 'helix_bundle_tails'  => '$',  #How many residues to leave before an after the helx bundle.
	 'debug'               => '$',  #Variable that facilites printing debug messages
	 'pdb'                 => 'R2::PDB::Parser'   #The PDB parser object
	}
       );



###########################################################################
#####        Override accessors that require default values           #####
#####  (This generates harmless warings, make sure to turn the off)   #####
###########################################################################


#==========================================================================
#Assign default output dir if not provided by the user.

sub outdir {

  my ($self, $dir) = @_;

  my $default_dir = "./helix_bundles";

  if ( $dir ) {
    system "mkdir -p $dir" unless (-d $dir);
    $self->{'R2::PDB::CutHelixBundle::outdir'} = $dir;
  }

  unless ($self->{'R2::PDB::CutHelixBundle::outdir'}) {
    system "mkdir $default_dir" unless (-d $default_dir);
    $self->{'R2::PDB::CutHelixBundle::outdir'} = $default_dir;
  }

  return $self->{'R2::PDB::CutHelixBundle::outdir'};
}



#==========================================================================
#Assign default helix bundle size to cut if not provided by the user.

sub helix_bundle_size {

  my ($self, $size) = @_;

  my $default_size = 3;

  if ($size) {

    #validate
    unless  ( $size > 1 ) {
      die "Error: helix_bundle_size only takes integer values greater than 1\n";
    }

    $self->{'R2::PDB::CutHelixBundle::helix_bundle_size'} = $size;
  }

  unless ($self->{'R2::PDB::CutHelixBundle::helix_bundle_size'}) {
    $self->{'R2::PDB::CutHelixBundle::helix_bundle_size'} = $default_size;
  }

  return $self->{'R2::PDB::CutHelixBundle::helix_bundle_size'};
}



#==========================================================================
#Assign default tail size for cutting helix bundles if not privided by
#he user. 

sub helix_bundle_tails {

  my ($self, $tail_size) = @_;

  my $default_tail_size = 1;

  if ($tail_size) {

    #Validate
    unless ($tail_size =~ /^\d+$/ && $tail_size > 0) {
      die "Error: helix_bundle_tails only takes integer values greater than zero.\n";
    }

    $self->{'R2::PDB::CutHelixBundle::helix_bundle_tails'} = $tail_size;
  }

  unless ($self->{'R2::PDB::CutHelixBundle::helix_bundle_tails'}){
    $self->{'R2::PDB::CutHelixBundle::helix_bundle_tails'} = $default_tail_size;
  }

  return $self->{'R2::PDB::CutHelixBundle::helix_bundle_tails'};
}






###########################################################################
#####                  The methods in this class                      #####
###########################################################################



#==========================================================================
#
#



sub cut_helix_bundles {
  my ($self, $chain) = @_;

  my $pdb_id = $self->pdb->id . "_$chain";


  #Determine whether the structure has helices annotated
  die "Error: No helices in structure $pdb_id\n" unless (%{ $self->pdb->helix });


  #Structures with no resolutions data have a 100 value assigned in the parser
#  return 3 if ($self->pdb->resolution > $self->resolution_cutoff);



  #There must be at least the number of helices that we want to cut
  unless (scalar keys %{ $self->pdb->stride->{$chain} } >= $self->helix_bundle_size) {
    print "Need at least ", $self->helix_bundle_size, " helices in strucure $pdb_id to make the cut\n";
    return;
  }




  #Identify the index of the first helix in the last bundle
  my @helices = sort {$a <=> $b} keys  %{ $self->pdb->stride->{$chain} };
  my $first_helix_idx = $helices[0];
  my $last_helix_idx  = $helices[$#helices]  - $self->helix_bundle_size + 1;



  #Cut helices according to the helix-bundle size provided by the user
  for (my $i = $first_helix_idx; $i <= $last_helix_idx ; $i++) {


    #Get the first and last helices in the bundle
    my $first_helix = $self->pdb->stride->{$chain}->{$i};
    my $last_helix  = $self->pdb->stride->{$chain}->{$i + $self->helix_bundle_size - 1};


    #Now get the first and last residues in the helix bundle (left and right tails)
    my $first_res = $first_helix->[1];
    my $last_res  = $last_helix->[3];


    #the left tail
    my ($lTail, $rTail) = undef;
    if ($first_res - $self->helix_bundle_tails < 1) {
      $lTail = 1;
    }
    else {
      $lTail = $first_res - $self->helix_bundle_tails;
    }

    #The right tail
    if ($last_res + $self->helix_bundle_tails > scalar keys %{ $self->pdb->atom_seq1->{$chain} }) {
      $rTail = scalar keys %{ $self->pdb->atom_seq1->{$chain} };
    }
    else {
      $rTail = $last_res + $self->helix_bundle_tails;
    }



#    print Data::Dumper->Dump([$first_helix, $last_helix, $first_res, $last_res],[qw(*first *last *first_res *last_res)]);
#    <STDIN>;


    #Create the PDB file with the cut helices
    my @helices_in_this_bundle = ($i .. $i + $self->helix_bundle_size - 1);
    my $helices_str = join("_", @helices_in_this_bundle);
    my $hbundle_name = "${pdb_id}_h$helices_str";

    $self->generate_helix_bundle_pdb_file($chain, $hbundle_name, \@helices_in_this_bundle, $lTail, $rTail);
    #print "Helix bundle PDB file generated: $hbundle_name\n";
  }

  return 1;

}


#==========================================================================
#Given the cordinates of the first and last amino acid in a helix bundle,
#generate a PDB file with the corresponding atom coordinates.

sub generate_helix_bundle_pdb_file {

  my ($self, $chain, $bundle_name, $bundle_helices, $firstRes, $lastRes) = @_;


  #The output directory
  my $outdir = $self->outdir;


  #The output file name
  my $pdb_file = "$outdir/${bundle_name}.pdb";


  #Do generate bundle file if it already exists
  return if (-f $pdb_file);


  open (my $fhandle, ">", $pdb_file) || die $!;


  #Print the header
  print $fhandle $self->pdb->header, "\n";


  #Resolution line
  if ($self->pdb->resolution_line) {
    print $fhandle $self->pdb->resolution_line, "\n";
  }


  #The DBREF section to extract SwissProt IDs
  print $fhandle join("\n", @{ $self->pdb->dbref_lines->{$chain} }), "\n";



  #The MODRES section
  if (exists $self->pdb->modres->{$chain} && @{ $self->pdb->modres->{$chain} }) {
    print $fhandle join ("\n", @{ $self->pdb->modres->{$chain} }), "\n";
  }


  #Print the helix section
  foreach my $helix (@{ $self->pdb->helix->{$chain} }) {
    print $fhandle "$helix\n";
  }


  #Print the atoms of the bundle
  foreach my $res_num ($firstRes .. $lastRes) {
    foreach my $atom (@{ $self->pdb->atom->{$chain}->{$res_num} }) {
      print $fhandle "$atom\n";
    }
  }


  #The CONNECT section if it exists
  if (@{ $self->pdb->connect }) {
    print $fhandle join("\n", @{ $self->pdb->connect }), "\n";
  }


  #The MASTER line if it exists
  print $fhandle $self->pdb->master, "\n" if ($self->pdb->master);



  print $fhandle "END\n";

  close $fhandle;


}


#==========================================================================
#Get the index of the first helix that overlaps the first TMS

sub get_idx_of_first_helix {

  my ($self, $chain) = @_;

  my $idx_last_helix = $#{ $self->helix_corrected->{$chain} };

 HELIX:foreach my $idx ( 0 .. $idx_last_helix ) {

    my $helix =  $self->helix_corrected->{$chain}->[$idx];

    my $lpos = $self->trim(substr($helix, 21, 4)); #chars 22-26
    my $rpos = $self->trim(substr($helix, 33, 4)); #chars 34-37


    #return index of first helix with overlap
    if ( $self->is_there_overlap([$lpos, $rpos], $self->hmmtop->{$chain}->[0]) ) {
      return $idx;
    }
  }

  #no overlap found
  return undef;
}



#==========================================================================
#Determine whether there 2 sets of amino acid positions overlap

sub is_there_overlap {

  my ($self, $coords1, $coords2) = @_;

  my ($l1, $r1) = undef;
  my ($l2, $r2) = undef;

  #assign l1 and r1 to the left most set of coordinates
  if ($coords1->[0] <= $coords2->[0]) {
    ($l1, $r1) = @$coords1;
    ($l2, $r2) = @$coords2;
  }
  else {
    ($l2, $r2) = @$coords1;
    ($l1, $r1) = @$coords2;
  }

  #coordinates must be grater than zero
  die "Error: coordinates must exist to estimate overlap" unless ($l1 && $r1 && $l2 && $r2);

  if ($r1 - $l2 >= 0) {
    return 1;
  }
  else {
    return 0;
  }
}






#==========================================================================
#Extract sequence in 1-letter code from the ATOM section

sub extract_atom_seq_1 {

  my ($self, $chain) = @_;

  my $seq = "";
  foreach my $idx (sort {$a <=> $b} keys %{ $self->atom_seq1->{$chain} }) {
    $seq .= $self->atom_seq1->{$chain}->{$idx};
  }

  die "Error: could not extract seq from ATOM section" unless ($seq);

  return $seq;

}



#==========================================================================
#Trim blanks from the begnning and end of a string

sub trim {

  my ($self, $string) = @_;

  die "Error: string passed to trim can't be empty\n" unless ($string);

  #triming
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;

  return $string;
}



#==========================================================================
#The help for this class

sub help {

  my $self = shift;

  my $help = << 'HELP';


HELP

  print $help;
  exit;
}




1;

