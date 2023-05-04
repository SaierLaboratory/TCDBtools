package TCDB::PDB::Parser;


use warnings;
no warnings;

use strict;
use Data::Dumper;

use Bio::DB::SwissProt;
use File::Type;

#use Bio::Structure::IO::pdb;
use Class::Struct;


#==========================================================================
# The purpose of this class is parse PDB files and retrieve individual
# Sections of files.
#
# The files and the different sections extracted as in the document, or
# the residue numbers in ATOM section can be reset to start at one. 
#
# -------
#
#
# GENERAL GUIDE FOR USAGE:
#
# Load class:
#
# use R2::PDB::parser;
#
#
# Constructor:
#
# my $myObj = new R2::PDB::parser();
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
# my $myObj = new R2::PDB::Parser (infile => '/path/to/pdb/file');
# my $myObj = R2::PDB::Parser->new(infile => '/path/to/pdb/file');
#
#
# Initializing object through accessors:
#
# my $myObj = new R2::PDB::parser();
# $myObj->infile('/path/to/pdb/file');
#
#
# Once the  object is initialized parse the PDB file as:
#
# $myObj->parse();
#
#
# List of available methods:
# parse()
#    Parse the PDB file
#
# extract_atom_seq_1($chain)
#    Extract the amino acid sequence of a chain in 1-letter code
#
# save_atom_sequence_in_fasta_format($chain, $outdir)
#    Save to a file the amino acid sequence of a chain in fasta format
#
# create_pdb_for_chain ($chain, $outdir)
#    Create a PDB files with the coordinates of only one chain
#
# -------
# Date created:  11/13/2015
# Programmer:    Arturo Medrano
# Institution:   UCSD, WLU
#==========================================================================



struct ("R2::PDB::Parser" =>
	{
	 'infile'              => '$',  #Input PDB file to cut
	 'tmpdir'              => '$',  #Directory to put temporal files (e.g. fasta files)
	 'renum_res_in_atom'   => '$',  #Indicate whether the residue number in ATOM section should 
	                                #be renumbered to start at 1.
	 'replace_helix'       => '$',  #Flag indicating if the HELIX coordinates will be corrected
	                                #and replaced by combining stride and hmmtop
	 'min_helix_size'      => '$',  #Define the minimum acceptable size for  helix
	 'max_helix_size'      => '$',  #Define the maximum acceptable size for  helix
	 'helix_separator'     => '$',  #The minimum number of residues between two helices

	 'header'              => '$',  #The header line of the file
	 'id'                  => '$',  #The PDB ID in file
	 'resolution'          => '$',  #The resolution of the structure in Angstroms
	 'resolution_line'     => '$',  #The resolution line that can be added to generated PDB files

         'dbref'               => '%',  #Hash containing the coordinates for each chain
	 'dbref_lines'         => '%',  #To include DBREF data when generating pdb files


	 'chain_ids'           => '@',  #List with the chains in the PDB file
	 'heterocomplex'       => '$',  #Flat indicating if structure is a hetercomplex
	 'chimeric'            => '$',  #Flag indicating if there are chimeric chains

	 'models'              => '%',  #Hash containing the different models IDs in the PDB file
	 'has_models'          => '$',  #Flag indicated whether the structure has more than 1 model

	 'seqres'              => '%',  #The seqres section
	 'sres_seq1'           => '%',  #Sequence per chain in 1 letter code
	 'sres_seq3'           => '%',  #Sequence per chain in 3 letter code


	 'modres'              => '%',  #The MODRES lines for each chain
	 'modres_corrected'    => '%',  #MODRES lines with renumbered residues
	 'modres_dic'          => '%',  #Mapings from modified residues to regular residues.
	 'has_red_modres'      => '$',  #Signal presence of redundant modified residues.

	 'helix'               => '%',  #Helices that will be used to cut the pdb file
	 'helix_corrected'     => '%',  #Helices with corrected positions according to atom_corrected


	 'atom'                => '%',  #Store here the information on atoms for each chain
	 'map_res_position'    => '%',  #hash containing the positions of the atoms starting from 1
	 'atom_corrected'      => '%',  #Store atoms with position of residues corrected for hmmtop comparisons
	 'hetatm'              => '%',  #Store here hereogenous atoms that are not modified residues.
	 'ter'                 => '%',  #The TER lines per chain

	 'atom_seq1'           => '%',  #Sequence in 1-letter code extracted from the atoms
	 'atom_seq3'           => '%',  #Sequence in 3-letter code extracted from the atoms
	 'atom_seq3_corrected' => '%',  #Corrected AA sequence in 3-letter code so it starts at amino acid 1

	 'connect'             => '@',  #If there is a connect section, save it here
	 'master'              => '$',  #save the master line here

	 'hmmtop'              => '%',  #In case it's necessary to run hmmtop for a specific chain.
	 'stride'              => '%',  #This will be used to replace the helix section of the structure
	 'debug'               => '$'   #To control the printing for debugging purposes
	}
       );





#==========================================================================
#Assign default temporal dir if not prived by the user.

sub tmpdir {
  my ($self, $dir) = @_;

  my $default_dir = "./tmp";

  if ( $dir ) {
    system "mkdir -p $dir" unless (-d $dir);
    $self->{'R2::PDB::Parser::tmpdir'} = $dir;
  }

  unless ($self->{'R2::PDB::Parser::tmpdir'}) {
    system "mkdir $default_dir" unless (-d $default_dir);
    $self->{'R2::PDB::Parser::tmpdir'} = $default_dir;
  }

  return $self->{'R2::PDB::Parser::tmpdir'};
}



#==========================================================================
#Override accessor for min_helix_size and set defaul value

sub min_helix_size {

  my ($self, $cutoff) = @_;

  my $default_cutoff = 15;

  if ($cutoff) {

    #validate
    unless ( $cutoff =~ /^\d+$/ && $cutoff >= $default_cutoff) {
      die "Error: min_helix_size should be at least $default_cutoff";
    }

    $self->{'R2::PDB::Parser::min_helix_size'} = $cutoff;
  }

  unless ($self->{'R2::PDB::Parser::min_helix_size'}) {
    $self->{'R2::PDB::Parser::min_helix_size'} = $default_cutoff;
  }

  return $self->{'R2::PDB::Parser::min_helix_size'};
}



#==========================================================================
#Override accessor for max_helix_sixe and set defaul value

sub max_helix_size {

  my ($self, $cutoff) = @_;

  my $default_cutoff = 40;

  if ($cutoff) {

    #validate
    unless ( $cutoff =~ /^\d+$/ && $cutoff >= $default_cutoff) {
      die "Error: max_helix_size should be at least $default_cutoff";
    }

    $self->{'R2::PDB::Parser::max_helix_size'} = $cutoff;
  }

  unless ($self->{'R2::PDB::Parser::max_helix_size'}) {
    $self->{'R2::PDB::Parser::max_helix_size'} = $default_cutoff;
  }

  return $self->{'R2::PDB::Parser::max_helix_size'};
}






#==========================================================================
#Parse PDB file


sub parse {

  my $self = shift;


  #Input PDB file must exist
  my $pdb_file = $self->infile;
  die "Input PDB file does not exist: $pdb_file --> $!\n" unless (-f $pdb_file);


  #Register here the current model
  my $model  = undef;

  my $fileType = File::Type->new->checktype_filename($pdb_file);
  my $infile = undef;

  if ($fileType =~ /gzip/) {
     open ($infile, "-|", "zcat $pdb_file") || die $!;
  }
  elsif ($fileType =~ /bzip2/) {
    open ($infile, "-|", "bzcat $pdb_file") || die $!;
  }
  elsif ($fileType =~ /octet-stream/) {
    open ($infile, "<", $pdb_file) || die $!;
  }
  else {
    die "Unknown file type for $pdb_file --> $fileType";
  }

  while(<$infile>) {

    chomp;


    #Get the PDB ID from the header section.
    if (/^HEADER\s+/) {
      my @headerContent = split(/\s+/, $_);

      #The PDB ID is the last element in the header
      $self->id($headerContent[$#headerContent]);
      $self->header($_);
    }



    #Get the resolution of the structure
    if (/^REMARK.+RESOLUTION\.\s+([0-9\.]+)\s+ANGSTROMS\./) {
      $self->resolution($1);
      $self->resolution_line($_);
    }
    elsif (/^REMARK.+RESOLUTION\.\s+NOT\s+APPLICABLE\./) {
      $self->resolution(100);  #Big number to indicate not applicable
    }



    #Identify if the PDB structure has models
    if (/^MODEL\s+(\d+)/) {
      $self->models->{$1} = 1;
      $model = $1;
    }



    #Get the info on subunits (DBREF record)q
    if (/^DBREF\s+/) {

      #Chain is char 13
      my $chain = $self->trim(substr($_, 12, 1));

      #Left position of sequence relative to PDB numbers (chars 15-18; right justified)
      my $lpdb =  $self->trim(substr($_, 14, 4));

      #Right position of sequence relative to PDB numbers (chars 21-24; right justified)
      my $rpdb =  $self->trim(substr($_, 20, 4));

      #The name of the external database (chars 27-32; left justified)
      my $xdb =   $self->trim(substr($_, 26, 6));

      #The external accession (chars 34-41; left justified)
      my $xacc =  $self->trim(substr($_, 33, 8));

      #The identificator code (chars 43-54; left justified)
      my $xid =   $self->trim(substr($_, 42, 12));

      #Left position relative to the external database protein (56-60; right justified)
      my $lxacc = $self->trim(substr($_, 55, 5));

      #Right position relative to the external database protein (63-67; right justified)
      my $rxacc = $self->trim(substr($_, 62, 5));


#      print "|$chain|,|$lpdb|,|$rpdb|,|$xdb|,|$xacc|,|$xid|,|$lxacc|,|$rxacc|\n\n";
#      <STDIN>;
#      next;


      #if segment has  coordinate 0 and length 0, ignore it!
      unless ($lpdb == 0 && $rpdb == 0) {
	push(@{ $self->dbref->{$chain} }, {lpdb=>$lpdb, rpdb=>$rpdb, xdb=>$xdb, xacc=>$xacc, xid=>$xid, lxacc=>$lxacc, rxacc=>$rxacc});
	push(@{ $self->dbref_lines->{$chain} }, $_);
      }

    }



    #Check for Modified residues
    if (/^MODRES/) {

      #Parse the info
      my $modres = substr($_, 12, 3);              #Cols 13-15
      my $chain  = substr($_, 16, 1);              #Col  17
      my $seqnum = $self->trim(substr($_, 18, 4)); #Cols 19-22
      my $mapres = substr($_, 24, 3);              #Cols 25-27


      #The same modified residue can map to more than one standard
      #amino acid
      $self->modres_dic->{$modres}->{$mapres} = 1;
      push (@{ $self->modres->{$chain} }, $_);
    }



    #Read the SEQRES section
#    if (/^SEQRES/) {
#      my @seqresLine = split (/\s+/, $_);
#      push (@{ $self->seqres->{$seqresLine[2]} }, \@seqresLine);
#
#    }



    #Read the coordinates of HELICES in this structure
    if (/^HELIX/) {

      #Extract the chain of the helix (char 20)
      my $chain = substr($_, 19, 1);
      push(@{ $self->helix->{$chain} }, $_);


      #Extract the helix parameters to simulate that
      #the reuslts of stride (useful for further processing).
      my $hNum = $self->trim(substr($_,  7, 3));  #helix number

      my $res1 = substr($_,  15, 3);              #C-terminal residue
      my $pos1 = $self->trim(substr($_,  21, 4)); #Position of C-terminal residue

      my $res2 = substr($_,  27, 3);           #N-terminal residue
      my $pos2 = $self->trim(substr($_,  33, 4)); #Position of C-terminal residue



      #Load helices in the same format they are parsed from stride
      $self->stride->{$chain}->{$hNum} = [$res1, $pos1, $res2, $pos2];

    }




    #Parse the ATOM section and extract the sequence
    if (/^ATOM/) {

      my $residue = substr($_, 17, 3);
      my $chain   = substr($_, 21, 1);
      my $resNum  = $self->trim(substr($_, 22, 4));


      #Store lines for the purpose of creating PDB files with a
      #specific chain
      push(@{ $self->atom->{$chain}->{$resNum} }, $_);


      #Store the sequence in 3-letter code
      $self->atom_seq3->{$chain}->{$resNum} = $residue;
    }



    #Parse the HETATM section
    if (/^HETATM/) {

      my $residue = substr($_, 17, 3);               #Char 18-20
      my $chain   = substr($_, 21, 1);               #Char 22
      my $resNum  = $self->trim(substr($_, 22, 4));  #Char 23-26

      #If residue is in the MODRES section, add the residue to the atom section.
      #
      #NOTE:
      #  It's ok to use $self->modres here, because by the time the ATOM and
      #  HETATM sections are reached, the MODRES section is already parsed.
      if (%{ $self->modres_dic } && exists $self->modres_dic->{$residue}) {

	#Store lines for the purpose of creating PDB files with the a
	#specific chain
	push(@{ $self->atom->{$chain}->{$resNum} }, $_);

	#Add the modified residue to the sequence extracted from ATOM
	$self->atom_seq3->{$chain}->{$resNum} = $residue;
      }
      else {

	#Ligands, waters, etc.  go in this section
	push(@{ $self->hetatm->{$chain}->{$resNum} }, $_);
      }
    }



    #parse the TER section
    if (/^TER/) {
      my $chain = substr($_, 21, 1);
      $self->ter->{$chain} = $_;
    }



    #parse the CONNECT section
    if (/^CONECT/) {
      push (@{ $self->connect }, $_);
    }



    #parse the MASTER line
    if (/^MASTER/) {
      $self->master($_);
    }

  }
  close $infile;



  #Verify that an ATOM section was found
  my $pdbFile = $self->infile;
  die "No ATOM section found in PDB file: $pdbFile " unless (%{ $self->atom });


  #Extract the different chain IDs in the structure
  @{ $self->chain_ids } = sort {$a cmp $b} keys %{ $self->atom };


  #Indicate if the structure has more than one model.
  #Note this parser does not deal with more than one model yet.
  if (%{ $self->models } && scalar( keys %{ $self->models }) > 1) {
    $self->has_models(1);
    return;
  }



  #Define if the PDB structure is chimeric chains
  $self->is_structure_chimeric();
  return if ($self->chimeric);



  #Identify if strcuture is a heterocomplex
  $self->is_structure_heterocomplex();


  #Notify if there are redundant modified residues that translate to more
  #than 1 regular amino acid
  $self->identify_redundant_modified_residues();



  #Get the amino acid sequences of all subunits based on the
  #SEQRES section
#  foreach my $chain (@{ $self->chain_ids }) {
#    $self->get_seqres_seq($chain);
#  }



  #If indicated by the user, renumber the residue positions in the
  #ATOM section for all chains.
  #
  #NOTE: This is done before correcting Helices so stride and
  #      hmmtop run on renumbered residues
  if ($self->renum_res_in_atom) {
    foreach my $chain (@{ $self->chain_ids }) {


      #renumber residue number in ATOM and fill atom_seq3_corrected
      $self->map_residue_to_new_position($chain);


      #Correct renumber positions of residues in the ATOM section.
      $self->correct_residue_positions($chain);


      #Correct Helix positions
      $self->correct_helix_coords($chain);


      #Correct MODRES if applicable
      $self->correct_modres_coords($chain);

    }
  }

#  print Data::Dumper->Dump([$self->helix, $self->helix_corrected],[qw(*helix *helix_corrected)]);
#  exit;




  #Get the amino acid sequence of the ATOM section in 1-letter code.
  #I do it here to take advantabe of (1) modified-residues translations
  #if available; and (2) the corrected atom positions if residues
  #were renumbered.
  $self->fill_hash_with_atom_seq1();

#  print Data::Dumper->Dump([$self->atom_seq1],[qw(*atom_seq1)]);
#  exit;




  #If necessary, run Stride and hmmtop to correct and replace the HELIX
  #section in the structure (To discard 2-residue helices, correct helices split
  #in two, etc.)
  if ($self->replace_helix) {

    #Clean the contents of output hash before
    $self->helix_corrected({});


    #Get HELIX assignments by stride
    $self->run_stride();

#    print Data::Dumper->Dump([$self->stride ], [qw( *stride )]);
#    exit;



    foreach my $chain (@{ $self->chain_ids }) {

      #Run hmmtop when helices will be replaced and residues renumbered.
      $self->run_hmmtop($chain) if ($self->renum_res_in_atom);

#      print Data::Dumper->Dump([ $self->hmmtop ], [qw( *hmmtop )]);
#      exit;


      #generate helix records with helices corrected
      $self->generate_helix_record($chain);
    }
  }


#  print $self->infile, "\n";
#  print Data::Dumper->Dump([$self->stride, $self->hmmtop],[qw(*stride *hmmtop)]);
#  print Data::Dumper->Dump([$self->helix],[qw(*helix)]);
#  exit;


}


#==========================================================================
#Generate new HELIX records y combining stride and hmmtop

sub generate_helix_record {

  my ($self, $chain) = @_;


  #
  #Before generating the helix record correct the helices by considering
  #hmmtop information. To this point the helix information already comes
  #from STRIDE.
  #
  #NOTE: This function is called when helix records will be replaced.
  #      So stride helices are waranteed to be available. But HMMTOP
  #      will only be used if the ATOM residues were renumbered.
  #

  #Solve intersections with TMS
  my @stride_hmmtop_helices = ();

  foreach my $hidx (sort {$a <=> $b} keys %{ $self->stride->{$chain} }) {

    my $thisHelix = $self->stride->{$chain}->{$hidx};

    $self->solve_overlap_w_hmmtop($chain, $thisHelix, \@stride_hmmtop_helices);


#    print Data::Dumper->Dump([$thisHelix, $self->hmmtop->{$chain}, \@stride_hmmtop_helices],
#			     [qw(*thisHelix *hmmtop *stride_hmmtop_helices)]);
#    <STDIN>;

  }


  #-----------------------------------------------------------------
  #To this point all HELICES are corrected but some may now be
  #redundant with other helices. So now keep all non overlapping helices.


  my %best_helices = ();

  #Sort helices by left position. This will allow solving overlaps easily
  my @sorted_helices = sort by_left_pos_and_size @stride_hmmtop_helices;

#  print Data::Dumper->Dump([$self->min_helix_size, $self->helix_separator],
#			   [qw(*min_helix_size *helix_separator)]);
#  print Data::Dumper->Dump([\@sorted_helices],[qw(*sorted_helices)]);
#  exit;




  my $helix_cnt = 1;
 HELIX:foreach my $idx (0 .. $#sorted_helices) {

    my $helix = $sorted_helices[$idx];


    #Ignore small helices
    next HELIX unless ($helix->{'len'} >= $self->min_helix_size);


    #First Helix with acceptable size is the reference
    if ($helix_cnt == 1 && ! exists $best_helices{$helix_cnt}) {
      $best_helices{$helix_cnt} = $helix;
      $helix_cnt++;
      next HELIX;
    }


    #Get the coordinates from the previous helix
    my $prev_helix = $best_helices{$helix_cnt - 1};

    if ($self->debug) {
      print "\n==================================================\n";
      print "Before:\n", Data::Dumper->Dump([$self->pdb->id, $prev_helix, $helix],[qw(*id *prevH *thisH)]), "\n";
    }

    #Ignore current helix if it has the same positions than the previous helix
    next HELIX if ($helix->{'lpos'} == $prev_helix->{'lpos'} &&
		   $helix->{'rpos'} == $prev_helix->{'rpos'});



    #If there is overlap on the right side, fix it here.
    if ($helix->{'rpos'} > $prev_helix->{'rpos'}  &&
	$prev_helix->{'rpos'} >= $helix->{'lpos'} &&
	$prev_helix->{'lpos'} <= $helix->{'lpos'}) {


      #
      #NOTE: If an overlap is solved by splitting a helix, leave
      #      'helix_separator' amino acids as linker between the
      #      two resulting helices.
      #


      my $diff  = $prev_helix->{'rpos'} - $helix->{'lpos'} + 1; #true overlap
      my $lDiff = $helix->{'lpos'} - $prev_helix->{'lpos'};     #left side
      my $rDiff = $helix->{'rpos'} - $prev_helix->{'rpos'};     #right side


      if ($self->debug) {
	print Data::Dumper->Dump([$diff, $lDiff, $rDiff ],[qw(*diff *lDiff *rDiff)]), "\n";
      }

      #Instersection is big enough to contain a TMS
      if ($diff >= $self->min_helix_size) {


	#If the overlap is greater than the non overlaping regions,
	#and the non-overlapping regions are smaller than min_helix_size,
	#then fuse the two helices.
	if (($diff >= $lDiff  && $diff >= $lDiff) &&
	    ($lDiff < $self->min_helix_size && $rDiff < $self->min_helix_size)) {

	  print "Large overlap. Fusing two helices:\n" if ($self->debug);

	  $self->fuse_helices(\%best_helices, $chain, $helix_cnt - 1, $helix);


	  if ($self->debug) {
	    print "After1:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1},
						   $best_helices{$helix_cnt}],
						  [qw(*id *prev *this)]);
	    <STDIN>;
	  }

	  next HELIX;
	}



	#If the left non-overlapping region is greater than min_helix_size,
	#cut the previous helix and add another one containing the overlap and
	#The rest of the right helix.
	if ($lDiff >= $self->min_helix_size) {

	  print "Large overlap. Cutting helix to the left:\n" if ($self->debug);

	  $self->split_previous_helix(\%best_helices, $chain, $helix_cnt, $diff, $helix);

	  if ($self->debug) {
	    print "After1:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],[qw(*id *prev *this)]);
	    <STDIN>;
	  }

	  $helix_cnt++;
	  next HELIX;
	}



	#If the right non-overlapping region is greater than min_helix_size,
	#cut the previous helix and add another one containing the overlap and
	#The rest of the right helix.
	if ($rDiff >= $self->min_helix_size) {

	  print "Large overlap. Cutting Helix to the right:\n" if ($self->debug);

	  $self->split_current_helix(\%best_helices, $chain, $helix_cnt, $helix);

	  if ($self->debug) {
	    print "After1:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],[qw(*id *prev *this)]);
	    <STDIN>;
	  }

	  $helix_cnt++;
	  next HELIX;
	}
      }


      #Overlapping region is smaller than a TMS
      else {

	#Conditions to fuse helices
	#1. If the non-overlapping regions on both sides are smaller than a TMS.
	#
	#2. If either non-overlapping region fits a TMS but not the other side.
	if (($lDiff <  $self->min_helix_size && $rDiff <  $self->min_helix_size) ||
	    ($lDiff >= $self->min_helix_size && $rDiff <  $self->min_helix_size) ||
	    ($lDiff <  $self->min_helix_size && $rDiff >= $self->min_helix_size)) {

	  print "Small overlap. Fusing two helices:\n" if ($self->debug);

	  $self->fuse_helices(\%best_helices, $chain, $helix_cnt - 1, $helix);


	  if ($self->debug) {
	    print "After:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1},
						   $best_helices{$helix_cnt}],
						  [qw(*id *prev *this)]);
	    <STDIN>;
	  }

	  next HELIX;
	}


	#If previous helix is larger, cut it
	elsif ($prev_helix->{'len'} > $helix->{'len'}) {

	  print "Small Overlap. Cutting helix to the left:\n" if ($self->debug);

	  $self->split_previous_helix(\%best_helices, $chain, $helix_cnt, $diff, $helix);

	  if ($self->debug) {
	    print "After1:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],[qw(*id *prev *this)]);
	    <STDIN>;
	  }

	  $helix_cnt++;
	  next HELIX;
	}

	else {
	  print "Small overlap. Cutting Helix to the right:\n" if ($self->debug);

	  $self->split_current_helix(\%best_helices, $chain, $helix_cnt, $helix);

	  if ($self->debug) {
	    print "After1:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],[qw(*id *prev *this)]);
	    <STDIN>;
	  }

	  $helix_cnt++;
	  next HELIX;
	}
      }
    }


    #Some times after cutting helices, the next helix in @sorted_helices may
    #start before than the most recent helix in %best_helices.
    #Deal with those cases here
    elsif ($helix->{'lpos'} < $prev_helix->{'lpos'}) {

      print "Atypical case: right helix starts before left helix.\n" if ($self->debug);

      #if current helix extends further to the right, the they right tail region
      #does not fit a TMS, fuse both helices.
      if ($helix->{'rpos'} > $prev_helix->{'rpos'}) {

	

	my $diff = $helix->{'rpos'} - $prev_helix->{'rpos'} - 1;

	#Fuse helices
	if ($diff < $self->min_helix_size) {

	  print "  Fuse right helix with previous helix:\n" if ($self->debug);

	  $self->fuse_helices(\%best_helices, $chain, $helix_cnt - 1, $helix);


	  if ($self->debug) {
	    print "After:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],[qw(*id *prev *this)]);
	    <STDIN>;
	  }
	}

	#Cut helice to the right
	else {

	  print "  Split right protruding helix:\n" if ($self->debug);

	  $self->split_current_helix(\%best_helices, $chain, $helix_cnt, $helix);


	  if ($self->debug) {
	    print "After:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],[qw(*id *prev *this)]);
	    <STDIN>;
	  }

	  $helix_cnt++;
	  next HELIX;
	}
      }
      else {
	if ($self->debug) {
	  print "  Do nothing!\n";
	  print "After:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],[qw(*id *prev *this)]);
	  <STDIN>;
	}
      }
    }  #atypical case


    #if the current helix starts less than 'helix_separator' amino acids after the previous helix,
    #shorten the longest helix such that there are 'helix_separator' amino acids separating the helices
    elsif (($helix->{'lpos'} - $prev_helix->{'rpos'}) > 0 &&
	   ($helix->{'lpos'} - $prev_helix->{'rpos'} - 1) < $self->helix_separator) {

      print "Contiguous helices:\n" if ($self->debug);

      $self->fix_contiguous_helices (\%best_helices, $chain, $helix_cnt, $helix);

      if ($self->debug) {
	print "After contiguous:\n", Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],
							[qw(*id *prev *this)]);
	<STDIN>;
      }

      $helix_cnt++;

    }


    #If the current helix is contained within the previos helix, do nothing.
    elsif ($helix->{'rpos'} <= $prev_helix->{'rpos'} && $helix->{'lpos'} >= $prev_helix->{'lpos'}) {

      if ($self->debug) {
	print "Helix inside: do nothing!\n",
	  Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],
			     [qw(*id *prev *this)]);
	<STDIN>;
      }
    }


    #Helix is not contained, nor overlaps, nor is continuous with the previous helix, so add it to results
    else {
      $best_helices{$helix_cnt} = $helix;

      if ($self->debug) {
	print "Helices is not contained, nor overlaps, nor contiguous:\n",
	  Data::Dumper->Dump([$self->pdb->id, $best_helices{$helix_cnt-1}, $best_helices{$helix_cnt}],[qw(*id *prev *this)]);
	<STDIN>;
      }

      $helix_cnt++;
    }
  }

#  print $self->id, "_$chain\n";
#  foreach my $h (sort {$a <=> $b} keys %best_helices) {
#    print Data::Dumper->Dump([$best_helices{$h}],[qw(*best_helix)]);
#  }
#  exit;





#  print Data::Dumper->Dump([$cleanHelices],[qw(*cleanHelices)]);
#  exit;



  #With short helices dealt with, proceed to generate the new Helix Records
  foreach my $helix_idx (sort {$a <=> $b} keys %best_helices) {

    #Get the coordinates of this helix
    my $stride_helix = $best_helices{$helix_idx};


    #Create initial empty line with HELIX tag
    my $line = sprintf ("%-80s", "HELIX");

    #Helix number
    my $hNum = sprintf ("%3s", $helix_idx);
    substr($line,  7, 3, $hNum);

    #Helix identifier
    my $hId  = sprintf ("%3s", $helix_idx);
    substr($line,  11, 3, $hId);



    #N-terminal residue
    substr($line,  15, 3, $stride_helix->{lres});

    #Chain of N-terminal residue
    substr($line,  19, 1, $chain);

    #Position of N-terminal residue
    my $pos1 = sprintf ("%4s", $stride_helix->{lpos});
    substr($line,  21, 4, $pos1);



    #C-terminal residue
    substr($line,  27, 3, $stride_helix->{rres});

    #Chain of c-terminal residue
    substr($line,  31, 1, $chain);

    #Position of c-terminal residtue
    my $pos2 = sprintf ("%4s", $stride_helix->{rpos});
    substr($line,  33, 4, $pos2);



    #Helix class: Right-handed alpha-helix
    substr($line,  38, 2, " 1");

    #Helix length
    my $hLen = sprintf ("%5s", $stride_helix->{len});
    substr($line,  71, 5, $hLen);

    push(@{ $self->helix_corrected->{$chain} }, $line);


#    my $len = length $line;
#    print "Length: $len\n|$line|\n";
#    <STDIN>;

  }

#  print "generate_helix_record: ", $self->id, "\n";
#  print Data::Dumper->Dump([$self->helix_corrected],[qw(*helix)]);
#  exit;
}






#==========================================================================
#when two helices do not overlap, but are contiguous, they look
#like a single helix. In this case the longest helix will be
#shortened to guarantee that $self->helix_separator  amino acids 
#separate the helices

sub fix_contiguous_helices {

  my ($self, $helices, $chain, $helix_idx, $helix) = @_;

  my $sep = $self->helix_separator;


  #Left helix is larger
  if ($helices->{$helix_idx - 1}->{'len'} >= $helix->{'len'}) {

    my $lpos_lHelix = $helix->{'lpos'} - $sep - 1;

    #fix left position of previous helix
    $helices->{$helix_idx - 1}->{'rpos'} = $lpos_lHelix;
    $helices->{$helix_idx - 1}->{'rres'} = $self->atom_seq3_corrected->{$chain}->{$lpos_lHelix};
    $helices->{$helix_idx - 1}->{'len'}  = $helices->{$helix_idx - 1}->{'rpos'} - $helices->{$helix_idx - 1}->{'lpos'} + 1;

    #Add the right helix as is to results
    $helices->{$helix_idx} = $helix;
  }


  #Right helix is larger
  else {

    my $lPos_rHelix = $helices->{$helix_idx - 1}{'rpos'} + $sep + 1;

    $helices->{$helix_idx}->{'lpos'} = $lPos_rHelix;
    $helices->{$helix_idx}->{'lres'} = $self->atom_seq3_corrected->{$chain}->{$lPos_rHelix};
    $helices->{$helix_idx}->{'rpos'} = $helix->{'rpos'};
    $helices->{$helix_idx}->{'rres'} = $helix->{'rres'};
    $helices->{$helix_idx}->{'len'}  = $helices->{$helix_idx}->{'rpos'} - $helices->{$helix_idx}->{'lpos'} + 1;
  }

}






#==========================================================================
#Cut the helix to the right at the end of the overlap with the left
#helix.


sub split_current_helix {

  my ($self, $helices, $chain, $helix_idx, $right_helix) = @_;


  my $left_helix = $helices->{$helix_idx - 1};


  #The position of the right helix that removes the overlap.
  my $lPos_rHelix  = $left_helix->{'rpos'} + $self->helix_separator + 1;


  #Correct left position from right helix and add it to the results.
  $helices->{$helix_idx}->{'lpos'} = $lPos_rHelix;
  $helices->{$helix_idx}->{'lres'} = $self->atom_seq3_corrected->{$chain}->{$lPos_rHelix};
  $helices->{$helix_idx}->{'rpos'} = $right_helix->{'rpos'};
  $helices->{$helix_idx}->{'rres'} = $right_helix->{'rres'};
  $helices->{$helix_idx}->{'len'}  = $helices->{$helix_idx}->{'rpos'} - $helices->{$helix_idx}->{'lpos'} + 1;


}





#==========================================================================
#Cut the helix to the left at the intersection with the
#helix to the right.

sub split_previous_helix {

  my ($self, $helices, $chain, $helix_idx, $overlap, $right_helix) = @_;


  my $left_helix = $helices->{$helix_idx - 1};


  #The positions that removes the overlap.
  my $prev_right_pos = $left_helix->{'rpos'}  - $overlap;
  my $this_left_pos  = $right_helix->{'lpos'} + $self->helix_separator + 1;


  #Remove overlap from previous helix
  $helices->{$helix_idx - 1}->{'rpos'} = $prev_right_pos;
  $helices->{$helix_idx - 1}->{'rres'} = $self->atom_seq3_corrected->{$chain}->{$prev_right_pos};
  $helices->{$helix_idx - 1}->{'len'}  = $helices->{$helix_idx - 1}->{'rpos'} - $helices->{$helix_idx - 1}->{'lpos'} + 1;


  #Correct left position from right helix and add it to the results.
  $helices->{$helix_idx}->{'lpos'} = $this_left_pos;
  $helices->{$helix_idx}->{'lres'} = $self->atom_seq3_corrected->{$chain}->{$this_left_pos};
  $helices->{$helix_idx}->{'rpos'} = $right_helix->{'rpos'};
  $helices->{$helix_idx}->{'rres'} = $right_helix->{'rres'};
  $helices->{$helix_idx}->{'len'}  = $helices->{$helix_idx}->{'rpos'} - $helices->{$helix_idx}->{'lpos'} + 1;

}





#==========================================================================
#Given an overlap to the right, this function fuses the helix to the right
#to the previous helix. In this case the index to the data of the previous
#helix is stored in variable $helix_idx.
#
#NOTE: This is an example of the info in $right_helix
#$right_helix = {
#                 'len'  => 22,
#                 'lpos' => 33,
#                 'rpos' => 54,
#                 'lres' => 'PRO',
#                 'rres' => 'MET'
#               );


sub fuse_helices {

  my ($self, $helices, $chain, $helix_idx, $right_helix) = @_;

  $helices->{$helix_idx}->{'rpos'} = $right_helix->{'rpos'};
  $helices->{$helix_idx}->{'rres'} = $self->atom_seq3_corrected->{$chain}->{$right_helix->{'rpos'}};
  $helices->{$helix_idx}->{'len'}  = $helices->{$helix_idx}->{'rpos'} - $helices->{$helix_idx}->{'lpos'} + 1;

}







#==========================================================================
#Determine whether the coordinates of a helix overlap with a stride-helix
#and or and inferred TMS.

sub solve_overlap_w_hmmtop {

  my ($self, $chain, $helix, $solved_helices) = @_;


  #Stride helix coordinates
  my $hLeft  = $helix->[1];
  my $hRight = $helix->[3];
  my $hLong  = $hRight - $hLeft + 1;


  #
  #For some unknown reason, STRIDE can produce large helices (e.g. 80+ residues long).
  #Got to solve that problem here
  #


  #if there is a very long helix, divide it by 2. The assumption is that by comparing to
  #TMS any overlap will be fixed.
  my @helices = ();
  if ($hLong >= $self->max_helix_size) {

    #first half of helix
    my $hl1   = $hLeft;
    my $hr1   = $hLeft + int($hLong / 2) - 1;
    my $hraa1 = $self->atom_seq3_corrected->{$chain}->{$hr1};
    my $hs1   = $hr1 - $hl1 + 1;     #size


    #second half of helix
    my $hl2 = $hr1 + 2;
    my $hlaa2 = $self->atom_seq3_corrected->{$chain}->{$hl2};
    my $hr2 = $hRight;
    my $hs2 = $hr2 - $hl2 + 1;


    push(@helices, {lres=>$helix->[0], lpos=>$hl1, rres=>$hraa1, rpos=>$hr1, len=>$hs1});
    push(@helices, {lres=>$hlaa2, lpos=>$hl2, rres=>$helix->[2], rpos=>$hr2, len=>$hs2});
  }
  else {
    push(@helices, {lres=>$helix->[0], lpos=>$hLeft, rres=>$helix->[2], rpos=>$hRight, len=>$hLong});
  }






  #----------------------------------------------------------------------
  #Identify Overlap with inferred TMS.

 HLX:foreach my $hlx (@helices) {

    my %newHelix = %{ $hlx };

  TMS:foreach my $tms (@{ $self->hmmtop->{$chain}->[1] }) {

      #TMS array organization:
      #0. Left position of TMS
      #1. Right position of TMS


#      print "Solve overlap with hmmtop:\n";
#      print Data::Dumper->Dump([$helix, $tms],[qw(*helix *tms)]);


      #Abort loop if TMS is completely beyond the helix
      last TMS if ($tms->[0] > $hRight);


      #Check if there is overlap between the helix and TMS
      if ($self->is_there_overlap([$hLeft, $hRight], $tms)) {

	#Detect if helix is contained in the TMS
	if ($tms->[0] < $hLeft && $tms->[1] > $hRight) {
	  $newHelix{'lres'} = $self->atom_seq3_corrected->{$chain}->{$tms->[0]};
	  $newHelix{'lpos'} = $tms->[0];
	  $newHelix{'rres'} = $self->atom_seq3_corrected->{$chain}->{$tms->[1]};
	  $newHelix{'rpos'} = $tms->[1];
	}
	else {

	  #correct helix position due to an overlap to the left
	  if ($tms->[0] < $hLeft) {
	    $newHelix{'lres'} = $self->atom_seq3_corrected->{$chain}->{$tms->[0]};
	    $newHelix{'lpos'} = $tms->[0];
	  }

	  #correct helix position due to an overlap to the right
	  if ($tms->[1] > $hRight) {
	    $newHelix{'rres'} = $self->atom_seq3_corrected->{$chain}->{$tms->[1]};
	    $newHelix{'rpos'} = $tms->[1];
	  }
	}

	$newHelix{'len'}  = $newHelix{'rpos'} - $newHelix{'lpos'} + 1;
      }
    }
    push(@{ $solved_helices }, \%newHelix);
  }
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

  if (($l1 < $l2) && ($r1 - $l2 > 0)) {
    return 1;
  }
  else {
    return 0;
  }
}





#==========================================================================
#The helix coordinates of many structures are horrible (i.e. short helices,
#contiguous helices, etc.) So run Stride to reassign 2D structure
#elements.


sub run_stride {

  my $self = shift;

  #Make sure that the pdb file of the chain exists
  my $infile = $self->infile;


  #Clean the contents of the output hash because if the ATOM section
  #is not renumbered, this hash will contain the original HELIX records
  #of the input PDB file.
  $self->stride({});


  #Run stride here
  my $outStrideFile = $self->tmpdir . "/" . $self->id . ".stride";
  unless (-f $outStrideFile) {
     system "stride -o -f$outStrideFile $infile";
   }
  die "No Stride output file found: $outStrideFile" unless (-f $outStrideFile);


  #Parse Stride Output
  my @aHelices = ();
  open (my $strideh, "<", $outStrideFile) || die $!;
  while(<$strideh>) {

    chomp;
    push (@aHelices, $_) if (/^LOC\s+AlphaHelix/);
  }
  close $strideh;

#  print Data::Dumper->Dump([\@aHelices ], [qw( *aHelices )]);
#  exit;


  #Save the alpha helices found by Stride
  my %hNum = ();
  foreach my $helix (@aHelices) {

    #Columns in stride ouput:
    # 0. Tag 'LOC'
    # 1. string 'AlphaHelix'
    # 2. three-letter code of first aa in the helix
    # 3. position of first amino acid in the helix
    # 4. Chain where the helix is located
    # 5. three-letter code of second aa in the helix
    # 6. position of last amino acid in the helix
    # 7. Chain where the helix is located
    my @helixData = split (/\s+/, $helix);

    $hNum{ $helixData[4] } = 1 unless (exists $hNum{ $helixData[4] });


    #If atoms where renumbered the positions of the helices must be corrected
    my ($left, $right) = undef;
    if ($self->renum_res_in_atom) {
      $left  = $self->map_res_position->{ $helixData[4] }->{ $helixData[3] };
      $right = $self->map_res_position->{ $helixData[4] }->{ $helixData[6] };
    }
    else {
      $left  = $helixData[3];
      $right = $helixData[6];
    }


    $self->stride->{ $helixData[4] }->{ $hNum{ $helixData[4] }} = [$helixData[2], $left, $helixData[5], $right];
    $hNum{ $helixData[4] }++;
  }


#  print "run_stride: ", $self->id, "\n";
#  print Data::Dumper->Dump([$self->stride],[qw(*stride)]);
#  exit;
}



#==========================================================================
#Verify that the sequences of the helices assigned by stride are correct.
#
#NOTE:
#  This function will give the mistaken results if the PDB file is original
#  and the user renumbers residues. This is because stride will report
#  helices based on the orginal residue numbers.
#

sub verify_stride_assignments {

  my ($self, $chain) = @_;


  #
  #Instead of using the sequence in a string to extract the helices
  #(effectively nullifying the residue number), which makes the first 
  #residue in the structure number 1, irrespectively of
  #it's position in the structure. Use the hash to take advantage of the
  #residue number in the atom sequence.
  #

  print $self->id, "_$chain\n";

  foreach my $idx (sort {$a <=> $b} keys %{ $self->stride->{$chain} }) {

    my $helix = $self->stride->{$chain}->{$idx};

    #extract the sequence of the helix here
    my $helix_seq = "";
    foreach my $pos ($helix->[1] .. $helix->[3]) {
      $helix_seq .= $self->atom_seq1->{$chain}->{$pos};
    }

    print "Helix [", join(", ", @$helix), "]: $helix_seq\n";
    <STDIN>;
  }
  exit;
}




#==========================================================================
#Download PDB files

sub download {

  my ($self, $pdb_ids, $overwrite) = @_;


  my $base_url = "http://files.rcsb.org/download";
  my $outdir = $self->tmpdir;

 PDB:foreach my $id (@{ $pdb_ids }) {

    my $name     = "${id}.pdb";
    my $out_file = "$outdir/$name";
    my $ftp_url  = "$base_url/$name";

    unless ($overwrite) {
      next PDB if (-f $out_file);
    }

    print "Downloading $name\n";

    system "wget -O $out_file $ftp_url";
    die "Could not download: $name" unless (-f $out_file);

  }
}


#==========================================================================
#Run HMMTOP in case the user needs it

sub run_hmmtop {

  my ($self, $chain, $cust_name) = @_;

  my $fasta_fname = undef;
  if ($cust_name) {
    $fasta_fname = "${cust_name}.faa";
  }

  #Get chain sequence in fasta format
  $self->save_atom_sequence_in_fasta_format($chain, $self->tmpdir, $fasta_fname);


  my $infile = "";
  if ($fasta_fname) {
    $infile =  $self->tmpdir . "/" . $fasta_fname;
  }
  else {
    $infile = $self->tmpdir . "/" . $self->id . "_${chain}.faa";
  }
  die "Error: fasta file does not exist -> $infile" unless (-f $infile);


  #run hmmtop
  my $output = `hmmtop -if=$infile -of=-- -sf=FAS -pi=spred -is=pseudo 2> /dev/null`;
  chomp($output);


  my @parts = split(/\s+/, $output);
  my @tms = ();


  #If no TMS are inferred return and indicate the error
  if ($parts[4] == 0) {
    $self->hmmtop->{$chain} = [$parts[4], \@tms];
    return;
  }


  #get the corrdinates of each TMS
  for (my $i=5; $i <= $#parts - 1; $i += 2) {

    die "Error: coult not retrieve both TMS cordinates:\n$output\nCurrent Idx=$i" unless ($parts[$i] && $parts[$i+1]);
    push (@tms, [$parts[$i], $parts[$i+1]]);
  }

  $self->hmmtop->{$chain} = [$parts[4], \@tms];

}




#==========================================================================
#Verify that the TMS is correct. Extract the sequence of the TMS and
#compare it to the output of running HMMTOP. This function is for
#debugging purposes only

sub verify_hmmtop_predictions {

  my ($self, $chain) = @_;

  my $atom_seq = $self->extract_atom_seq_1($chain);

  foreach my $tms (@{ $self->hmmtop->{$chain} }) {

    my $tms_length = $tms->[1] - $tms->[0] + 1;

    my $tms_seq = substr($atom_seq, $tms->[0] - 1, $tms_length);

    print "TMS [", join(", ", @$tms), "]: $tms_seq\n";
    <STDIN>;
  }
  exit;
}







#==========================================================================
#Correct the positions of helices based on the renumbering of residues
#in the ATOM section.


sub correct_modres_coords {

  my ($self, $chain) = @_;


  #Abort if there are no MODRES lines in structure
  return unless (exists $self->modres->{$chain} && @{ $self->modres->{$chain} });


  foreach my $modres_line (@{ $self->modres->{$chain} }) {


#    print "Old modres:\n|$modres_line|\n\n";


    #Position of modified residue
    my $pos  = $self->trim(substr($modres_line, 18, 4));
    my $npos = $self->map_res_position->{$chain}->{$pos};


#    next RES unless ($npos);
    die "Could not get new position for modres: |$modres_line|" unless ($npos);


    #format postion for instertiong into PDB file
    my $spos = sprintf("%4s", $npos);


    my $newModres = $modres_line;
    substr($newModres, 18, 4, $spos);

#    print "New MODRES:\n|$newModres|\n\n";
#    <STDIN>;


    push (@{ $self->modres_corrected->{$chain} }, $newModres);
  }
}





#==========================================================================
#Correct the positions of helices based on the renumbering of residues
#in the ATOM section.


sub correct_helix_coords {

  my ($self, $chain) = @_;


  foreach my $helix_line (@{ $self->helix->{$chain} }) {


#    print $self->id, "\n";
#    print "Old helix:\n|$helix_line|\n\n";


    #left position of initial residue based on ATOM sequence (chars 22-25)
    my $lpos  = $self->trim(substr($helix_line, 21, 4));
    my $nlpos = $self->map_res_position->{$chain}->{$lpos};
    my $slpos = sprintf("%4s", $nlpos);


    #right position of last residue based on ATOM sequence (chars 34-37)
    my $rpos  = $self->trim(substr($helix_line, 33, 4));
    my $nrpos = $self->map_res_position->{$chain}->{$rpos};
    my $srpos = sprintf("%4s", $nrpos);


    #Length of the helix (chars 72-76)
    my $hlen = $self->trim(substr($helix_line, 71, 5));


    #Verify that the length of the helix was not changed
    my $nhlen = $nrpos - $nlpos + 1;
    unless ($hlen == $nhlen) {
      die "Original helix length $hlen, is not equal to new length $nhlen";
    }

    my $newHelix = $helix_line;
    substr($newHelix, 21, 4, $slpos);
    substr($newHelix, 33, 4, $srpos);

#    print "New Helix:\n|$newHelix|\n\n";
#    <STDIN>;


    push (@{ $self->helix_corrected->{$chain} }, $newHelix);
  }

#  print "correct_helix_coords:\n";
#  print Data::Dumper->Dump([$self->helix, $self->helix_corrected],[qw(*helix *helix_corrected)]);
#  exit;

}




#==========================================================================
#Reset the residue numbering in the ATOM section, this will allow runing
#the ATOM seq against HMMTOP and correctly map the TMS to ATOMS in the
#structure.
#
#This function also generates an amino acid sequence in 3-letter code
#that starts in amino acid 1. This new residue numbers will also allow
#the correction of alpha-helices positions in the HELIX section.

sub map_residue_to_new_position {

  my ($self, $chain) = @_;

  my $new_position = 1;
  foreach my $res_idx (sort {$a <=> $b} keys %{ $self->atom_seq3->{$chain} }) {

    #New position
    $self->map_res_position->{$chain}->{$res_idx} = $new_position;
#    print "ATOM: $res_idx  -> SEQRES: $new_position\n";
#    <STDIN>;

    #Sequence with new residue numbering
    $self->atom_seq3_corrected->{$chain}->{$new_position} = $self->atom_seq3->{$chain}->{$res_idx};

    $new_position++;
  }
}




#==========================================================================
#Correct the atom positions based based on renumbered residues

sub correct_residue_positions {

  my ($self, $chain) = @_;


  foreach my $idx (sort {$a <=> $b} keys %{ $self->atom->{$chain} }) {

    my $new_idx = $self->map_res_position->{$chain}->{$idx};


    #Change the original amino acid number to the match in seqres
    my @atom_data = @{ $self->atom->{$chain}->{$idx} };


    #Format the residue number justified to the right and in 4 characters
    my $newRes = sprintf("%4s", $new_idx);

    foreach my $line (@atom_data) {

      my $len = length $line;
      my $newLine = $line;

      #The residue number is found in chars 23-26
      substr($newLine, 22, 4, $newRes);
      my $len2 = length $line;

      unless ($len == $len2) {
	die "Lengths of old and new ATOM lines are different\n|$line|\n|$newLine|\n";
      }

      push (@{ $self->atom_corrected->{$chain}->{$new_idx} }, $newLine);

    }
  }
}



#==========================================================================
#Save the ATOM section with renumbered amino acids into a PDB file for
#a given chain.

sub save_corrected_atoms {

  my ($self, $chain, $outdir) = @_;

#  print "save_corrected_atoms:\n";
#  print Data::Dumper->Dump([$self->helix],[qw(*helix)]);
#  exit;


  #If the chain is chimeric, thare could be no ATOM coordinates
  #for the chain, so abort if the structure is chimeric.
  return if ($self->chimeric);


  #get ouput directory
  $outdir = $self->tmpdir unless ($outdir);
  system "mkdir -p $outdir" unless (-d $outdir);


  my $name    = $self->id . "_$chain";
  my $outfile = "$outdir/${name}.pdb";


  #Do not overwrite output file if already exists
  return if (-f $outfile);


  open (my $pdb, ">", $outfile) || die $!;

  #THe header of the PDB file
  print $pdb $self->header, "\n";


  #the resolution of the structure
  if ($self->resolution_line) {
    print $pdb $self->resolution_line, "\n";
  }



  #The DBREF section to extract SwissProt IDs
  if (exists $self->dbref_lines->{$chain}) {
    print $pdb join("\n", @{ $self->dbref_lines->{$chain} }), "\n";
  }




  #Print MODRES section
  if (exists $self->modres_corrected->{$chain}) {
    print $pdb join("\n", @{ $self->modres_corrected->{$chain} }), "\n";
  }



  #print corrected helices
  if (exists $self->helix_corrected->{$chain}) {
    print $pdb join("\n", @{ $self->helix_corrected->{$chain} }), "\n";
  }



  #Print the residues in the chain
  foreach my $res_num (sort {$a <=> $b} keys %{ $self->atom_corrected->{$chain} }) {
    print $pdb join("\n", @{ $self->atom_corrected->{$chain}->{$res_num} }), "\n";
  }


  #The CONNECT section if it exists
  if (@{ $self->connect }) {
    print $pdb join("\n", @{ $self->connect }), "\n";
  }


  #The MASTER line if it exists
  print $pdb $self->master, "\n" if ($self->master);



  print $pdb "END\n";

  close $pdb;
}






#==========================================================================
#Fill the hash $self->atom_seq1

sub fill_hash_with_atom_seq1 {

  my $self = shift;

  #For all chains in structure
  foreach my $chain (@{ $self->chain_ids }) {


    #Residue positions depend on whether or not the ATOM section
    #was corrected.
    my $resPositions = undef;
    if ($self->renum_res_in_atom) {
      $resPositions = $self->atom_seq3_corrected->{$chain};
    }
    else {
      $resPositions = $self->atom_seq3->{$chain};
    }

#    print Data::Dumper->Dump([$resPositions],[qw(*resPositions)]);
#    exit;



    #for all amino acids in the chain
    foreach my $pos (keys %{ $resPositions }) {

      my $aa3 = $resPositions->{$pos};
      my $aa1 = "";

      #Search in standard amino acids first
      if ($self->aaDic3to1($aa3)) {
	$aa1 = $self->aaDic3to1($aa3);
      }

      #Search in modified residues
      elsif (exists $self->modres_dic->{$aa3}) {
	my @aa_tmp = keys %{ $self->modres_dic->{$aa3} };
	$aa1 = $self->aaDic3to1($aa_tmp[0]);
      }

      #Amino acid was not found, report an error.
      unless ($aa1) {
	print "No one letter code for aa $aa3 at pos $pos in: ", $self->infile, "\n";
	exit;
      }

      $self->atom_seq1->{$chain}->{$pos} = $aa1;
    }
  }
}




#==========================================================================
#Identify if there are redundant modified residues that map to more than
#one regular residues

sub identify_redundant_modified_residues {

  my $self = shift;

  #Default value
  my $redundant = 0;


  if (%{ $self->modres }) {

    #All modified residues in file
    my @modres = keys %{ $self->modres };


    #mappings to regular residues
  MRES:foreach my $res (@modres) {

      my @map = keys %{ $self->modres_dic->{$res} };

      #Identify if there is redundancy here
      if (scalar @map > 1) {
	my $redundant = 1;
	last MRES;
      }
    }
  }

  $self->has_red_modres($redundant);
}


#==========================================================================
#create individual PDB files for all the chains in this structure

sub split_chains {

  my ($self, $outdir) = @_;


  unless ($outdir) {
    $outdir = $self->tmpdir;
  }


  foreach my $chain (@{ $self->chain_ids }) {

    if ($self->renum_res_in_atom) {
      $self->save_corrected_atoms($chain, $outdir);
    }
    else {
      $self->create_pdb_for_chain($chain, $outdir);
    }
  }
}




#==========================================================================
#Create a PDB file with the atoms of only one chain

sub create_pdb_for_chain {

  my ($self, $chain, $outdir) = @_;


  $outdir = $self->tmpdir unless ($outdir);
  system "mkdir -p $outdir" unless (-d $outdir);


  system "mkdir -p $outdir" unless (-d $outdir);
  my $outfile = "$outdir/" . $self->id . "_${chain}.pdb";


  #Do not overwrite output file if it already exists
  return if (-f $outfile);


  open (my $pdb, '>', $outfile) || die $!;


  #The HEADER
  print $pdb $self->header, "\n" if ($self->header);


  #the resolution of the structure
  if ($self->resolution_line) {
    print $pdb $self->resolution_line, "\n";
  }


  #The DBREF section to extract SwissProt IDs
  print $pdb join("\n", @{ $self->dbref_lines->{$chain} }), "\n";



  #The MODRES section
  if (exists $self->modres->{$chain} && @{ $self->modres->{$chain} }) {
    print $pdb join ("\n", @{ $self->modres->{$chain} }), "\n";
  }



  #the HELIX section
  if ($self->replace_helix) {
    print $pdb join ("\n", @{ $self->helix_corrected->{$chain} }), "\n";
  }
  elsif (exists $self->helix->{$chain} && @{ $self->helix->{$chain} }) {
    print $pdb join ("\n", @{ $self->helix->{$chain} }), "\n";
  }


  #The ATOM section
  foreach my $res (sort {$a <=> $b} keys %{ $self->atom->{$chain} }) {
    print $pdb join("\n", @{ $self->atom->{$chain}->{$res} }), "\n";
  }


  #The CONNECT section if it exists
  if (@{ $self->connect }) {
    print $pdb join("\n", @{ $self->connect }), "\n";
  }


  #The MASTER line if it exists
  print $pdb $self->master, "\n" if ($self->master);


  #The END of pdb file
  print $pdb "END\n";


  close $pdb;

}



#==========================================================================
#Extract sequence in 1-letter code from the ATOM section

sub extract_atom_seq_1 {

  my ($self, $chain) = @_;

  my $seq = "";
  foreach my $idx (sort {$a <=> $b} keys %{ $self->atom_seq1->{$chain} }) {
    $seq .= $self->atom_seq1->{$chain}->{$idx};
  }

  my $msg = "PDB: " . $self->id .
    "_$chain. Error: Could not extract seq from ATOM section";
  die $msg unless ($seq);

  return $seq;

}


#==========================================================================
#Get the sequence of each chain in three- and  one-letter format

sub get_seqres_seq {

  my ($self, $chain) = @_;


  #return if no SEQRES section is present in the structure
  return  unless ($self->seqres && %{ $self->seqres });


  #The complete SEQRES record for this chain. This is a list of
  #ordered arrays, where each array contains a line of the SEQRES
  #section:
  #  Third column is chain.
  #  Fourth column is length of chain
  #  Fifth to end of line are amino acids in 3-letter codes
  foreach my $line (@{ $self->seqres($chain) }) {


    #Skip first two columns in line (SEQRES and sequence line number)
    foreach my $idx (4 .. $#{ $line }) {

      if ($self->aaDic3to1($line->[$idx])) {

	#The sequence in 3-letter code
	push (@{ $self->sres_seq3->{$chain} }, $line->[$idx]);

	#The sequence in 1-letter code
	$self->sres_seq1->{$chain} .= $self->aaDic3to1($line->[$idx]);	
      }


      #Get aa from the MODRES record
      if (exists $self->modres->{$line->[$idx]}) {

	my @aa_tmp = keys %{ $self->modres->{$line->[$idx]} };

#	if ($self->has_red_modres) {
#	  my $f = $self->infile;
#	  die "Error: ambiguous residures in MODRES --> $f";
#	}

	#the sequence in 3-letter code
	push (@{ $self->sres_seq3->{$chain} }, $aa_tmp[0]) if (@aa_tmp);


	#The sequence in 1-letter code
	$self->sres_seq1->{$chain} .= $self->aaDic3to1($aa_tmp[0]) if (@aa_tmp);
      }


      #Ignore Amino group at the end of the sequence
      elsif ($line->[$idx] eq "NH2" && $idx == $#{ $line }) { }


      #Ignore Acetyl group at the begining of the sequence
      elsif ($line->[$idx] eq "ACE" && $idx == 4 && $line->[1] == 1) { }


      #The unexpected cases
      else {
#	print "No one letter code for aa ", $line->[$idx], "\n";
#	print "In line: ", join (" ", @$line), "\n";
#	exit;
      }
    }
  }
}


#==========================================================================
#Determine if the structure is a heterocomplex

sub is_structure_heterocomplex {

  my $self = shift;

  #No need to look at protein IDs if structure is chimeric.
  #A chimeric protein is a heterocomplex by definition.
  if ($self->chimeric) {
    $self->heterocomplex(1);
    return;
  }


  #By default structure is not a heterocomplex
  my $heterocomplex = 0;


  #Return if no DBREF section was present
  return $heterocomplex unless ($self->dbref && %{ $self->dbref });


  my %xdb_ids = ();
  foreach my $chain (@{ $self->chain_ids }) {

    #extract the external DB id
    my $xid = $self->dbref->{$chain}->[0]->{xacc};

    if ($xid) {
      $xdb_ids{ $self->dbref->{$chain}->[0]->{xacc} } = 1;
    }
    else {

      #Some structures contain atoms for more chains than described in the
      #section DBREF (e.g. 3DTX.pdb). In these cases, just ignore the chain.

#      print "Could not access DBREF record for chain $chain in: ", $self->infile, "\n";
#      print Data::Dumper->Dump([$self->dbref],[qw(*dbref)]);
#      die "Error at ->";
    }
  }


#  print Data::Dumper->Dump([\%xdb_ids],[qw(*xdb_ids)]);
#  exit;


  #If there is more than one key in hash %xdb_ids, that means
  #that the structure is a heterocomplex.
  unless ($heterocomplex) {
    if (scalar(keys %xdb_ids) > 1) {
      $heterocomplex = 1;
    }
  }

  $self->heterocomplex($heterocomplex);
}


#==========================================================================
#After reading the DBREF section, determine whether the structure
#contains a chimeric arrangement of proteins.

sub is_structure_chimeric {

  my $self = shift;

  #By defult is no chimeric
  my $chimeric = 0;


  #return if the structure does not have DBREF section
  return $chimeric unless ($self->dbref && %{ $self->dbref });


  #If a chimeric subunit is found, flag it as chimeric and return.
  foreach my $chain (@{ $self->chain_ids }) {

    #sometimes not all chains in the atom section are annotated
    #in the DBREF section
    next unless (exists $self->dbref->{$chain});


    my $length = scalar @{ $self->dbref->{$chain} };

    if ($length > 1) {
      $chimeric = 1;
      last;
    }
  }
  $self->chimeric($chimeric);
}




#==========================================================================
#The dicctionary to translate amino acids from 3-letter to 1-letter code.

sub aaDic3to1 {

  my ($self, $threeLetterCode) = @_;

  my %dic = ( 'ALA' => 'A',
	      'ARG' => 'R',
	      'ASN' => 'N',
	      'ASP' => 'D',
	      'ASX' => 'B',
	      'CYS' => 'C',
	      'GLU' => 'E',
	      'GLN' => 'Q',
	      'GLX' => 'Z',
	      'GLY' => 'G',
	      'HIS' => 'H',
	      'ILE' => 'I',
	      'LEU' => 'L',
	      'LYS' => 'K',
	      'MET' => 'M',
	      'PHE' => 'F',
	      'PRO' => 'P',
	      'SER' => 'S',
	      'THR' => 'T',
	      'TRP' => 'W',
	      'TYR' => 'Y',
	      'VAL' => 'V'
	    );

  return $dic{$threeLetterCode};

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
#Save the protein sequence of a chain in fasta format.
#If output directory is not given it prints the sequence to the
#current directory.

sub save_atom_sequence_in_fasta_format {

  my ($self, $chain, $outdir, $fileName) = @_;

  $outdir = '.' unless ($outdir);

  unless ($chain && $outdir) {
    die "Error: Both chain ID and output directory are mandatory!";
  }

  system "mkdir $outdir" unless (-d $outdir);


  #Save chain sequence in FASTA format
  my ($fasta_id, $file_name, $seq_file) = "";
  if ($fileName) {
    my ($name, $ext) = split (/\./, $fileName);
    $fasta_id  = $name;
    $file_name = $fileName;
  }
  else {
    $fasta_id  = $self->id . "_${chain}";
    $file_name = "${fasta_id}.faa";
  }
  $seq_file  = "$outdir/$file_name";

  unless (-f $seq_file) {
    open (my $faaFile, ">", $seq_file) || die $!;
    print $faaFile ">$fasta_id\n";
    print $faaFile $self->extract_atom_seq_1($chain), "\n";
    close $faaFile;
  }

  #Verify that file was created successfully and is not empty
  die "Could not create fasta file: $file_name" unless (-f $seq_file && ! -z $seq_file);

}



#==========================================================================
#Download protein sequence from UNIPROT or other related database.
#Valid values so far for external databases: UNP, PDB, GB

sub download_prot_seq {

  my ($self, $xdb, $xid) = @_;

  #Temporary fix to proteins without uniprot ID in Kevin structures
  my %exceptions = ('R9G739' => 'A0A0H3NW57');

  #Validate external database
  unless ($xdb && (uc $xdb eq "UNP" || uc $xdb eq "PDB" || uc $xdb eq "GB")) {
    die "Missing or wrong external DB. Valid values: UNP, PDB, GB";
  }


  my $spObj  = Bio::DB::SwissProt->new();

  my $seqObj = undef;


  #------------------------------------------------------------
  #Process UNP

  if (uc $xdb eq "UNP") {

    #Just one accession. ref returns undef if argument is a scalar
    unless (ref $xid) {

      my $outfile = $self->tmpdir . "/${xid}.faa";

      unless (-f $outfile ) {
	my $tmpID = (exists $exceptions{$xid})? $exceptions{$xid} : $xid;

	my $seqObj = $spObj->get_Seq_by_acc($tmpID);

	if ($seqObj) {
	  my $seq    = ">$xid\n" . $seqObj->seq . "\n";
	  $self->save_sequence($xid, $seq);
	}
	else { die "No sequence found for: $xid" }
      }
    }


    #Multiple accessions
    elsif (ref $xid eq "ARRAY") {

      foreach my $id (@{ $xid }) {

	my $outfile = $self->tmpdir . "/${id}.faa";

	unless (-f $outfile) {
	  my $seqObj = $spObj->get_Seq_by_acc($id);

	  if ($seqObj) {
	    my $seq    = ">$id\n" . $seqObj->seq . "\n";
	    $self->save_sequence($id, $seq);
	  }
	  else { die "No sequence found for: $id" }
	}
      }
    }

    else { die "Only a scalar or array reference are accepted: (", ref $xid, ") --> $xid" }
  }



  #------------------------------------------------------------
  #Process GB as EMBL

  elsif (uc $xdb eq "GB") {

    #If ref $xid is undef it means it's a scalar value
    if (!(ref $xid) || ref $xid eq "SCALAR") {

      my $ids = $spObj->id_mapper(-from => "EMBL", -to => "ACC", -ids => $xid);

      if ($ids) {

	my @seqIDs = @{ $ids->{$xid} };

	#The GB id ($xid) maps to only one uniprot IDs
	if (scalar @seqIDs == 1) {

	  my $outfile = $self->tmpdir . "/${xid}.faa";

	  unless (-f $outfile) {
	    my $seqObj = $spObj->get_Seq_by_acc($seqIDs[0]);
	    my $seq    = ">$xid\n" . $seqObj->seq . "\n";
	    $self->save_sequence($xid, $seq);
	  }
	}

	#GB ID maps to more than one Uniprot ID
	else {
	  die $self->id, ": $xid is chimeric: $xid maps to muliple uniprot accessions";
	}
      }
    }

    #Still need when dealing with multiple GB's
    elsif (ref $xid eq "ARRAY") { die "Unsupported option"; }
  }


  #------------------------------------------------------------
  #Process PDB IDs

  elsif (uc $xdb eq "PDB") {

    my $ids = $spObj->id_mapper(-from => "PDB_ID", -to => "ACC", -ids => $xid);

    if ($ids) {

#      print Data::Dumper->Dump([$xid, $ids ], [qw(*xid *ids )]);
#      exit;

      my @seqIDs = @{ $ids->{$xid} };

      if (scalar @seqIDs == 1) {

	my $outfile = $self->tmpdir . "/${xid}.faa";

	unless (-f $outfile ) {

	  my $seqObj = $spObj->get_Seq_by_acc($seqIDs[0]);

	  if ($seqObj) {
	    my $seq    = ">$xid\n" . $seqObj->seq . "\n";
	    $self->save_sequence($xid, $seq);
	  }
	  else { die "No sequence found for: $xid" }
	}
      }
      else { die "PDB ID maps to multiple UniProt IDs" }
    }
    else  {die "Could not find Uniprot IDs for PDB ID: $xid"; }

#   foreach my $id (keys %{ $ids }) {
#
#     my @eqIDs = @{ $ids->{$id} };
#
#     foreach my $mapID (@eqIDs) {
#	my $seqObj = $spObj->get_Seq_by_acc($mapID);
#
#	if ($seqObj) {
#	  my $seq    = ">$mapIDq\n" . $seqObj->seq . "\n";
#	  $self->save_sequence($mapID, $seq);
#	}
#	else  { die "No sequence found for: $id" }
#     }
  }
}




#==========================================================================
#Print a sequence file

sub save_sequence {

  my ($self, $id, $seq) = @_;

  my $outfile = $self->tmpdir . "/${id}.faa";

  open (my $outh, ">", $outfile) || die $!;
  print $outh $seq;
  close $outh;
}




#==========================================================================
#Sort helices according to left position and size

sub by_left_pos_and_size {

  if ($a->{'lpos'} == $b->{'lpos'}) {
    $b->{'len'}  <=> $a->{'len'};
  }
  else {
    $a->{'lpos'} <=> $b->{'lpos'};
  }
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

