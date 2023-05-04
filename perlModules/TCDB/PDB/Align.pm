package TCDB::PDB::Align;


use warnings;
no warnings;

use strict;
use Data::Dumper;

use Bio::SeqIO;
use R2::PDB::Parser;

#use Bio::Structure::IO::pdb;
use Class::Struct;




#
#The functions in this Class align 2 structures. It is assumed that
#The structures contain only the chain that will be alignmed.
#
#If the structure has multiple chain, Then you should first
#cut it using the funcion:
#  R2::PDB::Parser->create_pdb_for_chain
#
#When structures are aligned with 'superpose' the coordinates
#of the aligned structure are saved for later visualization of
#the quality of the alignment.
#


struct ("R2::PDB::Align" =>
	{
	 'indir'          => '$',  #Directory where both input structures are located.
	 'outdir'         => '$',  #Directory where results of the alignment will be placed.
	 'rep_unit'       => '$',  #Number of TMS in the repeat unit, the alignment must cover at least these TMS
	 'max_resolution' => '$',  #Maximum acceptable resolution in both structures to perform comparisons
	 'min_aa_aln_hlx' => '$',  #Minimum number of aligned residues to considered that two helices where aligned
	 'hmmtop_seqs'    => '%',  #hmmtop inferences for the uniprot sequences associated with both structures.
	 'fam_seqs_file1' => '$',  #file with uniprot sequences of the family of pdb1 in fasta format
	 'fam_seqs_file2' => '$',  #file with uniprot sequences of the family of pdb2 in fasta format
	 'atom2prot_map'  => '%',  #map of residue equivalences between the ATOM seq and the uniprote seq
	 'rmsd'           => '%',  #The RMSD of the alignement
	 'len_cutoff'     => '$',  #Minimum acceptable length for the proteins being compared
	 'cov_cutoff'     => '$',  #Minimum acceptable coverage for the alignment
	 'nres_cutoff'    => '$',  #Minimum number of aligned residues
	 'algorithm'      => '$',  #Algorithm that will be used to align the structures (superpose, cealign, all)
	 'debug'          => '$',  #Variable that facilites printing debug messages
	 'id1'            => '$',  #ID to use for first structure  (overrides ID extracted from structure)
	 'id2'            => '$',  #ID to use for second structure (overrides ID extracted from structure)
	 'pdb1'           => 'R2::PDB::Parser',   #First  structure to be aligned
	 'pdb2'           => 'R2::PDB::Parser'    #Second structure to be aligned
	}
       );


#==========================================================================
#Override accessor for algorithms to validate possible values:
#  superpose, cealign or all

sub max_resolution {

  my ($self, $value) = @_;


  my $default_value = 4.0;


  if ($value) {

    #validate
    unless ( $value < 6.0 ) {
      die "Error: Max acceptable structure resolution is 6.0";
    }

    $self->{'R2::PDB::Align::max_resolution'} = $value;
  }

  unless ($self->{'R2::PDB::Align::max_resolution'}) {
    $self->{'R2::PDB::Align::max_resolution'} = $default_value;
  }

  return $self->{'R2::PDB::Align::max_resolution'};
}




#==========================================================================
#Override accessor for algorithms to validate possible values:
#  superpose, cealign or all

sub algorithm {

  my ($self, $value) = @_;


  my $default_prog = 'superpose';


  if ($value) {

    my $prog = lc $value;

    #validate
    unless ( $prog =~ /^(superpose|cealign|all)$/) {
      die "Error: algorithm takes only 3 values (superpose, cealign or all)";
    }

    $self->{'R2::PDB::Align::algorithm'} = $prog;
  }

  unless ($self->{'R2::PDB::Align::algorithm'}) {
    $self->{'R2::PDB::Align::algorithm'} = $default_prog;
  }

  return $self->{'R2::PDB::Align::algorithm'};
}




#==========================================================================
#Override accessor for cov_cutoff to set a default cutoff

sub min_aa_aln_hlx {

  my ($self, $cutoff) = @_;

  my $default_cutoff = 3;

  if ($cutoff) {

    #validate
    unless ( $cutoff =~ /^[0-9]+$/ && $cutoff >= $default_cutoff) {
      die "Error: min_aa_aln_hlx should be at least $default_cutoff\n";
    }

    $self->{'R2::PDB::Align::min_aa_aln_hlx'} = $cutoff;
  }

  unless ($self->{'R2::PDB::Align::min_aa_aln_hlx'}) {
    $self->{'R2::PDB::Align::min_aa_aln_hlx'} = $default_cutoff;
  }

  return $self->{'R2::PDB::Align::min_aa_aln_hlx'};
}




#==========================================================================
#Override accessor for cov_cutoff to set a default cutoff

sub cov_cutoff {

  my ($self, $cutoff) = @_;

  my $default_cutoff = 0.3;

  if ($cutoff) {

    #validate
    unless ( $cutoff =~ /^[0-9\.]+$/ && $cutoff > 0 && $cutoff <= 1) {
      die "Error: cov_cutoff should be an number between 0 and 1\n";
    }

    $self->{'R2::PDB::Align::cov_cutoff'} = $cutoff;
  }

  unless ($self->{'R2::PDB::Align::cov_cutoff'}) {
    $self->{'R2::PDB::Align::cov_cutoff'} = $default_cutoff;
  }

  return $self->{'R2::PDB::Align::cov_cutoff'};
}






#==========================================================================
#Override accessor for len_cutoff and set defaul value

sub len_cutoff {

  my ($self, $cutoff) = @_;

  my $default_cutoff = 20;

  if ($cutoff) {

    #validate
    unless ( $cutoff =~ /^\d+$/ && $cutoff >= 20) {
      die "Error: len_cutoff should be at least 20\n";
    }

    $self->{'R2::PDB::Align::len_cutoff'} = $cutoff;
  }

  unless ($self->{'R2::PDB::Align::len_cutoff'}) {
    $self->{'R2::PDB::Align::len_cutoff'} = $default_cutoff;
  }

  return $self->{'R2::PDB::Align::len_cutoff'};
}



#==========================================================================
#Override accessor for nres_cutoff (number of aligned residues) and set
#defaul value

sub nres_cutoff {

  my ($self, $cutoff) = @_;

  my $default_cutoff = 50;

  if ($cutoff) {

    #validate
    unless ( $cutoff =~ /^\d+$/ && $cutoff >= 20) {
      die "Error: nres_cutoff should be at least 20\n";
    }

    $self->{'R2::PDB::Align::nres_cutoff'} = $cutoff;
  }

  unless ($self->{'R2::PDB::Align::nres_cutoff'}) {
    $self->{'R2::PDB::Align::nres_cutoff'} = $default_cutoff;
  }

  return $self->{'R2::PDB::Align::nres_cutoff'};
}






###########################################################################
##                                                                       ##
##              Functions associated with this class                     ##
##                                                                       ##
###########################################################################


#==========================================================================
#Align two structures, assuming that both structures have only one chain.
#That is, if the strucutures have multiple chains, they should be split into
#chain before aligning.
#
#Function R2::PDB::Parser->split_chains() can be be used to cut a structure
#into its different chains.
#

sub align_structures {

  my $self = shift;


  #Both pdb files must exist
  return unless (-f $self->pdb1->infile && -f $self->pdb2->infile);

#  print Data::Dumper->Dump([ $self->pdb1->infile, $self->pdb2->infile],
#			   [qw(*infile1 *infile2)]);



  #Parse input structures if they have not been parsed yet
  $self->pdb1->parse() unless (%{ $self->pdb1->atom });
  $self->pdb2->parse() unless (%{ $self->pdb2->atom });

#  print Data::Dumper->Dump([$self->pdb1, $self->pdb2],[qw(*pdb1 *pdb2)]);
#  exit;



  #Do not analyze quimeric proteins
  return if ($self->pdb1->chimeric || $self->pdb2->chimeric);


  #Do not analyze if the resolution of the structures is poor
  return if ($self->pdb1->resolution > $self->max_resolution || $self->pdb2->resolution > $self->max_resolution);



  #Extract the chains that will be compared. It is assumed that there is
  #only one chain in both structures
  my ($chain1, $chain2) = undef;
  if (scalar keys %{ $self->pdb1->atom } == 1) {
    $chain1 = ( keys %{ $self->pdb1->atom })[0];
  }
  else {
    print "Structure: ", $self->pdb1->id, " has more than one chain!\n";
    exit;
  }

  if (scalar keys %{ $self->pdb2->atom } == 1) {
    $chain2 = ( keys %{ $self->pdb2->atom })[0];
  }
  else {
    print "Structure: ", $self->pdb2->id, " has more than one chain!\n";
    exit;
  }


  #Sequence length of both structures
  my $len1 = length $self->pdb1->extract_atom_seq_1($chain1);
  my $len2 = length $self->pdb2->extract_atom_seq_1($chain2);

#  print Data::Dumper->Dump([$len1, $len2 ], [qw(*len1 *len2 )]);
#  exit;



  #Do not align structures if the length of any of the structures is
  #smaller than $self->len_cutoff
  return if ($len1 < $self->len_cutoff || $len2 < $self->len_cutoff);



  #PDB IDs for both chains (some times the user preferes to name
  #chains differently than their official PDB IDs)
  my $id1 = $self->id1 ? $self->id1 : $self->pdb1->id . "_" . $chain1;
  my $id2 = $self->id2 ? $self->id2 : $self->pdb2->id . "_" . $chain2;

#  print "Structure Identification:\n";
#  print Data::Dumper->Dump([$id1, $chain1, $len1, $id2, $chain2, $len2],
#			   [qw(*id1 *chain1 *len1 *id2 *chain2 *len2)]);
#  exit;



  #Get the swissprot IDs of the chains in order to extract the coordinates
  #of TMS for the SwissProt sequences
  my $sp1   = $self->pdb1->dbref->{$chain1}->[0]->{xacc};
  my $sp1db = $self->pdb1->dbref->{$chain1}->[0]->{xdb};

  my $sp2   = $self->pdb2->dbref->{$chain2}->[0]->{xacc};
  my $sp2db = $self->pdb2->dbref->{$chain2}->[0]->{xdb};

#  print Data::Dumper->Dump([$sp1, $sp1db, $sp2, $sp2db], [qw(*sp1 *sp1db *sp2 *sp2db )]);
#  exit;


  #There should be a DBREF section in the structure and thus both $sp1 and $sp2 must have values
  return unless ($sp1 && $sp2);




  #There must be a protein sequence for both structures or return.
  #
  #NOTE:
  #   There are some cases where the uniprot ID in the structure does not correspond
  #   to the accession in the family. For example PDB:1E12 (uniprot P0DMH7) from
  #   family 3.E.1... P0DMH7 is not in the list of protein sequnces of the family.
  #
  #   To solve the problem, if the hmmtop resulst are not found for a given proteing id,
  #   the sequence of the identifier will be extracted from uniprot, hmmtop is run
  #   and the results added to hash hmmtop_seqs so other functions can use it.
  my @tms1 = @{ $self->get_tms_coordinates($sp1db, $sp1) };
  my @tms2 = @{ $self->get_tms_coordinates($sp2db, $sp2) };


  #If there are no TMS in sequences abort the alignment
  return unless ((@tms1 && @tms2) && ($tms1[0] > 0) && ($tms2[0] > 0));

#  print Data::Dumper->Dump([\@tms1, \@tms2],[qw(*tms1 *hmmtop2)]);
#  exit;

#  if ($self->pdb1->id eq '4PGS' && $self->pdb2->id eq '4QNC') {
#    print Data::Dumper->Dump([$self->pdb1->id, $sp1, \@tms1, $self->pdb2->id, $sp2, \@tms2 ],
#			     [qw(*pdb1 *sp1 *tms1 *pdb2 *sp2 *tms2 )]);
#    exit;
#  }


  #----------------------------------------------------------------------
  #The chains to be compared should also have TMS or abort the comparison
  #This is because some chains may have only soluble domains of a transporter).


  $self->pdb1->run_hmmtop($chain1, $self->id1);
  $self->pdb2->run_hmmtop($chain2, $self->id2);


  #There must be results to proceed with RMSD calculation, otherwise they are not transporters
  return unless (exists $self->pdb1->hmmtop->{$chain1} && exists $self->pdb2->hmmtop->{$chain2});
  return unless ($self->pdb1->hmmtop->{$chain1}->[0] > 0 && $self->pdb2->hmmtop->{$chain2}->[0]);

#  print Data::Dumper->Dump([ $self->pdb1->hmmtop, $self->pdb2->hmmtop ],[qw(*hmmtop1 *hmmtop2)]);
#  exit;




  #----------------------------------------------------------------------
  #Now we need to determine the regions of sequence alignment between
  #the ATOM sequence and the origimal uniprot sequence for both
  #input structures. This will help us  determine which TMS in the
  #structures correspond to which TMS in the sequence.

  $self->map_atom2prot_residues($chain1, $sp1db, $sp1, $id1, $chain2, $sp2db, $sp2, $id2);

#  print "Before map_atom2prot_residues($sp1, $id1, $sp2, $id2):\n";
#  print Data::Dumper->Dump([$sp1, $id1, $sp2, $id2],[qw(*sp1 *id1 *sp2 *id2)]);
#  print Data::Dumper->Dump([$self->atom2prot_map],[qw(*atom2prot_map)]);
#  <STDIN>;




  #
  #The largest structure will be fixed and the smallest structure will
  #be moving during the alignment.
  #

  #get the file names for the fixed and moving structures
  my (%fixed, %moving) = ();
  my ($aligned_pdb, $rmsd_file) = "";

  if ($len1 >= $len2) {
    %fixed  = ('id' => $id1, 'ch' => $chain1, 'sp' => $sp1, 'len' => $len1, 'pdb' => $self->pdb1->infile,
	       tms => $self->pdb1->hmmtop->{$chain1}, stride => $self->pdb1->stride->{$chain1});
    %moving = ('id' => $id2, 'ch' => $chain2, 'sp' => $sp2, 'len' => $len2, 'pdb' => $self->pdb2->infile,
	       tms => $self->pdb2->hmmtop->{$chain2}, stride => $self->pdb2->stride->{$chain2});
  }
  else {
    %fixed  = ('id' => $id2, 'ch' => $chain2, 'sp' => $sp2, 'len' => $len2, 'pdb' => $self->pdb2->infile,
	       tms => $self->pdb2->hmmtop->{$chain2}, stride => $self->pdb2->stride->{$chain2}),;
    %moving = ('id' => $id1, 'ch' => $chain1, 'sp' => $sp1, 'len' => $len1, 'pdb' => $self->pdb1->infile,
	       tms => $self->pdb1->hmmtop->{$chain1}, stride => $self->pdb1->stride->{$chain1});
  }





  #Estimate the RMSD values based on the algorithms specified by the user
  my $rmsd1 = {};
  my $rmsd2 = {};
  if ($self->algorithm eq "superpose" || $self->algorithm eq "all") {

    $rmsd1 = run_superpose(\%fixed, \%moving, $self->cov_cutoff, $self->nres_cutoff, $self->outdir);

#    print Data::Dumper->Dump([ $rmsd1 ], [qw( *rmsd )]);
#    <STDIN>;


    #There must be a valid RMSD to proceed.
    if (%{ $rmsd1 }) {

      #Make sure the TMS alignment make sense!
      $self->check_tms_alignment($rmsd1);

#      print "After checking TMS alignment:\n";
#      print Data::Dumper->Dump([$rmsd1],[qw(*rmsd1)]);
#      exit;



      if (exists $rmsd1->{tms_align}) {

	#Add here the number of TMS in each protein
	$rmsd1->{fsptms} = $self->hmmtop_seqs->{ $fixed{sp}  }->[0];
	$rmsd1->{msptms} = $self->hmmtop_seqs->{ $moving{sp} }->[0];

#	print Data::Dumper->Dump([$rmsd1->{fsptms}, $rmsd1->{msptms} ], [qw(*fsptms *msptms)]);
#	<STDIN>;


	#At least the number of helices in the repeat unit must be aligned
	unless (scalar keys (%{ $rmsd1->{tms_align} }) >= $self->{rep_unit}) {
	  $rmsd1 = undef;
	}
      }
      else {
	$rmsd1 = undef;
      }
    }
  }



  if($self->algorithm eq "cealign" || $self->algorithm eq "all") {
    $rmsd2 = run_cealign(\%fixed,   \%moving, $self->cov_cutoff, $self->nres_cutoff, $self->outdir);
  }





  #Do not consider this alignment if the minimum number of amino acids
  #is not aligned.
  if ($rmsd1 && %{ $rmsd1 }) {
    $self->rmsd->{'superpose'} = $rmsd1;
  }

  if ($rmsd2 && %{ $rmsd2 }) {
    $self->rmsd->{'cealign'}   = $rmsd2;
  }
}



#==========================================================================
#Checking the alignment.... Note: the analysis of 3HB is hard wired, it
#needs to change to be able to run with any helix bundle size!


sub check_tms_alignment_tmp {

  my ($self, $rmsd) = @_;

#  print Data::Dumper->Dump([$rmsd],[qw(*rmsd)]);
#  exit;



  my $fsp = $rmsd->{fsp};
  my $fch = $rmsd->{fch};

  my $msp = $rmsd->{msp};
  my $mch = $rmsd->{mch};



  #----------------------------------------------------------------------
  #Get the TMS coordintates for the uniprot sequences

  my $fsp_tms = $self->hmmtop_seqs->{$fsp};
  my $msp_tms = $self->hmmtop_seqs->{$msp};

#  print Data::Dumper->Dump([$fsp_tms, $msp_tms], [qw(*fsp_tms *msp_tms)]);
#  exit;

#  print Data::Dumper->Dump([ $rmsd->{aln} ], [qw( *stride_aln )]);
#  exit;




  #**********************************************************************
  #    Determine which helices are alignined between the structures
  #**********************************************************************

  #for both hashes the first key is always the fixed structure and
  #the second key is the moving structure
  my %alignedHelices = ();
  foreach my $coords (@{ $rmsd->{aln} }) {

    #the aligned residues in both structres
    my $fres = $coords->[0];
    my $mres = $coords->[1];



    #Determine in which helices the aligned residues are located. This
    #Is the structure of the returning hash:
    #
    #{nhelix=>helixNum, chres=>aaATOM, helix=>[start, end]};
    #
    my $fHlx = get_helix4res($fres, $rmsd->{fh});
    my $mHlx = get_helix4res($mres, $rmsd->{mh});

    next unless ($fHlx  && $mHlx);



    #Determine in which TMS from the uniprot sequence the aligned residues located.
    #This is the structure of the returning hash:
    #
    #{ntms=>tmsNumber, chres=>aaATOM, spres=>$aaSP, tms=>[start,  end]}
    #
    my $fspTMS = $self->get_prot_tms_for_res($fres, $fsp, $rmsd->{fix});
    my $mspTMS = $self->get_prot_tms_for_res($mres, $msp, $rmsd->{mov});


#    print "----------------------------------------------------------------------\n";
#    print "fPDB: ", $rmsd->{fix}, "   fRes: $fres   mPDB: ", $rmsd->{mov}, "   mRes:  $mres\n",
#      Data::Dumper->Dump([$fsp_tms, $fHlx, $fspTMS, $msp_tms, $mHlx, $mspTMS],
#			 [qw(*fTMS *fHlx *fspTMS *mTMS *mHlx *mspTMS)]), "\n";
#    <STDIN>;




    #
    #Both residues must be part of an alpha-helix and at least one of the residues
    #in the helix must be part of a TMS in the uniprot protein sequence.

    if (($fHlx && $fHlx->{nhelix}) && ($mHlx && $mHlx->{nhelix})) {

      my ($fTMS, $mTMS) = (0, 0);
      $fTMS = $fspTMS->{ntms} if ($fspTMS && exists $fspTMS->{ntms});
      $mTMS = $mspTMS->{ntms} if ($mspTMS && exists $mspTMS->{ntms});

      if ($fTMS && $mTMS) {
	push (@{ $alignedHelices{  $fHlx->{nhelix} }{ $mHlx->{nhelix} } },
	      [$fres, $mres, $fTMS, $mTMS]);
      }
    }
  }

#  print "Alined Helices:\n", Data::Dumper->Dump([\%alignedHelices],[qw(*alignedHelices)]);
#  <STDIN>;

  #For debugging purposes this is to see how helices are aligned before comparing to HMMTOP
#  $rmsd->{tms_align} = \%alignedHelices;
#  return;




  #Verify that the following conditions are met for all pairs of aligned helices:
  # (1) there are at least the minimum number of helices aligned.
  # (2) there are at least rep_unit contiguous helices in the
  #     alignment.
  # (3) the helices intersect contiguous TMS in the corresponding
  #     uniprot protein sequence.
  #

  my $cnt_contiguous_helices = 0;
  my $contiguous_tms = {};
  my $numAlignHelices = (sort {$b <=> $a} keys %alignedHelices)[0];


 F:foreach my $fhelix_num1 (1 .. $numAlignHelices) {

    #The TMS must exist
    next F unless (exists $alignedHelices{ $fhelix_num1 });


    #the scond helix in the fixed structure must be aligned
    my $fhelix_num2 = $fhelix_num1 + 1;
    next F unless (exists $alignedHelices{ $fhelix_num2 });



    #the third helix in the fixed structure must be aligned
    my $fhelix_num3 = $fhelix_num2 + 1;
    next F unless (exists $alignedHelices{ $fhelix_num3 });



    #Get the current TMS in the moving structure. Some times
    #more than one TMS in the moving structure align with
    #one TMS in the fixed structure.
    my @mhelix1 = keys %{ $alignedHelices{ $fhelix_num1 }};
    my @mhelix2 = keys %{ $alignedHelices{ $fhelix_num2 }};
    my @mhelix3 = keys %{ $alignedHelices{ $fhelix_num3 }};



  M1:foreach my $mhelix_num1 (@mhelix1) {

      #At least min_aa_aln_hlx residues should be aligned between both helices
      my $nResH1 = scalar (@{ $alignedHelices{ $fhelix_num1 }{ $mhelix_num1 }});
      next M1 if($nResH1 < $self->min_aa_aln_hlx);


    M2:foreach my $mhelix_num2 (@mhelix2) {

	#At least min_aa_aln_hlx residues should be aligned between both helices
	my $nResH2 = scalar (@{ $alignedHelices{ $fhelix_num2 }{ $mhelix_num2 }});
	next M2 if($nResH2 < $self->min_aa_aln_hlx);


      M3:foreach my $mhelix_num3 (@mhelix3) {

	  #At least min_aa_aln_hlx residues should be aligned between both helices
	  my $nResH3 = scalar (@{ $alignedHelices{ $fhelix_num3 }{ $mhelix_num3 }});
	  next M3 if($nResH2 < $self->min_aa_aln_hlx);



	  #Aligned helices on both structures are contiguous,
	  if (are_helices_sequential($mhelix_num1, $mhelix_num2, $mhelix_num3)) {

	    #get the information on TMS relative to the uniprot sequences.
	    get_tms_for_helices($alignedHelices{ $fhelix_num1 }{ $mhelix_num1 }, $contiguous_tms);
	    get_tms_for_helices($alignedHelices{ $fhelix_num2 }{ $mhelix_num2 }, $contiguous_tms);
	    get_tms_for_helices($alignedHelices{ $fhelix_num3 }{ $mhelix_num3 }, $contiguous_tms);
	    $cnt_contiguous_helices++;

#	    print Data::Dumper->Dump([$contiguous_tms],[qw(*alnTMS)]);
#	    <STDIN>;

#	    last M1;
	  }
	}
      }
    }
  }


  if ($cnt_contiguous_helices) {
    $rmsd->{hlx_align} =  \%alignedHelices;
    $rmsd->{tms_align} = $contiguous_tms;

#    print "check_tms_aligned:  ", $rmsd->{fix}, " vs ", $rmsd->{mov}, "\n";
#    print Data::Dumper->Dump([$cnt_contiguous_helices, \%alignedHelices, $rmsd],
#			     [qw(*cnt_seq_hlx *alignedHelices *rmsd)]);
#    exit;
  }
}






#==========================================================================
#Make sure the alignment makes sense in terms of the TMS aligned. This is
#the version of the function that works with any number of helices
#in a bundle


sub check_tms_alignment {

  my ($self, $rmsd) = @_;

#  print Data::Dumper->Dump([$rmsd],[qw(*rmsd)]);
#  exit;



  my $fsp = $rmsd->{fsp};
  my $fch = $rmsd->{fch};

  my $msp = $rmsd->{msp};
  my $mch = $rmsd->{mch};



  #----------------------------------------------------------------------
  #Get the TMS coordintates for the uniprot sequences

  my $fsp_tms = $self->hmmtop_seqs->{$fsp};
  my $msp_tms = $self->hmmtop_seqs->{$msp};

#  print Data::Dumper->Dump([$fsp_tms, $msp_tms], [qw(*fsp_tms *msp_tms)]);
#  exit;

#  print Data::Dumper->Dump([ $rmsd->{aln} ], [qw( *stride_aln )]);
#  exit;




  #**********************************************************************
  #    Determine which helices are alignined between the structures
  #**********************************************************************

  #for both hashes the first key is always the fixed structure and
  #the second key is the moving structure
  my %alignedHelices = ();
  foreach my $coords (@{ $rmsd->{aln} }) {

    #the aligned residues in both structres
    my $fres = $coords->[0];
    my $mres = $coords->[1];



    #Determine in which helices the aligned residues are located. This
    #Is the structure of the returning hash:
    #
    #{nhelix=>helixNum, chres=>aaATOM, helix=>[start, end]};
    #
    my $fHlx = get_helix4res($fres, $rmsd->{fh});
    my $mHlx = get_helix4res($mres, $rmsd->{mh});

#    print Data::Dumper->Dump([$fHlx, $mHlx ], [qw(*fHlx *mHlx )]);
#    <STDIN>;



    #Both residues must be part of helices
    next unless ($fHlx  && $mHlx);



    #Determine in which TMS from the uniprot sequence the aligned residues are located.
    #This is the structure of the returning hash:
    #
    #{ntms=>tmsNumber, chres=>aaATOM, spres=>$aaSP, tms=>[start,  end]}
    #
    my $fspTMS = $self->get_prot_tms_for_res($fres, $fsp, $rmsd->{fix});
    my $mspTMS = $self->get_prot_tms_for_res($mres, $msp, $rmsd->{mov});


#   print "\n----------------------------------------------------------------------\n";
#   print "fPDB: ", $rmsd->{fix}, "   fRes: $fres   mPDB: ", $rmsd->{mov}, "   mRes:  $mres\n",
#     Data::Dumper->Dump([$fsp_tms, $fHlx, $fspTMS, $msp_tms, $mHlx, $mspTMS],
#			 [qw(*fsp_tms *fHlx *fspTMS *msp_tms *mHlx *mspTMS)]), "\n";
#   <STDIN>;




    #
    #Both residues must be part of an alpha-helix and at least one of the residues
    #in the helix must be part of a TMS in the uniprot protein sequence.
    #
    if (($fHlx && $fHlx->{nhelix}) && ($mHlx && $mHlx->{nhelix})) {

      my ($fTMS, $mTMS) = (0, 0);
      $fTMS = $fspTMS->{ntms} if ($fspTMS && exists $fspTMS->{ntms});
      $mTMS = $mspTMS->{ntms} if ($mspTMS && exists $mspTMS->{ntms});

      if ($fTMS && $mTMS) {
	push (@{ $alignedHelices{  $fHlx->{nhelix} }{ $mHlx->{nhelix} } },
	      [$fres, $mres, $fTMS, $mTMS]);
      }
    }
  }

#  print "Alined Helices:\n", Data::Dumper->Dump([\%alignedHelices],[qw(*alignedHelices)]);
#  <STDIN>;

  #For debugging purposes this is to see how helices are aligned before comparing to HMMTOP
#  $rmsd->{tms_align} = \%alignedHelices;
#  return;



  #Verify that the following conditions are met for all pairs of aligned helices:
  # (1) there are at least the minimum number of helices aligned.
  # (2) there are at least rep_unit contiguous helices in the
  #     alignment.
  # (3) the helices intersect contiguous TMS in the corresponding
  #     uniprot protein sequence.
  #

  my $contiguous_moving_helices = 0;
  my $contiguous_tms = {};
  my $numAlignHelices = (sort {$b <=> $a} keys %{ $rmsd->{fh} })[0];


 F:foreach my $fhelix (1 .. ($numAlignHelices - $self->rep_unit + 1)) {



    #First check that at least  $self->rep_unit helices are aligned
    my @fixed_HlxBundle_aligned = ();
    foreach my $fHelixIdx ($fhelix .. ($fhelix + $self->rep_unit - 1)) {

      if (exists $alignedHelices{ $fHelixIdx }) {
	push (@fixed_HlxBundle_aligned,  $fHelixIdx);
      }
    }
    next F unless (scalar @fixed_HlxBundle_aligned == $self->rep_unit);

#    print Data::Dumper->Dump([\@fixed_HlxBundle_aligned], [qw(*fixed_HlxBundle_aligned )]);
#    <STDIN>;



    #
    #To this point I know that $self->rep_unit contiguous helices
    #of the fixed structure were aligned. I still need to verify
    #that the aligned helices in the moving structure are also
    #contiguous. But take into account that Some times
    #more than one TMS in the moving structure align with
    #one TMS in the fixed structure.
    #


    #Get the helices in the moving structure that were aligned
    #with at least $self->min_aa_aln_hlx residues aligned. This is the
    #structure of the hash that will be generated (%wellAlnHelices):
    #       {fixed_helix => [aligned moving_helices]}
    my %wellAlnHelices = ();
    foreach my $fHlxIdx ($fhelix .. ($fhelix + $self->rep_unit - 1)) {

      my @mHelices = keys %{ $alignedHelices{ $fHlxIdx }};

      foreach my $mHlxIdx (@mHelices) {

	my $nResAln = scalar @{ $alignedHelices{ $fHlxIdx }{ $mHlxIdx }};

#	print "Fixed: $fHlxIdx\t\tMoving:  $mHlxIdx\t\tAlnResidues: $nResAln\n";
#	<STDIN>;

	#keep only moving helices that have at least $self->min_aa_aln_hlx
	#aligned residues with the fixed helix
	if ($nResAln >= $self->min_aa_aln_hlx) {
	  push (@{ $wellAlnHelices{$fHlxIdx} }, $mHlxIdx);
	}
      }
    }


    #all the helices in the fixed structure should be well aligned.
    #That means the repeat unit was well aligned
    next F unless (scalar keys %wellAlnHelices == $self->rep_unit);


#    print Data::Dumper->Dump([\@fixed_HlxBundle_aligned, \%wellAlnHelices ], [qw(*fixed_HlxBundle_aligned *wellAlnHelices)]);
#    <STDIN>;



    #Get the sequence of aligned helices in the moving structure
    my %movHelicesPaths = ();
    my $helix_in_repeat = 1;
    foreach my $fHelixIdx (sort {$a <=> $b} keys %wellAlnHelices) {
      foreach my $mHelixIdx (sort {$a <=> $b} @{ $wellAlnHelices{$fHelixIdx} }) {

	if ($helix_in_repeat == 1) {
	  push (@{ $movHelicesPaths{$helix_in_repeat} }, [{fix=>$fHelixIdx, mov=>$mHelixIdx}]);
	}
	else {

	  foreach my $array (@{ $movHelicesPaths{$helix_in_repeat - 1} }) {

	    my @new_array = @{ $array };
	    push (@new_array, {fix=>$fHelixIdx, mov=>$mHelixIdx});

	    push (@{ $movHelicesPaths{$helix_in_repeat} }, \@new_array);
	  }
	}
      } #loop mHelixIdx

      $helix_in_repeat++;

    } #loop fHelixIdx

#    print Data::Dumper->Dump([\%movHelicesPaths ], [qw( *movHelicesPaths )]);
#    <STDIN>;


    #Now verify that the helices in the moving structure are contiguous and
    #determine wich TMS were aligned
  HlXALN:foreach my $alignment (@{ $movHelicesPaths{$self->rep_unit} }) {

      my @mov_hlx = map { $_->{mov} } @{ $alignment };

      if (are_helices_sequential(\@mov_hlx)) {

	foreach my $alnPair (@{ $alignment }) {
	  get_tms_for_helices($alignedHelices{ $alnPair->{fix} }{ $alnPair->{mov} }, $contiguous_tms);
	}

	$contiguous_moving_helices = 1;
	last HlXALN;
      }
    }
#    print Data::Dumper->Dump([ $contiguous_tms], [qw( *contiguous_tms)]);
#    <STDIN>;

  } # F loop


  if ($contiguous_moving_helices) {
    $rmsd->{hlx_align} =  \%alignedHelices;
    $rmsd->{tms_align} = $contiguous_tms;

#    print "check_tms_aligned:  ", $rmsd->{fix}, " vs ", $rmsd->{mov}, "\n";
#    print Data::Dumper->Dump([\%alignedHelices, $rmsd],
#			     [qw(*alignedHelices *rmsd)]);
#    exit;
  }
}








#==========================================================================
#Get the TMS coordinates of a uniprot ID, if the coordinates are not
#already available in $self->hmmtop_seqs, then the sequence is
#downloaded from uniprot, TMS are predicted with HMMTOP, and finally
#the TMS coordinates are added to $self->hmmtop_seqs.


sub get_tms_coordinates {

  my ($self, $xdb, $xid) = @_;


  #Accession is present in the TCDB family... great, return coordinates!
  return $self->hmmtop_seqs->{$xid} if (exists $self->hmmtop_seqs->{$xid});


  #Make sure the sequence file is there, even if the TMS coordinates are already
  #available in  $self->hmmtop_seqs->{$xid}. This is because other functions need
  #that sequence.
  my $seqfile = $self->outdir . "/${xid}.faa";
  $self->pdb1->download_prot_seq($xdb, $xid) unless (-f $seqfile && !( -z $seqfile ));


  #make sure sequence was downloaded
  die "Could not find file: $seqfile" unless (-f $seqfile && !( -z $seqfile ));





  #Run hmmtop
  my $output = `hmmtop -if=$seqfile -of=-- -sf=FAS -pi=spred -is=pseudo 2> /dev/null`;
  chomp($output);

  #parse HMMTOP output
  my ($header, $sep, $tmsData) = split(/\s+(IN|OUT)\s+/, $output);
  my @tmsHMMTOP = split(/\s+/, $tmsData);


  #Return null values If no TMS are inferred
  if ($tmsHMMTOP[0] == 0) {

    #Add data to main hash before returning
    $self->hmmtop_seqs->{$xid} = [0, []];
    return  $self->hmmtop_seqs->{$xid};
  }


  #get the corrdinates of each TMS
  my @tms = ();
  for (my $i=1; $i <= $#tmsHMMTOP - 1; $i += 2) {

    die "Error: coult not retrieve both TMS cordinates:\n$output\nCurrent Idx=$i" unless ($tmsHMMTOP[$i] && $tmsHMMTOP[$i+1]);
    push (@tms, [$tmsHMMTOP[$i], $tmsHMMTOP[$i+1]]);

  }


  #Add TMS coordinates to main hash
  $self->hmmtop_seqs->{$xid} = [$tmsHMMTOP[0], \@tms];


  #Return TMS coordinates
  return $self->hmmtop_seqs->{$xid};
}








#==========================================================================
#get the corresponding TMS in the uniprot sequnced based on the
#aligned helices


sub get_tms_for_helices {

  my ($data, $out) = @_;

  #
  # $data is an array of arrays with the following structure:
  #              [fres, mres, $ftms, mtms]
  #



  foreach my $data (@{ $data }) {
    $out->{ $data->[2] }->{ $data->[3] } = 1;
  }
}


#==========================================================================
#Determine if the helix numbers are contiguous in the structure

sub are_helices_sequential {

  my $helices = shift;

  my @sorted = sort {$a <=> $b} @{ $helices };

  my $sequencial_helices = 1;
  foreach my $idx (0 .. $#sorted - 1) {

    my $diff = $sorted[$idx + 1] - $sorted[$idx];

    unless ($diff == 1) {
      $sequencial_helices = 0;
      last;
    }
  }

  return $sequencial_helices;
}


#==========================================================================
#Get the corresponding TMS in the protein sequence for the aligned
#residue in the structure


sub get_prot_tms_for_res {

  my ($self, $res, $sp, $pdb_id) = @_;


  #There must be predicted TMS for uniprot sequence
  die "No TMS for uniprot sequence of $pdb_id: $sp" unless (exists $self->hmmtop_seqs->{$sp});


  #get the TMS inferred for this protein
  my $tms = $self->hmmtop_seqs->{$sp};


  #now get the corresponding residue in the protein sequence
  my $spRes = $self->atom2prot_map->{$pdb_id}->{$res};
  return unless ($spRes);

#  print "get_prot_tms_for_res()\n";
#  print "$res in $pdb_id maps to $spRes in $sp\n", Data::Dumper->Dump([$tms],[qw(*tms)]);
#  <STDIN>;


  my $tms_cnt = 1;
  foreach my $pair (@{ $tms->[1] }) {

    my ($start, $end) = @$pair;

#    print "$res in $pdb_id maps to $spRes in $sp\n";
#    print Data::Dumper->Dump([ $pair ], [qw( *current_tms )]);
#    <STDIN>;


    if ($spRes >= $start && $spRes <= $end) {
      return {ntms=>$tms_cnt, chres=>$res, spres=>$spRes, tms=>[$start,  $end]};
    }
    $tms_cnt++;
  }
}



#==========================================================================
#Given and aminoacid number and the coords of helices in the structure
#determine in which helix is the aminoacid located.

sub get_helix4res {

  my ($res, $helices) = @_;

  foreach my $h (sort {$a<=>$b} keys %{ $helices }) {

    my ($aa1, $hStart, $aa2, $hEnd) = @{ $helices->{$h} };

    if ($res >= $hStart && $res <= $hEnd) {
      return {nhelix=>$h, chres=>$res, helix=>[$hStart,  $hEnd]};
    }
  }
}




#==========================================================================
#Run sequence alignments between the ATOM sequence and the protein
#sequence for both structures. This function assumes that the
#fasta sequence for the ATOM regions of both structures already exist.

sub map_atom2prot_residues {

  my ($self, $chain1, $sp1db, $protAcc1, $chain_id1, $chain2,  $sp2db, $protAcc2, $chain_id2) = @_;

  #
  #Note:
  #  The chain1 and chain2 varaibles are used to track the position of the first
  #  residue in the structure. This is because when structures are cut in
  #  helix bundle, the first residue won't be always 1.
  #

#  print "map_atom2prot_residues($protAcc1, $chain_id1, $protAcc2, $chain_id2):\n";
#  print Data::Dumper->Dump([$self->fam_seqs_file1, $self->fam_seqs_file2],[qw(*file1 *file2)]);
#  exit;


  #First get the protein sequences of the two proteins
  $self->get_fasta_sequence($sp1db, $protAcc1, $self->fam_seqs_file1);
  $self->get_fasta_sequence($sp2db, $protAcc2, $self->fam_seqs_file2);


  #Verify that the sequences from the ATOM and uniprot exist;
  my $dir = $self->outdir;
  my $sp1 = "$dir/${protAcc1}.faa";
  my $sp2 = "$dir/${protAcc2}.faa";
  my $ch1 = "$dir/${chain_id1}.faa";
  my $ch2 = "$dir/${chain_id2}.faa";
  die "All sequences must exist:\n$sp1\n$ch1\n$sp2\n$ch2" unless (-f $sp1 && -f $ch1 && -f $sp2 && $ch2);


  #run the alignments between protein and ATOM sequence
  $self->align_atom_vs_prot_seqs($chain1, $protAcc1, $sp1, $chain_id1, $ch1);
  $self->align_atom_vs_prot_seqs($chain2, $protAcc2, $sp2, $chain_id2, $ch2);

  die "No alignment results for: $chain_id1" unless (exists $self->atom2prot_map->{$chain_id1});
  die "No alignment results for: $chain_id2" unless (exists $self->atom2prot_map->{$chain_id2});

}


#==========================================================================
#Align two sequences and extract their region of alignment

sub align_atom_vs_prot_seqs {

  my ($self, $chain, $sp, $spFile, $ch, $chFile) = @_;

  #
  #The output of Needleman & Wush algorithm is two fasta sequences
  #of the same size in one file
  #

  my $outfile = $self->outdir . "/${ch}.needle";
  my $atomSeq = "-asequence $chFile";
  my $protSeq = "-bsequence $spFile";
  my $params  = "-gapopen 10 -gapextend 0.5 -aformat fasta -outfile $outfile";


  #Run here Needleman & Wunsh
  unless (-f $outfile && ! -z $outfile) {
    system "needle $atomSeq $protSeq $params 2> /dev/null";
  }
  return unless (-f $outfile && ! -z $outfile);


  #Now parse the Output of needleman and Wunsh
  my $alnObj  = Bio::SeqIO->new(-file => $outfile , -format => "Fasta");
  my $atomObj = $alnObj->next_seq;
  my $protObj = $alnObj->next_seq;

  unless ($atomObj->length == $protObj->length) {
    die "Seqs in needle output are not of the same length: $outfile";
  }

  #
  #NOTE:
  #  Get the right numbering of the ATOM section. If this structure was
  #  cut, it many not start in the first residue. However, because my
  #  program cut the structure it is guarnateed that all atoms are
  #  contiguous. So, I just need the position of the first atom in
  #  the structure.
  #

  #determin which structure is bein analized.
  my $min_pdb_id = $self->pdb1->id;
  my $pdbObj = ($ch =~ /$min_pdb_id/)? $self->pdb1 : $self->pdb2;


  my $atom_idx = (sort {$a <=> $b} keys %{ $pdbObj->atom->{$chain} })[0];
#  print "Jenny:\nFirst Atom in structure $ch (chain $chain) is $atom_idx\n";
#  <STDIN>;


  my @atomSeq = split(//, $atomObj->seq);

  my $prot_idx = 1;
  my @protSeq = split(//, $protObj->seq);

  foreach my $i (1 .. $atomObj->length) {

    #Get residues in current position
    my $rAtom = $atomSeq[$i-1];
    my $rProt = $protSeq[$i-1];


    #Identify which positions are identical
    if ($rAtom =~ /[a-zA-Z]/ && $rProt =~ /[a-zA-Z]/) {
      $self->atom2prot_map->{$ch}->{$atom_idx} = $prot_idx;
      $atom_idx++;
      $prot_idx++;
    }
    elsif ($rAtom =~ /[a-zA-Z]/) {
      $atom_idx++;
    }
    elsif ($rProt =~ /[a-zA-Z]/) {
      $prot_idx++;
    }
  }


  #For debuggin purposes
#  print "$ch vs $sp\n";
#  foreach my $i (sort {$a <=> $b} keys %{ $self->atom2prot_map->{$ch} }) {
#    print "atom( $i ) --> prot ( ", $self->atom2prot_map->{$ch}->{$i}, " )\n";
#  }
#  <STDIN>;
}



#==========================================================================
#Retrieve the protein sequence from the file containing all the proteins
#in the corresponding  family


sub get_fasta_sequence {

  my ($self, $db, $acc, $seqsFile) = @_;


  #if sequence file already exists, return
  my $outfile = $self->outdir . "/${acc}.faa";
  return if (-f $outfile);


  #Sequence file must exist;
  die "Sequence file can't be open: $seqsFile" unless (-f $seqsFile);


  #open file with the family sequences
  my $famSeqsObj = Bio::SeqIO->new(-file => $seqsFile , -format => "Fasta");


  #Search for the corresponding protein now
  my $found_seq_flag = 0;
 SEQ:while (my $sobj = $famSeqsObj->next_seq) {

    if ($sobj->primary_id =~ /$acc/) {

      my $out = Bio::SeqIO->new(-file => ">$outfile", -format => "Fasta");
      $out->write_seq($sobj);

      $found_seq_flag = 1;
      last SEQ;
    }
  }

  #Sequence was not found in file with family sequences, try to download it from UNIPROT
  unless ($found_seq_flag) {
    $self->pdb1->download_prot_seq($db, $acc);
    die "For structure: Can't find $acc in file $seqsFile" unless (-f $outfile);
  }
}




#==========================================================================
#Run superpose

sub run_superpose {

  my ($fixed, $moving, $cov_cutoff, $nres_cutoff, $outdir) = @_;


  #Validate input structures and output directory
  unless (-f $fixed->{'pdb'}  || !(-z $fixed->{'pdb'})) {
    die "Reference structure not found: $fixed";
  }
  unless (-f $moving->{'pdb'} || !(-z $moving->{'pdb'})) {
    die "Moving structure not found: $moving";
  }
  system "mkdir -p $outdir" unless (-d $outdir);


  #input files
  my $fpdb = $fixed->{'pdb'};
  my $mpdb = $moving->{'pdb'};


  #output files
  my $root_name = $fixed->{'id'} . "_fix_" . $moving->{'id'} . "_mov";
  my $aln_file  = "$outdir/${root_name}.pdb";
  my $rmsd_file = "$outdir/${root_name}.rmsd";


  #run superpose here
  unless (-f $rmsd_file && !( -z $rmsd_file)) {
    system "superpose $mpdb -s -all $fpdb -s -all -o $aln_file > $rmsd_file";
  }


  #Validate output files
  unless (-f $aln_file && -f $rmsd_file) {
    system "touch $outdir/no_rmsd_for_$root_name";
    return {};
  }


  #parse RMSD output
  my ($rmsd, $aln_residues) = undef;
  my @alignedRes = ();

  open (my $superpose, "<", $rmsd_file) || die $!;
  while (<$superpose>) {
    chomp;


    #Trim trailing spaces
    s/^\s+//;
    s/\s+$//;


    #Read RMSD
    if (/r\.m\.s\.d:\s+([0-9\.\+\-]+)/) {
      $rmsd = $1;
    }


    #read number of residues aligned
    if (/Nalign:\s+(\d+)/) {
      $aln_residues = $1;
    }


    #Determine the residues that were aligned
    if (/^\|/) {

      #remove first and last pipes
      s/^\|//;
      s/\|$//;

      my ($q, $d, $t) = split(/\|/, $_);
      my $mov = trim($q);
      my $dis = trim($d);
      my $fix = trim($t);
      next unless ($mov && $dis && $fix);


      #get cordinates for the two structures
      my ($resFix, $resMov) = undef;
      if ($fix =~ /\s+(\d+)$/) {
	$resFix = $1;
      }
      if ($mov =~ /\s+(\d+)$/) {
	$resMov = $1;
      }
      next unless ($resFix && $resMov);

      push (@alignedRes, [$resFix, $resMov]);

    }
  }
  close $superpose;

  unless ($rmsd && $aln_residues) {
    die "Could not read RMSD and Nalign from: $rmsd_file";
  }

  #Get the length of the smaller structure
  my $lfix = $fixed->{'len'};
  my $lmov = $moving->{'len'};
  my $smaller = ($lmov <= $lfix)? $lmov : $lfix;


  #Now calculate the coverage of the structural alignment
  my $localCov = $aln_residues / $smaller;

#  print Data::Dumper->Dump([$lfix, $lmov, $smaller, $localCov ], [qw(*lfix *lmov *smaller *localCov )]);
#  <STDIN>;



  #Ignore RMSD if the coverage of the alignment does not meet
  #the minimum threshold.
  return {} unless ($localCov >= $cov_cutoff);


  #Ignore RMSD if the minimum number of amino acids
  #is not aligned
  return {} unless($aln_residues >= $nres_cutoff);


  #Format output hash
  my %out = (
	       fix  => $fixed->{id},  fch => $fixed->{ch},  fsp => $fixed->{sp},  lfix => $lfix, ftms => $fixed->{tms},  fh => $fixed->{stride},
	       mov  => $moving->{id}, mch => $moving->{ch}, msp => $moving->{sp}, lmov => $lmov, mtms => $moving->{tms}, mh => $moving->{stride},
	       rmsd => $rmsd, nres => $aln_residues,
	       cov  => $localCov, prog => 'superpose',
	       aln  => \@alignedRes
	    );


#  print Data::Dumper->Dump([\%out],[qw(*out)]);
#  <STDIN>;


  return \%out;
}





#==========================================================================
#Align structures with the pymol script cealign

sub run_cealign {

  my ($fixed, $moving, $cov_cutoff, $nres_cutoff, $outdir) = @_;


  #Validate input structures and output directory
  unless (-f $fixed->{'pdb'}  || !(-z $fixed->{'pdb'})) {
    die "Reference structure not found: $fixed";
  }
  unless (-f $moving->{'pdb'} || !(-z $moving->{'pdb'})) {
    die "Moving structure not found: $moving";
  }
  system "mkdir -p $outdir" unless (-d $outdir);


  #input files
  my $fpdb = $fixed->{'pdb'};
  my $mpdb = $moving->{'pdb'};


  #output files
  my $fID = $fixed->{'id'};
  my $mID = $moving->{'id'};


  #run pymol's cealign program
  my $cmd = "pymol -cq $fpdb $mpdb -d 'cealign $fID, $mID'";
  my $output = `$cmd`;


  #parse pymol's output and the the RMSD value
  my  ($rmsd, $aln_residues) = "";
  if ($output =~ /RMSD\s+([0-9\.]+)\s+over\s+(\d+)\s+residues/) {
    $rmsd = $1;
    $aln_residues = $2;
  }

#  print Data::Dumper->Dump([$rmsd, $aln_residues],[qw(*rmsd *aln_residues)]);
#  exit;


  unless ($rmsd && $aln_residues) {
    return {};
  }

  #Get the length of the smaller structure
  my $lfix = $fixed->{'len'};
  my $lmov = $moving->{'len'};
  my $smaller = ($lmov <= $lfix)? $lmov : $lfix;


  #Now calculate the coverage of the structural alignment
  my $localCov = $aln_residues / $smaller;


  #Ignore RMSD if the coverage of the alignment does not meet
  #the minimum threshold.
  return {} unless ($localCov >= $cov_cutoff);

  #Ignore RMSD if the minimum number of amino acids
  #is not aligned
  return {} unless($aln_residues >= $nres_cutoff);


  #Format output hash
  my %out = (
	       fix  => $fixed->{'id'},  lfix => $lfix,
	       mov  => $moving->{'id'}, lmov => $lmov,
	       rmsd => $rmsd, nres => $aln_residues,
	       cov  => $localCov, prog => 'cealign'
	    );

#  print Data::Dumper->Dump([\%out],[qw(*out)]);
#  <STDIN>;


  return \%out;
}



#==========================================================================
#run hmmtop with the ATOM sequence and get the position of the inferred TMSs

sub run_hmmtop {

  my ($self, $chain) = @_;

  #get chain sequence in fasta format
  $self->pdb->save_atom_sequence_in_fasta_format($chain, $self->tmpdir);

  my $infile = $self->tmpdir . "/" . $self->pdb->id . "_${chain}.faa";
  die "Error: fasta file does not exist -> $infile" unless (-f $infile);

  #run hmmtop
  my $output = `hmmtop -if=$infile -of=-- -sf=FAS -pi=spred -is=pseudo 2> /dev/null`;
  chomp($output);


  my @parts = split(/\s+/, $output);

  #If not TMS are inferred return and indicate the error
  return 0 unless ($parts[4] > 0);

  for (my $i=5; $i <= $#parts - 1; $i += 2) {

    die "Error: coult not retrieve both TMS cordinates:\n$output\nCurrent Idx=$i" unless ($parts[$i] && $parts[$i+1]);
    push (@{ $self->hmmtop->{$chain} }, [$parts[$i], $parts[$i+1]]);
  }

  return 1;

}



sub trim {

  my $string = shift;

  return $string unless ($string);

  #triming
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;

  return $string;
}



1;
