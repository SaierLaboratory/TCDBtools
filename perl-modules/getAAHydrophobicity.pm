#---------------------------------------------------------------------------
# Singapore, 21st October 2011

# ACADEMIC SOFTWARE LICENSE for TMSOC
# ***********************************

# Copyright:
# Wing-Cheong Wong, Sebastian Maurer-Stroh, Georg Schneider, Frank Eisenhaber
# Bioinformatics Institute (BII) A*STAR Singapore

# This software is subjected to license restrictions. The software is
# provided as is. In its present form and without written consent of 
# the copyright holder, the software is provided in accordance with GPL
# (www.gnu.org/licenses/gpl.html). Among other issues, this excludes any
# warranty and indemnity when using this software. For a commercial license,
# please approach the authors (via email to wongwc@bii.a-star.edu.sg).

# When publishing results with this software, please refer to:
# Wing-Cheong Wong, Sebastian Maurer-Stroh, Frank Eisenhaber, 2011, 
# "Not all transmembrane helices are born equal: Towards the extension of 
# the sequence homology concept to the membrane proteins", Biology Direct

# When observing bugs or strange behavior of the software, the user is 
# encouraged to contact the authors.

#---------------------------------------------------------------------------
# Hydrophobicity scales based on Stephen H. White laboratory at UC Irvine
#---------------------------------------------------------------------------
use strict;	# Forces variables to be declared
use POSIX;

sub getAAHydrophobicity {
	my($residue) = @_;
	
	my $hydrophobicityValue;
	my $octanolInterface;
	
	if ($residue eq 'A') {
		$octanolInterface = 0.33;
	}
	elsif ($residue eq 'R') {
		$octanolInterface = 1.00;
	}
	elsif ($residue eq 'N') {
		$octanolInterface = 0.43;
	}
	elsif ($residue eq 'D') {
		$octanolInterface = 2.41;
	}
	elsif ($residue eq 'C') {
		$octanolInterface = 0.22;
	}
	elsif ($residue eq 'E') {
		$octanolInterface = 1.61;
	}
	elsif ($residue eq 'Q') {
		$octanolInterface = 0.19;
	}
	elsif ($residue eq 'G') {
		$octanolInterface = 1.14;
	}
	elsif ($residue eq 'H') {
		$octanolInterface = 1.37;
	}
	elsif ($residue eq 'I') {
		$octanolInterface = -0.81;
	}
	elsif ($residue eq 'L') {
		$octanolInterface = -0.69;
	}
	elsif ($residue eq 'K') {
		$octanolInterface = 1.81;
	}
	elsif ($residue eq 'M') {
		$octanolInterface = -0.44;
	}
	elsif ($residue eq 'F') {
		$octanolInterface = -0.58;
	}
	elsif ($residue eq 'P') {
		$octanolInterface = -0.31;
	}
	elsif ($residue eq 'S') {
		$octanolInterface = 0.33;
	}
	elsif ($residue eq 'T') {
		$octanolInterface = 0.11;
	}
	elsif ($residue eq 'W') {
		$octanolInterface = -0.24;
	}
	elsif ($residue eq 'Y') {
		$octanolInterface = 0.23;
	}
	elsif ($residue eq 'V') {
		$octanolInterface = -0.53;
	}
	else {
		$octanolInterface = 'NaN';
	}
	
	$hydrophobicityValue = $octanolInterface;	
	return $hydrophobicityValue;
}
1;
