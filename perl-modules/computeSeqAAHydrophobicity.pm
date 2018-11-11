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

sub computeSeqAAHydrophobicity {
	my($sequence) = @_;

	my $i;
	my @seqArray;
	my $validlength;
#	my $validAA;

	my $interfaceSum;
	my $octanolSum;
	my $octanolInterfaceSum;
	my $interfaceScale;
	my $octanolScale;
	my $octanolInterfaceScale;

	@seqArray = split(//,uc($sequence));
	$validlength = scalar(@seqArray);

	if (scalar(@seqArray) > 0) {
		# Sum the interface and octanol scales for the sequence
		$interfaceSum = 0;
		$octanolSum = 0;
		$octanolInterfaceSum = 0;
		for ($i=0; $i<scalar(@seqArray); $i++) {
			if ($seqArray[$i] eq 'A') {
				$interfaceSum += 0.17;
				$octanolSum += 0.51;
				$octanolInterfaceSum += 0.33;
			}
			elsif ($seqArray[$i] eq 'R') {
				$interfaceSum += 0.81;
				$octanolSum += 1.81;
				$octanolInterfaceSum += 1.00;
			}
			elsif ($seqArray[$i] eq 'N') {
				$interfaceSum += 0.42;
				$octanolSum += 0.85;
				$octanolInterfaceSum += 0.43;
			}
			elsif ($seqArray[$i] eq 'D') {
				$interfaceSum += 1.23;
				$octanolSum += 3.64;
				$octanolInterfaceSum += 2.41;
			}
			elsif ($seqArray[$i] eq 'C') {
				$interfaceSum += -0.24;
				$octanolSum += -0.02;
				$octanolInterfaceSum += 0.22;
			}
			elsif ($seqArray[$i] eq 'E') {
				$interfaceSum += 2.02;
				$octanolSum += 3.63;
				$octanolInterfaceSum += 1.61;
			}
			elsif ($seqArray[$i] eq 'Q') {
				$interfaceSum += 0.58;
				$octanolSum += 0.77;
				$octanolInterfaceSum += 0.19;
			}
			elsif ($seqArray[$i] eq 'G') {
				$interfaceSum += 0.01;
				$octanolSum += 1.15;
				$octanolInterfaceSum += 1.14;
			}
			elsif ($seqArray[$i] eq 'H') {
				$interfaceSum += 0.96;
				$octanolSum += 2.33;
				$octanolInterfaceSum += 1.37;
			}
			elsif ($seqArray[$i] eq 'I') {
				$interfaceSum += -0.31;
				$octanolSum += -1.12;
				$octanolInterfaceSum += -0.81;
			}
			elsif ($seqArray[$i] eq 'L') {
				$interfaceSum += -0.56;
				$octanolSum += -1.25;
				$octanolInterfaceSum += -0.69;
			}
			elsif ($seqArray[$i] eq 'K') {
				$interfaceSum += 0.99;
				$octanolSum += 2.80;
				$octanolInterfaceSum += 1.81;
			}
			elsif ($seqArray[$i] eq 'M') {
				$interfaceSum += -0.23;
				$octanolSum += -0.67;
				$octanolInterfaceSum += -0.44;
			}
			elsif ($seqArray[$i] eq 'F') {
				$interfaceSum += -1.13;
				$octanolSum += -1.71;
				$octanolInterfaceSum += -0.58;
			}
			elsif ($seqArray[$i] eq 'P') {
				$interfaceSum += 0.45;
				$octanolSum += 0.14;
				$octanolInterfaceSum += -0.31;
			}
			elsif ($seqArray[$i] eq 'S') {
				$interfaceSum += 0.13;
				$octanolSum += 0.46;
				$octanolInterfaceSum += 0.33;
			}
			elsif ($seqArray[$i] eq 'T') {
				$interfaceSum += 0.14;
				$octanolSum += 0.25;
				$octanolInterfaceSum += 0.11;
			}
			elsif ($seqArray[$i] eq 'W') {
				$interfaceSum += -1.85;
				$octanolSum += -2.09;
				$octanolInterfaceSum += -0.24;
			}
			elsif ($seqArray[$i] eq 'Y') {
				$interfaceSum += -0.94;
				$octanolSum += -0.71;
				$octanolInterfaceSum += 0.23;
			}
			elsif ($seqArray[$i] eq 'V') {
				$interfaceSum += 0.07;
				$octanolSum += -0.46;
				$octanolInterfaceSum += -0.53;
			}
			else {
				$validlength -= 1;
			}
		}

		$interfaceScale = 0.0;
		$octanolScale = 0.0;
		$octanolInterfaceScale = 0.0;
		if ($validlength > 0) {
			$interfaceScale = $interfaceSum;
			$octanolScale = $octanolSum;
			$octanolInterfaceScale = $octanolInterfaceSum;
		}
	}
	return ($octanolScale, $interfaceScale, $octanolInterfaceScale);
}
1;