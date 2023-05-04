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
use strict;	# Forces variables to be declared
use POSIX;

sub computeSeqAAEntropyGrpIVL {
	my($sequence) = @_;

	my $i;
	my @seqArray;
	my $validlength;
	my $validAA;
	my @freqArray;
	my $entropy = -1;

	@seqArray = split(//,uc($sequence));
	$validlength = scalar(@seqArray);
	for ($i=0; $i<18; $i++) {
		$freqArray[$i] = 0;
	}

	if (scalar(@seqArray) > 0) {
		# Count the frequency of occurence for each character type
		for ($i=0; $i<scalar(@seqArray); $i++) {
			if ($seqArray[$i] eq 'I') {	# hydrophobic group (simple TM)
				$freqArray[0] += 1;
			}
			elsif ($seqArray[$i] eq 'L') {
				$freqArray[0] += 1;
			}
			elsif ($seqArray[$i] eq 'V') {
				$freqArray[0] += 1;
			}
			elsif ($seqArray[$i] eq 'C') {	# hydrophobic 
				$freqArray[1] += 1;
			}
			elsif ($seqArray[$i] eq 'A') {	
				$freqArray[2] += 1;
			}
			elsif ($seqArray[$i] eq 'M') {
				$freqArray[3] += 1;
			}
			elsif ($seqArray[$i] eq 'W') {
				$freqArray[4] += 1;
			}
			elsif ($seqArray[$i] eq 'F') {
				$freqArray[5] += 1;
			}
			elsif ($seqArray[$i] eq 'P') {	# Structural residues
				$freqArray[6] += 1;
			}
			elsif ($seqArray[$i] eq 'G') {
				$freqArray[7] += 1;
			}
			elsif ($seqArray[$i] eq 'R') {	# other individual residues
				$freqArray[8] += 1;
			}
			elsif ($seqArray[$i] eq 'N') {
				$freqArray[9] += 1;
			}
			elsif ($seqArray[$i] eq 'D') {
				$freqArray[10] += 1;
			}
			elsif ($seqArray[$i] eq 'E') {
				$freqArray[11] += 1;
			}
			elsif ($seqArray[$i] eq 'Q') {
				$freqArray[12] += 1;
			}
			elsif ($seqArray[$i] eq 'H') {
				$freqArray[13] += 1;
			}
			elsif ($seqArray[$i] eq 'K') {
				$freqArray[14] += 1;
			}
			elsif ($seqArray[$i] eq 'S') {
				$freqArray[15] += 1;
			}
			elsif ($seqArray[$i] eq 'T') {
				$freqArray[16] += 1;
			}
			elsif ($seqArray[$i] eq 'Y') {
				$freqArray[17] += 1;
			}
			else {
				$validlength -= 1;
			}
		}

		# Convert frequency to probability
		for ($i=0; $i<scalar(@freqArray); $i++) {
			if ($validlength > 0) {
				$freqArray[$i] /= $validlength;
			}
			else {
				$freqArray[$i] = 0.0;
			}
#			print " ".$freqArray[$i];
		}

		for ($i=0,$entropy=0.0; $i<scalar(@freqArray); $i++) {
			if ($freqArray[$i] > 0) {
				$entropy += $freqArray[$i]*log($freqArray[$i])/log(2);
			}
		}
		$entropy = $entropy*(-1);

#		print " ".$entropy."\n";
	}
	return $entropy;
}
1;