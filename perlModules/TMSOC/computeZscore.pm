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

sub computeZscore {
	my($entropy, $hydroscale, $entropycenter, $entropyspread, $hydroscalecenter, $hydroscalespread, $correlation, $correlationflag) = @_;

	my $term1;
	my $term2;
	my $term3;
	my $zscore;

	$term1 = (($entropy-$entropycenter)*($entropy-$entropycenter))/($entropyspread*$entropyspread);
	$term2 = (2*$correlation*($entropy-$entropycenter)*($hydroscale-$hydroscalecenter))/($entropyspread*$hydroscalespread);
	$term3 = (($hydroscale-$hydroscalecenter)*($hydroscale-$hydroscalecenter))/($hydroscalespread*$hydroscalespread);
	
	if ($correlationflag == 1) {
		$zscore = $term1 - $term2 + $term3;
	} else { $zscore = $term1 + $term3; }
#	print "\n".$term1." ".$term3." ".$zscore."\n";
#	print $zscore."\n";

	my $p = 0.436;	# Slope of the regressed line (note that high hydrophobicity is -ve)

#	if ($entropy<$entropycenter && $hydroscale<$hydroscalecenter) {	
	if ($p*$entropyspread*($hydroscale-$hydroscalecenter) < -1*$hydroscalespread*($entropy-$entropycenter)) {
#		print -1*$zscore."\n";
		return -1*$zscore;
	}
	else {
#		print $zscore."\n";
		return $zscore;
	}
}
1;
