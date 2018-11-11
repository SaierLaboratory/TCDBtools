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

sub computeSummaryStatistics {
	my($vector) = @_;	# Note that the vector must be comma-separated

	my @numarray;
	@numarray = split(/,/,$vector);


	my $i;
	my $totalsum;
	my $squaresum;
	my $count;
	my $center;
	my $spread;

	my @sortednumarray;
	my @absdeviation;
	my $median;
	my $MAD;
	if (scalar(@numarray) >= 1) {
		for ($i=0,$totalsum=0,$squaresum=0,$count=0; $i<scalar(@numarray); $i++) {
			$totalsum = $totalsum + $numarray[$i];
			$squaresum = $squaresum + ($numarray[$i]*$numarray[$i]);
			$count = $count + 1;
		}
		if ($count > 0) {
			$center = $totalsum/$count;
			$spread = ($squaresum/$count)-($center*$center);
			if ($spread < 0) {		# To handle rounding errors
				$spread = 0.0;
			}
			else {
				$spread = sqrt($spread);
			}
		}
		else {
			$center = -1;
			$spread = -1;
		}

		# Calculates the median of the numbers
		@sortednumarray = sort { $a <=> $b } @numarray;	# numerical sort
		if ( scalar(@sortednumarray)%2 ) {	# odd length
			$median = $sortednumarray[scalar(@sortednumarray)/2];
		}
		else {	# even length
			$median = 0.5*($sortednumarray[floor(scalar(@sortednumarray)/2)-1] + $sortednumarray[floor(scalar(@sortednumarray)/2)]);
		}

		# Calculates the median absolute deviation of the numbers
		@absdeviation = ();
		for ($i=0; $i<scalar(@numarray); $i++) {
			push(@absdeviation, abs($numarray[$i]-$median));
		}
		@sortednumarray = sort { $a <=> $b } @absdeviation;	# numerical sort
		if ( scalar(@sortednumarray)%2 ) {	# odd length
			$MAD = $sortednumarray[scalar(@sortednumarray)/2];
		}
		else {	# even length
			$MAD = 0.5*($sortednumarray[floor(scalar(@sortednumarray)/2)-1] + $sortednumarray[floor(scalar(@sortednumarray)/2)]);
		}
		if ($MAD < 0) {	# To handle rounding errors
			$MAD = 0.0;
		}

#		$center = $median;
#		$spread = $MAD;
#		print "Median/MAD: ".$median.",".$MAD."\n";
	}
	else {
		$center = -1;
		$spread = -1;
	}
	return ($center, $spread);
}
1;