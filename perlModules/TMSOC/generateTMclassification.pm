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

use TMSOC::getAAHydrophobicity;
use TMSOC::computeSummaryStatistics;

use TMSOC::computeSeqAAEntropyGrpIVL;
use TMSOC::computeSeqAAHydrophobicity;
use TMSOC::computeZscore;

sub generateTMclassification {
	my($TMregions, $FASTAseq) = @_;

	my @results = ();

	my $winSizeX = 12;
	my $winSizeY = 19;

	my $entropy;		# Complexity measure variable based on hydrophobicity & structural grouping
	my $aveEntropy;
	my $stdEntropy;
	my $entropyVector;

	my $octanolScale;	# Hydrophobicity/N-octanol measure
	my $aveOctanol;
	my $stdOctanol;
	my $octanolVector;
	my @octanolstats;

	my $interfaceScale;	# Hydrophobicity/interface measure
	my $aveinterface;
	my $stdinterface;
	my $interfaceVector;
	my @interfacestats;

	my $oct2IntfScale;	# Hydrophobicity/N-octanol-Interface measure
	my $aveOct2Intf;
	my $stdOct2Intf;
	my $oct2IntfVector;
	my @oct2Intfstats;

	my $zscore;		# Result section variables
	my $type;
	my $maskedFASTAseq;

	my $i;
	my $j;
	my $k;
	my $l;
	my $subseq;
	my @segments;
	my @seqAA;
	@segments = split(/\s+/, $TMregions);
	@seqAA = split(//,$FASTAseq);

	my @position;
	my @helixAA;
	my $helixseq;
	for ($i=0; $i<scalar(@segments); $i++) {
		@position = split(/,/, $segments[$i]);
		$helixseq = "";
		for ($j=$position[0]; $j<=$position[1]; $j++) {
			$helixseq = $helixseq.$seqAA[$j];
		}
		@helixAA = split(//, $helixseq);

#		print $segments[$i]." ".$helixseq." ".scalar(@helixAA)."\n";

		# Check if the TM helix sequence is smaller than either of the window sizes
		my $extendLength;
		my $leftpos;
		my $rightpos;
		my $rightResidueHydro;
		my $leftResidueHydro;
		
# 		if (scalar(@helixAA)<$winSizeX || scalar(@helixAA)<$winSizeY) {
# 			if ($winSizeX > $winSizeY) {
# 				$extendLength = $winSizeX-scalar(@helixAA);
# 			}
# 			else {
# 				$extendLength = $winSizeY-scalar(@helixAA);
# 			}
# 			$leftpos = $position[0];
# 			$rightpos = $position[1];
# 			for ($k=0; $k<$extendLength; $k++) {
# 				$leftpos = $leftpos - 1;	# Extend left position
# 				$rightpos = $rightpos + 1;	# Extend right position
# 				if ( $leftpos >= 0 && $rightpos<(scalar(@seqAA)) ) {
# 					$leftResidueHydro = getAAHydrophobicity($seqAA[$leftpos]);
# 					$rightResidueHydro = getAAHydrophobicity($seqAA[$rightpos]);
# 					if ($leftResidueHydro < $rightResidueHydro) {
# 						$rightpos = $rightpos - 1;	# Restore right position to previous
# 					}
# 					elsif ($leftResidueHydro > $rightResidueHydro) {
# 						$leftpos = $leftpos + 1;	# Restore left position to previous
# 					}
# 					else {
# 						# do nothing i.e. extend both left and right position
# 					}
# 				}
# 				elsif ($leftpos < 0 && $rightpos<(scalar(@seqAA)) ) {	# left position out of valid index
# 					$leftpos = $leftpos + 1;	# Restore left position to previous
# 				}
# 				elsif ($leftpos >= 0 && $rightpos>=(scalar(@seqAA)) ) {	# left position out of valid index
# 					$rightpos = $rightpos - 1;	# Restore right position to previous
# 				}
# 				else {	# do nothing
# 				}
# 				if (($rightpos-$leftpos+1) >= (scalar(@helixAA)+$extendLength)) {
# 					last;
# 				}
# 			}

# 			# Create new extended TM-helix sequence
# 			$helixseq = "";
# 			for ($k=$leftpos; $k<=$rightpos; $k++) {
# 				$helixseq = $helixseq.$seqAA[$k];
# 			}
# 			@helixAA = split(//, $helixseq);

# #			print "After extension : ".$segments[$i]." ".$helixseq." ".scalar(@helixAA)."\n";
# 		}

		# Calculate the entropy for each predicted segment
		$entropyVector = "";
		for ($l=0; $l<scalar(@helixAA)-$winSizeX+1; $l++) {
			$subseq = "";
			for ($k=$l; $k<$l+$winSizeX; $k++) {
				$subseq = $subseq.$helixAA[$k];
			}

			# Compute the entropy of the subsequence
			$entropy = computeSeqAAEntropyGrpIVL($subseq);

			if ($l==0) {
				$entropyVector = $entropy;	# Based on hydrophobic IVL group as one, other residues as individual
			}
			else {
				$entropyVector = $entropyVector.",".$entropy;	# Based on hydrophobic & structural residues(G/P) group
			}
		}

		# Calculate the hydrophobicity for each predicted segment
		$octanolVector = "";
		$interfaceVector = "";
		$oct2IntfVector = "";
		for ($l=0; $l<scalar(@helixAA)-$winSizeY+1; $l++) {
			$subseq = "";
			for ($k=$l; $k<$l+$winSizeY; $k++) {
				$subseq = $subseq.$helixAA[$k];
			}

			# Compute the hydrophobicity of the subsequence in terms of octanol, interface, octanol-interface scales
			($octanolScale, $interfaceScale, $oct2IntfScale) = computeSeqAAHydrophobicity($subseq);

			if ($l==0) {
				$octanolVector = $octanolScale;
				$interfaceVector = $interfaceScale;
				$oct2IntfVector = $oct2IntfScale;
			}
			else {
				$octanolVector = $octanolVector.",".$octanolScale;
				$interfaceVector = $interfaceVector.",".$interfaceScale;
				$oct2IntfVector = $oct2IntfVector.",".$oct2IntfScale;
			}
		}

		# Compute the summary statistics of the entropy & hydrophobicity vectors
		($aveEntropy, $stdEntropy) = computeSummaryStatistics($entropyVector);
		($aveOctanol, $stdOctanol) = computeSummaryStatistics($octanolVector);
		($aveinterface, $stdinterface) = computeSummaryStatistics($interfaceVector);
		($aveOct2Intf, $stdOct2Intf) = computeSummaryStatistics($oct2IntfVector);

		# Compute the zscore of the TM segment
		$zscore = computeZscore($aveEntropy, $aveOct2Intf, 2.4, 0.30, -0.64, 2.85, 0, 0);	# against functional-TM (UniProt)
#		$zscore = computeZscore($aveEntropy, $aveOct2Intf, 2.42, 0.295, -0.41, 2.91, 0, 0);	# against functional-TM (SCOP)

		# Determine the type of TM segment whether complex(z>=-3.29), simple(z<=-5.41) or twilight(-5.41<z<-3.29)
		if ($zscore >= -3.29) {
			$type = "complex";
		}
		elsif ($zscore <= -5.41) {
			$type = "simple";
		}
		else {
			$type = "twilight";
		}

		# Format decimal places to 2
		$aveEntropy = sprintf "%.2f", $aveEntropy; $stdEntropy = sprintf "%.2f", $stdEntropy;
		$aveOctanol = sprintf "%.2f", $aveOctanol; $stdOctanol = sprintf "%.2f", $stdOctanol;
		$aveinterface = sprintf "%.2f", $aveinterface; $stdinterface = sprintf "%.2f", $stdinterface;
		$aveOct2Intf = sprintf "%.2f", $aveOct2Intf; $stdOct2Intf = sprintf "%.2f", $stdOct2Intf;
		$zscore = sprintf "%.2f", $zscore;

#		print "(".$aveEntropy.",".$stdEntropy.")";
#		print "(".$aveOctanol.",".$stdOctanol.")";
#		print "(".$aveinterface.",".$stdinterface.")";
#		print "(".$aveOct2Intf.",".$stdOct2Intf.")"."\n";

		push(@results, $helixseq.";".$segments[$i].";".$aveEntropy.";".-1*$aveOct2Intf.";".$zscore.";".$type);
	}

	# Create the masked sequence
	my @maskedseqAA = split(//,lc($FASTAseq));
	my @tmp;
	for ($i=0; $i<scalar(@results); $i++) {
		@tmp = split(/\;/, $results[$i]);
		if ($tmp[5] eq "simple") {
			@position = split(/\,/, $tmp[1]);
			for ($j=$position[0]; $j<=$position[1]; $j++) {
				$maskedseqAA[$j] = 'X';
			}
		}
		elsif ($tmp[5] eq "complex" || $tmp[5] eq "twilight") {
			@position = split(/\,/, $tmp[1]);
			for ($j=$position[0]; $j<=$position[1]; $j++) {
				$maskedseqAA[$j] = uc($maskedseqAA[$j]);
			}
		}
		else {}
	}
	$maskedFASTAseq = join('', @maskedseqAA);
#	print $maskedFASTAseq."\n";

	return (\@results, $maskedFASTAseq);

}
1;
