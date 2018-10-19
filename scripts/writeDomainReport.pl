#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

use TCDB::Assorted;


#Generte a report of all the domains. Example:
#    writeDomainReport ../../negative_control.txt superfamily_dc50_e01


my $infile = $ARGV[0];
my @sDirs  = split(/,/, $ARGV[1]);


#==========================================================================
#THe main domains in the supperfamily clan

my @mfs = ("PF07690", "PF00083", "PF12832", "PF00854", "PF13000", "PF01306",
	   "PF01733", "PF01770", "PF03092", "PF03137", "PF03209", "PF03219",
	   "PF03239", "PF03825", "PF05631", "PF05977", "PF05978", "PF06609",
	   "PF06779", "PF06813", "PF06963", "PF07672", "PF11700", "PF13347");
my %refDomains = map {$_ => 1} @mfs;


my %domainStats = ();
my %headerDomains = ();



#==========================================================================
#The families for which the report will be built

open (my $fh, "<", $infile) || die $!;
chomp(my @rFams = <$fh>);
close $fh;

#print Data::Dumper->Dump([\@rFams ], [qw(*rFams)]);
#exit;



#==========================================================================
#get the Counts for each family in each input directory. This function
#Works for the families in the Negative control



getCountsForNegativeControlFamilies(\@rFams, \@sDirs, \%domainStats, \%headerDomains);

#print Data::Dumper->Dump([\%headerDomains ], [qw(*headerDomains )]);
#exit;




#==========================================================================
#Get the counts for each of the families in the positive control


#getCountsForPositiveControlFamilies(\@rFams, \@sDirs, \%domainStats, \%headerDomains);

#print Data::Dumper->Dump([ \%domainStats], [qw(*domainStats )]);





#Generate the repport now
createReport(\%domainStats);





###########################################################################
##                               Functions                               ##
###########################################################################



sub getCountsForPositiveControlFamilies {


  my ($fams, $indirs, $outStats, $hDomains) = @_;

  foreach my $dir (@{ $indirs }) {

    foreach my $fam (@{ $fams }) {


      #-----------------------------------------------------------------
      #Read file with rescued domains

      my $infile = "$dir/2.A.1/reports/2.A.1_rescuedDomains.tsv";
      my $famEx = get_regex_for_tcid($fam);
      my $perlGrep = qq(perl -ne 'print if (/$famEx/);' $infile);

      open (my $fileh, "-|", $perlGrep) || die $!;
      chomp(my @matches = <$fileh>);
      close $fileh;

      die "Error: No pfam matches for family $fam!" unless (@matches);

#      print Data::Dumper->Dump([$fam, \@matches ], [qw(*fam *matches)]);
#      exit;



      #-----------------------------------------------------------------
      #Count the total number of rescues


      #First initalize counters
      my ($num, $id, @dom) = split (/\s+/, $matches[0]);
      foreach my $d (@dom) {
	my ($pfamid, @kk) = split (/\|/, $d);

	next unless (exists $refDomains{$pfamid});

	#Collect domain for final report
	$hDomains->{$dir}->{$pfamid} = 1;

	#initialize counters
	$outStats->{$dir}->{$fam}->{prots} = 0;
	$outStats->{$dir}->{$fam}->{$pfamid} = {DirectHit => 0, Rescued1 => 0, Rescued2 => 0, Nohit => 0};
      }

#      print Data::Dumper->Dump([$hDomains ], [qw(*header )]);
#      exit;


      #Now make the counts
      foreach my $sys (@matches) {

	my ($n, $tc, @domains) = split(/\s+/, $sys);

	foreach my $dom (@domains) {

	  my ($pfam, @rest) = split (/\|/, $dom);
	  my $status = $rest[$#rest];

	  #Work only with domains in the reference supperfamily
	  next unless (exists $refDomains{$pfam});

	  #Count proteins with this domain in the family
	  $outStats->{$dir}->{$fam}->{$pfam}->{$status}++
	}

	#Count each protein per system
	$outStats->{$dir}->{$fam}->{prots}++;
      }
    }
  }
}







#==========================================================================
#Create report for the families in the negative control. Each family
#in the negative control was compared to all the proteins in the
#positive control to determine if the domains of the reference
#superfamily (MFS in this case) can be found in the each individual
#family of the negative control.




sub createReport {

  my $data = shift;



 DIR:foreach my $dir (keys %{ $data }) {

    #first sort the reference domains based on the number of hits,
    #This is the order in which they'll appear in the report
    my @header = sortDomainsByPfamHits($dir);

#    print Data::Dumper->Dump([\@header ], [qw(*header )]);
#    exit;

    my $outFile = "$dir/domainReport.tbl";
    open (my $outh, ">", $outFile) || die $!;

  FAM:foreach my $fam (TCDB::Assorted::sort_by_family(keys %{ $data->{$dir} })) {

      my $prots = $data->{$dir}->{$fam}->{prots};
      die "Number of proteins can't be zero! -> Dir: $dir, Fam: $fam" unless ($prots);


      print $outh "$fam\t$prots";
    DOM:foreach my $domain (@header) {

	next DOM if ($domain eq 'prots');

	my $counts = (exists $data->{$dir}->{$fam}->{$domain})? $data->{$dir}->{$fam}->{$domain} : {};

	if (%{ $counts }) {
	  my $dh  = $counts->{DirectHit};
	  my $dhp = sprintf ("%.1f", $dh/$prots * 100);

	  my $r1 = $counts->{Rescued1};
	  my $r1p = sprintf ("%.1f", $r1/$prots * 100);

	  my $r2 = $counts->{Rescued2};
	  my $r2p = sprintf ("%.1f", $r2/$prots * 100);

	  my $nh = $counts->{Nohit};
	  my $nhp = sprintf ("%.1f", $nh/$prots * 100);

	  print $outh qq(\t$domain|$dh(${dhp}%),$r1(${r1p}%),$r2(${r2p}%),$nh(${nhp}%));
	}
	else {
	  print $outh qq(\t$domain|0(0.0%),0(0%),0(0%),0(0%));
	}
      }
      print $outh "\n";
    }
  }
}



sub sortDomainsByPfamHits {

  my $dir = shift;


  my %stats = ();


  #First initialize counters
  foreach my $dom (keys %{ $headerDomains{$dir} }) {
    $stats{$dom} = {DirectHit => 0, Rescued1 => 0, Rescued2 => 0, Nohit => 0};
  }



  #Get global counts of matches per domain
 FAM:foreach my $fam (keys %{ $domainStats{$dir} }) {
  DOM:foreach my $dom (keys %{ $domainStats{$dir}{$fam} }) {

      next DOM if ($dom eq 'prots');

      $stats{$dom}{DirectHit} += $domainStats{$dir}{$fam}{$dom}{DirectHit};
      $stats{$dom}{Rescued1}  += $domainStats{$dir}{$fam}{$dom}{Rescued1};
      $stats{$dom}{Rescued2}  += $domainStats{$dir}{$fam}{$dom}{Rescued2};
      $stats{$dom}{Nohit}     += $domainStats{$dir}{$fam}{$dom}{Nohit};
    }
  }


  my @out = sort {
                   if ($stats{$a}{DirectHit} == $stats{$b}{DirectHit}) {
		     if ($stats{$a}{Rescued1} == $stats{$b}{Rescued1}) {
		       if ($stats{$a}{Rescued2} == $stats{$b}{Rescued2}) {
			 $stats{$a}{Nohit} <=> $stats{$b}{Nohit};
		       }
		       else {
			 $stats{$b}{Rescued2} <=> $stats{$a}{Rescued2};
		       }
		     }
		     else {
		       $stats{$b}{Rescued1} <=> $stats{$a}{Rescued1};
		     }
		   }
		   else {
		     $stats{$b}{DirectHit} <=> $stats{$a}{DirectHit};
		   }
                 } keys %stats;

  return @out;
}





sub getCountsForNegativeControlFamilies {

  my ($fams, $dirs, $outStats, $hDomains) = @_;

  foreach my $fam (@{ $fams }) {

    foreach my $dir (@{ $dirs }) {


      #-----------------------------------------------------------------
      #Read file with rescued domains

      my $infile = "$dir/$fam/2.A.1/reports/2.A.1_rescuedDomains.tsv";
      my $famEx = get_regex_for_tcid($fam);
      my $perlGrep = qq(perl -ne 'print if (/$famEx/);' $infile);

      open (my $fileh, "-|", $perlGrep) || die $!;
      chomp(my @matches = <$fileh>);
      close $fileh;

      die "Error: No pfam matches for family $fam!" unless (@matches);

#      print Data::Dumper->Dump([$fam, \@matches ], [qw(*fam *matches)]);
#      exit;



      #-----------------------------------------------------------------
      #Count the total number of rescues


      #First initalize counters
      my ($num, $id, @dom) = split (/\s+/, $matches[0]);
      foreach my $d (@dom) {
	my ($pfamid, @kk) = split (/\|/, $d);

	next unless (exists $refDomains{$pfamid});

	#Collect domain for final report
	$hDomains->{$dir}->{$pfamid} = 1;

	#initialize counters
	$outStats->{$dir}->{$fam}->{prots} = 0;
	$outStats->{$dir}->{$fam}->{$pfamid} = {DirectHit => 0, Rescued1 => 0, Rescued2 => 0, Nohit => 0};
      }


      #Now make the counts
      foreach my $sys (@matches) {

	my ($n, $tc, @domains) = split(/\s+/, $sys);

	foreach my $dom (@domains) {

	  my ($pfam, @rest) = split (/\|/, $dom);
	  my $status = $rest[$#rest];

	  #Work only with domains in the reference supperfamily
	  next unless (exists $refDomains{$pfam});

	  #Count proteins with this domain in the family
	  $outStats->{$dir}->{$fam}->{$pfam}->{$status}++
	}

	#Count each protein per system
	$outStats->{$dir}->{$fam}->{prots}++;
      }
    }
  }
}



sub get_regex_for_tcid {

  my $tcid = shift;

  die "Function must receive a tcid: $tcid" unless ($tcid);


  my $regex = "";
  my $lenID = scalar split(/\./, $tcid);

  if ($lenID <= 4) {
    $regex = qq(\\s$tcid\\.);
  }
  else {
    $regex = qq(\\s$tcid\\s);
  }

  return $regex;

}
