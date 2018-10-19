#!/usr/bin/env perl -w

use strict;
use warnings;
use Data::Dumper;

use TCDB::Assorted;


#Based on several directories determine the final report:
#Example:
#writeDomainStats superfamily_dc80_e-5,superfamily_dc80_e001,superfamily_dc80_e01 PF07690 > report_dc80.txt


my @indirs  = split(/,/, $ARGV[0]);                 #input directories
my %domains = map {$_ => 1} split (/,/, $ARGV[1]);  #Domains for which statistis will be computed
my $outdir = ".";


my %stats = ();
my %famSize = ();


#Get the counts for each domain of interest in each file
foreach my $dir (@indirs) {

  getStatsForFile($dir, \%domains, \%stats, \%famSize);

#  print Data::Dumper->Dump([\%stats ], [qw(*stats )]);
#  exit;
}



#Generate the report


print "Family\tSize\t\t", join ("\t\t", @indirs), "\n";

foreach my $fam (sort_by_family(keys %{ $stats{ $indirs[0] }})) {

  my @results = ();

  print $fam, "\t", $famSize{$fam};

  foreach my $dir (@indirs) {

    foreach my $dom (sort keys %domains) {

      my $dstats = $stats{$dir}{$fam}{$dom};

#      print Data::Dumper->Dump([$dstats ], [qw(*dstats )]);
#      <STDIN>;

      my $dirHitsN = $dstats->{dirHit}->{nhits};
      my $dirHitsP = $dstats->{dirHit}->{perc};

      my $oneN     = $dstats->{'1st'}->{nhits};
      my $oneP     = $dstats->{'1st'}->{perc};

      my $twoN     = $dstats->{'2nd'}->{nhits};
      my $twoP     = $dstats->{'2nd'}->{perc};

      my $missN    = $dstats->{noHit}->{nhits};
      my $missP    = $dstats->{noHit}->{perc};

      print "\t\t$dirHitsN ($dirHitsP\%)|$oneN ($oneP\%)|$twoN ($twoP\%)|$missN ($missP\%)";

    }
  }

  print "\n";
}







sub getStatsForFile {

  my ($indir, $rDomains, $outStats, $fMembers) = @_;


  my $infile = "$indir/domainReport.tbl";
  die "Report file not found or empty: $infile" unless (-f $infile && !(-z $infile));

  my %tags = (0=>'dirHit', 1=>'1st', 2=>'2nd', 3=>'noHit');

  open (my $fileh, "<", $infile);
  while(<$fileh>) {
    chomp;
    my ($fam, $fSize, @doms) = split(/\t/);
    $fMembers->{$fam} = $fSize;

  D:foreach my $dom (@doms) {

      my ($pfam, $counts_str) = split(/\|/, $dom);


      #Only work with the user specified domains.
      next D unless (exists $rDomains->{$pfam});


      #Now process the counts in the domain of interest
      my @counts = split(/,/, $counts_str);
      my $cnt = 0;
      foreach my $c (@counts) {

	my ($freq, $perc) = undef;
	if ($c =~ /^(\d+)\(([0-9\.]+)\%\)$/) {
	  $freq = $1;
	  $perc = $2;
	}
	else {
	  die "Could not extract counts: $c";
	}

	$outStats->{$indir}->{$fam}->{$pfam}->{ $tags{$cnt} }->{nhits} = $freq;
	$outStats->{$indir}->{$fam}->{$pfam}->{ $tags{$cnt} }->{perc} = $perc;
	$cnt++;
      }

#      print Data::Dumper->Dump([$outStats ], [qw(*outStats )]);
#      <STDIN>;

    }

  }
  close $fileh;


}
