#!/usr/bin/env perl

my $infile = $ARGV[0]
    or die "\tI need a fasta file and an output file\n"
    ."$0 fastafile.faa fastafile.clm\n\n";
my $outfile = $ARGV[1]
    or die "\tI need an output file\n$0 fastafile.faa fastafile.clm\n\n";

if( -f "$infile" ) {
    if( my($refSeq,$refNCBI) = readFasta("$infile") ) {
        my $printed = printClms($outfile,$refSeq,$refNCBI);
        if( $printed > 1 ) {
            print "\t$printed sequences printed\n\n";
        }
        elsif( $printed == 1 ) {
            print "\t$printed sequence printed\n\n";
        }
        else {
            print "\tno sequences printed\n\n";
        }
    }
    else {
        die "\tapparently no sequences in $infile\n\n";
    }
}
else {
    die "\t$infile is not here\n\n";
}

sub readFasta {
    my $infile = $_[0];
    if( -f "$infile" ) {
        my %seq = ();
        my %ncbi = ();
        my $id = "";
        open( my $INFAA,"<","$infile" );
        while(<$INFAA>) {
            if( m{^>(\S+)} ) {
                $id = $1;
                my @items = split;
                if( $items[1] =~ m{gi\|\d+} ) {
                    $ncbi{"$id"} = $items[1];
                }
                $seq{"$id"} = "";
            }
            else {
                chomp;
                $seq{"$id"} .= $_;
            }
        }
        close($INFAA);
        my @ids = sort keys %seq;
        my $count = @ids;
        if( $count > 0 ) {
            return(\%seq,\%ncbi);
        }
        else {
            return();
        }
    }
    else {
        return();
    }
}

sub printClms {
    my($outfile,$refSeqs,$refNCBI) = @_;
    my @ids   = sort keys %$refSeqs;
    my @wncbi =  sort keys %$refNCBI;
    my $countwNCBI = @wncbi;
    open( my $CLMN,">","$outfile" );
    if( $countwNCBI > 0 ) {
        print {$CLMN} join("\t","Id","NCBI","Seq"),"\n";
    }
    else {
        print {$CLMN} join("\t","Id","Seq"),"\n";
    }
    my $printed = 0;
    for my $id ( @ids ) {
        my $seq = $refSeqs->{$id};
        if( $countwNCBI > 0 ) {
            my $ncbi = exists $refNCBI->{$id} ? $refNCBI->{$id} : "NA";
            print {$CLMN} join("\t",$id,$ncbi,$seq),"\n";
        }
        else {
            print {$CLMN} join("\t",$id,$seq),"\n";
        }
        $printed++;
    }
    close($CLMN);
    return($printed);
}

sub bringTCDB {
    system("wget -N http://www.tcdb.org/public/tcdb -O $results_dir/tcdb");
    open( my $TCDBI,"<","$results_dir/tcdb" );
    open( my $TCDBO,">","$results_dir/tcdb.faa" );
    while(<$TCDBI>) {
        if( m{^>} ) {
            s{>gnl\|TC-DB\|(\S+?)\|(\S+)}{>lcl\|$2-$1};
            print {$TCDBO} $_;
        }
        else {
            print {$TCDBO} $_;
        }
    }
    close($TCDBI);
    close($TCDBO);
}
