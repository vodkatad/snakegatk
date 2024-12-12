#!/usr/bin/env perl
use warnings;
use strict;
if (scalar(@ARGV) != 2) {
    die "Usage $0 SNV|indel|both nomulti|multi < headerlessvcf"
}
my $what = $ARGV[0];
my $multi = $ARGV[1];
my $nmulti = 0;
my $header = <STDIN>;
chomp $header;
my @line = split("\t", $header);
print "ID\t" . $line[9]; # assume always at least one sample
for (my $i = 10; $i < scalar(@line); $i++) {
	print "\t" . $line[$i];
}
print "\n";

while (<STDIN>) {
    chomp;
    my @line = split("\t", $_);
    my @l = split(':', $line[8]);
    #die "Wrong vcf FORMAT" if ($line[8] ne 'GT:AD:AF:DP:F1R2:F2R1:SB' && $line[8] ne 'GT:AD:DP:GQ:PL' && $line[8] ne 'GT:AD:AF:F1R2:F2R1:DP:SB:MB');
    #   TODO inefficient split only here...
    die "Wrong vcf FORMAT" if ($l[4] ne 'NR' || $l[5] ne 'NV');
    if ($line[4] =~ /,/) {
        if ($multi eq 'multi') {
            die "Sorry still to be implemented"; # probably will need to use a library for this
        }
        $nmulti++;
        next;     
    } else {
	my $first = 1;
	for (my $i = 9; $i < scalar(@line); $i++) {
	       &manage_entry($line[2], $line[3], $line[4], $line[$i], $what, $line[0], $line[1], $first);
	       $first = 0;
	}
	print "\n";
    }
}

sub manage_entry {
    my $id = shift;
    my $ref = shift;
    my $alt = shift;
    my $g = shift;
    my $what = shift;
    my $chr = shift;
    my $b = shift;
    my $first = shift;
    $b = $b-1; #switch to zero based
    my $e = $b + length($ref); # end escluded considering length of ref, TODO FIXME for long indels
    if ($what eq 'SNV') {
        return if (length($ref) != 1 || length($alt) != 1);
    } elsif ($what eq 'indel') { 
        return if (length($ref) == 1 && length($alt) == 1);
    } # we do not check for both
    my @afs = split(':',$g);
    my $af = 0;
    if ($afs[4] != 0) {
        $af = $afs[5]/$afs[4];
    }
    if ($first) {
        print $id . "\t" . $af;
    } else {
        print "\t" . $af;
    }
}


print STDERR "multiallelic\t$nmulti";
