#!/usr/bin/perl
use strict;
use warnings;
#use PerlIO::gzip;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
##
## Usage ./sampleFastq.pl <fastq r1> <fastq r2> <outFastq r1> <outFastq r2> <prob of keeping reads>
## modified from https://www.biostars.org/p/110107/
#
if (scalar(@ARGV) != 5) {
	die "Usage: ./sampleFastq.pl <fastq r1> <fastq r2> <outFastq r1> <outFastq r2> <prob of keeping reads>";
}

srand(42);
my $fqf = new IO::Uncompress::Gunzip $ARGV[0], MultiStream => 1 or die "could not read $ARGV[0]: $!";
my $fqr = new IO::Uncompress::Gunzip $ARGV[1], MultiStream => 1 or die "could not read $ARGV[1]: $!";
my $fqoutf = new IO::Compress::Gzip $ARGV[2] or die "could not write $ARGV[2]: $!";
my $fqoutr = new IO::Compress::Gzip $ARGV[3] or die "could not write $ARGV[3]: $!";
my $proba = $ARGV[4];

print STDERR "Keeping $proba reads\n";

my $nbLines = 1;
my $fqRecordf = '';
my $fqRecordr = '';
while (my $line1=<$fqf> ){
    #print STDERR "$nbLines $line1\n";
    my $line2 = <$fqr>;
    $fqRecordf .= $line1;
    $fqRecordr .= $line2;
    if ($nbLines % 4 == 0) {
        my $random = rand(1);
        if ($random <= $proba) {
            print $fqoutf $fqRecordf;
            print $fqoutr $fqRecordr;
        }
        $fqRecordf = '';
        $fqRecordr = '';
    }
    $nbLines++;
}
close $fqf;
close $fqr;
close $fqoutf;
close $fqoutr;
