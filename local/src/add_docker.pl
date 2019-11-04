#!/usr/bin/perl

use strict;
use warnings;
my $docker = "FIXME";
if (scalar(@ARGV) == 1) {
    $docker = $ARGV[0];
}
my $next = "";
while (<STDIN>) {
    if ($next ne "") {
        print $next;
        $next = "";
    }
    if ($_ =~ /^(\s+)singularity:/) {
        $next = $1 . 'docker: ' . $docker . "\n";
    } 
    print $_;
}
# singularity cannot be last line so no last loop issue
