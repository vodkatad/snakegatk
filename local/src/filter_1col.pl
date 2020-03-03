#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

$,="\t";
$\="\n";

$SIG{__WARN__} = sub {die @_};

my $usage="$0 [-i] [-v] [-w] [-s SEPARATOR] COL_NUM filter_file < input_file\n
        filter_file must be a single column file
        -i ignore case
        -v invert the filter: in output are printed lines that do not have in COL_NUM an element present in filter_file
        -w COL_NUM may have more than one field separed by ';', input rows are reported if at least one word is present in filter_file
        -s set the separator
        -f each row starting by '>' is always passed to the standard output
";

my $help=0;
my $invert=0;
my $ignore_case=0;
my $multiple_word=0;
my $separator = ';';
my $pass_fasta_headers=0;
GetOptions (
        'h|help' => \$help,
        'invert|v' => \$invert,
        'ignore-case|i' => \$ignore_case,
        'multiple_word|w' => \$multiple_word,
        'separator|s=s' => \$separator,
        'fasta|f' => \$pass_fasta_headers,
) or die($usage);

if($help){
        print $usage;
        exit(0);
}

my $col_1=shift;
die($usage) if !defined($col_1) or $col_1 !~ /\d+/;
$col_1--;
die($usage) if $col_1<0;

my $filter_filename = shift;
open FH, $filter_filename or die("can't open file ($filter_filename)");

my %filter=();
while(<FH>){
        chomp;
        my @F= split /\t/;
        die("Error in ($filter_filename) file; 1! column allowed") if scalar(@F)!=1;
        my $key = $F[0];
        if ($ignore_case) {
                $key = lc($key);
        }
        $filter{$key}=1;
}

while(<>){
        chomp;
        if($pass_fasta_headers and m/^>/){
                print;
                next;
        }
        my @F=split /\t/;
        die("Insufficent column number in STDIN") if !defined($F[$col_1]);
        my $k=$F[$col_1];
        if ($ignore_case) {
                $k = lc($k);
        }
        my $found=$filter{$k};
        if($multiple_word and not defined $found){
                for(split /$separator/, $k){
                        $found = $filter{$_};
                        last if defined $found;
                }
        }

        print if (($invert and !defined($found)) or (!$invert and defined($found)));
}

