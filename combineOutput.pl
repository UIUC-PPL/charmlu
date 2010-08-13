#!/usr/bin/perl

use strict;
use warnings;

my ($BLKSIZE, $numBlks) = @ARGV;

my @data;

for (my $i = 0; $i < $numBlks; $i++) {
    for (my $j = 0; $j < $numBlks; $j++) {
        open FILE, "<", "input-generated-LU-$i-$j";
        my $offi = 0;
        for (<FILE>) {
            my @nums = split(/ /, $_);
            my $offj = 0;
            foreach (@nums) {
                my $index = ($i + $offi) * $BLKSIZE + $j + $offj;
                $data[$index] = $_;
                $offj++;
            }
            $offi++;
            #print "@nums";
        }
        close FILE;
    }
}

print "@data";


