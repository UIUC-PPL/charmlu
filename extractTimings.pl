#!/usr/bin/perl

use strict;
use warnings;

my %blockTimingStart;
my %blockTimingStop;

for (<>) {
    if (/^Block (\d+) queueing local LU at internalStep \d+, start time = \d+\.\d+, time = (\d+\.\d+)$/) {
        $blockTimingStart{$1} = $2;
    }
    if (/^Block (\d+) finished local LU at internalStep \d+, time = (\d+\.\d+)$/) {
        $blockTimingStop{$1} = $2;
    }
}

for my $key (sort {$a <=> $b} (keys %blockTimingStart)) {
    if (exists $blockTimingStop{$key} and exists $blockTimingStart{$key}) {
        my $timeStop = $blockTimingStop{$key};
        my $timeStart = $blockTimingStart{$key};
        my $timeStopPrev = $blockTimingStop{$key-1};
        my $delay = 0;

        my $diff = $timeStop - $timeStart;
        if (exists $blockTimingStop{$key-1}) {
            $delay = $timeStart - $timeStopPrev;
        }
        if ($diff >= 0 and $delay >= 0) {
            print "$key, $diff, $delay\n";
        }
    }
}
