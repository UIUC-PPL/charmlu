#!/usr/bin/perl

use strict;
use warnings;

my @files = @ARGV;

print "Run number, job id, total time, percent peak\n";

foreach my $file (@files) {
    open FILE, "<$file";
    my $runnum = -1;
    my $jobid = -1;
    my $time = -1;
    if ($file =~ /[A-Za-z\-]*([0-9]+)\.o([0-9]+).*/) {
        $runnum = $1;
        $jobid = $2;
    }
    for (<FILE>) {
        if (/Time\(s\):\s+([0-9.]+)/) {
            $time = $1;
        }
        # Assumption that time always comes before peak: should hold because
        # main chare prints this in order
        if (/Jaguar/) {
            my @lst = split ' ', $_;
            my $percent = substr $lst[9], 0, (length $lst[9]) - 1;
            print "$runnum, $jobid, $time, $percent\n";
        }
    }
    close FILE;
}
