#!/usr/bin/perl

use strict;
use warnings;

if ($#ARGV < 5) {
    print "usage: [iterations] [vary_size] @[matrix size] @[block size] @[min] @[max]\n";
    exit 1;
} 

my ($iterations, $vary_size, $matrixsize, $blocksize, $min, $max) = @ARGV;

print "matrixsize, blocksize, min, max, timings, percentGPU, GPUsize, totalMsgs\n";

for (my $i = 0; $i < $iterations; $i++) {

#`make clean && make OPTS1=\'-DMIN_SIZE=$min -DMAX_SIZE=$max\' -j4 2>/dev/null`;
    
    my @timings;
    my @percentGPU;
    my @GPUsize;
    my @totalMsgs;

    my $trials = 5;

    for (my $j = 0; $j < $trials; $j++) {
	`charmrun +p2 ./lu-proj $matrixsize 1000 1 $blocksize 2>/dev/null > file1`;

	open FILE, "<", "file1";

	for (<FILE>) {
	    if (/total wall time ([\d\.]*)/) {
		push @timings, $1;
	    }

	    if (/percent GPU work = (.*)/) {
		push @percentGPU, $1;
	    }

	    if (/average GPU msg size = (.*)/) {
		push @GPUsize, $1;
	    }

	    if (/message instances = (.*)/) {
		push @totalMsgs, $1;
	    }
	}

	`rm file1`;
    }

    print "$matrixsize, $blocksize, $min, $max, ";

    for (@timings) {
	print "$_, ";
    }

    for (@percentGPU) {
	print "$_, ";
    }

    for (@GPUsize) {
	print "$_, ";
    }
    
    for (@totalMsgs) {
	print "$_, ";
    }

    print "end\n";

    #do the vary thing
    $blocksize *= 2;
    #$n += $vary_size;
}
