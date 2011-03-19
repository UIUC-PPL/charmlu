#!/usr/bin/perl

use strict;
use warnings;

my $N = 64;
my $mode = "vn";
my $runnum = "167";
my $BLK = 500;
my $matrix = 110000;
my $thresh = 390;
my $agglom = 50;

`mkdir lu-$N-$BLK-$runnum`;

print "Starting on $N nodes, " . ($mode eq "vn" ? $N * 4 : $N) . " processors\n";
my $dcmf = 16*1024*1024;
`qsub --env "BG_MAXALIGNEXP=1000:DCMF_RECFIFO=$dcmf" -A CharmRTS -q prod-devel --mode $mode -n $N -t 55 -o lu-$N-$BLK-$runnum.out -e lu-$N-$BLK-$runnum.error ./lu.trace $matrix $BLK $thresh $agglom 3 1\
6 16 +gz-trace +traceroot lu-$N-$BLK-$runnum +logsize 10000`;
