#!/usr/bin/perl

use strict;
use warnings;

my $N =  64;
my $mode = "vn";
my $runnum = "285";
my $BLK = 300;
my $matrix = 75000;
my $thresh = 200;
my $agglom = 100;
my $mapargs = "3 64 4";
my $peRotate = 8;


print "Starting job $runnum on $N nodes, " . ($mode eq "vn" ? $N * 4 : $N) . " processors\n";
print "inputs: $matrix $BLK $thresh $agglom $mapargs $peRotate\n";

`mkdir run-$runnum-$N-$BLK`;

print "Starting on $N nodes, " . ($mode eq "vn" ? $N * 4 : $N) . " processors\n";
my $dcmf = 16*1024*1024;
`qsub --env "BG_MAXALIGNEXP=1000:DCMF_RECFIFO=$dcmf" -A CharmRTS -q prod-devel --mode $mode -n $N -t 55 -o lu-$N-$BLK-$runnum.out -e lu-$N-$BLK-$runnum.error ./lu.trace $matrix $BLK $thresh $agglom 3 1\
6 16 +gz-trace +traceroot lu-$N-$BLK-$runnum +logsize 10000`;
`qsub --env "BG_MAXALIGNEXP=1000:DCMF_RECFIFO=$dcmf" -A PARTS --mode $mode -n $N -t 55 -o run-$runnum-$N-$BLK.out -e run-$runnum-$N-$BLK.error ./lu.trace $matrix $BLK $thresh $agglom $mapargs $peRotate +gz-trace +traceroot run-$runnum-$N-$BLK +logsize 10000 +LBOff`;
