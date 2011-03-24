#!/usr/bin/perl

use strict;
use warnings;

my $N =  256;
my $mode = "vn";
my $runnum = "262";
my $BLK = 300;
my $matrix = 96000;
my $thresh = 310;
my $agglom = 100;
my $mapargs = "3 128 8";
my $peRotate = 16;


print "Starting job $runnum on $N nodes, " . ($mode eq "vn" ? $N * 4 : $N) . " processors\n";
print "inputs: $matrix $BLK $thresh $agglom $mapargs $peRotate\n";

my $dcmf = 16*1024*1024;
`qsub --env "BG_MAXALIGNEXP=1000:DCMF_RECFIFO=$dcmf" -A PARTS --mode $mode -n $N -t 55 -o run-$runnum-$N-$BLK.out -e run-$runnum-$N-$BLK.error ./lu.prod $matrix $BLK $thresh $agglom $mapargs $peRotate +LBOff`;
