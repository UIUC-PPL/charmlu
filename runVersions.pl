#!/usr/bin/perl


foreach $m (500,1000) {
    foreach $s (32) {
	$tracedir = "traces-${s}--memconstrained-${m}MB";
	print `rm -fr $tracedir`;
	print `mkdir $tracedir`;
	$size = 1024 * $s;
	$strat = 6;
	print `qsub -n 64 --mode smp -t 10 ./lu-proj $size $strat $m +traceroot $tracedir`;
#       	print `qsub -n 64 --mode smp -t 10 ./lu $size $strat $m`;

    }
}

