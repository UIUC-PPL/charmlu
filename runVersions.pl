#!/usr/bin/perl


foreach $s (-1..1) {
    $tracedir = "traces-strat-$s";
    print `rm -fr $tracedir`;
    print `mkdir $tracedir`;
    print `qsub -n 64 --mode smp -t 10 ./lu-proj 16384 $s +traceroot $tracedir`;

}
