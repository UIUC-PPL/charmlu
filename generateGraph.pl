#!/opt/local/bin/perl

use strict;
use warnings;

use Text::CSV;

my $parser = Text::CSV->new();

my $temp = "tmp-graph.dat";
chomp $temp;
open DATFILE, ">", "$temp" or die $!;

my %xmap;

for (<>) {
    if ($parser->parse($_)) {
        my @cols = $parser->fields();
        my ($runnum, $commit, $machine, $type, $status, $jobid, $hour, $min,
            $procs, $matrix, $block, $memory, $pivot, $mapping, $tilex, $tiley,
            $rotate, $stride,$t1, $send, $time, $peak) = @cols;
        my $gflopsMax = 10.39;
        my $gflops = $peak * $gflopsMax * $procs;
        push @{$xmap{$procs}}, $gflops;
    }
}

for my $proc (keys %xmap) {
    print DATFILE "$proc @{$xmap{$proc}}\n"
}

close DATFILE;

my $plotScript = <<END;
set terminal postscript color solid
set title 'test'
set ylabel 'GFlops'
set xlabel 'Number of processors'
#set pointsize 0.5

plot '$temp' using 1:2 title 'baseline' with points lt 1 pt 5, \\
     '$temp' using 1:3 title 'noisolateactive' with points lt 2 pt 9, \\
     '$temp' using 1:4 title 'isolateu' with points lt 3 pt 6
END

my $opFile = "tmp-graph.ps";
open (OPIPE, "| gnuplot > $opFile && ps2pdf $opFile") or die $!;
print OPIPE "$plotScript";
