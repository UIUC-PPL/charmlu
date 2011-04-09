#!/opt/local/bin/perl

use strict;
use warnings;

use Text::CSV;

my $parser = Text::CSV->new();
my %xmap;

for (<>) {
    if ($parser->parse($_)) {
        my @cols = $parser->fields();
        my ($runnum, $commit, $machine, $type, $status, $jobid, $hour, $min,
            $procs, $matrix, $block, $memory, $pivot, $mapping, $tilex, $tiley,
            $rotate, $stride,$t1, $send, $time, $peak) = @cols;
        my $gflopsMax = 10.39;
        my $gflops = ($peak / 100) * $gflopsMax;
        push @{$xmap{$procs}}, $gflops;
    }
}

my $temp = "tmp-graph.dat";
chomp $temp;
open DATFILE, ">", "$temp" or die $!;

for my $proc (sort {$a <=> $b} (keys %xmap)) {
    print DATFILE "$proc @{$xmap{$proc}}\n"
}

close DATFILE;

my $plotScript = <<END;
set terminal postscript color
#set title 'test'
set border 1+2+4+8 lw 1.3
set xlabel 'Number of processors'
set logscale x
set xrange [100:2300]


set ylabel 'GFlops'

set xtics axis in nomirror 132,4,2112
set ytics autofreq nomirror
set mxtics 20
set mytics 5

set style line 1 lt 2 lc rgb "red" lw 3
set style line 2 lt 1 lc rgb "blue" lw 3
set style line 3 lt 4 lc rgb "green" lw 3
set style line 4 lt 8 lc rgb "orange" lw 3

# Exclusive prio classes
# plot '$temp' using 1:5 title 'redncb' with linespoints ls 1 pt 7, \\
#      '$temp' using 1:2 title 'baseline' with linespoints ls 2 pt 5, \\
#      '$temp' using 1:4 title 'isolateu' with linespoints ls 3 pt 13, \\
#      '$temp' using 1:3 title 'noisolateactive' with linespoints ls 4 pt 9
# set yrange [5:7.5]

plot '$temp' using 1:2 title 'send1' with linespoints ls 1 pt 7, \\
     '$temp' using 1:3 title 'send2' with linespoints ls 2 pt 5, \\
     '$temp' using 1:4 title 'send3' with linespoints ls 3 pt 13, \\
     '$temp' using 1:5 title 'send4' with linespoints ls 4 pt 9, \\
     '$temp' using 1:6 title 'send3' with linespoints ls 5 pt 9
set yrange [6.5:7.5]
END

my $opFile = "tmp-graph.ps";
open (OPIPE, "| gnuplot > $opFile && ps2pdf $opFile") or die $!;
print OPIPE "$plotScript";
