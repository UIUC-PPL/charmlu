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
#unset key
set terminal postscript color solid
#set title 'test'
set border 1+2+4+8 lw 1.3
set xlabel 'Number of processors'
set logscale x
set xrange [100:2300]
set yrange [5:7.5]

set ylabel 'GFlops'

#set xtics (132, "" 264, 528, "" 1056, 2112)
set xtics axis in nomirror 132,4,2112
set ytics autofreq nomirror
set mxtics 20

plot '$temp' using 1:5 title 'redncb' with linespoints lt 5 pt 7, \\
     '$temp' using 1:2 title 'baseline' with linespoints lt 1 pt 5, \\
     '$temp' using 1:4 title 'isolateu' with linespoints lt 3 pt 6, \\
     '$temp' using 1:3 title 'noisolateactive' with linespoints lt 2 pt 9

# bm = 0.10
# lm = 0.10
# rm = 0.90
# gap = 0.03
# angle = gap / 2.0
# xoffset = gap / 1.6
# yoffset = angle / 2.0

# numSizes = 3
# size1 = 0.33 - bm / numSizes
# size2 = 0.33 - bm / numSizes
# size3 = 0.33 - bm / numSizes

# y1 = 700; y2 = 950;
# y3 = 2500; y4 = 4000;
# y5 = 10000; y6 = 15000;

# set multiplot
# set border 1+2+8
# set xtics nomirror
# set ytics nomirror
# set lmargin at screen lm
# set rmargin at screen rm
# set bmargin at screen bm
# set tmargin at screen bm + size1
# set yrange [y1:y2]

# plot '$temp' using 1:2 title 'baseline' with points lt 1 pt 5, \\
#      '$temp' using 1:3 title 'noisolateactive' with points lt 2 pt 9, \\
#      '$temp' using 1:4 title 'isolateu' with points lt 3 pt 6, \\
#      '$temp' using 1:5 title 'redncb' with points lt 5 pt 7
# #set pointsize 0.5

# prevSize = size1

# unset xtics
# unset xlabel
# set border 2+8
# set lmargin at screen lm
# set rmargin at screen rm
# set bmargin at screen bm + prevSize + gap
# set tmargin at screen bm + prevSize + size2 + gap
# set yrange [y3:y4]

# lx = lm + xoffset;
# ly = bm + prevSize + yoffset;
# rx = rm + xoffset;
# ry = bm + prevSize + yoffset;

# set arrow from screen lx,       ly                 to screen \\
#                       lx - gap, ly - angle         nohead
# set arrow from screen lx,       ly + gap           to screen \\
#                       lx - gap, ly + gap - angle   nohead
# set arrow from screen rx,       ry                 to screen \\
#                       rx - gap, ry - angle         nohead
# set arrow from screen rx,       ry + gap           to screen \\
#                       rx - gap, ry + gap - angle   nohead

# set ylabel 'GFlops'
# plot '$temp' using 1:2 title 'baseline' with points lt 1 pt 5, \\
#      '$temp' using 1:3 title 'noisolateactive' with points lt 2 pt 9, \\
#      '$temp' using 1:4 title 'isolateu' with points lt 3 pt 6, \\
#      '$temp' using 1:5 title 'redncb' with points lt 5 pt 7
# unset ylabel

# prevSize = size1 + size2 + gap

# unset xtics
# unset xlabel
# set border 2+8+4
# set lmargin at screen lm
# set rmargin at screen rm
# set bmargin at screen bm + prevSize + gap
# set tmargin at screen bm + prevSize + size3 + gap
# set yrange [y5:y6]

# lx = lm + xoffset;
# ly = bm + prevSize + yoffset;
# rx = rm + xoffset;
# ry = bm + prevSize + yoffset;

# set arrow from screen lx,       ly                 to screen \\
#                       lx - gap, ly - angle         nohead
# set arrow from screen lx,       ly + gap           to screen \\
#                       lx - gap, ly + gap - angle   nohead
# set arrow from screen rx,       ry                 to screen \\
#                       rx - gap, ry - angle         nohead
# set arrow from screen rx,       ry + gap           to screen \\
#                       rx - gap, ry + gap - angle   nohead

# set key at graph 0.3,0.9

# plot '$temp' using 1:2 title 'baseline' with points lt 1 pt 5, \\
#      '$temp' using 1:3 title 'noisolateactive' with points lt 2 pt 9, \\
#      '$temp' using 1:4 title 'isolateu' with points lt 3 pt 6, \\
#      '$temp' using 1:5 title 'redncb' with points lt 5 pt 7
END

my $opFile = "tmp-graph.ps";
open (OPIPE, "| gnuplot > $opFile && ps2pdf $opFile") or die $!;
print OPIPE "$plotScript";
