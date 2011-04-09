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
        my $aspect = $tiley / $tilex;
        my $upar;
        if (int($rotate) == 0) {
            $upar = $tiley;
        } else {
            $upar = $tilex / $rotate * $tiley;
        }
        my $percentU = $upar / $procs * 100.0;

        ### Replace first hash key with x to grah
        $xmap{$send}{$procs} = $gflops;

        ### Exclusive
        #$xmap{$procs}{$type} = $gflops;
    }
}

my $temp = "tmp-graph.dat";
chomp $temp;
open DATFILE, ">", "$temp" or die $!;

my @tlist = (132, 528, 2112);
my @exclusive = ("baseline", "noisolateactive", "isolateu", "redncb");

for my $key (sort {$a <=> $b} (keys %xmap)) {
    print DATFILE "$key ";
    foreach my $p (@tlist) {
        if (exists $xmap{$key}{$p}) {
            print DATFILE "$xmap{$key}{$p} ";
        } else {
            print DATFILE "nan ";
        }
    }
    print DATFILE "\n";
}

close DATFILE;

my $plotScript = <<END;
set terminal postscript color
#set title 'test'
set border 1+2+4+8 lw 1.3

set ylabel 'GFlops'

set style line 1 lt 2 lc rgb "red" lw 3
set style line 2 lt 1 lc rgb "blue" lw 3
set style line 3 lt 4 lc rgb "green" lw 3
set style line 4 lt 8 lc rgb "orange" lw 3

### Exclusive prio classes
# set logscale x
# set xlabel 'Number of processors'
# set xrange [100:2300]
# set yrange [5:7.5]
# set xtics axis in nomirror 132,4,2112
# set ytics autofreq nomirror
# set mxtics 20
# set mytics 5
# plot '$temp' using 1:5 title 'redncb' with linespoints ls 1 pt 7, \\
#      '$temp' using 1:2 title 'baseline' with linespoints ls 2 pt 5, \\
#      '$temp' using 1:4 title 'isolateu' with linespoints ls 3 pt 13, \\
#      '$temp' using 1:3 title 'noisolateactive' with linespoints ls 4 pt 9

### Send limit
set xlabel 'Locally incomplete message sends'
set xtics autofreq nomirror 1,1,5
set ytics autofreq nomirror
set yrange [6:7.5]
set mytics 5
plot '$temp' using 1:2 title '132 processors' with linespoints ls 1 pt 7, \\
     '$temp' using 1:3 title '528 processors' with linespoints ls 2 pt 5, \\
     '$temp' using 1:4 title '2112 processors' with linespoints ls 3 pt 13

### Aspect ratio
# set xlabel 'Aspect ratio'
# set xtics autofreq nomirror
# set ytics autofreq nomirror
# set xrange [0:5]
# set mxtics 5
# set mytics 5
# plot '$temp' using 1:2 title '132 processors' with linespoints ls 1 pt 7, \\
#      '$temp' using 1:3 title '528 processors' with linespoints ls 2 pt 5, \\
#      '$temp' using 1:4 title '2112 processors' with linespoints ls 3 pt 13

### Stride
# set xlabel 'Active panel PEs/node'
# set xtics autofreq nomirror
# set ytics autofreq nomirror
# set mxtics 5
# set mytics 5
# plot '$temp' using 1:2 title '132 processors' with linespoints ls 1 pt 7, \\
#      '$temp' using 1:3 title '528 processors' with linespoints ls 2 pt 5, \\
#      '$temp' using 1:4 title '2112 processors' with linespoints ls 3 pt 13

### Rotation
# set xlabel 'Percent of PEs performing triangluar solves'
# set xtics autofreq nomirror
# set ytics autofreq nomirror
# #set xrange [0:500]
# set mxtics 5
# set mytics 5
# plot '$temp' using 1:2 title '132 processors' with linespoints ls 1 pt 7, \\
#      '$temp' using 1:3 title '528 processors' with linespoints ls 2 pt 5, \\
#      '$temp' using 1:4 title '2112 processors' with linespoints ls 3 pt 13

END

my $opFile = "tmp-graph.ps";
open (OPIPE, "| gnuplot > $opFile && ps2pdf $opFile") or die $!;
print OPIPE "$plotScript";
