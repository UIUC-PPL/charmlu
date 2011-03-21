#!/usr/bin/perl

use strict;
use warnings;

my %blockTimingStart;
my %blockTimingStop;

die "usage: <filename> <x-axis> <y-axis>" if ($#ARGV + 1 < 1);

my $file = $ARGV[0];
shift @ARGV;

my $xaxis = 0;
my $yaxis = 0;

if ($#ARGV + 1 > 0) {
    ($xaxis, $yaxis) = @ARGV;
}

open FILE, ">", "temp" or die $!;
open RFILE, "<", "$file";

my $chareArray;
my $PEs;

for (<RFILE>) {
    if (/Chare Array size: (\d+) X/) {
        $chareArray = $1;
    }
    if (/RESULT procs:\s+(\d+)/) {
        $PEs = $1;
    }
    if (/^Block (\d+) queueing local LU at internalStep \d+, start time = \d+\.\d+, time = (\d+\.\d+)$/) {
        $blockTimingStart{$1} = $2;
    }
    if (/^Block (\d+) finished local LU at internalStep \d+, time = (\d+\.\d+)$/) {
        $blockTimingStop{$1} = $2;
    }
}

for my $key (sort {$a <=> $b} (keys %blockTimingStart)) {
    if (exists $blockTimingStop{$key} and exists $blockTimingStart{$key}) {
        my $timeStop = $blockTimingStop{$key};
        my $timeStart = $blockTimingStart{$key};
        my $timeStopPrev = $blockTimingStop{$key-1};
        my $delay = 0;

        my $diff = $timeStop - $timeStart;
        if (exists $blockTimingStop{$key-1}) {
            $delay = $timeStart - $timeStopPrev;
        }
        if ($diff >= 0 and $delay >= 0) {
            my $DGEMMsecs = 0.03;
            my $TUmodel = (($chareArray-$key-1)*($chareArray-$key-1)/$PEs)*$DGEMMsecs;
            print FILE "$key $diff $delay $TUmodel\n";
        }
    }
}

close FILE;
close RFILE;

# Now output for gnuplot
my $plotCmds = <<END;
set terminal pdf
END

if ($xaxis == 0) {
  $plotCmds .= "set autoscale x\n";
} else {
  $plotCmds .= "set xrange [0:$xaxis]\n";
}

if ($yaxis == 0) {
  $plotCmds .= "set autoscale y\n";
} else {
  $plotCmds .= "set yrange [0:$yaxis]\n";
}

$plotCmds .= <<END;
set title '$file'
set xlabel 'Execution Time (seconds)'
set ylabel 'LU Step'
set pointsize 0.5
plot 'temp' using 1:2 title 'Active panel time' with points lt 1 pt 5, \\
     'temp' using 1:3 title 'Time outside active panel' with points lt 2 pt 9, \\
     'temp' using 1:4 title 'Trailing update time per step' with lines lw 2
END

print "$plotCmds";
