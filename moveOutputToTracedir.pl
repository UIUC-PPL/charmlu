#!/usr/bin/perl


foreach $filename (@ARGV) {
    open(FILE, $filename);

    my $strategy = "unknown";
    my $result = "";

    $jobnum = -1;
    $tracedir = "";

    if($filename =~ m/([0-9]*)\.output/ ){
	$jobnum = $1;

	while($line = <FILE>){
	    
	    if($line =~ m/traceroot: (.*)\/lu-proj/) {
		$tracedir = $1;
	    }
	    
	    
	}
	
	
    }

    print "running: mv $jobnum.* $tracedir\n";
    print `mv $jobnum.* $tracedir`;

    close FILE;
}
