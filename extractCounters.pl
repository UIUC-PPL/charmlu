#!/usr/bin/perl


foreach $filename (@ARGV) {
    open(FILE, $filename);

    my %chunkCounts = ();

    while($line = <FILE>){
	
	if($line =~ m/BGP_UPC_Read_Counter_Value returned torus ([xyzmp]*) chunks = (.*)/){
#	    print "$1  $2\n";
	    $chunkCounts{$1} += $2;
	}
    }


    print "=========================\n";
    while ( my ($link, $chunks) = each(%chunkCounts) ) {
	$mb = $chunks * 32 / 1000000;
        print "$link => $mb MB\n";
    }

}
