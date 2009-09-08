#!/usr/bin/perl


foreach $filename (@ARGV) {
    open(FILE, $filename);

    my %chunkCounts = ();
    my $strategy = "unknown";
    my $result = "";

    while($line = <FILE>){
	
	if($line =~ m/BGP_UPC_Read_Counter_Value returned torus ([xyzmp]*) chunks = (.*)/){
#	    print "$1  $2\n";
	    $chunkCounts{$1} += $2;
	}

	if($line =~ m/Using (.*)/){
	    $strategy = $1;
	}

	if($line =~ m/RESULT (.*)/){
	    $result = $1;
	}

    }


 #   print "==============  $strategy  ===========\n";
  #  print "$result\n";
    $totalmb = 0;
    while ( my ($link, $chunks) = each(%chunkCounts) ) {
	$mb = $chunks * 32 / 1000000;
#        print "$link => $mb MB\n";
	$totalmb += $mb;
    }

    print "$strategy\t$filename\t$totalmb\t$result\n";


}
