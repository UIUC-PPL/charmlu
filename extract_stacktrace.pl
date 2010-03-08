#!/usr/bin/perl

while ($line = <>){
    if($line =~ m/(0x[0-9a-fA-F]+)/ ){
#	print "$1\n";
	print `addr2line $1 -e lu-proj`;
    }
}
