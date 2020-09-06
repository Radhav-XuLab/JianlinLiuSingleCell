#!/usr/bin/perl

use strict;
use warnings;


my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
open(INFO1, $file1);
open(INFO2, $file2);

my %hash;
foreach my $lines2(<INFO2>){
	my @array2 = split /\t/, $lines2;
	$hash{$array2[3]} += $array2[5];
}


#foreach my $k (sort keys %hash) {
    #print "$k => $hash{$k}\n";
#}




foreach my $lines1(<INFO1>){
	chomp $lines1;
	#print $lines1;
	my @array1 = split /\t/, $lines1;
	if (exists $hash{$array1[3]} ){
		print join("\t",@array1,$hash{$array1[3]},"\n");
	}

}

close(INFO1);
close(INFO2);
