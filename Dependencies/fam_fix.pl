#!/usr/bin/perl

use strict;
use warnings;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

unless ($infile and $outfile) {
	die "\n\nERROR: You need to specify a fam file and an output file name.\n\n";
}

unless (-f $infile) {
	die "\n\nERROR: Cannot find $infile. Check file path";
}

open IN, "$infile";
open OUT, ">$outfile";

while (<IN>) {
	my $line = $_;
	chomp $line;
	
	my @l = split/\s+/, $line;
	
	my $fid = $l[0];
	my $iid = $l[1];
	
	my $joined = $fid."_".$iid;
	
	unless ($fid =~ /msp\d/ and $iid =~ /msp\d/) {
		$fid = $joined;
		$iid = $joined;
	}
	
	print OUT "$fid\t$iid\t0\t0\t0\t-9\n";

}

close IN;
close OUT;

exit;
