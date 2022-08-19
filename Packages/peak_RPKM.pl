#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my $col = 7;
my $depth = 0;
my($output,$input);
GetOptions ("col=i" => \$col,"output=s" => \$output,"input=s" => \$input,"depth=f" => \$depth) or printCMD();
printCMD() if $depth == 0;

sub printCMD {
	select STDERR;
	print "\tUseage: peak_RPKM.pl \n";
	print "\tOptions:\n";
	print "\t\t--depth <reads depth/million>\n";
	print "\t\t--col <reads column number> (default: 7)\n";
	print "\t\t--input <input file>\n";
	print "\t\t--output <output file>\n";
	exit;
}

open IN,'<', $input or die "Can't open file $input: $!\n";
open OUT,'>', $output or die "Can't open file $output: $!\n";
while(<IN>){
	chomp;
	my $reads = (split /\t/)[$col - 1];
	my $start = (split /\t/)[1];
	my $end = (split /\t/)[2];
	my $length = ($end - $start)/1000;
	my $rpkm = $reads/($length * $depth);
	print OUT $_."\t".$rpkm."\n";
}

close IN;
close OUT;
