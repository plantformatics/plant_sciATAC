#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 [file.bam] [threads]\n" unless @ARGV == 2;

my $version = $ARGV[1];

open F, "samtools view -h -@ $ARGV[1] $ARGV[0] |" or die;
while(<F>){
	chomp;
	if($_ =~ /^@/){
		print "$_\n";
		next;
	}
	my @col = split("\t",$_);
	if($col[4] == 0){
		next;
	}
	my $num_mapped;
	my $lib;
	my $cb;
	my $ub;
	foreach(@col[11..$#col]){
		if($_ =~ /^NH/){
			my @val = split(":",$_);
			$num_mapped = $val[2];
		}elsif($_ =~ /^RG/){
			my @val = split(":",$_);
			$lib = $val[2];
		}elsif($_ =~ /^CB/){
			my @val = split(":",$_);
			$cb = $val[2];
			my $len = length($cb) - 2;
			$cb = substr($cb, 0, $len);
		}elsif($_ =~ /^UB/){
			my @val = split(":",$_);
			$ub = $val[2];
		}
	}
	if($num_mapped && $lib && $cb && $ub){
		if($num_mapped < 3){
			my $um = "UM:Z:" . $cb . $ub;
			print "$_\t$um\n";
		}
	}
}
close F;
