#!/usr/bin/perl
use strict;
use warnings;

my %mtx;
my $total = 0;

open F, "samtools view $ARGV[0] | " or die;
while(<F>){
	$total++;
	if(($total % 1000000) == 0){
		print STDERR " - iterated over $total BAM records ...\n";
	}
	chomp;
	my $geneID;
	my $cellID;
	my @col = split("\t",$_);
	foreach(@col){
		if($_ =~ /^CB:Z/){
			my @cid = split(":",$_);
			my @id = split("-", $cid[2]);
			$cellID = $id[0] . "-" . $ARGV[1];
		}elsif($_ =~ /^GX:Z/){
			my @gid = split(":", $_);
			my @id = split(/\./, $gid[2]);
			$geneID = $id[0];
		}

	}
	if($geneID && $cellID){
		$mtx{$cellID}{$geneID}++;
	}
	
}
close F;

my @cells = keys %mtx;
for (my $i = 0; $i < @cells; $i++){
	my @genes = keys %{$mtx{$cells[$i]}};
	for (my $j = 0; $j < @genes; $j++){
		print "$genes[$j]\t$cells[$i]\t$mtx{$cells[$i]}{$genes[$j]}\n";
	}
}
