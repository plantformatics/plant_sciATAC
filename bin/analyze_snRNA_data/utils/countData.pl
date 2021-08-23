#!/usr/bin/perl
use strict;
use warnings;

my %nuclear;
my %pt;
my %mt;
my %ids;

open F, "samtools view $ARGV[0] | " or die;
while(<F>){
	chomp;
	my @col = split("\t",$_);
	my $bc;
	foreach(@col){
		if($_ =~ /^CB:Z:/){
			my @cid = split(":",$_);
                        my @id = split("-", $cid[2]);
                        $bc = $id[0] . "-" . $ARGV[1];
			last;
		}
	}
	if($col[2] eq 'chrC'){
		$pt{$bc}++;
	}elsif($col[2] eq 'chrM'){
		$mt{$bc}++;
	}else{
		$nuclear{$bc}++;
	}
	$ids{$bc}++;
}
close F;

print "cellID\ttotal\tnuclear\tpt\tmt\tlibrary\n";

my @keys = keys %ids;
for (my $i = 0; $i < @keys; $i++){
	my $total = $ids{$keys[$i]};
	my $npt = 0;
	my $nmt = 0;
	my $nnc = 0;
	my $lib = $ARGV[1];
	if(exists $pt{$keys[$i]}){
		$npt = $pt{$keys[$i]};
	}
	if(exists $mt{$keys[$i]}){
		$nmt = $mt{$keys[$i]};
	}
	if(exists $nuclear{$keys[$i]}){
		$nnc = $nuclear{$keys[$i]};
	}
	print "$keys[$i]\t$keys[$i]\t$total\t$nnc\t$npt\t$nmt\t$lib\n";
}
