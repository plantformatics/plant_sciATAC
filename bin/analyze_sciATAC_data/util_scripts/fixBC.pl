#!/usr/bin/perl
use strict;
use warnings;

my %bcs;
my %chl;
my %mit;
my $tag = $ARGV[1];

open F, "samtools view -h $ARGV[0] | " or die;
while(<F>){
	chomp;
	if($_ =~ /^@/){
		print "$_\n";
		next;
	}
	my @col = split("\t",$_);
	if($_ =~ /XA:Z:/ && $col[4] < 30){
		my @tag = grep(/NM:i:/, @col);
		my @xa = grep(/XA:Z:/, @col);
		my @edits = split(":",$tag[0]);
		my $ed = $edits[2];
		my @mapped = split(";",$xa[0]);
		my $near = 0;
		foreach my $aln (@mapped){
			my @vals = split(",",$aln);
			my $dif = $vals[$#vals] - $ed;
			if($dif < 3){
				$near++;
			}
		}
		if($near > 1){
			next;
		}
	}
	my $bc;
	my $bcindex;
	my $hasbc = 0;
	for (my $i = 11; $i < @col; $i++){
		if($col[$i] =~ /BC:Z:/){
			$bcindex = $i;
			$bc = $col[$i];
			$hasbc++;
		}
	}
	if($hasbc > 0){
		$bc = "$bc" . "-" . "$tag";

		# initialize Pt/Mt hashes
		if(!$chl{$bc}){
			$chl{$bc} = 0;
		}
		if(!$mit{$bc}){
			$mit{$bc} = 0;
		}

		# correct array index
		$col[$bcindex] = $bc;
		my $line = join("\t",@col);
		$bcs{$bc}++;

		# skip Pt/Mt on print
		if($col[2] =~ /Pt/){
			$chl{$bc}++;
		}elsif($col[2] =~ /Mt/){
			$mit{$bc}++;
		}else{
			print "$line\n";
		}
	}else{
		next;
	}
	
}
close F;

my $temp = $tag . "_bc_counts.txt";
open (my $t1, '>', $temp) or die;
print $t1 "cellID\ttotal\tnuclear\tPt\tMt\tlibrary\ttissue\n";

my @ids = sort {$bcs{$b} <=> $bcs{$a}} keys %bcs;
for (my $j = 0; $j < @ids; $j++){
	my $nuc = $bcs{$ids[$j]} - $chl{$ids[$j]} - $mit{$ids[$j]};
	my $tissue = "leaf";
	print $t1 "$ids[$j]\t$ids[$j]\t$bcs{$ids[$j]}\t$nuc\t$chl{$ids[$j]}\t$mit{$ids[$j]}\t$tag\t$tissue\n";
}
close $t1;
