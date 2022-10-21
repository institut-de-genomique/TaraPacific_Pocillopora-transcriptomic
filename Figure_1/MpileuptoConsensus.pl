#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);
my($mpileup,$out,$Mincov);
GetOptions(
	'-mpileup=s' => \$mpileup, #FILE in mpileup format from bcftools or samtools
	'-Mincov=s' => \$Mincov, #INTEGER Minimum coverage per position to call the consensus base
	'-out=s' => \$out #FILE to write the consensus sequence
);

my %h;
my %maxcov;
my %Mis;
my ($Nb_read,$List_reads);
open(MP,"$mpileup") or die "$mpileup not found";
<MP>;	
while (<MP>){
	chomp;
	my @cols=split("\t",$_);
	unless (exists $maxcov{$cols[0]}){$maxcov{$cols[0]}=0;}
	unless (exists $Mis{$cols[0]}){$Mis{$cols[0]}=0;}
	if($cols[3]>$maxcov{$cols[0]}){$maxcov{$cols[0]}=$cols[3];}
	if($cols[3]<$Mincov){
		$h{$cols[0]}.="N";
	}
	else{
		my %count=();
		while ($cols[4] =~/([-\+])(\d+)(\D+?)[,.\*\^]/){
			if($1 eq "+"){
				$count{uc($3)}++;
			}
			$cols[4]=~s/$2$3//;
		}
		foreach my $k (keys %count){
			if($count{$k}>0.5*$cols[3]){$h{$cols[0]}.=$k;$Mis{$cols[0]}++;}
		}
		foreach my $i ("A","T","C","G",",",".","a","t","c","g","*"){
			$count{$i}=0;
		}
		$count{$1}++ while ($cols[4] =~ /(.)/g);
		if($count{"."}+$count{","}>0.5*$cols[3]){$h{$cols[0]}.=$cols[2];}
		elsif($count{"A"}+$count{"a"}>0.5*$cols[3]){$h{$cols[0]}.="a";$Mis{$cols[0]}++}
		elsif($count{"T"}+$count{"t"}>0.5*$cols[3]){$h{$cols[0]}.="t";$Mis{$cols[0]}++}
		elsif($count{"C"}+$count{"c"}>0.5*$cols[3]){$h{$cols[0]}.="c";$Mis{$cols[0]}++}
		elsif($count{"G"}+$count{"g"}>0.5*$cols[3]){$h{$cols[0]}.="g";$Mis{$cols[0]}++}
		else {
			$h{$cols[0]}.="N";
		}
	}
}
close MP;

my $i=0;
my %tri;
foreach my $name (keys %maxcov) {
	$i+=0.1;
	$tri{$maxcov{$name}+$i}=$name;
}
my @list;
foreach my $qtt (sort {$b <=> $a} keys %tri) {
	push @list,$tri{$qtt};
}
my $file=substr($mpileup,0,4);
open (OUT, ">$out");	
foreach my $name (@list) {
	if (length($h{$name})>10){
		while($h{$name}=~/^N/){$h{$name}=substr($h{$name},1)}
		while($h{$name}=~/N$/){$h{$name}=substr($h{$name},0,-1)}
		if (length($h{$name})>10){
			my $Ncount = ($h{$name} =~ tr/N/n/);
			my $alternativcount = ($h{$name} =~ tr/atcg/atcg/);
			print OUT ">$file"."_$name".":$maxcov{$name}:$Ncount:$alternativcount:$Mis{$name}\n$h{$name}\n";	
		}
	}
}
close OUT;
