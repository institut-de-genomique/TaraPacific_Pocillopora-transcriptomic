#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw( min max sum);
my ($in,$out);
my $maxMis=100;
my $minPCaligned=0;
my $minPCidentity=0;
my $minReadSize=0;
my $uniqOnly="NO";
GetOptions(
	'-in=s' => \$in, #FILE bam format
	'-out=s' => \$out, #FILE output path for the filtered bam file
	'-maxMis:i' => \$maxMis, #[1-100] Max number of mismatch on one read, default 100
	'-minPCaligned:i' => \$minPCaligned, #[0-100] Minimum percentage of read length aligned, default 0
	'-minPCidentity:i' => \$minPCidentity, #[0-100] Minimum percentage of read length aligned, default 0
	'-minReadSize:i' => \$minReadSize, #[0-200] Minimum read size, default 0
	'-uniqBestMatch=s' => \$uniqOnly #[YES|NO]  NO = random best match (default)
);

my %h;
open(BAM,"samtools view -h $in|") or die "$in not found";

unless($out =~/\.bam$/){$out=$out.".bam";}
open (OUT, "| samtools view -b - -o $out");
while (<BAM>){
	chomp;
	my $line=$_;
	my @cols=split(/\t/,$line);
	#header print	
	if($cols[0] =~/^@/){
		print OUT join("\t",@cols)."\n";
	}
	#Information about the match
	else{
		my $Alignlength=0;
		my $Readlength=0;
		my $TotalErrors=0;
		my $SecondaryAlignment="NO";
		if($cols[1]>=256){$SecondaryAlignment="YES";}

		$cols[5] =~ s/(\d+)[MX=DN]/$Alignlength+=$1/eg;

		$Readlength = length($cols[9]);

		($TotalErrors) = $line =~/NM:i:(\d+)/;
		
		my $uniqBestMatch="YES";
		if($line =~/XA/ and $cols[4]==0){$uniqBestMatch="NO";}

		#Comparison with selected criteria
		if ($SecondaryAlignment eq "YES"){
			$h{"Supplementary or secondary Alignement"}++;
		}
		if ($Readlength<$minReadSize){
			$h{"Too short reads"}++;
		}
		elsif ($Alignlength/$Readlength*100<$minPCaligned){
			$h{"Too short alignments"}++;	
		}
		elsif ($TotalErrors>$maxMis or $TotalErrors/$Alignlength*100>100-$minPCidentity){
			$h{"Too many mismatches or indels"}++;
		}
		else {
			if ($uniqBestMatch eq "YES"){
				$h{"uniqBestMatch"}++;
				print OUT "$line\n";
			}
			else{
				$h{"MultipleBestMatch"}++;			
				if ($uniqOnly eq "NO") {
				 	print OUT "$line\n";
				}
			}
		}
	}
}
close OUT;

my @report_list=("Supplementary or secondary Alignement", "Too short reads", "Too short alignments", "Too many mismatches or indels", "uniqBestMatch", "MultipleBestMatch");

open (REPORT, ">".substr($out,0,-3)."report");
foreach my $k (@report_list){
	if(exists $h{$k}){
		print REPORT substr($out,0,-4)."\t$k\t$h{$k}\n";
	}
	else{print REPORT substr($out,0,-4)."\t$k\t0\n"}
}
close REPORT;



