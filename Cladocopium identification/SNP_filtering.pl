#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Long;
my($in,$out,$cover);
GetOptions(
	'-in=s' => \$in, #vcf file. gz format accepted
	'-MinCover=s' => \$cover, #minimum coverage of a SNP in a sample to be kept
	'-out=s' => \$out #4 files: table with frequencies;table with coverages; combined table;vcf filtered
);


#4 files
open (FREQ, ">$out.freq.tab");
open (COV, ">$out.cov.tab");
open (COMB, ">$out.combined.tab");
open (VCF, ">$out.vcf");

if($in=~/.*\.gz$/){
	open(IN,"gunzip -c $in|") or die "$in not found";
}
else {
	open(IN,"$in") or die "$in not found";
}

my @samples;
while (<IN>){	
	chomp;	
	my @cols=split("\t",$_);
	if ($cols[0] =~ /^#/){
		print VCF "$_\n";
		if ($cols[0] eq "#CHROM"){
			my $header=substr($cols[0],1)."\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\tNbSampleOver4reads\tNbSampleWithAltAllele\tNbSampleNA\t";
			foreach(my $i=9;$i<=$#cols;$i++){
				my ($SampleName)=$cols[$i]=~/\/*(.+)_/;
				push (@samples, $SampleName);
			}
			$header.=join("\t",@samples);
			#Ecriture du header
			print FREQ "$header\n";
			print COV "$header\n";
			print COMB "CHROM\tPOS\tSAMPLE\tCOV\tFREQ\n";
		}
	}
	#For each biallelic SNP
	elsif ($cols[0]!~/^#/ and length($cols[3])==1 and length($cols[4])==1){
		#order of GT AD PL DP...
		my @list=split(":",$cols[8]);
		my $i=0;
		my $DP="NA";
		my $AD="NA";
		foreach my $info (@list){
			if($info eq "DP"){$DP=$i;}
			if($info eq "AD"){$AD=$i;}
			$i++;		
		}
		if($AD eq "NA"){die "SNP info not found"}
		my @frequences;
		my @coverages;
		my $NbSampleOver4reads=0;
		my $NbSampleWithAltAllele=0;
		my $NbSampleNA=0;

		#For each sample in the vcf file
		foreach(my $i=9;$i<=$#cols;$i++){	
			my $sample=$samples[$i-9];		
			my @SNPInfo=split(":",$cols[$i]);
			if($SNPInfo[$AD] eq "."){$SNPInfo[$AD]=".,.";}	
			my ($reference,$alternative)=$SNPInfo[$AD]=~/(.*),(.*)/;
			if($reference eq "."){$reference=0};
			if($alternative eq "."){$alternative=0};
			if($reference+$alternative<$cover){
				push (@frequences, "NA");
				push (@coverages, $reference+$alternative);
				$NbSampleNA++;
			}
			else{
				$NbSampleOver4reads++;
				print COMB "$cols[0]\t$cols[1]\t$sample\t".($reference+$alternative)."\t".($alternative/($reference+$alternative))."\n";
				my $frequence=sprintf("%.4g", $alternative/($reference+$alternative) );
				push (@frequences, $frequence);
				push (@coverages, $reference+$alternative);
				#Nombre d'échantillons avec allèle alternatif >= 0.05
				if($frequence>=0.05){
					$NbSampleWithAltAllele++;
				}	
			}
		}
		#Print results in FREQ and COV files
		if ($NbSampleOver4reads/($NbSampleOver4reads+$NbSampleNA+$NbSampleWithAltAllele)>=0.9){
			print VCF "$_\n";
		}
		print FREQ "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$NbSampleOver4reads\t$NbSampleWithAltAllele\t$NbSampleNA\t".join("\t",@frequences)."\n";
		print COV "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$NbSampleOver4reads\t$NbSampleWithAltAllele\t$NbSampleNA\t".join("\t",@coverages)."\n";
	}
}
close FREQ;
close COV;
