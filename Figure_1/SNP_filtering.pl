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


#Ouverture de 2 fichiers pour écrire les résultats
open (FREQ, ">$out.freq.tab");
open (COV, ">$out.cov.tab");
open (COMB, ">$out.combined.tab");
open (VCF, ">$out.vcf");

#Lecture du fichier $in nommé IN
if($in=~/.*\.gz$/){
	open(IN,"gunzip -c $in|") or die "$in not found";
}
else {
	open(IN,"$in") or die "$in not found";
}

my @samples;
while (<IN>){	
	#supprime le dernier caractère de chaque ligne (\n)
	chomp;	
	#construit un tableau @cols avec comme séaparateur la tabulation pour chaque ligne du fichier
	my @cols=split("\t",$_);
	#Enregistrement des noms des échantillons dans la variable $header
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
	#Traitement des SNPs biallelic
	elsif ($cols[0]!~/^#/ and length($cols[3])==1 and length($cols[4])==1){
		#tableau pour noter l'ordre des GT AD PL DP
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

		#Création d'un tableau pour enregistrer les valeurs de fréquence et et différentes infos
		my @frequences;
		my @coverages;
		my $NbSampleOver4reads=0;
		my $NbSampleWithAltAllele=0;
		my $NbSampleNA=0;

		#Boucle pour chaque sample
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
		#Ecriture des résultats
		if ($NbSampleOver4reads/($NbSampleOver4reads+$NbSampleNA+$NbSampleWithAltAllele)>=0.75){
			print VCF "$_\n";
		}
		print FREQ "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$NbSampleOver4reads\t$NbSampleWithAltAllele\t$NbSampleNA\t".join("\t",@frequences)."\n";
		print COV "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t$NbSampleOver4reads\t$NbSampleWithAltAllele\t$NbSampleNA\t".join("\t",@coverages)."\n";
	}
}
close FREQ;
close COV;
