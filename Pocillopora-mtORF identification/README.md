
# Pocillo mt-ORF analysis

## Identification of mtORF sequences in _Pocillopora_ _Tara_ Pacific samples.
  
## Phylogeny of Tara Pocillopora mtORF with Johnston et al 2022 sequences

Johnston et al 2022 https://doi.org/10.1007/s00338-021-02107-9 
Accession IDs and Taxonomy: Sample-haplotype_johnston2022_NicheDifferences.tab (675 sequences including 8a and 11) 
Sequences : mtORF_johnston2022.fa 

Fasta modification to add sequences names of Johnston et al 2017
```bash
perl -e 'open(REF,"Sample-haplotype_johnston2022_NicheDifferences.tab");while(<REF>){chomp;@t=split("\t",$_);$t[1]=~s/Haplotype /type/; if($t[1]=~/^P\. (\D+)$/){$t[1]="type1a-$1"}elsif($t[1]=~/P\. verrucosa (.+)$/){$t[1]="type$1-verr
ucosa"}; $h{$t[0]}=$t[1]}; while(<>){if($_=~/^>(.+)/){$flag=0;if(exists $h{$1}){$flag=1;print ">$h{$1}-$1-johnston\n"}}elsif($flag==1){print "$_"}}' mtORF_johnston2022.fa > mtORF_johnston2022_NicheDifferences_annotate.fa
```

Concatenation of both files
```bash
cat Pocillopora_meandrina_v3_mtORF_11Islands.fa mtORF_johnston2022_NicheDifferences_annotate.fa > TP-Johnston_NichesDifferences_mtORF-Pocillo.fa
```

Multiple alignment with mafft v 7.490.
```bash
mafft TP-Johnston_NichesDifferences_mtORF-Pocillo.fa > TP-Johnston_NichesDifferences_mtORF-Pocillo.mafft
```

Neighbor-joining phylogeny with Mega version 10.1.5. 
Default parameters and 100 bootstrap 
File export in newick : TP-Johnston_mtORF-Pocillo.mafft.ML.nwk 

## mtORF Phylogeny representation in R

```r
setwd("/env/cns/proj/TaraPacifique/scratch/METAG_WORK/Quentin/Pocillo-mtORF/")
# 2.4 - Load associated metaData

CladeTable <- read.table("Variables11Islands.txt",h=T,sep="\t")

#Load tree
library(ape)
library(ggtree)
library(dendextend)

mtORFTree<-read.tree(file="TP-Johnston_NichesDifferences_mtORF-Pocillo.mafft.NJ-consensus.nwk")

#Groups tips by Pocillopora SVD group
mtORFTree <- groupOTU(mtORFTree, split(CladeTable$Sample,CladeTable$PocilloGG))

#Simplify labels.
mtORFTree$tip.label<-sub("-johnston","",mtORFTree$tip.label)
mtORFTree$tip.label<-sub("POC","",mtORFTree$tip.label)

#Plot the figure.
pdf(file="mtORF_MLtree_TP-Johnston.pdf",width=11,height = 11)
ggtree(mtORFTree,aes(x,y),layout="circular",branch.length = "none")+
  geom_tippoint(aes(color=group),size=3)+
  geom_tiplab2(size=3,hjust=-0.05)+
  scale_color_manual(values=c("black","tomato2","goldenrod","turquoise4", "palegreen", "darkorchid"))+
  geom_nodelab(aes(label=as.numeric(label)*100),size=3,hjust = 1.5)
dev.off()
```
