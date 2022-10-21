# Identifications of Cladocopium lineages 
We used three independant methodes to identify Cladocopium lineages.  

## 1. Using ITS2 sequences 
### A. Barplot of ITS2 profiles (Supplementary Figure X) 

ITS2 sequences are available under the PRJEB52458 project and were obtained following the procedure detailed here : https://arxiv.org/abs/2207.02475  
ITS2 sequences were processed with the Symportal pipeline and the results are available here : https://doi.org/10.5281/zenodo.4061796  
ITS2 profil abundances in each colony is in : TARA_PACIFIC_METAB_ITS2_coral_its2_type_profiles_absolute_abund_and_meta_v1.csv  

The R code to produce Supplementary Figure 4 is below  
```r
#Required library
library(ggplot2)

ITS2<-read.table("TARA_PACIFIC_METAB_ITS2_coral_its2_type_profiles_absolute_abund_and_meta_v1.csv",sep=",",check.names = F,h=T)
#columns and lignes selection
colnames(ITS2)<-ITS2[6,]
ITS2<-ITS2[7:nrow(ITS2),]
colnames(ITS2)[2]<-"TaraSampleName"

#Addition of concatenated sample Names
Corres<-read.table("Variables11Islands.txt",sep="\t",h=T)
Corres<-Corres[,c("Samples","TaraSampleName")]
tab<-merge(Corres,ITS2,by="TaraSampleName")

#conversion in column format with column names
tab2<-melt(tab,id.vars = colnames(tab[,1:3]), measure.vars = colnames(tab[,4:ncol(tab)]))
colnames(tab2)<-c("TaraSampleName","SampleName","ProfileNumber","ProfileName","Value")
tab2$Value<-as.numeric(tab2$Value)
#Table filtering
tab3<-tab2[tab2$Value>0,2:ncol(tab2)]
#Parsing of Profile name
tab3$Clade<-substr(tab3$ProfileName,1,1)
tab3$Type<-sub('-.*','',tab3$ProfileName,perl = T)
#Parsing of Sample name
tab3$Island<-substr(tab3$SampleName,1,3)
tab3$SiteColo<-substr(tab3$SampleName,4,10)
#Island ordering
tab3$Island<-factor(tab3$Island,levels=sort(unique(tab3$Island),decreasing = T))

#Specific color selection
set<-c("#3a867c","#3a9a7c",
       "#fa04b3",
       "#3a86af",
       "#c4ffba",
       "#B47229",
       "#658E67",
       "#91569B","#a5569b",
       "#042955",
       "#3a727c","#3a867c","#000000",
       "#3aae7c",
       "#e96b13","#e97f13","#e99313",
       "#3a86b9","#3a86cd","#3a86e1","#3a86f5",
       "#fa04c7",
       "#d7d40e",
       "#4e9af5",
       "#E41A1C","#ee1a1c","#f81a1c","#ff1a1c")

pdf(file="ITS2_type2_Pocillopora.pdf",width=12)
ggplot(tab3,aes(x=SiteColo,y=Value,fill=Type))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=set)+
  facet_grid(.~Island,scale="free",space="free")+  theme(axis.text.x=element_text(angle=90),panel.grid=element_blank(),axis.ticks.x=element_blank(),panel.background=element_blank(),axis.line=element_blank(),legend.key.size=unit(2,"mm"))+
  guides(fill=guide_legend(ncol=1))
dev.off()
```

### Hierarchical clustering of Samples ITS2 distances (Figure 1B)

#Genetic clades for the host and the symbiont based on the SNP.

```r
Variable<-read.table("Variables11Islands.txt",sep="\t",h=T)

#Distance between samples based on ITS2 sequences.
library(vegan)
library(RColorBrewer)
library(ape)
library(ggrepel)
colors<-brewer.pal(8,"Set1")
pal<-colorRampPalette(colors)
set<-sample(pal(10),replace = F)

#Distance table
tab<-read.table("TARA_PACIFIC_METAB_ITS2_coral_between_sample_distances_C_braycurtis_distances_sqrt_v1.dist",sep="\t",h=F)
#Selection of Pocillopora samples from the 11 first islands.
tab2<-tab[tab$V1%in%Variable$TaraSampleName,c(TRUE,TRUE,tab$V1%in%Variable$TaraSampleName)]
colnames(tab2)<-c("TaraSampleName","Id",tab2$V1)
#Addition of SampleNames
tab2<-merge(Variable[,c("Samples","TaraSampleName")],tab2,by.x="TaraSampleName",by.y="TaraSampleName")
rownames(tab2)<-tab2$TaraSampleName

#Conversion in matrix 
mat<-as.matrix(tab2[,4:ncol(tab2)])
mat<-mat[,rownames(mat)]

#Addition of Sample names
rownames(mat)<-tab2$Samples
colnames(mat)<-tab2$Samples

#Elimination of samples containing 2 Cladocopium ITS2 profiles.
GoodSamples<-Variable[!Variable$Sample%in%c("I10S01C006POC","I05S02C002POC"),1]
mat2<-mat[rownames(mat)%in%GoodSamples,rownames(mat)%in%GoodSamples]

#Hierarchical clustering of Samples ITS2 distances.

dendITS<-hclust(as.dist(mat2),method = "average") 
dendITS2 <- as.dendrogram(dendITS) %>%
  set("labels_cex", c(.5)) %>%
  set("labels_cex", c(.5)) %>% 
  set("branches_lwd",c(0.5)) %>%
  set("labels_col",c("#117733", "#cc6677" ,"#44aa99", "#332288","#6699cc")[as.factor(unlist(lapply(dend$labels[dend$order],function(x){Variable[Variable$Samples==x,"SymbioGG"]})))]) 

ggdendITS2 <- as.ggdend(dendITS2)


#ggplot(ggdendITS2, horiz = TRUE, theme=NULL) + 
  theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.x=element_line(color="black"),panel.background = element_blank())+
  labs(y="",x="ITS2 dendrogram")
```
## 2. Using SNPs identified on transcriptomic reads (Figure 1A).  

```r
tab<-read.table("SNP-Filtering-Coverage-Biallelic-Population-Structure-DualT-Cladocopium-C1-CDS-All-Samples-Perform-Genotyping-Removed-D2-Biallelic-SQ-30-Frequency-NA-Filtering.txt",h=T,check.names = F)
colnames(tab)<-sub("_.*","",colnames(tab))
#Selection of 82 samples containing a single Cladocopium species.
tab2<-tab[,colnames(tab)%in%Variable[Variable$Symbio_Genus=="C","Samples"]]
dendSNP<-hclust(dist(t(tab2)),method = "ward.D2")

#Optimal number of cluster
library(factoextra)
OptimalKGapSNP<-fviz_nbclust(tab2, hcut, method = "gap_stat",k.max=12,nboot=50) +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "gap_stat method")

ggplot(OptimalKGapSNP$data,aes(x=clusters,y=gap,group = 1))+
  geom_point(size=2,alpha=0.5)+
  geom_line()+
  geom_errorbar(aes(ymin=ymin,ymax=ymax),width=0.2)+
  geom_vline(xintercept = 5,linetype=2,color="red")+
  theme_bw()

dendSNP2 <- as.dendrogram(dendSNP) %>%
  set("labels_cex", c(.5)) %>%
  set("labels_cex", c(.5)) %>% 
  set("branches_lwd",c(0.5)) %>%
  set("labels_col",c("#117733", "#cc6677" ,"#44aa99", "#332288","#6699cc")[as.factor(unlist(lapply(dendSNP$labels[dendSNP$order],function(x){Variable[Variable$Samples==x,"SymbioGG"]})))]) 

ggdendSNP2 <- as.ggdend(dendSNP2)

#ggplot(ggdendSNP2, horiz = TRUE, theme=NULL) + 
  theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.x=element_line(color="black"),panel.background = element_blank())+
  labs(y="",x="SNPs dendrogram")


#Ordered dendrogram according to Clado Lineages
dendITS4 <-rotate(dendITS2,Variable[order(Variable$SymbioGG),"Samples"][Variable[order(Variable$SymbioGG),"Samples"]%in%labels(dendITS2)])
dendSNP4 <-rotate(dendSNP2,Variable[order(Variable$SymbioGG),"Samples"][Variable[order(Variable$SymbioGG),"Samples"]%in%labels(dendSNP2)])
  
dendCombined<-dendlist(dendITS4, dendSNP4)
pdf(file="combined_dendrogram.pdf",height=10,width=12)
dend_diff(dendCombined)
dev.off()
```

## 3. Analysis of psbA<sup>ncr</sup> sequences in _Pocillopora_ samples

We analyzed the non-coding region of the psbA gene (psbA<sup>ncr</sup>) to assign _Cladocopium_ lineages to known _Cladocopium_ species following the method used in Johnston et al 2022 and Turnham et al 2021. The code below describes : (A) The mapping of metagenomic reads of the 82 _Pocillopora_ corals containing Cladocopium symbionts on 3 psbAncr sequences (C. latusorum MW819767.1, C. pacificum MW861717, and C. goreaui KF572161.1) the build of consensus sequences for each sample. (B) The alignement then construction of Bayesian phylogeny using 2 psbA<sup>ncr</sup> sequences for each Cladocopium clade identified in Johnston et al 2022 and and the 82 consensus sequences. (C) The representation of the phylogeny with R.


### A. Identification of psbA<sup>ncr</sup> sequences in each _Pocillopora_ colonies


#1) Extraction of psbAncr reads from 82 Pocillopra samples
#require bwa-mem2/2.2.1 and samtools/1.15.1
#

bwa-mem2 index psbA_Cladocopium.fa
cat ReadsetsPOC_MetaG.list | while read a b c ;do jobify -b -q normal -c 4 -t 4:00:00 "bwa-mem2 mem -t 4 -M psbA_Cladocopium.fa ${b} ${c} | samtools view -b -@ 4 -F 4 - -o ${a}_PsBA.aln.bam";done

for i in *_PsBA.aln.bam ;do jobify -b -q normal -c 1 -t 1:00:00 perl /env/cns/proj/projet_CNM/script/Quentin/BamFiltration.pl -in $i -out ${i%.bam}.filtered100-90.bam -minPCaligned 100 -minPCidentity 90; d
one

for i in *_PsBA.aln.filtered.bam ;do echo -ne "$i\t" ;samtools view $i | awk '{print $3}' | sort | uniq -c | sort -nrk1,1 |head -1 | awk '{print $2"\t"$1}' ;done >

Bilan_mapping_Cladocopium100-90.tab

#Manual curation of Bilan file (to remove not analyzed samples).
cat Bilan_mapping_Cladocopium100-90.tab | while read a b c ;do jobify -b -q normal -c 1 -t 1:00:00 samtools sort -o ${a%.bam}.sort.bam $a; done
for i in *_PsBA.aln.filtered100-90.sort.bam ;do jobify -b -q normal -c 1 -t 1:00:00 samtools index $i; done

#Manual consensus
cat Bilan_mapping_Cladocopium100-90.tab | while read a b c ;do jobify -b -q normal -c 1 -t 00:05:00 samtools mpileup -Aa -f psbA_Cladocopium.fa -r $b ${a%.bam}.sort.bam -o ${a%.bam}.sort.mpileup;done
for i in *PsBA.aln.filtered100-90.sort.mpileup ;do jobify -b -q normal -c 1 -t 00:15:00 perl /env/cns/proj/TaraPacifique/bin/Quentin/MpileuptoConsensus3.pl -mpileup $i -Mincov 4 -out ${i}.consensus.fa;done

#concat
for i in *_PsBA.aln.filtered100-90.sort.mpileup.consensus.fa;do echo -e ">${i%_*}"; tail -n +2 $i ;done > All_psbA-Cladocopium100-90.aln.mpileup.consensus.fa

#OPTIONEL
#Second mapping sur la ref sélectionnée.
bwa-mem2 index C_goreauii_rt113_psbA.fa
bwa-mem2 index C_latusorum_MW819767.1_Pal09_620_psbA.fa
bwa-mem2 index C_pacificum_MW861717_psbA.fa

cat ReadsetsPOC_MetaG.list | while read a b c ;do REF=`grep ${a} Bilan_mapping_CladocopiumBis2.tab | awk '{print $2}'`; jobify -b -q normal -c 4 -t 4:00:00 "bwa-mem2 mem -t 4 ${REF}.fa ${b} ${c} | samtools
 view -b -@ 4 -F 4 - -o ${a}_${REF}_2ndMapping.aln.bam";done


for i in *_psbA_2ndMapping.aln.bam ;do jobify -b -q normal -c 1 -t 1:00:00 samtools sort -o ${i%.bam}.sort.bam $i; done
for i in *_psbA_2ndMapping.aln.sort.bam ;do jobify -b -q normal -c 1 -t 00:05:00 samtools index $i; done
for i in *_psbA_2ndMapping.aln.sort.bam ;do jobify -b -q normal -c 1 -t 00:05:00 samtools mpileup -Aa -f psbA_Cladocopium.fa $i -o ${i%.bam}.mpileup;done
for i in *_psbA_2ndMapping.aln.sort.bam ;do jobify -b -q normal -c 1 -t 00:15:00 perl /env/cns/proj/TaraPacifique/bin/Quentin/MpileuptoConsensus3.pl -mpileup ${i%.bam}.mpileup -Mincov 2 -out ${i%.bam}.mpil
eup.consensus.fa;done

#concat
for i in *_psbA_2ndMapping.aln.sort.mpileup.consensus.fa;do echo -e ">${i%_*}"; tail -n +2 $i ;done > All_psbA-CladocopiumV3.aln.mpileup.consensus.fa

#module load bcftools/1.15
for i in *_PsBA.aln.filtered.sort.bam;do jobify -b -q normal -c 1 "bcftools mpileup -ABOu -f psbA_C.pacificum.fa $i | bcftools call --ploidy 1 -mv -Oz -o ${i%.bam}.call.vcf.gz;bcftools norm -f psbA_C.pacif
icum.fa ${i%.bam}.call.vcf.gz -Ob -o ${i%.bam}.call.norm.bcf;bcftools filter --IndelGap 0 ${i%.bam}.call.norm.bcf -Ob -o ${i%.bam}.calls.norm.flt-indels.bcf;bcftools index ${i%.bam}.calls.norm.flt-indels.b
cf;cat psbA_C.pacificum.fa | bcftools consensus -M n -a N ${i%.bam}.calls.norm.flt-indels.bcf > ${i%.bam}.consensus.fa";done



for i in *_PsBA.aln.filtered.sort.bam;do jobify -b -q normal -c 1 "bcftools mpileup -A -f psbA_Cladocopium.fa $i | bcftools call -c | vcfutils.pl vcf2fq -Q 0 -l 0 -d 0 > ${i%.sort.bam}.bcftools-consensus.f
a";done

#####################
#Phylogeny Revision##
#####################
#MrByes
#installé en local : make --preset "current directory"
#puis ./mb
#ne fonctionne pas en sshell...
#Fichier de cofig MrbayesCommands.nex
begin mrbayes;
   set autoclose=yes nowarn=yes;
   execute All_psbA.aln.mpileup.consensus_woN_AllJohnston2022_muscle.nexus;
   lset nst=6 rates=invgamma;
   prset brlenspr=clock:birthdeath;
   mcmc nruns=1 ngen=1000000 samplefreq=100 printfreq=100 diagnfreq=1000 file=All_psbA.aln.mpileup.consensus_woN_AllJohnston2022_muscle.nexus1;
   sump;
   sumt;
end;

jobify -b -q normal -c 1 --mem="10G" ./mb -i MrbayesCommands.nex


### B. Bayesian phylogeny of psbA<sup>ncr</sup> sequences

### C. Tree representation on R (Figure 1C)

################################
##SymbioSpecies identification##
################################
#Johnston et al 2022: https://zenodo.org/record/6710608#.YxCxB7TP1aQ from https://pubmed.ncbi.nlm.nih.gov/35960256/
#Tuhrnam et al 2021: Dataset 2 from https://www.nature.com/articles/s41396-021-01007-8
#directory: /env/cns/proj/TaraPacific/scratch/METAT_WORK/Quentin/SymbioSpecies/

setwd("/env/cns/proj/TaraPacifique/scratch/METAT_WORK/Quentin/SymbioSpecies/")
library(ape)
library(ggtree)

#Test 1 Johnston tree

Tree<-ape::read.nexus("psbA_Johnston2022.tree.nexus")

Tree2<-root(Tree, node = 179)

test<-as_tibble(Tree2)

ggtree(Tree2,size=0.5,aes(x,y)) + 
  geom_treescale()+
  geom_text(aes(label=node), hjust=-.3,size=3)+
  geom_nodepoint()+
  geom_tiplab(size=3)+
  theme_tree2()+
  xlim(0,0.15)
  

theme(plot.margin = unit(c(10,10,10,10), "mm"))
#as_ylab=
############################
#Test 2 Jonhnston alignment#
############################
setwd("/env/cns/proj/TaraPacifique/scratch/METAT_WORK/Quentin/SymbioSpecies/bin/")

Tree<-ape::read.nexus("Figure5a-Cladocopium_psbA.nex.con.tre")

pdf(file="Johnston-tree.pdf",height=12,width=9)
ggtree(Tree,size=0.5,aes(x,y)) + 
  geom_treescale()+
  geom_text(aes(label=node), hjust=-.3,size=3)+
  geom_nodepoint()+
  geom_tiplab(size=2.5)+
  theme_tree2()+
  xlim(0,0.20)
dev.off()

####################
library(ape)
library(ggtree)

clade<-read.table("/env/cns/proj/TaraPacifique/ARTICLES/MetaTArticle/ITS2-analysis/Variables11Islands.txt",h=T,sep="\t")
clade<-clade[!(clade$SymbioGG%in%c("CD","D")),]

psbA<-read.tree("All_psbA.aln.mpileup.consensus_woN_muscle.newick")
plot(psbA)

psbA<-root(psbA, node = 209)

psbA2 <- groupOTU(psbA, split(clade$Sample,clade$Islands))


ggtree(psbA2,size=0.5,aes(x,y))+ geom_treescale()+
  geom_text(aes(label=node), hjust=-.3,size=3)+
  geom_nodepoint() + 
  geom_tiplab(size=3,aes(color=group))+
  theme_tree2()+
  xlim(0,20)


############
#johnston+TP
############
setwd("/env/cns/proj/TaraPacifique/scratch/METAT_WORK/Quentin/SymbioSpecies/bin/")

Tree<-ape::read.nexus("All_psbA-Cladcopium100-90.aln.mpileup.consensus_woN_SelectedJohnston2022_clustal_simple.nexus1.con.tre")
clade<-read.table("/env/cns/proj/TaraPacifique/ARTICLES/MetaTArticle/ITS2-analysis/Variables11Islands.txt",h=T,sep="\t")
clade<-clade[!(clade$SymbioGG%in%c("CD","D")),]
clade<-clade[!(clade$SymbioGG%in%c("CD","D")),]


Tree2 <- groupOTU(Tree$con_50_majrule, split(clade$Sample,clade$SymbioGG))

pdf(file="SelectedJohnston-tree_TP_Cladocopium100-90_clustal_1000000_linear.pdf",height=9,width=7)
ggtree(Tree2,size=0.5,aes(x,y)) + 
  geom_treescale()+
  geom_tiplab(size=2,aes(color=group))+
  scale_color_manual(values=c("red","#117733", "#cc6677" ,"#44aa99", "#332288","#6699cc"))+
  geom_nodelab(size=2,aes(label=as.numeric(label)*100,x=x+0.20/100))+
  theme_tree2()+
  xlim(0,0.20)
dev.off()

pdf(file="SelectedJohnston-tree_TP_Cladocopium100-90_clustal_1000000_circular.pdf",height=10,width=10)
ggtree(Tree2,size=0.5,aes(x,y),layout="circular") + 
  geom_treescale()+
  geom_tiplab(size=2.5,aes(color=group))+
  scale_color_manual(values=c("red","#117733", "#cc6677" ,"#44aa99", "#332288","#6699cc"))+
  geom_nodelab(size=2,aes(label=as.numeric(label)*100,x=x+0.20/100))+
  theme_tree2()+
  xlim(0,0.20)
dev.off()


install.packages("remotes")
remotes::install_github("fmichonneau/phyloch")
Tree<-phyloch::read.mrbayes("All_psbA-Cladocopium100-90.aln.mpileup.consensus_woN_SelectedJohnston2022_clustal.nexus1_simple.con.tre")
clade<-read.table("/env/cns/proj/TaraPacifique/ARTICLES/MetaTArticle/ITS2-analysis/Variables11Islands.txt",h=T,sep="\t")
clade<-clade[!(clade$SymbioGG%in%c("CD","D")),]
clade<-clade[!(clade$SymbioGG%in%c("CD","D")),]
Tree2 <- groupOTU(Tree, split(clade$Sample,clade$SymbioGG))

ggtree(Tree,size=0.5,aes(x,y)) + 
  geom_treescale()+
  geom_tiplab(size=2.5)+
  geom_text(aes(label=prob))+
  theme_tree2()+
  xlim(0,0.20)


