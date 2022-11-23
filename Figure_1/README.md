# Identification of _Cladocopium_ lineages 
We used three independant methods to identify _Cladocopium_ lineages.  

## 1. Identification of _Cladocopium_ lineages with ITS2 sequences.
ITS2 sequences are available under the PRJEB52458 project and were obtained following the procedure detailed here : https://arxiv.org/abs/2207.02475  
ITS2 sequences were processed with the Symportal pipeline and the results are available here : https://doi.org/10.5281/zenodo.4061796  
ITS2 profil abundances in each colony is in : TARA_PACIFIC_METAB_ITS2_coral_its2_type_profiles_absolute_abund_and_meta_v1.csv  

### A. Barplot of ITS2 profiles
The R code to produce Supplementary Fig. 2 is below  

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

#Supplementary Fig. 2
ggplot(tab3,aes(x=SiteColo,y=Value,fill=Type))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=set)+
  facet_grid(.~Island,scale="free",space="free")+  theme(axis.text.x=element_text(angle=90),panel.grid=element_blank(),axis.ticks.x=element_blank(),panel.background=element_blank(),axis.line=element_blank(),legend.key.size=unit(2,"mm"))+
  guides(fill=guide_legend(ncol=1))
```

### B. Hierarchical clustering of ITS2 profile Bray Curtis distances (Supplementary Fig. 3b)
ITS2 sequences were processed with the Symportal pipeline and the results are available here : https://doi.org/10.5281/zenodo.4061796  
ITS2 profil Braycurtis distances is available on Zenodo (https://doi.org/10.5281/zenodo.4061796) : TARA_PACIFIC_METAB_ITS2_coral_between_sample_distances_C_braycurtis_distances_sqrt_v1.dist

```r
#File of metadata
Variable<-read.table("Variables11Islands.txt",sep="\t",h=T)

colors<-brewer.pal(8,"Set1")
pal<-colorRampPalette(colors)
set<-sample(pal(10),replace = F)

#BrayCurtis distance table
tab<-read.table("TARA_PACIFIC_METAB_ITS2_coral_between_sample_distances_C_braycurtis_distances_sqrt_v1.dist",sep="\t",h=F)

#Selection of _Pocillopora_ samples from the 11 first islands.
tab2<-tab[tab$V1%in%Variable$TaraSampleName,c(TRUE,TRUE,tab$V1%in%Variable$TaraSampleName)]
colnames(tab2)<-c("TaraSampleName","Id",tab2$V1)

#Addition of SampleNames to tab2
tab2<-merge(Variable[,c("Samples","TaraSampleName")],tab2,by.x="TaraSampleName",by.y="TaraSampleName")
rownames(tab2)<-tab2$TaraSampleName

#Conversion in matrix 
mat<-as.matrix(tab2[,4:ncol(tab2)])
mat<-mat[,rownames(mat)]

#Addition of Sample names to mat
rownames(mat)<-tab2$Samples
colnames(mat)<-tab2$Samples

#Elimination of samples containing 2 _Cladocopium_ ITS2 profiles.
GoodSamples<-Variable[!Variable$Sample%in%c("I10S01C006POC","I05S02C002POC"),1]
mat2<-mat[rownames(mat)%in%GoodSamples,rownames(mat)%in%GoodSamples]

#Hierarchical clustering of Samples ITS2 distances.
dendITS<-hclust(as.dist(mat2),method = "average") 
dendITS2 <- as.dendrogram(dendITS) %>%
  set("labels_cex", c(.5)) %>%
  set("labels_cex", c(.5)) %>% 
  set("branches_lwd",c(0.5)) %>%
  set("labels_col",c("#117733", "#cc6677" ,"#44aa99", "#332288","#6699cc")[as.factor(unlist(lapply(dend$labels[dend$order],function(x){Variable[Variable$Samples==x,"SymbioGG"]})))]) 

#Conversion into ggplot tree
ggdendITS2 <- as.ggdend(dendITS2)

#Supplementary Fig. 3b
ggplot(ggdendITS2, horiz = TRUE, theme=NULL) + 
  theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.x=element_line(color="black"),panel.background = element_blank())+
  labs(y="",x="ITS2 dendrogram")
 
#Write dendrogram in newick/nexus format
library(ctc)
write(file="Cladocopium_ITS2.dendrogram.newick",hc2Newick(dendITS,flat=T))
library(ape)
write.nexus(as.phylo(dendITS), file="Cladocopium_ITS2.dendrogram.nex")

```

## 2. Identification of _Cladocopium_ lineages with SNPs called on transcriptomic reads.  
### A. SNP calling

### B. Hierarchical clustering of SNP frequencies (Supplementary Fig. 3a)

```r
library(dendextend)

#SNP frequencies loading
tab<-read.table("Cladocopium-selected-SNPs.tab",h=T,check.names = F)
colnames(tab)<-sub("_.*","",colnames(tab))

#Selection of 82 samples containing a single Cladocopium species.
tab2<-tab[,colnames(tab)%in%Variable[Variable$Symbio_Genus=="C","Samples"]]

#Hierarchical clustering
dendSNP<-hclust(dist(t(tab2)),method = "complete")

#Write dendrogram in newick/nexus format
library(ctc)
library(ape)
write(file="Cladocopium_SNP.dendrogram.newick",hc2Newick(dendSNP,flat=T))
write.nexus(as.phylo(dendSNP), file="Cladocopium_SNP.dendrogram.nex")

#Calculation of the optimal number of cluster
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

#Addition of tree features
dendSNP2 <- as.dendrogram(dendSNP) %>%
  set("labels_cex", c(.5)) %>%
  set("labels_cex", c(.5)) %>% 
  set("branches_lwd",c(0.5)) %>%
  set("labels_col",c("#117733", "#cc6677" ,"#44aa99", "#332288","#6699cc")[as.factor(unlist(lapply(dendSNP$labels[dendSNP$order],function(x){Variable[Variable$Samples==x,"SymbioGG"]})))]) 

#Conversion in ggplot format
ggdendSNP2 <- as.ggdend(dendSNP2)

#Supplementary Fig. 3a
ggplot(ggdendSNP2, horiz = TRUE, theme=NULL) + 
  theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.x=element_line(color="black"),panel.background = element_blank())+
  labs(y="",x="SNPs dendrogram")


#Ordered dendrogram according to Clado Lineages
dendITS4 <-rotate(dendITS2,Variable[order(Variable$SymbioGG),"Samples"][Variable[order(Variable$SymbioGG),"Samples"]%in%labels(dendITS2)])
dendSNP4 <-rotate(dendSNP2,Variable[order(Variable$SymbioGG),"Samples"][Variable[order(Variable$SymbioGG),"Samples"]%in%labels(dendSNP2)])
  
dendCombined<-dendlist(dendITS4, dendSNP4)
#Supplementary Fig. 3a and b
dend_diff(dendCombined)
```

## 3. Analysis of psbA<sup>ncr</sup> sequences in _Pocillopora_ samples  (Supplementary Fig. 3c)

We analyzed the non-coding region of the psbA gene (psbA<sup>ncr</sup>) to assign _Cladocopium_ lineages to known _Cladocopium_ species following the method used in Johnston et al 2022 (https://pubmed.ncbi.nlm.nih.gov/35960256/) and Turnham et al 2021 (https://www.nature.com/articles/s41396-021-01007-8). The code below describes : (A) The mapping of metagenomic reads of the 82 _Pocillopora_ corals containing Cladocopium symbionts on 3 psbAncr sequences (C. latusorum MW819767.1, C. pacificum MW861717, and C. goreaui KF572161.1) the build of consensus sequences for each sample. (B) The alignement then construction of Bayesian phylogeny using 2 psbA<sup>ncr</sup> sequences for each Cladocopium clade identified in Johnston et al 2022 and and the 82 consensus sequences. (C) The representation of the phylogeny with R (Supplementary Fig. 3c).

### A. Identification of psbA<sup>ncr</sup> sequences in each _Pocillopora_ colonies  

Alignement of meatgenomic reads from 82 _Pocillopra_ samples on 3 psbA<sup>ncr</sup> sequences 
Metagenomic fastq files are available at the ENA under project PRJEB47249  

Required tools: bwa-mem2/2.2.1, samtools/1.15.1  

The following command should be executed on each metagenomic readset.  
```bash
bwa-mem2 index psbA_Cladocopium.fa
bwa-mem2 mem -t 4 -M psbA_Cladocopium.fa read1.fastq read2.fastq | samtools view -b -@ 4 -F 4 - -o SAMPLE_PsBA.aln.bam"
perl /env/cns/proj/projet_CNM/script/Quentin/BamFiltration.pl -in SAMPLE_PsBA.aln.bam -out SAMPLE_PsBA.aln.filtered100-90.bam -minPCaligned 100 -minPCidentity 90
```
To count the number of reads aligned on the most covered psbA<sup>ncr</sup> sequence:  
```bash
for i in *_PsBA.aln.filtered100-90.bam ;do echo -ne "$i\t" ;samtools view $i | awk '{print $3}' | sort | uniq -c | sort -nrk1,1 |head -1 | awk '{print $2"\t"$1}' ;done > Summary_mapping_Cladocopium100-90.tab
```

Sort and index of all bam files
```bash
cat Summary_mapping_Cladocopium100-90.tab | while read a b c ;do samtools sort -o ${a%.bam}.sort.bam $a; done
for i in *_PsBA.aln.filtered100-90.sort.bam ;do samtools index $i; done
```

Build of consensus sequences for each sample
```bash
cat Bilan_mapping_Cladocopium100-90.tab | while read a b c ;do samtools mpileup -Aa -f psbA_Cladocopium.fa -r $b ${a%.bam}.sort.bam -o ${a%.bam}.sort.mpileup;done
for i in *PsBA.aln.filtered100-90.sort.mpileup ;do perl /env/cns/proj/TaraPacifique/bin/Quentin/MpileuptoConsensus.pl -mpileup $i -Mincov 4 -out ${i}.consensus.fa;done
for i in *_PsBA.aln.filtered100-90.sort.mpileup.consensus.fa;do echo -e ">${i%_*}"; tail -n +2 $i ;done > All_psbA-Cladocopium100-90.aln.mpileup.consensus.fa
```
### B. Multiple alignement and Bayesian phylogeny of psbA<sup>ncr</sup> sequences

The multiple alignement is done with MegaX software (https://www.megasoftware.net/) using ClustalW method.  
The result is a nexus alignment file All_psbA.aln.mpileup.consensus_woN_AllJohnston2022_muscle.nexus.
The bayesian phylogeny is computed with MrByes (http://nbisweden.github.io/MrBayes/)

Config file : MrbayesCommands.nex  
```bash
begin mrbayes;
   set autoclose=yes nowarn=yes;
   execute psbA-Cladocopium_Johnston2022.clustal.nexus;
   lset nst=6 rates=invgamma;
   prset brlenspr=clock:birthdeath;
   mcmc nruns=1 ngen=1000000 samplefreq=100 printfreq=100 diagnfreq=1000 file=psbA-Cladocopium_Johnston2022.clustal.phylogeny.nexus;
   sump;
   sumt;
end;
```
```bash
./mb -i MrbayesCommands.nex
```

### C. Tree representation on R (Supplementary Fig. 3c)

```r
library(ape)
library(ggtree)

#Phylogeny loading
Tree<-ape::read.nexus("psbA-Cladocopium_Johnston2022.clustal.phylogeny.nexus")
#Metadata on _Cladocopium_ lineages
clade<-read.table("Variables11Islands.txt",h=T,sep="\t")

#Elimination of samples containing _Durusdinium_
clade<-clade[!(clade$SymbioGG%in%c("CD","D")),]

#To color leaves according to _Cladocopium_ lineages
Tree2 <- groupOTU(Tree$con_50_majrule, split(clade$Sample,clade$SymbioGG))

#Supplementary Fig. 3c
ggtree(Tree2,size=0.5,aes(x,y)) + 
  geom_treescale()+
  geom_tiplab(size=2,aes(color=group))+
  scale_color_manual(values=c("red","#117733", "#cc6677" ,"#44aa99", "#332288","#6699cc"))+
  geom_nodelab(size=2,aes(label=as.numeric(label)*100,x=x+0.20/100))+
  theme_tree2()+
  xlim(0,0.20)
dev.off()
```
