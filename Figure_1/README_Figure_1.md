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
### A. Identification of psbA<sup>ncr</sup> sequences in each _Pocillopora_ colonies

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


