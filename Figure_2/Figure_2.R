#Symbio ITS2/SNPs
#Comparison between ITS2 profiles and SNP clustering
#-->SNP clustering = Julie Lê-Hoang's cahier

##ITS2 analysis

#Analysis of ITS2 from Ben Hume Symportal results for the first 11 islands.
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colors<-brewer.pal(8,"Set1")
pal<-colorRampPalette(colors)
set<-sample(pal(13),replace = F)

#ITS2.txt table is converted from TARA_PACIFIC_METAB_ITS2_coral_its2_type_profiles_absolute_abund_and_meta_v1.csv
#On google drive (from Ben Hume) but soon accessible on zenodo: https://docs.google.com/document/d/1mDSttur-4suuO-p5etcRSt1YK4FibeaH/edit#bookmark=id.gjdgxs
tab<-read.table(file = "ITS2.txt",sep="\t",h=T,check.names = F)
tab2<-melt(tab,id.vars = colnames(tab[,1:4]), measure.vars = colnames(tab[,5:ncol(tab)]))
colnames(tab2)<-c("ProfileNumber","SampleName","Taxon","RefCollab","ProfileName","Value")
tab3<-tab2[tab2$Value>0,2:ncol(tab2)]
tab3[order(tab3$SampleName),][1:10,]
SampleList<-read.table("SampleList.txt",h=F,sep="\t")
tab3$Experiment<-apply(tab3,1,function(x){if(length(SampleList[SampleList$V1==x[1],])>0){"METAT"}else{"Other"}})
tab3$Clade<-substr(tab3$ProfileName,1,1)
tab3$Island<-substr(tab3$SampleName,1,3)
tab3$Type<-sub('-.*','',tab3$ProfileName,perl = T)
tab3$Type2<-sub('-.*|\\/.*','',tab3$ProfileName,perl = T)
tab3$SiteColo<-substr(tab3$SampleName,5,12)
tab4<-tab3[tab3$Experiment=="METAT"&tab3$Taxon=="Pocillopora",]
tab4$SampleName<-sub(".*?_","",sub("_OA.*","",tab4$SampleName))
tab4$Island<-factor(tab4$Island,levels=sort(unique(tab4$Island),decreasing = T))

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
#Supplementary figure 2
pdf(file="ITS2_type2_Pocillo.pdf",width=12)
ggplot(tab4)+
  geom_bar(aes(x=SampleName,y=Value,fill=Type),stat="identity",position="fill")+
  scale_fill_manual(values=set)+
  facet_grid(.~Island,scale="free",space="free")+
  theme(axis.text.x=element_text(angle=90),panel.grid=element_blank(),axis.ticks.x=element_blank(),panel.background=element_blank(),axis.line=element_blank(),legend.key.size=unit(2,"mm"))+guides(fill=guide_legend(ncol=1))
dev.off()

#Idem for all colonies in each island separately (not published)
pdf(file="ITS2_All_poc_colonies.pdf",width=15)
for (i in sort(unique(tab3$Island))){
  tab4<-tab3[tab3$Taxon=="Pocillopora"&tab3$Island==i,]
  if(nrow(tab4)>0){
    plot(ggplot(tab4)+geom_bar(aes(x=SampleName,y=Value,fill=ProfileName),stat="identity",position="fill")+facet_grid(.~Island,scale="free",space="free")+theme(axis.text.x=element_text(angle=90),legend.key.size=unit(2,"mm"))+guides(fill=guide_legend(ncol=1)))
  }
}
dev.off()

#####################
#Comparison with SNP#
#####################
#Genetic clades for the host and the symbiont based on the SNP.
CladoLineage<-read.table("Variables11Islands.txt",sep="\t",h=T)
head(tab3)
tab3$Sample<-paste(gsub("_","",sub("_OA.*","",tab3$SampleName),),"POC",sep="")
tab3$Island<-factor(tab3$Island,levels=sort(unique(tab4$Island),decreasing = F))
Merging<-merge(CladoLineage,tab3[,c(5,8,9,12)],by="Sample",all.X=T)

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
#Elimination of I64 ??
list<-unique(tab4[tab4$Island!="I64",3])
tab2<-tab[tab$V1%in%list,c(TRUE,TRUE,tab$V1%in%list)]
#Conversion in matrix 
mat<-as.matrix(tab2[,3:ncol(tab2)])
rownames(mat)<-tab2$V1
#2 dimensional reduction analyses: nmds and PCoA. 
nmds<-metaMDS(as.dist(mat))
PCOA<-pcoa(mat)
#Inclusion of results in a dataframe.
dat<-data.frame(RefCollab=rownames(nmds$points),nmds$points[,1:2],PCOA$vectors[,1:4])
dat$Island<-apply(dat,1,function(x){tab4[tab4$RefCollab==x[1],8][1]})
dat$Type<-apply(dat,1,function(x){tab4[tab4$RefCollab==x[1],10][1]})
dat$SiteColo<-apply(dat,1,function(x){tab4[tab4$RefCollab==x[1],11][1]})
dat$clustkmeans<-kmeans(as.dist(mat),5)$cluster
dat$clusthclust<-cutree(hclust(as.dist(mat)),5)
dat$Sample<-sub("_","",paste(dat$Island,dat$SiteColo,"POC",sep=""))

#PCoA or nmds representation of unifrac distances.
pdf(file="ITS2_profiles_unifrac_distance_sqrt_pcoa_k5.pdf")
ggplot(dat)+geom_text(aes(x=MDS1,y=MDS2,color=Island,label=SiteColo),size=3)+scale_color_manual(values=pal(length(unique(dat$Island))))+theme(panel.background = element_rect(fill="transparent",color="black"),panel.grid=element_blank())+labs(title="Pocillo_ITS2Profiles_BrayCurtis_nmds")
ggplot(dat)+geom_text(aes(x=MDS1,y=MDS2,color=factor(clusthclust),label=SiteColo),size=3)+scale_color_manual(values=pal(length(unique(dat$Island))))+theme(panel.background = element_rect(fill="transparent",color="black"),panel.grid=element_blank())+labs(title="Pocillo_ITS2Profiles_BrayCurtis_nmds_hclust")
ggplot(dat)+geom_text(aes(x=MDS1,y=MDS2,color=factor(clustkmeans),label=SiteColo),size=3)+scale_color_manual(values=pal(length(unique(dat$Island))))+theme(panel.background = element_rect(fill="transparent",color="black"),panel.grid=element_blank())+labs(title="Pocillo_ITS2Profiles_BrayCurtis_nmds_kmeans")
ggplot(dat)+geom_text(aes(x=Axis.1,y=Axis.2,color=Island,label=Island),size=3)+scale_color_manual(values=pal(length(unique(dat$Island))))+theme(panel.background = element_rect(fill="transparent",color="black"),panel.grid=element_blank())+labs(title="Pocillo_ITS2Profiles_unifrac_pcoa")
ggplot(dat)+geom_text(aes(x=Axis.1,y=Axis.2,color=factor(clusthclust),label=Island),size=3)+scale_color_manual(values=pal(length(unique(dat$Island))))+theme(panel.background = element_rect(fill="transparent",color="black"),panel.grid=element_blank())+labs(title="Pocillo_ITS2Profiles_unifrac_pcoa_hclust")
ggplot(dat)+geom_text(aes(x=Axis.1,y=Axis.2,color=factor(clustkmeans),label=Island),size=3)+scale_color_manual(values=pal(length(unique(dat$Island))))+theme(panel.background = element_rect(fill="transparent",color="black"),panel.grid=element_blank())+labs(title="Pocillo_ITS2Profiles_unifrac_pcoa_kmeans")
ggplot(dat)+geom_text(aes(x=Axis.1,y=Axis.2,color=factor(clustkmeans),label=Island),size=3)+scale_color_manual(values=pal(length(unique(dat$Island))))+theme(panel.background = element_rect(fill="transparent",color="black"),panel.grid=element_blank())+labs(title="Pocillo_ITS2Profiles_unifrac_pcoa_kmeans")
dev.off()

#merging of ITS2 profile distances and SNP cluster
library(RColorBrewer)
library(rcartocolor) 
mycarto_cols <- carto_pal(12, "Safe")
mycols3 <- c(mycarto_cols[1], mycarto_cols[4],
             mycarto_cols[2], mycarto_cols[7], mycarto_cols[5],
             mycarto_cols[11], "white")

dat2<-merge(dat,CladoLineage,by="Sample",all=T)
dat2$SymbioGG<-factor(dat2$SymbioGG,levels=c("CD","D","Group1","Group2","Group3","Group4","Group5",NA))
#Supplementary Figure 4
pdf(file="PcoAITS2-ColorClade.pdf",width=8)
ggplot(dat2[dat2$SymbioGG!="D",])+
  geom_point(aes(x=Axis.1,y=Axis.2,color=SymbioGG),size=2)+
  geom_text_repel(aes(x=Axis.1,y=Axis.2,label=Island),size=3, max.overlaps=100)+
  scale_color_manual(values=mycols3)+
  theme(panel.background = element_rect(fill="transparent",color="black"),panel.grid=element_blank())+labs(title="Pocillo_ITS2Profiles_unifrac_pcoa",x=paste("Axis 1, ",round(PCOA$values$Relative_eig[1]*100),"%",sep=""),y=paste("Axis 2, ",round(PCOA$values$Relative_eig[2]*100),"%",sep=""))
dev.off()

#Hierarchical clustering of Samples ITS2 distances.
library(dendextend)
temp<-dat2
mat2<-mat
temp$RefCollab<-factor(temp$RefCollab,levels=unique(rownames(mat2)))
temp<-temp%>%arrange(RefCollab)
rownames(mat2)<-temp[is.na(temp$RefCollab)==F,1]
#Elimination of samples containing 2 Symbiodinium clades.
GoodSamples<-CladoLineage[CladoLineage$Symbio_Genus=="C"&!CladoLineage$Sample%in%c("I10S01C006POC","I05S02C002POC"),1]
mat2<-mat2[rownames(mat2)%in%GoodSamples,rownames(mat2)%in%GoodSamples]
hclu<-hclust(as.dist(mat2))


pdf(file="Hclus_ITS2-distance.pdf",width=17,height = 11)
plot(hclu)
dev.off()
