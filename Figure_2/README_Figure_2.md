# Identifications of Cladocopium lineages 
We used three independant methodes to identify Cladocopium lineages.  

## 1. Using ITS2 sequences  

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
## 2. Using psbA<sup>ncr</sup> sequences  

## 3. Using SNVs identified on transcriptomic reads.  
