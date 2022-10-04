# Variance partitioning

To obtain Figure 3 of the article required files are :  
-the associations between host, symbiont and island : Variables11Islands.txt in this directory
-normalized gene expression for the host: Pocillopora_MetaT_TPM.tab available at https://doi.org/10.5281/zenodo.6341761 
-normalized gene expression for the symbiont: CladocopiumC1_MetaT_TPM.tab available at https://doi.org/10.5281/zenodo.6341761

The R script VariancePartition_subsampling.R execute variance partition 5 times with a random subsampling of 2 samples.  
To get 100 different subsampling, the following command was performed 22 times.

```bash
for i in seq `1 22`;do Rscript --vanilla VariancePartition_subsampling.R --variables=Variables11Islands.txt --HostExpression=Pocillopora_MetaT_TPM.tab --OutHost=Varpart_Host_rep10_${i}.tab --SymbiontExpression=CladocopiumC1_MetaT_TPM.tab --OutSymbiont=Varpart_Symbiont_rep10_${i}.tab --CPU=6;sleep 30;done
```

Concatenation of output files

```bash
cat Varpart_Symbiont_rep10_* | awk '$4>0&&$1!="Gene"{print $0}' > Varpart_Symbiont-combined.tab
cat Varpart_Host_rep10_* | awk '$4>0&&$1!="Gene"{print $0}' > Varpart_Host-combined.tab
```


Conversion of Varpart_Symbiont-combined.tab file into table with R (lot of memory required).   
Similar commands for the host.

```r
library(data.table)
library(reshape2)
tab<-fread("Varpart_Symbiont-combined.tab",h=F,sep="\t")
tab2<-acast(tab,V1+V2~V4,fill=0,fun.aggregate=mean,value.var="V3")
tab3<-tab2[grep("Residuals",rownames(tab2),invert=T),]
write.table(tab3,file="Varpart_Symbiont-combined_formated.tab",sep="\t",quote=F)
```
Similar commands for the host.

```r
library(data.table)
library(reshape2)
tab<-fread("Varpart_Host-combined.tab",h=F,sep="\t")
tab2<-acast(tab,V1+V2~V4,fill=0,fun.aggregate=mean,value.var="V3")
tab3<-tab2[grep("Residuals",rownames(tab2),invert=T),]
write.table(tab3,file="Varpart_Host-combined_formated.tab",sep="\t",quote=F)
```

Generation of Figure 3 for the host and the Symbiont using Varpart_Host-combined_formated.tab file.

```r
Genet<-read.table("Variables11Islands.txt",sep="\t",h=T)
rownames(Genet)<-Genet$Samples
library(RColorBrewer)
set<-brewer.pal(8,"Set2")

#For the Symbiont.
library(ggplot2)
Varpart<-read.table(file="Varpart_Symbiont-combined_formated.tab",sep="\t",h=T)
Varpart2<-Varpart[,2:101]
Varpart3<-data.frame(Gene=sub("_(Islands|PocilloGG|SymbioGG)","",rownames(Varpart2)),Variable=sub(".*_","",rownames(Varpart2)),VarianceMedian=apply(Varpart2,1,median),sd=apply(Varpart2,1,sd),Q1=apply(Varpart2,1,function(x){quantile(x,probs = 0.25)}),Q3=apply(Varpart2,1,function(x){quantile(x,probs = 0.75)}),row.names = NULL)
Varpart3<-Varpart3[order(Varpart3$VarianceMedian,decreasing = T),]
Varpart3$Gene<-factor(Varpart3$Gene,levels=unique(Varpart3$Gene))

head(Varpart3)

pdf(file="VariancePartitionSymbiont_batchcorrected_I04corrected.pdf",width=9)
ggplot(Varpart3[Varpart3$VarianceMedian>=0.5,])+ 
  geom_violin(scale = "count",size=0.1,aes(x=Variable,y=VarianceMedian*100,fill=Variable))+ 
  geom_boxplot(width=0.05,fill="grey60",size=0.1,aes(x=Variable,y=VarianceMedian*100))+
  scale_fill_manual(values=set[1:3])+
  scale_x_discrete(expand=c(0,0))+
  geom_text(data=data.frame(table(Varpart3[Varpart3$VarianceMedian>=0.5,2])),aes(x=Var1,y=105,label=paste(Freq,"genes")))+
  theme_bw()+ 
  theme(panel.grid.major.x=element_blank(),legend.position="none")
dev.off()

#For the Host.
library(ggplot2)
Varpart<-read.table(file="Varpart_Host-combined_formated.tab",sep="\t",h=T)
Varpart2<-Varpart[,2:101]
Varpart3<-data.frame(Gene=sub("_(Islands|PocilloGG|SymbioGG)","",rownames(Varpart2)),Variable=sub(".*_","",rownames(Varpart2)),VarianceMedian=apply(Varpart2,1,median),sd=apply(Varpart2,1,sd),Q1=apply(Varpart2,1,function(x){quantile(x,probs = 0.25)}),Q3=apply(Varpart2,1,function(x){quantile(x,probs = 0.75)}),row.names = NULL)
Varpart3<-Varpart3[order(Varpart3$VarianceMedian,decreasing = T),]
Varpart3$Gene<-factor(Varpart3$Gene,levels=unique(Varpart3$Gene))

pdf(file="VariancePartitionHost_batchcorrected_I04corrected.pdf",width=9)
ggplot(Varpart3[Varpart3$VarianceMedian>=0.5,])+ 
  geom_violin(scale = "count",size=0.1,aes(x=Variable,y=VarianceMedian*100,fill=Variable))+ 
  geom_boxplot(width=0.05,fill="grey60",size=0.1,aes(x=Variable,y=VarianceMedian*100))+
  scale_fill_manual(values=set[1:3])+
  scale_x_discrete(expand=c(0,0))+
  geom_text(data=data.frame(table(Varpart3[Varpart3$VarianceMedian>=0.5,2])),aes(x=Var1,y=105,label=paste(Freq,"genes")))+
  theme_bw()+ 
  theme(panel.grid.major.x=element_blank(),legend.position="none")
dev.off()

```
