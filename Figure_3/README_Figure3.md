# Variance partitioning

The R script execute variance partition 5 times with a random subsampling of 2 samples.
To get 100 different subsampling, the following command was performed 22 times.

``for i in seq `1 22`;do Rscript --vanilla VariancePartition_subsampling.R --variables=Variables11Islands.txt --HostExpression=Pocillopora_MetaT_TPM.tab --OutHost=Varpart_Host_rep10_${i}.tab --SymbiontExpression=CladocopiumC1_MetaT_TPM.tab --OutSymbiont=Varpart_Symbiont_rep10_${i}.tab --CPU=6;sleep 30;done``

Concatenation of output files

`cat Varpart_Symbiont_rep10_* | awk '$4>0&&$1!="Gene"{print $0}' > Varpart_Symbiont-combined.tab
cat Varpart_Host_rep10_* | awk '$4>0&&$1!="Gene"{print $0}' > Varpart_Host-combined.tab`


Conversion of Varpart_Symbiont-combined.tab file into table with R (lot of memory required)
Similar commands for the host.

`library(data.table)
library(reshape2)
tab<-fread("Varpart_Symbiont-combined.tab",h=F,sep="\t")
tab2<-acast(tab,V1+V2~V4,fill=0,fun.aggregate=mean,value.var="V3")
tab3<-tab2[grep("Residuals",rownames(tab2),invert=T),]
write.table(tab3,file="Varpart_Symbiont-combined_formated.tab",sep="\t",quote=F)`
