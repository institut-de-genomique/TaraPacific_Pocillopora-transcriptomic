#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
	--variables FILE samples in lignes (rownames) variable in columns (header)
	--HostExpression FILE Table normalized for the host
	--SymbiontExpression FILE Table normalized for the symbiont
	--CPU INT number of available CPU
	--OutSymbiont file name for the network graph
	--OutHost file name for the network graph
\n\n")
  q(save="no")
}
 

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
## mclresults
if(is.null(argsL$variables)) {writeLines("Variables not loaded See --help"); q(save="no")}
## hostexpr
if(is.null(argsL$HostExpression)) {writeLines("HostExpression not loaded See --help"); q(save="no")}
## symbiont expr
if(is.null(argsL$SymbiontExpression)) {writeLines("SymbiontExpression not loaded , See --help"); q(save="no")}
## OutExpressionBarplot
if(is.null(argsL$CPU)) {writeLines("OutExpressionBarplot not loaded , See --help"); q(save="no")}
## OutSymbiont
if(is.null(argsL$OutSymbiont)) {writeLines("Out not loaded , See --help"); q(save="no")}
## OutHost
if(is.null(argsL$OutHost)) {writeLines("Out not loaded , See --help"); q(save="no")}

writeLines("\n\n***Parameters***")
writeLines(paste(names(unlist(argsL)),as.data.frame(unlist(argsL))[,1],sep=":"))
writeLines("***\n\n")

##########################
###########CORE###########
##########################


#article
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1323-z
#r pipeline
#https://bioconductor.org/packages/release/bioc/html/variancePartition.html
#library(devtools)
#install_github("GabrielHoffman/variancePartition")
library(variancePartition)
library(reshape2)
param = SnowParam(as.numeric(argsL$CPU), "SOCK")
#Samples informations (Islands and genetic groups for the host and symbiodinium)
Genet<-read.table(argsL$variables,sep="\t",h=T)
rownames(Genet)<-Genet$Samples
#removing NAs 
#Genet<-Genet[rowSums(cbind(is.na(Genet$PocilloGG),is.na(Genet$SymbioGG)))==0,]
#Genet<-Genet[is.na(Genet$SymbioGG)==F,]


#
#Analysis of Cladocopium
#


Total<-data.frame()
#subsampling test
tab <- read.table(argsL$SymbiontExpression, header=TRUE,stringsAsFactors = F)
rownames(tab)<-tab$Gene
GenetC<-Genet[Genet$Islands!="I04",]
GenetC<-GenetC[is.na(GenetC$SymbioGG)==F,]
GenetC<-GenetC[GenetC$Samples%in%colnames(tab),]
#Subsamples<-"All"
Subsamples<-strsplit(unique(paste(as.vector(sample(GenetC$Samples)),as.vector(sample(GenetC$Samples))))[1:10]," ")
for (i in 1:5){
  print(paste("Symbiont Run ",i,"/100",sep=""))
#Samples with genetic group information
Subsample<-Subsamples[[i]]
GenetC2<-GenetC[!(GenetC$Samples%in%Subsample),]
datC<-tab[,colnames(tab)%in%GenetC2$Samples]
#Expressed in 3 samples at least : 23474 genes and 82 samples
datC<-datC[rowSums(datC>0)>=3,]
#VariancePartition
#removing batch effect first

residListC <- fitVarPartModel(datC, ~ (1|Batch), GenetC2,fxn=residuals,BPPARAM = param)
residMatrixC <- do.call(rbind, residListC)
# fit model on residuals
form <- ~ (1|Islands) + (1|PocilloGG) + (1|SymbioGG)
varPartC <- fitExtractVarPartModel(residMatrixC, form, GenetC2,BPPARAM = param)
rm(residListC)
rm(residMatrixC)

#pdf(file="VarPart_batchcorrected_3Samples_woI04_Symbio.pdf")
#plotVarPart(varPartC)
#dev.off()


df<-cbind(melt(as.matrix(varPartC)),SubSampling=paste(sort(Subsample),collapse="-"))
colnames(df)<-c("Gene","Variable","Variance","SubSampling")
df$Variance<-round(df$Variance,digits=4)
Total<-rbind(Total,df)
}
write.table(file=argsL$OutSymbiont,row.names=F,Total,sep="\t",quote=F)

####
#Pocillo
####
#
#Analysis of Pocillopora
#

tab<-read.table(file = argsL$HostExpression,h=T,sep="\t")
GenetP<-Genet[is.na(Genet$PocilloGG)==F,]
GenetP<-GenetP[GenetP$Samples%in%colnames(tab),]
Total<-data.frame()
#Subsamples<-"All"
Subsamples<-strsplit(unique(paste(as.vector(sample(GenetP$Samples)),as.vector(sample(GenetP$Samples))))[1:10]," ")
#subsampling test
for (i in 1:5){
  print(paste("Pocillopora Run ",i,"/100",sep=""))
  Subsample<-Subsamples[[i]]
  GenetP2<-GenetP[!(GenetP$Samples%in%Subsample),]
	GenetP2
  datP<-tab[,colnames(tab)%in%GenetP2$Samples]

#Expressed in 5 samples at least.
datP<-datP[rowSums(datP>0)>=3,]
#VariancePartition
#removing batch effect first

residListP <- fitVarPartModel(datP, ~ (1|Batch), GenetP2,fxn=residuals,BPPARAM=param)
residMatrixP <- do.call(rbind, residListP)

# fit model on residuals
form <- ~ (1|Islands) + (1|PocilloGG) + (1|SymbioGG)
varPartP <- fitExtractVarPartModel(residMatrixP, form, GenetP2,BPPARAM=param)
rm(residListP)
rm(residMatrixP)

#random subsampling
df<-cbind(melt(as.matrix(varPartP)),SubSampling=paste(sort(Subsample),collapse="-"))
colnames(df)<-c("Gene","Variable","Variance","SubSampling")
df$Variance<-round(df$Variance,digits=4)
Total<-rbind(Total,df)
}
write.table(file=argsL$OutHost,Total,row.names=F,sep="\t",quote=F)
