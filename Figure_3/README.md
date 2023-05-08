# Variance partitioning

To obtain Figure 3 of the article required files are :  
-the associations between host, symbiont and island : Variables11Islands.txt in this directory <br>
-normalized gene expression for the host: Pocillopora_MetaT_TPM.tab available at https://doi.org/10.5281/zenodo.6341761 <br>
-normalized gene expression for the symbiont: CladocopiumC1_MetaT_TPM.tab available at https://doi.org/10.5281/zenodo.6341761 <br>
-the summary table of top variance genes for the host: SummaryTable_upsetR_Host.tab in this directory <br>
-the summary table of top variance geens for the symbiont: SummaryTable_upsetR_Symbiont.tab in this directory <br>
-the gene lengths for the host: Pocillopora_meandrina_v3.1.annot.mrna_SeqLengths.txt in this directory <br>
-the read counts for the host: Pocillopora_MetaT_ReadCount.tab in this directory in this directory <br>
-the functional annotations for the host: Pocillopora_meandrina_v3.2_annot.xlsx in the Figure_4 directory <br>

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
Varpart<-read.table(file="Varpart_Symbiont-combined_formated.tab",sep="\t",h=T)
Varpart2<-Varpart[,2:101]
Varpart3<-data.frame(Gene=sub("_(Islands|PocilloGG|SymbioGG)","",rownames(Varpart2)),Variable=sub(".*_","",rownames(Varpart2)),VarianceMedian=apply(Varpart2,1,median),sd=apply(Varpart2,1,sd),Q1=apply(Varpart2,1,function(x){quantile(x,probs = 0.25)}),Q3=apply(Varpart2,1,function(x){quantile(x,probs = 0.75)}),row.names = NULL)
Varpart3<-Varpart3[order(Varpart3$VarianceMedian,decreasing = T),]
Varpart3$Gene<-factor(Varpart3$Gene,levels=unique(Varpart3$Gene))

#Addition of residuals
library(data.table)
Varpart4<-data.frame(dcast(setDT(Varpart3),Gene~Variable,fill = 0,value.var = c("VarianceMedian","sd","Q1","Q3")))
Varpart4$Residuals<-1-rowSums(Varpart4[,grep("VarianceMedian",colnames(Varpart4))])
#Elimination of genes <0.5
Varpart5<-Varpart4[apply(Varpart4[,grep("VarianceMedian",colnames(Varpart4))],1,max)>0.5,]
write.table(Varpart5,sep="\t",quote=F,row.names = F,file="VarpartSubsampled_Symbiont_selected.tsv")

#Figure 3a
library(ggplot2)
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

#Addition of residuals
library(data.table)
Varpart4<-data.frame(dcast(setDT(Varpart3),Gene~Variable,fill = 0,value.var = c("VarianceMedian","sd","Q1","Q3")))
Varpart4$Residuals<-1-rowSums(Varpart4[,grep("VarianceMedian",colnames(Varpart4))])
#Elimination of genes <0.5
Varpart5<-Varpart4[apply(Varpart4[,grep("VarianceMedian",colnames(Varpart4))],1,max)>0.5,]
write.table(Varpart5,sep="\t",quote=F,row.names = F,file="VarpartSubsampled_Host_selected.tsv")

#Figure 3b
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


# Figure 3c
## 1. INITIALIZATION ------------------------------------------------------------------------------------------------------------------------------------

  # Set the prefix for each output file name
    outputPrefix <- "TaraPacific_Pocillopora_VarPart_GitHub_2022-11"

    
### 2. Load Data and Dependencies -----------------------------------------------------------------------------------------------------------------------------
    
    varPart_data <- read.table("SummaryTable_upsetR_Host.tab", header =T)
    varPart_data_sym <- read.table("SummaryTable_upsetR_Symbiont.tab", header =T)
    
    topVarGenes_Ile <- subset(varPart_data, Variable %in% c("Islands") & VarianceMedian >= 0.5)
    topVarGenes_Clade <- subset(varPart_data, Variable %in% c("PocilloGG") & VarianceMedian >= 0.5)
    topVarGenes_SymClade <- subset(varPart_data, Variable %in% c("SymbioGG") & VarianceMedian >= 0.5)
    
    # topVarGenes_Ile_sym <- subset(varPart_data_sym, Variable %in% c("Islands") & VarianceMedian >= 0.5)
    # topVarGenes_Clade_sym <- subset(varPart_data_sym, Variable %in% c("PocilloGG") & VarianceMedian >= 0.5)
    # topVarGenes_SymClade_sym <- subset(varPart_data_sym, Variable %in% c("SymbioGG") & VarianceMedian >= 0.5)

  
##### III - GO Enrichment Analysis #####
########################################
    # source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R");
    # installAnRichment() 
    #       
    # options(stringsAsFactors = FALSE)
    # library("anRichment")
    library("goseq")
    #BiocManager::install("GO.db")
    library("GO.db")
    library("geneLenDataBase")
    #BiocManager::install("qvalue")
    library("qvalue")
    
    
    # Step 1- Get gene lengths
      gene_lengths <- read.table("Pocillopora_meandrina_v3.1.annot.mrna_SeqLengths.txt", header=T, row.names=1, com='')
      gene_lengths <- as.matrix(gene_lengths[,1,drop=F])
      
    # Step 2 - Get host background gene list (All genes)
      cts_background <- read.table("Pocillopora_MetaT_ReadCount.tab", header=TRUE, com='', row.names=1, check.names=FALSE)
      Allbackground.gene_ids <- c(gsub(";.*","",rownames(cts_background)))
      
    # Step 3 - Parse GO assignments
      library(openxlsx)
      GO_dat <- read.xlsx("Pocillopora_meandrina_v3.2_annot.xlsx", sheet = "Poc_v3.2_annot.ipr", rowNames = F, colNames = F)
      
      library(plyr)
      GO_dat <- ddply(GO_dat, .(X1), summarise, 
                      paste0(unique(unlist(strsplit(X14, split=","))), collapse=","))
      colnames(GO_dat) <- c("GeneID", "GOID")
      GO_dat$GOID <- gsub("NA,","",GO_dat$GOID)
      GO_dat$GOID <- gsub(",NA","",GO_dat$GOID)
      GO_dat$GOID <- gsub("NA","",GO_dat$GOID)
      GO_info <- data.frame("X2"=GO_dat$GOID)
      rownames(GO_info) <- GO_dat$GeneID
      
      GO_info_listed <- apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
      names(GO_info_listed) <- rownames(GO_info)
      
      get_GO_term_descr =  function(x) {
        d = 'none';
        go_info = GOTERM[[x]];
        if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
        return(d);}
    
    
    topVarGenes <- rbind(topVarGenes_Ile,
                         topVarGenes_Clade,
                         topVarGenes_SymClade)
    
    module <- unique(topVarGenes$Variable)
    
    for (j in 1:length(module)) {
      # GoSeq Enrichment Analysis
      infile <- paste("Pocillopora_varPart_",module[j],"_",sep='')
      factor_labeling <- subset(topVarGenes, Variable %in% module[j])
      rownames(factor_labeling) <- gsub(";.*","", factor_labeling$Gene)
      factor_labeling[,1] <- rep('custom_list', dim(factor_labeling)[1])
      factor_labeling <- factor_labeling[,1,drop=F]
      colnames(factor_labeling) <- c('type')
      factor_list <- unique(factor_labeling[,1])
      select <- grep("Pmea_",rownames(factor_labeling))
      DE_genes <- rownames(factor_labeling)[select]
      #NOTE: Background genes MUST also include DEGs!
      background.gene_ids <- unique(c(Allbackground.gene_ids, DE_genes))
      sample_set_gene_ids = background.gene_ids
      # Step 5 - Organize go_id -> list of genes
      GO_to_gene_list = list()
      for (gene_id in intersect(names(GO_info_listed), sample_set_gene_ids)) {
        go_list = GO_info_listed[[gene_id]]
        for (go_id in go_list) {
          GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)
        }}
      # Step 6 - GO-Seq protocol: build pwf based on ALL DE features
      sample_set_gene_lengths <- gene_lengths[sample_set_gene_ids,]
      GO_info_listed <- GO_info_listed[ names(GO_info_listed) %in% sample_set_gene_ids ]
      cat_genes_vec <- as.integer(sample_set_gene_ids %in% rownames(factor_labeling))
      pwf <- nullp(cat_genes_vec, bias.data=sample_set_gene_lengths)
      rownames(pwf) <- sample_set_gene_ids
      # perform functional enrichment testing for each category.
      for (feature_cat in factor_list) {
        message('Processing category: ', feature_cat)
        gene_ids_in_feature_cat = rownames(factor_labeling)[factor_labeling$type == feature_cat]
        cat_genes_vec = as.integer(sample_set_gene_ids %in% gene_ids_in_feature_cat)
        pwf$DEgenes = cat_genes_vec
        res = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=TRUE)
        ## over-represented categories:
        pvals = res$over_represented_pvalue
        pvals[pvals > 1 - 1e-10] = 1 - 1e-10
        q = qvalue(pvals, lambda=0)
        res$over_represented_FDR = q$qvalues
        go_enrich_filename = paste("./varPart_GOSeq/",infile,"GOseq.enriched.txt", sep='')
        result_table = res[res$over_represented_pvalue<=0.05,]
        descr = unlist(lapply(result_table$category, get_GO_term_descr))
        result_table$go_term = descr;
        result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) {
          gene_list = GO_to_gene_list[[x]]
          gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
          paste(gene_list, collapse=', ');
        }) )
        write.table(result_table[order(result_table$over_represented_pvalue),], file=go_enrich_filename, sep=';', quote=F, row.names=F)
        ## under-represented categories:
        pvals = res$under_represented_pvalue
        pvals[pvals>1-1e-10] = 1 - 1e-10
        q = qvalue(pvals, lambda=0)
        res$under_represented_FDR = q$qvalues
        go_depleted_filename = paste("./varPart_GOSeq/",infile,"GOseq.depleted.txt", sep='')
        result_table = res[res$under_represented_pvalue<=0.05,]
        descr = unlist(lapply(result_table$category, get_GO_term_descr))
        result_table$go_term = descr;
        result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) {
          gene_list = GO_to_gene_list[[x]]
          gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
          paste(gene_list, collapse=', ');
          write.table(result_table[order(result_table$under_represented_pvalue),], file=go_depleted_filename, sep=';', quote=F, row.names=F)
        }))}
    }    

    
## 3. VISUALIZATION OF DAPC GO ENRICHMENTS (Dotplots) ---------------------------------------------------------------------------------------------------------------
  
  library(ggplot2)
      
  # Step 1 - Merge the GOSeq results of interest and create an input table for dotplot
    GOSeqfile_Ile <- "Pocillopora_varPart_Islands_GOseq.enriched.txt"
    EnrichTab_Ile <- read.table(paste0("./varPart_GOSeq/",GOSeqfile_Ile,sep=''),
                           header = TRUE,
                           quote="\"",
                           sep=";")
    EnrichTab_Ile$variable <- "Island"
    
    GOSeqfile_Clade <- "Pocillopora_varPart_PocilloGG_GOseq.enriched.txt"
    EnrichTab_Clade <- read.table(paste0("./varPart_GOSeq/",GOSeqfile_Clade,sep=''),
                           header = TRUE,
                           quote="\"",
                           sep=";")
    EnrichTab_Clade$variable <- "PocilloGG"
    
    GOSeqfile_SymClade <- "Pocillopora_varPart_SymbioGG_GOseq.enriched.txt"
    EnrichTab_SymClade <- read.table(paste0("./varPart_GOSeq/",GOSeqfile_SymClade,sep=''),
                           header = TRUE,
                           quote="\"",
                           sep=";")
    EnrichTab_SymClade$variable <- "SymbioGG"
    
    EnrichTab <- rbind(EnrichTab_Ile, EnrichTab_Clade, EnrichTab_SymClade)
    EnrichTab$ratio <- EnrichTab$numDEInCat/EnrichTab$numInCat
    EnrichTab$pval <- EnrichTab$over_represented_pvalue
    
  # Step 2 - Create function to make Dotplot of top 10000 (showCategory = 10000) enriched Biological Process categories  
    library(grid)
    dotplot_goseq <- function(df, showCategory=10000){
      df <- df[with(df, order(ratio, pval, decreasing = c(TRUE, FALSE))),]
      df <- head(df, n=showCategory)
      d_plot <- ggplot(subset(df, ontology %in% c("BP")), aes_string(x="term", 
                                                                     y="ratio", 
                                                                     fill="pval",
                                                                     size="numDEInCat")) + 
        geom_point(color="black",
                   pch=21) +
        scale_fill_gradient(low="#00dbde",
                            high="#FFF94C",
                            name = "p-value") +
        scale_size_continuous(range = c(5, 20),
                              breaks=c(1,5,15,25,45,60)) +
        facet_grid(.~variable,
                   drop = TRUE,
                   scales = "free_x", space = "free_x") +
        xlab("") +
        ylab("Count Ratio") +
        scale_x_discrete(expand = c(0.01, 0.2)) +
        #coord_flip() +
        theme_bw(base_size=9) + 
        theme(text = element_text(size = 20, face = 'bold'),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        guides(size=guide_legend(title="N"))
      return(d_plot)
    }
    
    dotplot_goseq(EnrichTab) 
    
  # Step 3 - Plot results

    pdf(file = paste0(outputPrefix,"_AllVariables_Dotplot.pdf"), h = 15, w = 30)
    dotplot_goseq(EnrichTab)
    dev.off()
    
    png(file = paste0(outputPrefix, "_AllVaribales_Dotplot.png"), h = 15, w = 30, units = "in", res = 300)
    dotplot_goseq(EnrichTab)
    dev.off()
 ```
