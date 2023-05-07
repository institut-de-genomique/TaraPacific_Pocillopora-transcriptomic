## TARA Pacific Project - Pocillopora Holobiont MetaT Analyses ##
#################################################################

# Author: ARMSTRONG Eric 
# Created: 21 August 2018
# Last Edited: 28 November 2022

# Modified from a script by SW Davies obtained here: 
# https://github.com/daviessw/Cladocopium_Micronesia/blob/master/Cladocopium_Community_Analyses.R
# https://rstudio-pubs-static.s3.amazonaws.com/706490_2e17c2cf656a42cd96252b43ccc98755.html#3_Discriminant_Analysis_of_Principal_Components_(DAPC)

## 1. INITIALIZATION ------------------------------------------------------------------------------------------------------------------------------------

  # Step 1 - Set Defaults

    # set directory
      setwd("C:/Users/uax75/OneDrive/Documents/R/TaraCoral_2020/Final_Results/Pocillopora/ManuscriptDocs/Scripts/Figure4_DAPC")

    # Set the prefix for each output file name
      outputPrefix <- "TaraPacific_Pocillopora_DAPC_GitHub_2022-11"
    
## 2. LOAD DATA & METADATA -------------------------------------------------------------------------------------------------------------------------------

  # Step 1 - Load Pocillopora Host and Photosymiont raw count data (reads mapped to the predicted coding sequences)
    cts_raw_host <- read.table("../Pocillopora_MetaT_ReadCount.tab", header=TRUE, com='', sep='\t', row.names=1, check.names=FALSE)
    cts_raw_sym <- read.table("../CladocopiumC1_MetaT_ReadCount.tab", header=TRUE, com='', sep='\t', row.names=1, check.names=FALSE)


  # Step 2 - Load table with lineage assignations
    CladeTable <- read.table("../Variables11Islands.txt", header = T)

  # Step 3 - Load historical and in situ environmental data
    library(xlsx)
    envData <- read.xlsx("../TaraPacific_Nutrents_SST_timeseries_mean_products-20220317_11Islands.xlsx", sheetName = "Sheet1")
        
      # Remove extraneous information
        envData <- envData[,-c(1:3)]
    
  # Step 4 - Create Metadata tables with all relevant sampling information
    library(stringr)
    library(expss)
    
    metaData_host <- data.frame("Colony" = colnames(cts_raw_host),
                                "Ile" = substr(colnames(cts_raw_host), 1, 3),
                                "Site" = substr(colnames(cts_raw_host), 4, 6), 
                                "IleSite" = substr(colnames(cts_raw_host), 1, 6),
                                "Clade" = vlookup(lookup_value = colnames(cts_raw_host), 
                                                  dict = CladeTable, 
                                                  result_column = "PocilloGG", 
                                                  lookup_column = "Samples"),
                                "SymClade" = vlookup(lookup_value = colnames(cts_raw_host), 
                                                  dict = CladeTable, 
                                                  result_column = "SymbioGG", 
                                                  lookup_column = "Samples"))
    library(dplyr)
    metaData_host <- left_join(metaData_host, envData, by="Colony")
    rownames(metaData_host) <- colnames(cts_raw_host)
    
    metaData_sym <- data.frame("Colony" = colnames(cts_raw_sym),
                                "Ile" = substr(colnames(cts_raw_sym), 1, 3),
                                "Site" = substr(colnames(cts_raw_sym), 4, 6),
                                "IleSite" = substr(colnames(cts_raw_sym), 1, 6),
                                "Clade" = vlookup(lookup_value = colnames(cts_raw_sym), 
                                                  dict = CladeTable, 
                                                  result_column = "PocilloGG", 
                                                  lookup_column = "Samples"),
                                "SymClade" = vlookup(lookup_value = colnames(cts_raw_sym), 
                                                  dict = CladeTable, 
                                                  result_column = "SymbioGG", 
                                                  lookup_column = "Samples"))
    metaData_sym <- left_join(metaData_sym, envData, by="Colony")
    rownames(metaData_sym) <- colnames(cts_raw_sym)
    
## 3. DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS (DAPC) -------------------------------------------------------------------

  # Step 1 - Merge host and photosymbiont read count data frames by colony
    cols_keep <- intersect(colnames(cts_raw_host), colnames(cts_raw_sym))
    cols_keep <- cols_keep[!cols_keep %in% c("I09S03C010POC", "I09S02C001POC")] # remove the hybrid and unassigned host colonies
    cts_raw_common <- rbind(cts_raw_host[,cols_keep], cts_raw_sym[,cols_keep])
    
  # Step 2 - Split host and photosymbiont reads so that each are vsd normalized separately
    cts_raw_common_host <- cts_raw_common[grep("Pmea", row.names(cts_raw_common)),]
    cts_raw_common_sym <- cts_raw_common[grep("SymbC1", row.names(cts_raw_common)),]
    
  # Step 3 - Prefilter data to remove low count genes
    
    # Remove genes with fewer than 10 counts summed across all samples
      # making a vector counting number of samples with counts <=10 acros 90% of colonies
        lowcountgenes_host <- apply(cts_raw_common_host[,1:ncol(cts_raw_common_host)],1,function(x){sum(x<=10)})
        lowcountgenes_sym <- apply(cts_raw_common_sym[,1:ncol(cts_raw_common_sym)],1,function(x){sum(x<=10)})
      # removing genes with counts < 10 in more than 90% of samples
        highcts_raw_common_host <- cts_raw_common_host[-which(lowcountgenes_host>(0.9*ncol(cts_raw_common_host))),]
        highcts_raw_common_sym <- cts_raw_common_sym[-which(lowcountgenes_sym>(0.9*ncol(cts_raw_common_sym))),]

      nrow(highcts_raw_common_host) #28972 high count genes pass filter for the host
      nrow(highcts_raw_common_sym)  #20517 high count genes pass filter for the photosymbiont
      
  # Step 4 - Normalize the expression data
    library(DESeq2)
    dds_host <- DESeqDataSetFromMatrix(countData = highcts_raw_common_host, 
                                      colData = subset(metaData_host, Colony %in% colnames(highcts_raw_common_host)), 
                                      design = ~ Ile)
    
    dds_sym <- DESeqDataSetFromMatrix(countData = highcts_raw_common_sym,
                                      colData = subset(metaData_host, Colony %in% colnames(highcts_raw_common_sym)), 
                                      design = ~ Ile)

    # DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
    # variance stabilization is very good for heatmaps, etc.
    
      vsd_mat_host <- varianceStabilizingTransformation(dds_host, blind=FALSE)
      cts_common_vsd_host <- assay(vsd_mat_host)

      vsd_mat_sym <- varianceStabilizingTransformation(dds_sym, blind=FALSE)
      cts_common_vsd_sym <- assay(vsd_mat_sym)

  # Step 5 - flip count tables
    cts_common_vsd_host <- t(cts_common_vsd_host)
    cts_common_vsd_sym <- t(cts_common_vsd_sym)
      
  # Step 6 - Perform DAPC with a large number of PCs
    library(adegenet)
    set.seed(999) 

    # Three DAPC models for the host (grouped by island, host genotype, and phtosymbiont genotype)
      HostDAPC_ile <- dapc(x = cts_common_vsd_host, 
                       grp = factor(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Ile), 
                       var.contrib = TRUE,
                       scale = FALSE, 
                       n.pca = 100,
                       n.da = length(unique(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Ile)) - 1)
        
      HostDAPC_clade <- dapc(x = cts_common_vsd_host, 
                         grp = factor(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade), 
                         var.contrib = TRUE, 
                         scale = FALSE, 
                         n.pca = 100, 
                         n.da = length(unique(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade)) - 1)
      
      HostDAPC_symclade <- dapc(x = cts_common_vsd_host, 
                                grp = factor(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$SymClade), 
                                var.contrib = TRUE,
                                scale = FALSE,
                                n.pca = 100,
                                n.da = length(unique(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$SymClade)) - 1)

    # Three DAPC models for the photosymbiont (grouped by island, host genotype, and phtosymbiont genotype)
      SymDAPC_ile <- dapc(x = cts_common_vsd_sym, 
                          grp = factor(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Ile), 
                          var.contrib = TRUE,
                          scale = FALSE,
                          n.pca = 100,
                          n.da = length(unique(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Ile)) - 1)
      
      SymDAPC_clade <- dapc(x = cts_common_vsd_sym, 
                          grp = factor(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Clade), 
                          var.contrib = TRUE,
                          scale = FALSE,
                          n.pca = 100,
                          n.da = length(unique(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Clade)) - 1)
      
      SymDAPC_symclade <- dapc(x = cts_common_vsd_sym, 
                          grp = factor(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade), 
                          var.contrib = TRUE,
                          scale = FALSE,
                          n.pca = 100,
                          n.da = length(unique(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade)) - 1)
      
  # Step 7 - Use a-scores to select optimal number of PCs to retain
    optim.a.score(HostDAPC_ile) #11 PCs optimal
    optim.a.score(HostDAPC_clade) #12 PCs optimal
    optim.a.score(HostDAPC_symclade) #11 PCs optimal
    # Kept 11 PCs for all host DAPCs
    
    optim.a.score(SymDAPC_ile) #11 PCs optimal
    optim.a.score(SymDAPC_clade) #12 PCs optimal
    optim.a.score(SymDAPC_symclade) #11 PCs optimal
    # Kept 11 PCs for all symbiont DAPCs as no large difference between the three maxima
    
  # Step 8 - Rerun DAPCs using optimal number of PCs
    set.seed(999) 

    # Three DAPC models for the host (grouped by island, host genotype, and phtosymbiont genotype)
      HostDAPC_ile <- dapc(x = cts_common_vsd_host, 
                       grp = factor(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Ile), 
                       var.contrib = TRUE,
                       scale = FALSE, 
                       n.pca = 11,
                       n.da = length(unique(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Ile)) - 1)
        
      HostDAPC_clade <- dapc(x = cts_common_vsd_host, 
                         grp = factor(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade), 
                         var.contrib = TRUE, 
                         scale = FALSE, 
                         n.pca = 11, 
                         n.da = length(unique(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade)) - 1)
      
      HostDAPC_symclade <- dapc(x = cts_common_vsd_host, 
                                grp = factor(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$SymClade), 
                                var.contrib = TRUE,
                                scale = FALSE,
                                n.pca = 11,
                                n.da = length(unique(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$SymClade)) - 1)

    # Three DAPC models for the photosymbiont (grouped by island, host genotype, and phtosymbiont genotype)
      SymDAPC_ile <- dapc(x = cts_common_vsd_sym, 
                          grp = factor(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Ile), 
                          var.contrib = TRUE,
                          scale = FALSE,
                          n.pca = 11,
                          n.da = length(unique(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Ile)) - 1)
      
      SymDAPC_clade <- dapc(x = cts_common_vsd_sym, 
                          grp = factor(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Clade), 
                          var.contrib = TRUE,
                          scale = FALSE,
                          n.pca = 11,
                          n.da = length(unique(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Clade)) - 1)
      
      SymDAPC_symclade <- dapc(x = cts_common_vsd_sym, 
                          grp = factor(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade), 
                          var.contrib = TRUE,
                          scale = FALSE,
                          n.pca = 11,
                          n.da = length(unique(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade)) - 1)
      
  # Step 9  Mean and STDEV of proportion assigned
    VarExp_Host_mean_ile <- round(mean(round(summary(HostDAPC_ile)$assign.per.pop, 2)), 2)
    VarExp_Host_sd_ile <- round(sd(round(summary(HostDAPC_ile)$assign.per.pop, 2)), 2)
    
    VarExp_Host_mean_clade <- round(mean(round(summary(HostDAPC_clade)$assign.per.pop, 2)), 2)
    VarExp_Host_sd_clade <- round(sd(round(summary(HostDAPC_clade)$assign.per.pop, 2)), 2)
    
    VarExp_Host_mean_symclade <- round(mean(round(summary(HostDAPC_symclade)$assign.per.pop, 2)), 2)
    VarExp_Host_sd_symclade <- round(sd(round(summary(HostDAPC_symclade)$assign.per.pop, 2)), 2)
    
    VarExp_Sym_mean_ile <- round(mean(round(summary(SymDAPC_ile)$assign.per.pop, 2)), 2)
    VarExp_Sym_sd_ile <- round(sd(round(summary(SymDAPC_ile)$assign.per.pop, 2)), 2)
    
    VarExp_Sym_mean_clade <- round(mean(round(summary(SymDAPC_clade)$assign.per.pop, 2)), 2)
    VarExp_Sym_sd_clade <- round(sd(round(summary(SymDAPC_clade)$assign.per.pop, 2)), 2)
    
    VarExp_Sym_mean_symclade <- round(mean(round(summary(SymDAPC_clade)$assign.per.pop, 2)), 2)
    VarExp_Sym_sd_symclade <- round(sd(round(summary(SymDAPC_clade)$assign.per.pop, 2)), 2)

  # Step 10 - Write discriminant genes
    
  ##### Host Discriminant Genes #####
      # Write list of Island discriminant genes (upper quantile loading scores) for DF1 and DF2: 
        DF1_Ile_TopLoadGenes <- data.frame(HostDAPC_ile$var.contr[which(HostDAPC_ile$var.contr[,1]>=quantile(HostDAPC_ile$var.contr,0.75)),1])
        DF2_Ile_TopLoadGenes <- data.frame(HostDAPC_ile$var.contr[which(HostDAPC_ile$var.contr[,2]>=quantile(HostDAPC_ile$var.contr,0.75)),2])
        HostDAPC_Ile_TopLoadGenes <- data.frame(geneLabels = c(gsub(";.*","",rownames(DF1_Ile_TopLoadGenes)),gsub(";.*","",rownames(DF2_Ile_TopLoadGenes))),
                                            LDScore = c(DF1_Ile_TopLoadGenes[,1],DF2_Ile_TopLoadGenes[,1]),
                                            Loading = c(rep("LD1",length(DF1_Ile_TopLoadGenes[,1])),rep("LD2",length(DF2_Ile_TopLoadGenes[,1]))),
                                            Group = "Ile")
        write.table(HostDAPC_Ile_TopLoadGenes, 
                  file = paste0(outputPrefix, "_HostDAPC_Ile_Both_LoadScores.csv"), sep = ',', row.names = F)
        
      # Write list of Host Clade discriminant genes (upper quantile loading scores) for DF1 and DF2: 
        DF1_Clade_TopLoadGenes <- data.frame(HostDAPC_clade$var.contr[which(HostDAPC_clade$var.contr[,1]>=quantile(HostDAPC_clade$var.contr,0.75)),1])
        DF2_Clade_TopLoadGenes <- data.frame(HostDAPC_clade$var.contr[which(HostDAPC_clade$var.contr[,2]>=quantile(HostDAPC_clade$var.contr,0.75)),2])
        HostDAPC_Clade_TopLoadGenes <- data.frame(geneLabels = c(gsub(";.*","",rownames(DF1_Clade_TopLoadGenes)),gsub(";.*","",rownames(DF2_Clade_TopLoadGenes))),
                                            LDScore = c(DF1_Clade_TopLoadGenes[,1],DF2_Clade_TopLoadGenes[,1]),
                                            Loading = c(rep("LD1",length(DF1_Clade_TopLoadGenes[,1])),rep("LD2",length(DF2_Clade_TopLoadGenes[,1]))),
                                            Group = "Clade")
        write.table(HostDAPC_Clade_TopLoadGenes, 
                  file = paste0(outputPrefix, "_HostDAPC_Clade_Both_LoadScores.csv"), sep = ',', row.names = F)
        
      # Write list of Photosymbiont Clade discriminant genes (upper quantile loading scores) for DF1 and DF2: 
        DF1_SymClade_TopLoadGenes <- data.frame(HostDAPC_symclade$var.contr[which(HostDAPC_symclade$var.contr[,1]>=quantile(HostDAPC_symclade$var.contr,0.75)),1])
        DF2_SymClade_TopLoadGenes <- data.frame(HostDAPC_symclade$var.contr[which(HostDAPC_symclade$var.contr[,2]>=quantile(HostDAPC_symclade$var.contr,0.75)),2])
        HostDAPC_SymClade_TopLoadGenes <- data.frame(geneLabels = c(gsub(";.*","",rownames(DF1_SymClade_TopLoadGenes)),gsub(";.*","",rownames(DF2_SymClade_TopLoadGenes))),
                                            LDScore = c(DF1_SymClade_TopLoadGenes[,1],DF2_SymClade_TopLoadGenes[,1]),
                                            Loading = c(rep("LD1",length(DF1_SymClade_TopLoadGenes[,1])),rep("LD2",length(DF2_SymClade_TopLoadGenes[,1]))),
                                            Group = "SymClade")
        write.table(HostDAPC_SymClade_TopLoadGenes, 
                  file = paste0(outputPrefix, "_HostDAPC_SymClade_Both_LoadScores.csv"), sep = ',', row.names = F)
        
  ##### Photosymbiont Discriminant Genes #####
    # Write list of Island discriminant genes (upper quantile loading scores) for DF1 and DF2: 
      DF1_Sym_Ile_TopLoadGenes <- data.frame(SymDAPC_ile$var.contr[which(SymDAPC_ile$var.contr[,1]>=quantile(SymDAPC_ile$var.contr,0.75)),1])
      DF2_Sym_Ile_TopLoadGenes <- data.frame(SymDAPC_ile$var.contr[which(SymDAPC_ile$var.contr[,2]>=quantile(SymDAPC_ile$var.contr,0.75)),2])
      SymDAPC_Ile_TopLoadGenes <- data.frame(geneLabels = c(gsub(";.*","",rownames(DF1_Sym_Ile_TopLoadGenes)),gsub(";.*","",rownames(DF2_Sym_Ile_TopLoadGenes))),
                                          LDScore = c(DF1_Sym_Ile_TopLoadGenes[,1],DF2_Sym_Ile_TopLoadGenes[,1]),
                                          Loading = c(rep("LD1",length(DF1_Sym_Ile_TopLoadGenes[,1])),rep("LD2",length(DF2_Sym_Ile_TopLoadGenes[,1]))),
                                          Group = "Ile")
      write.table(SymDAPC_Ile_TopLoadGenes, 
                file = paste0(outputPrefix, "_SymDAPC_Ile_Both_LoadScores.csv"), sep = ',', row.names = F)
      
    # Write list of Host Clade discriminant genes (upper quantile loading scores) for DF1 and DF2: 
      DF1_Sym_Clade_TopLoadGenes <- data.frame(SymDAPC_clade$var.contr[which(SymDAPC_clade$var.contr[,1]>=quantile(SymDAPC_clade$var.contr,0.75)),1])
      DF2_Sym_Clade_TopLoadGenes <- data.frame(SymDAPC_clade$var.contr[which(SymDAPC_clade$var.contr[,2]>=quantile(SymDAPC_clade$var.contr,0.75)),2])
      SymDAPC_Clade_TopLoadGenes <- data.frame(geneLabels = c(gsub(";.*","",rownames(DF1_Sym_Clade_TopLoadGenes)),gsub(";.*","",rownames(DF2_Sym_Clade_TopLoadGenes))),
                                          LDScore = c(DF1_Sym_Clade_TopLoadGenes[,1],DF2_Sym_Clade_TopLoadGenes[,1]),
                                          Loading = c(rep("LD1",length(DF1_Sym_Clade_TopLoadGenes[,1])),rep("LD2",length(DF2_Sym_Clade_TopLoadGenes[,1]))),
                                          Group = "Clade")
      write.table(SymDAPC_Clade_TopLoadGenes, 
                file = paste0(outputPrefix, "_SymDAPC_Clade_Both_LoadScores.csv"), sep = ',', row.names = F)
      
    # Write list of Photosymbiont Clade discriminant genes (upper quantile loading scores) for DF1 and DF2: 
      DF1_Sym_SymClade_TopLoadGenes <- data.frame(SymDAPC_symclade$var.contr[which(SymDAPC_symclade$var.contr[,1]>=quantile(SymDAPC_symclade$var.contr,0.75)),1])
      DF2_Sym_SymClade_TopLoadGenes <- data.frame(SymDAPC_symclade$var.contr[which(SymDAPC_symclade$var.contr[,2]>=quantile(SymDAPC_symclade$var.contr,0.75)),2])
      SymDAPC_SymClade_TopLoadGenes <- data.frame(geneLabels = c(gsub(";.*","",rownames(DF1_Sym_SymClade_TopLoadGenes)),gsub(";.*","",rownames(DF2_Sym_SymClade_TopLoadGenes))),
                                          LDScore = c(DF1_Sym_SymClade_TopLoadGenes[,1],DF2_Sym_SymClade_TopLoadGenes[,1]),
                                          Loading = c(rep("LD1",length(DF1_Sym_SymClade_TopLoadGenes[,1])),rep("LD2",length(DF2_Sym_SymClade_TopLoadGenes[,1]))),
                                          Group = "SymClade")
      write.table(SymDAPC_SymClade_TopLoadGenes, 
                file = paste0(outputPrefix, "_SymDAPC_SymClade_Both_LoadScores.csv"), sep = ',', row.names = F)

## 4. VISUALIZATIONS (Violin and Star Plots) -------------------------------------------------------------------------------------------------------------------------------------------------------------
  ##### Figure 4a - Violin Plots of Discriminant Genes and Loading Scores #####
    library(ggplot2)
    Ile_DiscrimGenes <- setdiff(HostDAPC_Ile_TopLoadGenes$geneLabels,c(HostDAPC_SymClade_TopLoadGenes$geneLabels, HostDAPC_Clade_TopLoadGenes$geneLabels))
    Clade_DiscrimGenes <- setdiff(HostDAPC_Clade_TopLoadGenes$geneLabels,c(HostDAPC_Ile_TopLoadGenes$geneLabels, HostDAPC_SymClade_TopLoadGenes$geneLabels))
    SymClade_DiscrimGenes <- setdiff(HostDAPC_SymClade_TopLoadGenes$geneLabels,c(HostDAPC_Ile_TopLoadGenes$geneLabels, HostDAPC_Clade_TopLoadGenes$geneLabels))
                                   
    Host_DiscrimGenes <- rbind(HostDAPC_Ile_TopLoadGenes,
                               HostDAPC_Clade_TopLoadGenes,
                               subset(HostDAPC_SymClade_TopLoadGenes, geneLabels %in% SymClade_DiscrimGenes))
    Host_DiscrimGenes$Group2 <- gsub("SymClade|Clade","Lineage",Host_DiscrimGenes$Group)
    
      
    Ile_Sym_DiscrimGenes <- setdiff(SymDAPC_Ile_TopLoadGenes$geneLabels,c(SymDAPC_SymClade_TopLoadGenes$geneLabels, SymDAPC_Clade_TopLoadGenes$geneLabels))
    Clade_Sym_DiscrimGenes <- setdiff(SymDAPC_Clade_TopLoadGenes$geneLabels,c(SymDAPC_Ile_TopLoadGenes$geneLabels, SymDAPC_SymClade_TopLoadGenes$geneLabels))
    SymClade_Sym_DiscrimGenes <- setdiff(SymDAPC_SymClade_TopLoadGenes$geneLabels,c(SymDAPC_Ile_TopLoadGenes$geneLabels, SymDAPC_Clade_TopLoadGenes$geneLabels))
    
    Symb_DiscrimGenes <- rbind(SymDAPC_Ile_TopLoadGenes,
                                SymDAPC_Clade_TopLoadGenes,
                                subset(SymDAPC_SymClade_TopLoadGenes, geneLabels %in% SymClade_Sym_DiscrimGenes))
    Symb_DiscrimGenes$Group2 <- gsub("SymClade|Clade","Lineage",Symb_DiscrimGenes$Group)
    
    
    Symb_DiscrimGenes$Model <- "Symbiont"
    Host_DiscrimGenes$Model <- "Pocillopora"
    Both_DiscrimGenes <- rbind(Host_DiscrimGenes, Symb_DiscrimGenes)
    
    write.table(file = paste0(outputPrefix,"_Figure4a_Data.tab"), Both_DiscrimGenes, quote = F, sep = '\t', row.names = F)
    
    discrim_Host_violplot <- ggplot() +
      geom_violin(data=subset(Both_DiscrimGenes, Model %in% c("Pocillopora")), aes(x = Group2, y = LDScore*1000, fill = Group2)) +
      annotate(geom = "text", x = 1, y = 0.015, label = "4,386 genes", size = 10) +
      annotate(geom = "text", x = 2, y = 0.015, label = "1,396 genes", size = 10) +
      annotate(geom = "text", x = 2, y = 0.01, label = "687 host-dependent/n 709 symbiont-dependent", size = 6) +
      scale_fill_manual(values=c("#66C2A5","#8a128a")) +
      coord_trans(y = "log10") +
      #scale_y_log10() + 
      #scale_y_continuous(trans='log10') +
      scale_y_continuous(limits = c(0.01,10), breaks = c(0,0.03,0.3,1,3,6,9)) +
      scale_x_discrete(breaks=c("Ile","Lineage"), labels=c("Environment", "Lineage/n(Host + Photosymbiont)")) +
      #facet_wrap(~Model) +
      xlab("Variables") +
      #ylab("Loading Score") +
      ylab(expression(bold(paste("Loading Score (x", 10^-3,")")))) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            text = element_text(size=20, face="bold"),
            legend.position ='none')
    discrim_Host_violplot
    
    pdf(file = paste0(outputPrefix, "_DiscrimGenes_Host.pdf"), w = 8.5, h=11)
    discrim_Host_violplot
    dev.off()
    
    png(file = paste0(outputPrefix, "_DiscrimGenes_Host.png"),  w = 8.5, h=11, units="in", res=300)
    discrim_Host_violplot
    dev.off()
    
    discrim_Symbiont_violplot <- ggplot() +
      geom_violin(data=subset(Both_DiscrimGenes, Model %in% c("Symbiont")), aes(x = Group2, y = LDScore*1000, fill = Group2)) +
      annotate(geom = "text", x = 1, y = 0.015, label = "2,495 genes", size = 10) +
      annotate(geom = "text", x = 2, y = 0.015, label = "2,041 genes", size = 10) +
      annotate(geom = "text", x = 2, y = 0.01, label = "1,544 host-dependet genes/n 497 symbiont-dependent genes", size = 6) +
      scale_fill_manual(values=c("#66C2A5","#8a128a")) +
      coord_trans(y = "log10") +
      #scale_y_log10() + 
      #scale_y_continuous(trans='log10') +
      scale_y_continuous(limits = c(0.01,10), breaks = c(0,0.03,0.3,1,3,6,9)) +
      scale_x_discrete(breaks=c("Ile","Lineage"), labels=c("Environment", "Lineage/n(Symbiont + Photosymbiont)")) +
      #facet_wrap(~Model) +
      xlab("Variables") +
      #ylab("Loading Score") +
      ylab(expression(bold(paste("Loading Score (x", 10^-3,")")))) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            text = element_text(size=20, face="bold"),
            legend.position ='none')
    discrim_Symbiont_violplot
    
    pdf(file = paste0(outputPrefix, "_DiscrimGenes_Symbiont.pdf"), w = 8.5, h=11)
    discrim_Symbiont_violplot
    dev.off()
    
    png(file = paste0(outputPrefix, "_DiscrimGenes_Symbiont.png"),  w = 8.5, h=11, units="in", res=300)
    discrim_Symbiont_violplot
    dev.off()

  ##### Figure 4b - Clustered Star Plots of Discriminant Genes #####
          
    # Step 1 - Assign Aesthetics (Colors, etc.)
      library(RColorBrewer)
      library(rcartocolor)
      
      # Island Colors
        mycols <- brewer.pal(11, "Spectral")
        mycols[6] <- "grey40"
      
      # Pocillopora Host SVD Colors
        mybrew_cols <- brewer.pal(12, "Paired")
        mycols_clade <- c(mybrew_cols[5], mybrew_cols[1], mybrew_cols[3],
                          mybrew_cols[2], mybrew_cols[4])
        # these correspond (in order) to: SVD1, SVD2, SVD3, SVD4, and SVD5
      
      # Symbiont Lineage Colors
        mycarto_cols <- carto_pal(12, "Safe")
        mycols_symclade <- c(mycarto_cols[1], mycarto_cols[3], mycarto_cols[4],
                             mycarto_cols[2], mycarto_cols[7], mycarto_cols[5],
                             mycarto_cols[11], "white")
        # these correspond (in order) to: Cladocopium C1 (undefined), D1, Lineage 1, Lineage 2, Lineage 3, Lineage 4, Lineage 5, Unassigned
      
    # Step 2 - Build ggplot data frame with points (x,y) and corresponding groups
      # For the Host:
        gg_ile <- as.data.frame(HostDAPC_ile$ind.coord)
        gg_ile$Ile <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Ile
        gg_ile$Clade <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade
        gg_ile$SymClade <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$SymClade
        
        gg_clade <- as.data.frame(HostDAPC_clade$ind.coord)
        gg_clade$Ile <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Ile
        gg_clade$Clade <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade
        gg_clade$SymClade <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$SymClade
        
        gg_symclade <- as.data.frame(HostDAPC_symclade$ind.coord)
        gg_symclade$Ile <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Ile
        gg_symclade$Clade <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade
        gg_symclade$SymClade <- subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$SymClade
      
      # For the Photosymbiont:
        gg_ile_symb <- as.data.frame(SymDAPC_ile$ind.coord)
        gg_ile_symb$Ile <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Ile
        gg_ile_symb$Clade <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Clade
        gg_ile_symb$SymClade <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade
        
        gg_clade_symb <- as.data.frame(SymDAPC_clade$ind.coord)
        gg_clade_symb$Ile <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Ile
        gg_clade_symb$Clade <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Clade
        gg_clade_symb$SymClade <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade
        
        gg_symclade_symb <- as.data.frame(SymDAPC_symclade$ind.coord)
        gg_symclade_symb$Ile <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Ile
        gg_symclade_symb$Clade <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Clade
        gg_symclade_symb$SymClade <- subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade

  # Step 3 - calculate group centroid locations
    library(Rmisc)
    ileSST_hist <- summarySE(data = metaData_host,
                        measurevar = "SST_mean_DegC",
                        groupvars = c("Ile"),
                        na.rm = TRUE)
    ileSST_insitu <- summarySE(data = metaData_host,
                        measurevar = "SST_snapshot_sampling_day_DegC",
                        groupvars = c("Ile"),
                        na.rm = TRUE)
    TSA_hist <- summarySE(data = metaData_host,
                        measurevar = "TSA_heat_mean_DegC",
                        groupvars = c("Ile"),
                        na.rm = TRUE)
    
    ilematch_df <- data.frame("Ile" = unique(metaData_host$Ile),
                              "Name" = c("Las Perlas", "Coiba", "Malpelo",
                                    "Rapa Nui", "Ducie", "Gambier",
                                    "Moorea","Aitutaki","Niue","Upolu","Guam"))
    linmatch_df <- data.frame("Clade" = unique(subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade),
                              "Name" = c("P. grandis (SVD4)","SSH5_Pver (SVD5)","P. effusa (SVD1)","P. meandrina (SVD2)","P. verrucosa (SVD3)"))
    symlinmatch_df <- data.frame("SymClade" = unique(subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade),
                              "Name" = c("C.pacificum (L5)","C. goreaui (L1)","C. pacificum (L4)","C. latusorum (L2)","C. latusorum (L3)"))
    
    centroids_ile <- aggregate(cbind(LD1,LD2)~subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Ile,data=gg_ile,mean)
    colnames(centroids_ile) <- c("Ile","LD1","LD2")
    centroids_ile$Label <- paste(vlookup(lookup_value = centroids_ile$Ile,
                                         dict = ilematch_df,
                                         lookup_column = "Ile",
                                         result_column = "Name"),
                                 formatC(round(summary(HostDAPC_ile)$assign.per.pop, 3), format='f', digits=2), sep=": ")
    centroids_ile$SST_hist <- vlookup(lookup_value = centroids_ile$Ile,
                                 dict = ileSST_hist,
                                 lookup_column = "Ile",
                                 result_column = "SST_mean_DegC")
    centroids_ile$SST_insitu <- vlookup(lookup_value = centroids_ile$Ile,
                                 dict = ileSST_insitu,
                                 lookup_column = "Ile",
                                 result_column = "SST_snapshot_sampling_day_DegC")
    centroids_ile$TSA_hist <- vlookup(lookup_value = centroids_ile$Ile,
                                 dict = TSA_hist,
                                 lookup_column = "Ile",
                                 result_column = "TSA_heat_mean_DegC")
    
    centroids_clade <- aggregate(cbind(LD1,LD2)~subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$Clade,data=gg_clade,mean)
    colnames(centroids_clade) <- c("Clade","LD1","LD2")
    centroids_clade$Label <- paste(vlookup(lookup_value = centroids_clade$Clade,
                                         dict = linmatch_df,
                                         lookup_column = "Clade",
                                         result_column = "Name"),
                                 formatC(round(summary(HostDAPC_clade)$assign.per.pop, 3), format='f', digits=2), sep=": ")
    
    centroids_symclade <- aggregate(cbind(LD1,LD2)~subset(metaData_host, Colony %in% colnames(highcts_raw_common_host))$SymClade,data=gg_symclade,mean)
    colnames(centroids_symclade) <- c("SymClade","LD1","LD2")
    centroids_symclade$Label <- paste(vlookup(lookup_value = centroids_symclade$SymClade,
                                         dict = symlinmatch_df,
                                         lookup_column = "SymClade",
                                         result_column = "Name"),
                                 formatC(round(summary(HostDAPC_symclade)$assign.per.pop, 3), format='f', digits=2), sep=": ")
    
    centroids_ile_symb <- aggregate(cbind(LD1,LD2)~subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Ile,data=gg_ile_symb,mean)
    colnames(centroids_ile_symb) <- c("Ile","LD1","LD2")
    centroids_ile_symb$Label <- paste(vlookup(lookup_value = centroids_ile_symb$Ile,
                                         dict = ilematch_df,
                                         lookup_column = "Ile",
                                         result_column = "Name"),
                                 formatC(round(summary(SymDAPC_ile)$assign.per.pop, 3), format='f', digits=2), sep=": ")
    
    centroids_clade_symb <- aggregate(cbind(LD1,LD2)~subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$Clade,data=gg_clade_symb,mean)
    colnames(centroids_clade_symb) <- c("Clade","LD1","LD2")
    centroids_clade_symb$Label <- paste(vlookup(lookup_value = centroids_clade_symb$Clade,
                                         dict = linmatch_df,
                                         lookup_column = "Clade",
                                         result_column = "Name"),
                                 formatC(round(summary(SymDAPC_clade)$assign.per.pop, 3), format='f', digits=2), sep=": ")
    
    centroids_symclade_symb <- aggregate(cbind(LD1,LD2)~subset(metaData_sym, Colony %in% colnames(highcts_raw_common_sym))$SymClade,data=gg_symclade_symb,mean)
    colnames(centroids_symclade_symb) <- c("SymClade","LD1","LD2")
    centroids_symclade_symb$Label <-  paste(vlookup(lookup_value = centroids_symclade_symb$SymClade,
                                         dict = symlinmatch_df,
                                         lookup_column = "SymClade",
                                         result_column = "Name"),
                                 formatC(round(summary(SymDAPC_symclade)$assign.per.pop, 3), format='f', digits=2), sep=": ")
    
  # Step 4 - merge centroid locations into ggplot dataframe
    gg_ile <- merge(gg_ile, centroids_ile, by="Ile", suffixes=c("",".centroid"))
    gg_clade <- merge(gg_clade, centroids_clade, by="Clade", suffixes=c("",".centroid"))
    gg_symclade <- merge(gg_symclade, centroids_symclade, by="SymClade", suffixes=c("",".centroid"))
    
    gg_ile_symb <- merge(gg_ile_symb, centroids_ile_symb, by="Ile", suffixes=c("",".centroid"))
    gg_clade_symb <- merge(gg_clade_symb, centroids_clade_symb, by="Clade", suffixes=c("",".centroid"))
    gg_symclade_symb <- merge(gg_symclade_symb, centroids_symclade_symb, by="SymClade", suffixes=c("",".centroid"))

      
  # Step 5 - Generate Star Plots
    library(ggplot2)
    library(ggrepel)
    library(ggnewscale)
    library(cowplot)
    
    ##### Host DAPC Star Plots #####
      
      dapc_ile_DFeig <- data.frame("DF" = seq(from=1, to=length(HostDAPC_ile$eig), by=1),
                                   "DFeigval" = HostDAPC_ile$eig,
                                   "DFkept" = c(rep("Kept",2),rep("Not",length(HostDAPC_ile$eig)-2)))
      dapc_ile_PCeig <- data.frame("PC" = seq(from=1, to=length(HostDAPC_ile$pca.eig), by=1),
                                 "PCcumvar" = 100 * cumsum(HostDAPC_ile$pca.eig)/sum(HostDAPC_ile$pca.eig),
                                 "PCkept" = c(rep("Kept",HostDAPC_ile$n.pca),rep("Not",length(HostDAPC_ile$pca.eig)-HostDAPC_ile$n.pca)))
      
      write.table(file = paste0(outputPrefix,"_FigureS4a_Data.tab"), dapc_ile_PCeig, quote = F, sep = '\t', row.names = F)
      write.table(file = paste0(outputPrefix,"_FigureS4b_Data.tab"), dapc_ile_DFeig, quote = F, sep = '\t', row.names = F)
  
      DFbox_ile <- ggplot() +
        geom_bar(data = dapc_ile_DFeig, aes(x=DF,y=DFeigval, fill=DFkept), 
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0.5,10.5), breaks = seq(from=1, to=10, by=1)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("Discriminant Function") +
        ylab("DF eigenvalues") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      DFbox_ile
      
      PCbox_ile <- ggplot() +
        geom_bar(data = dapc_ile_PCeig, aes(x=PC,y=PCcumvar, fill=PCkept), 
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0,length(HostDAPC_ile$pca.eig)+0.5), breaks = seq(from=0, to=length(HostDAPC_ile$pca.eig), by=20)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("PCA Axis") +
        ylab("Cumulative Variance Explained (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      PCbox_ile
      
      write.table(file = paste0(outputPrefix,"_Figure4b_Data.tab"), gg_ile, quote = F, sep = '\t', row.names = F)
      
      options(ggrepel.max.overlaps = Inf)
      plot_dapc_ile <- ggplot() +
        stat_ellipse(data=gg_ile, geom="polygon", aes(x=LD1,y=LD2,color=Ile, fill=TSA_hist), 
                     type="t", level=0.95, alpha=0.5, show.legend=F) +
        geom_segment(data=gg_ile, show.legend=F,
                     aes(x=LD1.centroid, y=LD2.centroid, xend=LD1, yend=LD2, color=Ile)) +
        annotate(geom="text", x=7, y=-10, 
                 label=paste("Proportion Reassigned /n", VarExp_Host_mean_ile, "?", VarExp_Host_sd_ile, sep=' '),
                 color="black", size = 14) +
        scale_fill_distiller(palette = "RdBu") +
        scale_color_manual(values=mycols,
                           name = "Island",
                           labels = c("Las Perlas", "Coiba", "Malpelo",
                                      "Rapa Nui", "Ducie", "Gambier",
                                      "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        new_scale("fill") +
        geom_point(data=gg_ile, aes(x=LD1,y=LD2,fill=Clade), color="black",
                   size=10, pch=21) +
        scale_fill_manual(values=mycols_clade,
                          name = expression(bold(paste(bolditalic("Pocillopora"),"SVD Clade", sep=''))),
                          labels = c("SVD 1", "SVD 2", "SVD 3",
                                     "SVD 4", "SVD 5", "Unassigned")) +
        new_scale("fill") +
        geom_text_repel(data=centroids_ile, aes(x=LD1, y=LD2, label=Label),
                         alpha=0.5, color="black", size=20, show.legend = FALSE) +
        xlab("Discriminant Function 1") +
        ylab("Discriminant Function 2") +
        scale_x_continuous(limits = c(-15,15), breaks = seq(from=-15, to=15, by=5)) +
        scale_y_continuous(limits = c(-12,10), breaks = seq(from=-10, to=10, by=5)) +
        scale_fill_manual(values=mycols,
                          name = "Island",
                          labels = c("Las Perlas", "Coiba", "Malpelo",
                                     "Rapa Nui", "Ducie", "Gambier",
                                     "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=40, face="bold"), legend.box ='horizontal',
              legend.position = "top")
      plot_dapc_ile
      
      pdf(file = paste0(outputPrefix, "_Host_Ile_TSAhistmean.pdf"), w = 20, h=20)
      plot_dapc_ile
      dev.off()
        
      png(file = paste0(outputPrefix, "_Host_Ile_TSAhistmean.png"),  w = 20, h=20, units="in", res=300)
      plot_dapc_ile
      dev.off()
      
      # Plot thermal stress anomaly data
      library(reshape2)
      plot_TSA <- ggplot() +
        geom_point(data=TSA_hist, aes(x=Ile, y = TSA_heat_mean_DegC, fill=TSA_heat_mean_DegC), shape=21) +
        scale_fill_distiller(palette = "RdBu",
                             name = "Climatological Mean TSA (?C)")
      plot_TSA
        
      plot_TSA_leg <- get_legend(plot_TSA)
      plot_TSA_leg <- plot_grid(plot_TSA_leg)
        
      pdf(file = paste0(outputPrefix, "_TSALegend.pdf"), w = 2, h=2)
      plot_TSA_leg
      dev.off()
        
      png(file = paste0(outputPrefix, "_TSALegend.png"),  w = 2, h=2, units="in", res=300)
      plot_TSA_leg
      dev.off()
        
      plot_TSA_gg <- ggplot_build(plot_TSA)
      mycols_TSA <- as.character(plot_TSA_gg$data[[1]]$fill)
  
      library(ggplotify)
      plot_HostDAPC_ile_DF2_TSAfill <- as.ggplot(~scatter(HostDAPC_ile, 2, 2, col=mycols_TSA, bg="white", scree.da=T, legend=FALSE, solid=0.4, ylim = c(-12, 10)))
      plot_HostDAPC_ile_DF2_TSAfill
        
      pdf(file = paste0(outputPrefix, "_Host_Ile_DF2_Density.pdf"), w = 11, h=8.5)
      plot_HostDAPC_ile_DF2_TSAfill
      dev.off()
      
      png(file = paste0(outputPrefix, "_Host_Ile_DF2_Density.png"), w = 11, h=8.5, units = "in", res = 300)
      plot_HostDAPC_ile_DF2_TSAfill
      dev.off()
      
      
      dapc_clade_DFeig <- data.frame("DF" = seq(from=1, to=length(HostDAPC_clade$eig), by=1),
                                   "DFeigval" = HostDAPC_clade$eig,
                                   "DFkept" = c(rep("Kept",2),rep("Not",length(HostDAPC_clade$eig)-2)))
      dapc_clade_PCeig <- data.frame("PC" = seq(from=1, to=length(HostDAPC_clade$pca.eig), by=1),
                                 "PCcumvar" = 100 * cumsum(HostDAPC_clade$pca.eig)/sum(HostDAPC_clade$pca.eig),
                                 "PCkept" = c(rep("Kept", HostDAPC_clade$n.pca),rep("Not",length(HostDAPC_clade$pca.eig)-HostDAPC_clade$n.pca)))
  
      write.table(file = paste0(outputPrefix,"_FigureS4c_Data.tab"), dapc_clade_DFeig, quote = F, sep = '\t', row.names = F)
      
      DFbox_clade <- ggplot() +
        geom_bar(data = dapc_clade_DFeig, aes(x=DF,y=DFeigval, fill=DFkept), 
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0.5,length(HostDAPC_clade$eig)+0.5), breaks = seq(from=1, to=length(HostDAPC_clade$eig), by=1)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("Discriminant Function") +
        ylab("DF eigenvalues") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      DFbox_clade
      
      PCbox_clade <- ggplot() +
        geom_bar(data = dapc_clade_PCeig, aes(x=PC,y=PCcumvar, fill=PCkept), 
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0,length(HostDAPC_clade$pca.eig)+0.5), breaks = seq(from=0, to=length(HostDAPC_clade$pca.eig), by=20)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("PCA Axis") +
        ylab("Cumulative Variance Explained (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      PCbox_clade
      
      write.table(file = paste0(outputPrefix,"_FigureS6b_Data.tab"), gg_clade, quote = F, sep = '\t', row.names = F)
      
      plot_dapc_clade <- ggplot() +
        stat_ellipse(data=gg_clade, geom="polygon", aes(x=LD1,y=LD2,color=Clade, fill=Clade), 
                     type="t", level=0.95, alpha=0.5, show.legend=F,) +
        geom_segment(data=gg_clade, show.legend=F,
                     aes(x=LD1.centroid, y=LD2.centroid, xend=LD1, yend=LD2, color=Clade)) +
        annotate(geom="text", x=-27.5, y=-16,
                 label=paste("Proportion Reassigned /n", VarExp_Host_mean_clade, "?", VarExp_Host_sd_clade, sep=' '),
                 color="black", size = 14) +
        scale_color_manual(values=mycols_clade,
                           name = "SVD Clade",
                           labels = c("SVD 1", "SVD 2", "SVD 3",
                                      "SVD 4", "SVD 5", "Unassigned")) +
        scale_fill_manual(values=mycols_clade,
                          name = "SVD Clade",
                          labels = c("SVD 1", "SVD 2", "SVD 3",
                                     "SVD 4", "SVD 5", "Unassigned")) +
        new_scale("fill") +
        geom_point(data=gg_clade, aes(x=LD1,y=LD2,fill=Ile), color="black", 
                   size=10, pch=21) +
        scale_fill_manual(values=mycols,
                          name = "Island",
                          labels = c("Las Perlas", "Coiba", "Malpelo",
                                     "Rapa Nui", "Ducie", "Gambier",
                                     "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        new_scale("fill") +
        geom_text_repel(data=centroids_clade, aes(x=LD1, y=LD2, label=Label),
                         alpha=0.5, color="black", size=20, show.legend = FALSE) +
        scale_fill_manual(values=mycols_clade,
                          name = "SVD Clade",
                          labels = c("SVD 1", "SVD 2", "SVD 3",
                                     "SVD 4", "SVD 5", "Unassigned")) +
        xlab("Discriminant Function 1") +
        ylab("Discriminant Function 2") +
        scale_x_continuous(limits = c(-35,11), breaks = seq(from=-35, to=10, by=5)) +
        scale_y_continuous(limits = c(-17,15), breaks = seq(from=-15, to=15, by=5)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=40, face="bold"),
              legend.box ='horizontal', legend.position = "top")
      plot_dapc_clade
      
      pdf(file = paste0(outputPrefix, "_Host_Clade.pdf"), w = 20, h=20)
      plot_dapc_clade
      dev.off()
      
      png(file = paste0(outputPrefix, "_Host_Clade.png"),  w = 20, h=20, units="in", res=300)
      plot_dapc_clade
      dev.off()
      
      dapc_symclade_DFeig <- data.frame("DF" = seq(from=1, to=length(HostDAPC_symclade$eig), by=1),
                                   "DFeigval" = HostDAPC_symclade$eig,
                                   "DFkept" = c(rep("Kept",2),rep("Not",length(HostDAPC_symclade$eig)-2)))
      dapc_symclade_PCeig <- data.frame("PC" = seq(from=1, to=length(HostDAPC_symclade$pca.eig), by=1),
                                 "PCcumvar" = 100 * cumsum(HostDAPC_symclade$pca.eig)/sum(HostDAPC_symclade$pca.eig),
                                 "PCkept" = c(rep("Kept", HostDAPC_symclade$n.pca),rep("Not",length(HostDAPC_symclade$pca.eig)-HostDAPC_symclade$n.pca)))
  
      write.table(file = paste0(outputPrefix,"_FigureS4d_Data.tab"), dapc_symclade_DFeig, quote = F, sep = '\t', row.names = F)
      
      DFbox_symclade <- ggplot() +
        geom_bar(data = dapc_symclade_DFeig, aes(x=DF,y=DFeigval, fill=DFkept), 
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0.5,length(HostDAPC_symclade$eig)+0.5), breaks = seq(from=1, to=length(HostDAPC_symclade$eig), by=1)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("Discriminant Function") +
        ylab("DF eigenvalues") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      DFbox_symclade
      
      PCbox_symclade <- ggplot() +
        geom_bar(data = dapc_symclade_PCeig, aes(x=PC,y=PCcumvar, fill=PCkept), 
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0,length(HostDAPC_symclade$pca.eig)+0.5), breaks = seq(from=0, to=length(HostDAPC_symclade$pca.eig), by=20)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("PCA Axis") +
        ylab("Cumulative Variance Explained (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      PCbox_symclade
      
      write.table(file = paste0(outputPrefix,"_FigureS7b_Data.tab"), gg_symclade, quote = F, sep = '\t', row.names = F)
      
      plot_dapc_symclade <- ggplot() +
        stat_ellipse(data=gg_symclade, geom="polygon", aes(x=LD1,y=LD2,color=Clade, fill=Clade), 
                     type="t", level=0.95, alpha=0.5, show.legend=F,) +
        geom_segment(data=gg_symclade, show.legend=F,
                     aes(x=LD1.centroid, y=LD2.centroid, xend=LD1, yend=LD2, color=Clade)) +
        annotate(geom="text", x=8, y=-16,
                 label=paste("Proportion Reassigned /n", VarExp_Host_mean_symclade, "?", VarExp_Host_sd_symclade, sep=' '),
                 color="black", size = 14) +
        scale_color_manual(values=mycols_symclade[c(3:7)],
                             name = "Cladocopium Lineage",
                             labels = c("L1", "L2", "L3",
                                        "L4", "L5", "Unassigned")) +
          scale_fill_manual(values=mycols_symclade[c(3:7)],
                            name = "Cladocopium Lineage",
                            labels = c("L1", "L2", "L3",
                                       "L4", "L5", "Unassigned")) +
        new_scale("fill") +
        geom_point(data=gg_symclade, aes(x=LD1,y=LD2,fill=Ile), color="black", 
                   size=10, pch=21) +
        scale_fill_manual(values=mycols,
                          name = "Island",
                          labels = c("Las Perlas", "Coiba", "Malpelo",
                                     "Rapa Nui", "Ducie", "Gambier",
                                     "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        new_scale("fill") +
        geom_text_repel(data=centroids_symclade, aes(x=LD1, y=LD2, label=Label),
                         alpha=0.5, color="black", size=20, show.legend = FALSE) +
        scale_fill_manual(values=mycols_symclade[c(3:7)],
                            name = "SVD SymClade",
                            labels = c("SVD 1", "SVD 2", "SVD 3",
                                       "SVD 4", "SVD 5", "Unassigned")) +
        xlab("Discriminant Function 1") +
        ylab("Discriminant Function 2") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=40, face="bold"),
              legend.box ='horizontal', legend.position = "top")
      plot_dapc_symclade
      
      pdf(file = paste0(outputPrefix, "_Host_SymClade.pdf"), w = 20, h=20)
      plot_dapc_symclade
      dev.off()
      
      png(file = paste0(outputPrefix, "_Host_SymClade.png"),  w = 20, h=20, units="in", res=300)
      plot_dapc_symclade
      dev.off()
    
    ##### Symbiont DAPC Plots #####
      
      dapc_ile_symb_DFeig <- data.frame("DF" = seq(from=1, to=length(SymDAPC_ile$eig), by=1),
                                   "DFeigval" = SymDAPC_ile$eig,
                                   "DFkept" = c(rep("Kept",2),rep("Not",length(SymDAPC_ile$eig)-2)))
      dapc_ile_symb_PCeig <- data.frame("PC" = seq(from=1, to=length(SymDAPC_ile$pca.eig), by=1),
                                 "PCcumvar" = 100 * cumsum(SymDAPC_ile$pca.eig)/sum(SymDAPC_ile$pca.eig),
                                 "PCkept" = c(rep("Kept",SymDAPC_ile$n.pca),rep("Not",length(SymDAPC_ile$pca.eig)-SymDAPC_ile$n.pca)))
  
      write.table(file = paste0(outputPrefix,"_FigureS4e_Data.tab"), dapc_ile_symb_PCeig, quote = F, sep = '\t', row.names = F)
      write.table(file = paste0(outputPrefix,"_FigureS4f_Data.tab"), dapc_ile_symb_DFeig, quote = F, sep = '\t', row.names = F)
      
      DFbox_ile_symb <- ggplot() +
        geom_bar(data = dapc_ile_symb_DFeig, aes(x=DF,y=DFeigval, fill=DFkept), 
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0.5,length(SymDAPC_ile$eig)+0.5), breaks = seq(from=1, to=length(SymDAPC_ile$eig), by=1)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("Discriminant Function") +
        ylab("DF eigenvalues") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      DFbox_ile_symb
      
      PCbox_ile_symb <- ggplot() +
        geom_bar(data = dapc_ile_symb_PCeig, aes(x=PC,y=PCcumvar, fill=PCkept), 
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0,length(SymDAPC_ile$pca.eig)+0.5), breaks = seq(from=0, to=length(SymDAPC_ile$pca.eig), by=20)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("PCA Axis") +
        ylab("Cumulative Variance Explained (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      PCbox_ile_symb
      
      write.table(file = paste0(outputPrefix,"_FigureS5d_Data.tab"), gg_ile_symb, quote = F, sep = '\t', row.names = F)
      
      plot_dapc_ile_symb <- ggplot() +
        stat_ellipse(data=gg_ile_symb, geom="polygon", aes(x=LD1,y=LD2,color=Ile, fill=Ile), 
                     type="t", level=0.95, alpha=0.5, show.legend=F,) +
        geom_segment(data=gg_ile_symb, show.legend=F,
                     aes(x=LD1.centroid, y=LD2.centroid, xend=LD1, yend=LD2, color=Ile)) +
        annotate(geom="text", x=70, y=10,
                 label=paste("Proportion Reassigned /n", VarExp_Sym_mean_ile, "?", VarExp_Sym_sd_ile, sep=' '),
                 color="black", size = 14) +
        scale_color_manual(values=mycols,
                           name = "Island",
                           labels = c("Las Perlas", "Coiba", "Malpelo",
                                      "Rapa Nui", "Ducie", "Gambier",
                                      "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        scale_fill_manual(values=mycols,
                          name = "Island",
                          labels = c("Las Perlas", "Coiba", "Malpelo",
                                     "Rapa Nui", "Ducie", "Gambier",
                                     "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        new_scale("fill") +
        geom_point(data=gg_ile_symb, aes(x=LD1,y=LD2,fill=SymClade), color="black",
                   size=10, pch=21) +
        scale_fill_manual(values=mycols_symclade[c(3:7)],
                          #name = expression(bold(atop(bolditalic("Cladocopium"), "Genetic Group",))),
                          name = expression(bold(paste(bolditalic("Cladocopium"), " Lineage",))),
                          labels = c("L1", "L2", "L3",
                                     "L4", "L5")) +
        new_scale("fill") +
        geom_text_repel(data=centroids_ile_symb, aes(x=LD1, y=LD2, label=Label),
                         alpha=0.5, color="black", size=20, show.legend = FALSE) +
        xlab("Discriminant Function 1") +
        ylab("Discriminant Function 2") +
        scale_fill_manual(values=mycols,
                          name = "Island",
                          labels = c("Las Perlas", "Coiba", "Malpelo",
                                     "Rapa Nui", "Ducie", "Gambier",
                                     "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=40, face="bold"),
              legend.box ='horizontal', legend.position = "top")
      plot_dapc_ile_symb
      
      pdf(file = paste0(outputPrefix, "_Symbiont_Ile.pdf"), w = 20, h=20)
      plot_dapc_ile_symb
      dev.off()
      
      png(file = paste0(outputPrefix, "_Symbiont_Ile.png"),  w = 20, h=20, units="in", res=300)
      plot_dapc_ile_symb
      dev.off()
      
      dapc_clade_symb_DFeig <- data.frame("DF" = seq(from=1, to=length(SymDAPC_clade$eig), by=1),
                                   "DFeigval" = SymDAPC_clade$eig,
                                   "DFkept" = c(rep("Kept",2),rep("Not",length(SymDAPC_clade$eig)-2)))
      dapc_clade_symb_PCeig <- data.frame("PC" = seq(from=1, to=length(SymDAPC_clade$pca.eig), by=1),
                                 "PCcumvar" = 100 * cumsum(SymDAPC_clade$pca.eig)/sum(SymDAPC_clade$pca.eig),
                                 "PCkept" = c(rep("Kept",SymDAPC_clade$n.pca),rep("Not",length(SymDAPC_clade$pca.eig)-SymDAPC_clade$n.pca)))

      write.table(file = paste0(outputPrefix,"_FigureS4g_Data.tab"), dapc_clade_symb_DFeig, quote = F, sep = '\t', row.names = F)
      
      DFbox_clade_symb <- ggplot() +
        geom_bar(data = dapc_clade_symb_DFeig, aes(x=DF,y=DFeigval, fill=DFkept),
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0.5,length(SymDAPC_clade$eig)+0.5), breaks = seq(from=1, to=length(SymDAPC_clade$eig), by=1)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("Discriminant Function") +
        ylab("DF eigenvalues") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      DFbox_clade_symb
      
      PCbox_clade_symb <- ggplot() +
        geom_bar(data = dapc_clade_symb_PCeig, aes(x=PC,y=PCcumvar, fill=PCkept),
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0,length(SymDAPC_clade$pca.eig)+0.5), breaks = seq(from=0, to=length(SymDAPC_clade$pca.eig), by=20)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("PCA Axis") +
        ylab("Cumulative Variance Explained (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      PCbox_clade_symb
      
      write.table(file = paste0(outputPrefix,"_FigureS5c_Data.tab"), gg_clade_symb, quote = F, sep = '\t', row.names = F)

      plot_dapc_clade_symb <- ggplot() +
        stat_ellipse(data=gg_clade_symb, geom="polygon", aes(x=LD1,y=LD2,color=Clade, fill=Clade),
                     type="t", level=0.95, alpha=0.5, show.legend=F,) +
        geom_segment(data=gg_clade_symb, show.legend=F,
                     aes(x=LD1.centroid, y=LD2.centroid, xend=LD1, yend=LD2, color=Clade)) +
        annotate(geom="text", x=-35, y=10,
                 label=paste("Proportion Reassigned /n", VarExp_Sym_mean_clade,  " ? ", VarExp_Sym_sd_clade, sep =' '),
                 color="black", size = 14) +
        scale_color_manual(values=mycols_clade,
                           name = "SVD Clade",
                           labels = c("SVD 1", "SVD 2", "SVD 3",
                                      "SVD 4", "SVD 5", "Unassigned")) +
        scale_fill_manual(values=mycols_clade,
                          name = "SVD Clade",
                          labels = c("SVD 1", "SVD 2", "SVD 3",
                                     "SVD 4", "SVD 5", "Unassigned")) +
        new_scale("fill") +
        geom_point(data=gg_clade_symb, aes(x=LD1,y=LD2,fill=Ile), color="black",
                   size=10, pch=21) +
        scale_fill_manual(values=mycols,
                          name = "Island",
                          labels = c("Las Perlas", "Coiba", "Malpelo",
                                     "Rapa Nui", "Ducie", "Gambier",
                                     "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        new_scale("fill") +
        geom_text_repel(data=centroids_clade_symb, aes(x=LD1, y=LD2, label=Label),
                         alpha=0.5, color="black", size=20, show.legend = FALSE) +
        scale_fill_manual(values=mycols_clade,
                          name = "SVD Clade",
                          labels = c("SVD 1", "SVD 2", "SVD 3",
                                     "SVD 4", "SVD 5", "Unassigned")) +
        xlab("Discriminant Function 1") +
        ylab("Discriminant Function 2") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=40, face="bold"),
              legend.box ='horizontal', legend.position = "top")
      plot_dapc_clade_symb

      pdf(file = paste0(outputPrefix, "_Symbiont_Clade.pdf"), w = 20, h =20)
      plot_dapc_clade_symb
      dev.off()

      png(file = paste0(outputPrefix, "_Symbiont_Clade.png"),  w = 20, h =20, units="in", res=300)
      plot_dapc_clade_symb
      dev.off()

      dapc_symclade_symb_DFeig <- data.frame("DF" = seq(from=1, to=length(SymDAPC_symclade$eig), by=1),
                                   "DFeigval" = SymDAPC_symclade$eig,
                                   "DFkept" = c(rep("Kept",2),rep("Not",length(SymDAPC_symclade$eig)-2)))
      dapc_symclade_symb_PCeig <- data.frame("PC" = seq(from=1, to=length(SymDAPC_symclade$pca.eig), by=1),
                                 "PCcumvar" = 100 * cumsum(SymDAPC_symclade$pca.eig)/sum(SymDAPC_symclade$pca.eig),
                                 "PCkept" = c(rep("Kept",SymDAPC_symclade$n.pca),rep("Not",length(SymDAPC_symclade$pca.eig)-SymDAPC_symclade$n.pca)))

      write.table(file = paste0(outputPrefix,"_FigureS4h_Data.tab"), dapc_symclade_symb_DFeig, quote = F, sep = '\t', row.names = F)
      
      DFbox_symclade_symb <- ggplot() +
        geom_bar(data = dapc_symclade_symb_DFeig, aes(x=DF,y=DFeigval, fill=DFkept),
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0.5,length(SymDAPC_symclade$eig)+0.5), breaks = seq(from=1, to=length(SymDAPC_symclade$eig), by=1)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("Discriminant Function") +
        ylab("DF eigenvalues") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      DFbox_symclade_symb

      PCbox_symclade_symb <- ggplot() +
        geom_bar(data = dapc_symclade_symb_PCeig, aes(x=PC,y=PCcumvar, fill=PCkept),
                 stat="identity", color = "black", show.legend = F) +
        scale_x_continuous(limits = c(0,length(SymDAPC_symclade$pca.eig)+0.5), breaks = seq(from=0, to=length(SymDAPC_symclade$pca.eig), by=20)) +
        scale_fill_manual(values=c("grey30","white")) +
        xlab("PCA Axis") +
        ylab("Cumulative Variance Explained (%)") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=30, face="bold"),
              legend.box ='horizontal')
      PCbox_symclade_symb
      
      write.table(file = paste0(outputPrefix,"_FigureS5b_Data.tab"), gg_symclade_symb, quote = F, sep = '\t', row.names = F)

      plot_dapc_symclade_symb <- ggplot() +
        stat_ellipse(data=gg_symclade_symb, geom="polygon", aes(x=LD1,y=LD2,color=SymClade, fill=SymClade),
                     type="t", level=0.95, alpha=0.5, show.legend=F,) +
        geom_segment(data=gg_symclade_symb, show.legend=F,
                     aes(x=LD1.centroid, y=LD2.centroid, xend=LD1, yend=LD2, color=SymClade)) +
        annotate(geom="text", x=-35, y=10,
                 label=paste("Proportion Reassigned /n", VarExp_Sym_mean_symclade,  " ? ", VarExp_Sym_sd_symclade, sep =' '),
                 color="black", size = 14) +
        scale_color_manual(values=mycols_symclade[c(3:7)],
                           name = "Cladocopium Group",
                           labels = c("L1", "L2", "L3",
                                      "L4", "L5")) +
        scale_fill_manual(values=mycols_symclade[c(3:7)],
                          name = "Cladocopium Group",
                          labels = c("L1", "L2", "L3",
                                     "L4", "L5")) +
        new_scale("fill") +
        geom_point(data=gg_symclade_symb, aes(x=LD1,y=LD2,fill=Ile), color="black",
                   size=10, pch=21) +
        scale_fill_manual(values=mycols,
                          name = "Island",
                          labels = c("Las Perlas", "Coiba", "Malpelo",
                                     "Rapa Nui", "Ducie", "Gambier",
                                     "Moorea","Aitutaki","Niue","Upolu","Guam")) +
        new_scale("fill") +
        geom_text_repel(data=centroids_symclade_symb, aes(x=LD1, y=LD2, label=Label),
                        alpha=0.5, color="black", size=20, show.legend = FALSE) +
        scale_fill_manual(values=mycols_symclade[c(3:7)],
                          name = "Cladocopium Group",
                          labels = c("L1", "L2", "L3",
                                     "L4", "L5")) +
        xlab("Discriminant Function 1") +
        ylab("Discriminant Function 2") +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              text = element_text(size=40, face="bold"),
              legend.box ='horizontal', legend.position = "top")
      plot_dapc_symclade_symb

      pdf(file = paste0(outputPrefix, "_Symbiont_SymClade.pdf"), w = 20, h =20)
      plot_dapc_symclade_symb
      dev.off()

      png(file = paste0(outputPrefix, "_Symbiont_SymClade.png"),  w = 20, h=20, units="in", res=300)
      plot_dapc_symclade_symb
      dev.off()

## 5. HOST GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS OF DAPC DISCRIMINANT GENES -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  library(goseq)
  library(GO.db)
  library(geneLenDataBase)
  library(qvalue)    
    
  # Step 1- Get gene lengths
    gene_lengths <- read.table("../Pocillopora_meandrina_v3.1.annot.mrna_SeqLengths.txt", header=T, row.names=1, com='')
    gene_lengths <- as.matrix(gene_lengths[,1,drop=F])
    
  # Step 2 - Get host background gene list (All genes)
    # cts_background <- read.table("../DiffExpress/Pocillopora_meandrina_v3.1.annot.mrna_RawCounts.tab", header=TRUE, com='', row.names=1, check.names=FALSE)
    Allbackground.gene_ids <- c(gsub(";.*","",rownames(cts_raw_host)))
    
  # Step 3 - Parse GO assignments
    detach("package:xlsx", unload=T)
    library(openxlsx)
    
    GO_dat <- read.xlsx("../Pocillopora_meandrina_v3.2_annot.xlsx", sheet = "Poc_v3.2_annot.ipr", rowNames = F, colNames = F)
      
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
    
  # Step 4 - Perform GOSeq Enrichment Analysis
    library(goseq)
      
    # Create a module that contains all discriminant functions (DF1 and DF2)
      module <- c("LD1","LD2")
     
    # Select the set of DAPC results to analyze
      DAPCGOI <- "HostDAPC_Ile_TopLoadGenes"
      # Options are: "HostDAPC_Ile_TopLoadGenes", "HostDAPC_Clade_TopLoadGenes", "HostDAPC_SymClade_TopLoadGenes"
    
    # Give a file prefix for the output results
      DAPCRes <- "Pocillopora_DAPC_Ile_"
      # Options are: "Pocillopora_DAPC_Ile_", "Pocillopora_DAPC_Clade_", "Pocillopora_DAPC_SymClade_"
      
    # Create a directory for the GOSeq output results
      dir.create("./DAPC_GOSeq/", 
                 showWarnings = TRUE, recursive = FALSE, mode = "0777")
      
    # Identify enriched/depleted GO Terms among genes of interest for each discriminant function  
      for (j in 1:length(module)) {
        # GoSeq Enrichment Analysis
        infile <- paste(DAPCRes, module[j], "_",sep='')
        factor_labeling <- subset(get(DAPCGOI), Loading %in% module[j]) #Replace the dataframe with your genes of interest
        rownames(factor_labeling) <- gsub(";.*","", factor_labeling$geneLabels)
        factor_labeling[,1] <- rep('custom_list', dim(factor_labeling)[1])
        factor_labeling <- factor_labeling[,1,drop=F]
        colnames(factor_labeling) <- c('type')
        factor_list <- unique(factor_labeling[,1])
        select <- grep("Pmea_",rownames(factor_labeling))
        DE_genes <- rownames(factor_labeling)[select]
        #NOTE: Background genes MUST also include DEGs!
        background.gene_ids <- unique(c(Allbackground.gene_ids, DE_genes))
        sample_set_gene_ids <- background.gene_ids
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
          go_enrich_filename = paste("./DAPC_GOSeq/",infile,"GOseq.enriched.txt", sep='')
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
          go_depleted_filename = paste("./DAPC_GOSeq/",infile,"GOseq.depleted.txt", sep='')
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
      
## 6. VISUALIZATION OF DAPC GO ENRICHMENTS (Dotplots) ---------------------------------------------------------------------------------------------------------------
    library(ggplot2)
      
  # Step 1 - Select the GOSeq results of interest and create an input table for dotplot
    GOSeqfile1 <- "Pocillopora_DAPC_Ile_LD1_GOseq.enriched.txt" #Replace with your GO Enrichment Results from above
    EnrichTab1 <- read.table(paste0("./DAPC_GOSeq/",GOSeqfile1,sep=''),
                           header = TRUE, quote="\"",
                           sep=";")
    EnrichTab1$ratio <- EnrichTab1$numDEInCat/EnrichTab1$numInCat
    EnrichTab1$pval <- EnrichTab1$over_represented_pvalue
    EnrichTab1$Loading <- "Discriminant  Function 1"
    
    GOSeqfile2 <- "Pocillopora_DAPC_Ile_LD2_GOseq.enriched.txt" #Replace with your GO Enrichment Results from above
    EnrichTab2 <- read.table(paste0("./DAPC_GOSeq/",GOSeqfile2,sep=''),
                           header = TRUE, quote="\"",
                           sep=";")
    EnrichTab2$ratio <- EnrichTab2$numDEInCat/EnrichTab2$numInCat
    EnrichTab2$pval <- EnrichTab2$over_represented_pvalue
    EnrichTab2$Loading <- "Discriminant Function 2"
    
    EnrichTab <- rbind(EnrichTab1, EnrichTab2)
    
    write.table(file = paste0(outputPrefix,"_Figure4c_Data.tab"), EnrichTab, quote = F, sep = '\t', row.names = F)
    
  # Step 2 - Create function to make Dotplot of top 100 (showCategory = 100) enriched Biological Process categories
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
        scale_size_continuous(range = c(5, 10)) +
        facet_grid(.~Loading,
                   drop = TRUE,
                   scales = "free_x") +
        xlab("") +
        ylab("Count Ratio") +
        scale_x_discrete(expand = c(0.01, 0.2)) +
        #coord_flip() +
        theme_bw(base_size=9) + 
        theme(text = element_text(size = 40, face = 'bold'),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        guides(size = guide_legend(title="DEGs", order = 1))
    }
    
  # Step 3 - Plot results
    
    plotDot <- dotplot_goseq(EnrichTab) 
    plotDot
    
    library(stringr)
    
    # pdf(file = paste0(outputPrefix, "_", substr(paste0(GOSeqfile),1,nchar(paste0(GOSeqfile))-4),"_Dotplot.pdf", sep=''), h = 8.5, w = 20)
    # plotDot
    # dev.off()
    # 
    # png(file = paste0(outputPrefix, "_", substr(paste0(GOSeqfile),1,nchar(paste0(GOSeqfile))-4), "_Dotplot.png", sep=''), h = 8.5, w = 20, units = "in", res = 300)
    # plotDot
    # dev.off()
    
      # Step 3 - Plot results

    pdf(file = paste0(outputPrefix,"_Host_Ile-LD1-LD2_Dotplot.pdf"), h = 20, w = 30)
    plotDot
    dev.off()
    
    png(file = paste0(outputPrefix, "_Host_Ile-LD1-LD2_Dotplot.png"), h = 20, w = 30, units = "in", res = 300)
    plotDot
    dev.off()