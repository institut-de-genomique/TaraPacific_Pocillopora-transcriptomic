## TARA Pacific Project - Pocillopora Holobiont MetaT Analyses ##
#################################################################

# Author: ARMSTRONG Eric
# Created: 21 August 2018
# Last Edited: 9 December 2022

## 1. INITIALIZATION ------------------------------------------------------------------------------------------------------------------------------------

  # Step 1 - Set Defaults

    # Set working directory
      setwd("C:/Users/uax75/OneDrive/Documents/R/TaraCoral_2020/Final_Results/Pocillopora/ManuscriptDocs/Scripts/Figure5_CCA/")
    
    # Set the prefix for each output file name
      outputPrefix <- "TaraPacific_Holobiont_CCA_GitHub_2022-12"


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
    
    
  
  # Step 5 - Subset count data to include only colonies with species assignations
    colonies_keep <- subset(metaData_host, !Clade %in% c("Unassigned"))
    cts_raw_host_keep <- cts_raw_host[,colonies_keep$Colony]

    colonies_keep_sym <- subset(metaData_sym, !Clade %in% c("Unassigned") & !SymClade %in% c("Unassigned") &  !Ile %in% c("I11"))
    cts_raw_sym_keep <- cts_raw_sym[,colonies_keep_sym$Colony]
    
    
  # Step 6 - Add select historical and environmental data to sample info table
    library(dplyr)
    metaData_env <- left_join(metaData_host, envData, by="Colony")
    
    metaData_env_sym <- left_join(metaData_sym, envData, by="Colony")
    
    rownames(metaData_env) <- metaData_env$Colony
    rownames(metaData_env_sym) <- metaData_env_sym$Colony
   
    
## 3. Transform Count Data -------------------------------------------------------------------

  # Step 1 - Prefilter data to remove low count genes
    library("DESeq2")
    # Remove genes with fewer than 10 counts summed across all samples
      # making a vector counting number of samples with counts <=10 within each colony 
    lowcountgenes_host <- apply(cts_raw_host_keep[,1:ncol(cts_raw_host_keep)],1,function(x){sum(x<=10)}) 
    # get rid of genes with counts < 10 in more than 90% of samples
    highcts_raw_host_keep <- cts_raw_host_keep[-which(lowcountgenes_host>(0.9*ncol(cts_raw_host_keep))),]
    
    lowcountgenes_sym <- apply(cts_raw_sym_keep[,1:ncol(cts_raw_sym_keep)],1,function(x){sum(x<=10)}) 
    # get rid of genes with counts < 10 in more than 90% of samples
    highcts_raw_sym_keep <- cts_raw_sym_keep[-which(lowcountgenes_sym>(0.9*ncol(cts_raw_sym_keep))),]
    
    
    nrow(highcts_raw_host_keep) #28730 genes pass filter for host 
    nrow(highcts_raw_sym_keep) #20390 genes pass filer for symbiont
    
  # Step 2 - Normalize expression data
    
    dds_host <- DESeqDataSetFromMatrix(countData = highcts_raw_host_keep, 
                                      colData = metaData_env, 
                                      design = ~ Ile)
    
    dds_sym <- DESeqDataSetFromMatrix(countData = highcts_raw_sym_keep, 
                                      colData = metaData_env_sym, 
                                      design = ~ Ile)
    
    
    # DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
    
    # using vsd transform:
      vsd_mat_host <- varianceStabilizingTransformation(dds_host, blind=TRUE)
      vsd_mat_sym <- varianceStabilizingTransformation(dds_sym, blind=FALSE)
    
    # Extract variance-stabilization-transformed read counts from dds object 
      cts_vsd_host <- assay(vsd_mat_host)
      cts_vsd_sym <- assay(vsd_mat_sym)
      

## 4. Ordination Analysis -------------------------------------------------------------------------------------------------------------------------------------------------------------

  # Step 1 - Identify environmental variables that are highly Correlated with on another
    # Load environmental data
      data_env <- metaData_env[,c(14,16:ncol(metaData_env))]
      data_env_sym <- metaData_env_sym[,c(14,16:ncol(metaData_env_sym))]

    # Remove colonies that still have missing environmental data
      data_env_cca <- data_env[complete.cases(data_env),]
      data_env_sym_cca <- data_env_sym[complete.cases(data_env_sym), ]
        
       # Flip count data frame for ordination analysis 
        data_cts_vsd_env <- data.frame(t(cts_vsd_host))
        data_cts_sym_vsd_env <- data.frame(t(cts_vsd_sym))
        
      CommonColonies <- intersect(rownames(data_cts_vsd_env), rownames(data_env_cca))
      data_env_cca <- data_env_cca[CommonColonies,]
      data_cts_vsd_env <- data_cts_vsd_env[CommonColonies,]
      
      CommonColonies_sym <- intersect(rownames(data_cts_sym_vsd_env), rownames(data_env_sym_cca))
      data_env_sym_cca <- data_env_sym_cca[CommonColonies_sym,]
      data_cts_sym_vsd_env <- data_cts_sym_vsd_env[CommonColonies_sym,]

      
    # Remove environmental variables missing data
      data_env_cca <- data.frame(sapply(data_env_cca, as.numeric))
      data_env_sym_cca <- data.frame(sapply(data_env_sym_cca, as.numeric))
      data_env_cca[data_env_cca == "NaN"] <- NA
      data_env_cca <- data_env_cca[ , colSums(is.na(data_env_cca))==0]
      data_env_sym_cca[data_env_sym_cca == "NaN"] <- NA
      data_env_sym_cca <- data_env_sym_cca[ , colSums(is.na(data_env_sym_cca))==0]
      
    # Remove environmental variables that show no variation across islands (uninformative)
      data_env_cca <- data_env_cca[, !sapply(data_env_cca, function(x) { sd(x) == 0} )]
      data_env_sym_cca <- data_env_sym_cca[, !sapply(data_env_sym_cca, function(x) { sd(x) == 0} )]
      rownames(data_env_cca) <- rownames(data_cts_vsd_env)
      rownames(data_env_sym_cca) <- rownames(data_cts_sym_vsd_env)
      
    
    # Add genotype data to environmental variables
      data_env_cca$Clade <- as.factor(vlookup(lookup_value = rownames(data_env_cca),
                                                      dict = CladeTable,
                                                      lookup_column = "Samples",
                                                      result_column = "PocilloGG"))  
      data_env_sym_cca$SymClade <- as.factor(vlookup(lookup_value = rownames(data_env_sym_cca),
                                                             dict = CladeTable,
                                                             lookup_column = "Samples",
                                                             result_column = "SymbioGG"))
      
    # Identify clusters of varibles that are co-correlated at R2 > 0.7
      library(klaR)
      ccres <- corclust(x = data_env[,intersect(colnames(data_env_cca),colnames(data_env_sym_cca))])
      
      pdf(file = paste0(outputPrefix, "_Both-CorrEnvVariables.pdf"), h = 8, w = 16)
      plot(ccres, mincor = 0.7)
      dev.off()
      
      png(file = paste0(outputPrefix, "_Both-CorrEnvVariables.png"), h = 8, w = 16, units = "in", res = 300)
      plot(ccres, mincor = 0.7)
      dev.off()
      
    # Choose one representative variable to keep for each correlation cluster
      envkept <- c("TSA_DCW_lastrecovery_day","TSA_DCW_minlenght_day",
                        "TSA_DCW_mean_DegC","SST_anomaly_max_DegC",
                        "TSA_DCW_meanrecovery_day","PAR_Sat",
                        "TSA_DCW_lastduration_day","merged_pH",
                        "mergedPO4_mole_L","TSA_DHW_maxrecovery_day",
                        "TSA_DHW_meanrecovery_day","TSA_DCW_freq_snapshot_sampling_day_DegC",
                        "TSA_DHW_lastmax_DegC","TSA_DHW_mean_DegC",
                        "TSA_heat_max_DegC","TSA_heat_freq_snapshot_sampling_day_day",
                        "TSA_DHW_max_DegC","TSA_cold_freq_mean_day",
                        "merged_chl","SST_anomaly_mean_DegC",
                        "SST_max_DegC","TSA_DCW_stdlenght_day",
                        "SST_min_DegC","TSA_DCW_snapshot_sampling_day_DegC",
                        "SST_anomaly_snapshot_sampling_day_DegC")
      
  # Step 2 - Identify environmental variables that are correlated with expression data
    set.seed(1)
    library(vegan)
    
    Poc.cca <- cca(data_cts_vsd_env ~ ., data_env_cca[,envkept], na.action= "na.omit")
    host_env_fit <- envfit(Poc.cca, data_env_cca[,envkept], permutations = 10000, na.rm = TRUE)
    
    Sym.cca <- cca(data_cts_sym_vsd_env ~ ., data_env_sym_cca[,envkept], na.action = "na.omit")
    sym_env_fit <- envfit(Sym.cca, data_env_sym_cca[,envkept], permutations = 10000, na.rm = TRUE)

    # Select all significant environmental variables at p-adj <= 0.05
      host_top_env <- names(host_env_fit$vectors$pvals[host_env_fit$vectors$pvals <= 0.05])
      sym_top_env <- names(sym_env_fit$vectors$pvals[sym_env_fit$vectors$pvals <= 0.05])

  # Step 3 - Automatic step-wise selection of variables that contribute significantly and non-redundantly to CCA 
    # Intercept only model for the host
      Poc.cca.mod0 <- cca(data_cts_vsd_env ~ 1, data_env_cca[,host_top_env], na.action= "na.omit")
    # Full model with all variables for the host
      Poc.cca.mod1 <- cca(data_cts_vsd_env ~ ., data_env_cca[,host_top_env], na.action= "na.omit")
    
    # Intercept only model for the photosymbiont
      Sym.cca.mod0 <- cca(data_cts_sym_vsd_env ~ 1, data_env_sym_cca[,sym_top_env], na.action = "na.omit")
    # Full model with all variable for the photosymbiont
      Sym.cca.mod1 <- cca(data_cts_sym_vsd_env ~ ., data_env_sym_cca[,sym_top_env], na.action = "na.omit")

    # Perform automatic step-wise analysis
      Poc.step.res <- ordiR2step(Poc.cca.mod0, Poc.cca.mod1, perm.max = 1000)
      Poc.step.res$anova  # Summary table
      Poc.step.res_adj <- Poc.step.res
      Poc.step.res_adj$anova$`Pr(>F)` <- p.adjust(Poc.step.res$anova$`Pr(>F)`, method = 'holm', n = ncol (data_env_cca[,host_top_env]))
      write.table(file="POC-Step-Results.tab", Poc.step.res_adj$anova, sep="\t", quote=F)
      
      Sym.step.res <- ordiR2step(Sym.cca.mod0, Sym.cca.mod1, perm.max = 1000)
      Sym.step.res$anova  # Summary table
      Sym.step.res_adj <- Sym.step.res
      Sym.step.res_adj$anova$`Pr(>F)` <- p.adjust(Sym.step.res$anova$`Pr(>F)`, method = 'holm', n = ncol (data_env_sym_cca[,sym_top_env]))
      write.table(file="SYM-Step-Results.tab", Sym.step.res_adj$anova, sep="\t", quote=F)


  # Step 4 - Ordination with top environmental variables (CCA)

      # Find Top Env Variables w/Ordination
        library("vegan")
        library("ggvegan")
        library("MASS")
      
      # Select all significant env variables at p-adj <= 0.05
      host_top_env2 <- c("SST_anomaly_mean_DegC","SST_min_DegC","TSA_DCW_lastrecovery_day",
                         "TSA_DCW_freq_snapshot_sampling_day_DegC","PAR_Sat","merged_chl",
                         "SST_max_DegC","TSA_DHW_mean_DegC","TSA_DCW_meanrecovery_day",
                         "TSA_DHW_lastmax_DegC","TSA_DHW_maxrecovery_day")
      sym_top_env2 <- c("SST_max_DegC","TSA_heat_freq_snapshot_sampling_day_day","merged_SiOH_mole_L",
                        "TSA_DCW_freq_snapshot_sampling_day_DegC","SST_anomaly_snapshot_sampling_day_DegC",
                        "PAR_Sat","TSA_cold_freq_snapshot_sampling_day_day","TSA_heat_max_DegC")
        
  # Step 2 - CCA using reduced set of environmental data
      data_env_cca2 <- data_env_cca[,host_top_env2]
      data_env_sym_cca2 <- data_env_sym_cca[,sym_top_env2]
      
          
          Poc.cca2 <- cca(data_cts_vsd_env, data_env_cca2, na.action= "na.omit")
          Sym.cca2 <- cca(data_cts_sym_vsd_env, data_env_sym_cca2, na.action = "na.omit")
          
          
        
      
      
##### V - Plot #####   
####################
      

  # Step 2 - Visualization     
    library("RColorBrewer")
    library("ggbiplot")
    
    ##### Note - Modification to ggbiplot function #####
    # Modified base ggbiplot function so that it only plots the top 10 loadings (arrows) on the PCA.
    # fix(ggbiplot)
    # somewhere around line 89
      # g <- g + geom_segment(data = df.v[1:10,]
    # somewhere around lline 127
      # g <- g + geom_segment(data = df.v[1:10,]
    ##### End Note #####
    ####################

  
      ford_POC <- fortify(Poc.cca2, axes= 1:2) # fortify the ordination
      ford_SYM <- fortify(Sym.cca2, axes= 1:2) # fortify the ordination
      
      take_POC <- c('CCA1', 'CCA2')  # which columns contain the scores we want
      take_SYM <- c('CCA1', 'CCA2')  # which columns contain the scores we want
      
      arrows_POC <- subset(ford_POC, Score == 'biplot')  # take only biplot arrow scores
      arrows_SYM <- subset(ford_SYM, Score == 'biplot')  # take only biplot arrow scores
      
        ## multiplier for arrows to scale them to the plot range
        mul <- ggvegan:::arrowMul(arrows_POC[, take_POC],
                                  subset(ford_POC, select = take_POC, Score == 'sites'))
        arrows_POC[, take_POC] <- arrows_POC[, take_POC] * mul  # scale biplot arrows
        
        mul_sym <- ggvegan:::arrowMul(arrows_SYM[, take_SYM],
                                  subset(ford_SYM, select = take_SYM, Score == 'sites'))
        arrows_SYM[, take_SYM] <- arrows_SYM[, take_SYM] * mul_sym  # scale biplot arrows
        
        
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
      
      labelnames <- c(arrows_POC$Label,arrows_SYM$Label)
      arrowlabels <- data.frame("Label" = c(unique(labelnames)),
                                "Number" = c(1:length(unique(labelnames))))
        
      library("ggnewscale")
      library("expss")
        
        plot3 <- ggplot() +
          stat_ellipse(data=subset(ford_POC, Score == 'sites'), geom="polygon", 
                       aes(x=CCA1,y=CCA2,
                           color=factor(vlookup(lookup_value = rownames(data_cts_vsd_env),
                                                dict = metaData_host,
                                                lookup_column = "Colony",
                                                result_column = "Clade")), 
                           fill=factor(vlookup(lookup_value = rownames(data_cts_vsd_env),
                                                dict = metaData_host,
                                                lookup_column = "Colony",
                                                result_column = "Clade"))), 
                   type="t", level=0.68, alpha=0.5, show.legend=F,) +
          scale_fill_manual(name='', values=mycols_clade,
                            labels=c("SVD 1", "SVD 2","SVD 3",
                                    "SVD 4", "SVD 5")) +
          scale_color_manual(name='', values=mycols_clade,
                            labels=c("SVD 1", "SVD 2","SVD 3", 
                                    "SVD 4", "SVD 5")) +
          new_scale("fill") +
          new_scale("color") +
          geom_point(data = subset(ford_POC, Score == 'sites'),
                     mapping = aes(x = CCA1, y = CCA2,
                                   shape=factor(vlookup(lookup_value = rownames(data_cts_vsd_env),
                                                dict = metaData_host,
                                                lookup_column = "Colony",
                                                result_column = "Clade")), 
                                   fill=factor(vlookup(lookup_value = rownames(data_cts_vsd_env),
                                                dict = metaData_host,
                                                lookup_column = "Colony",
                                                result_column = "Ile")), size=3)) +
          geom_segment(data = arrows_POC,
                       mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
                       arrow = arrow(length = unit(0.01, "npc"))) +
          geom_text(data = arrows_POC, # crudely push labels away arrow heads
                    mapping = aes(label = vlookup(lookup_value = arrows_POC$Label,
                                                  dict = arrowlabels,
                                                  lookup_column = "Label",
                                                  result_column = "Number"), x = CCA1 * 1.1, y = CCA2 * 1.1)) +
          scale_fill_manual(name='', values = mycols,
                            labels = c("Las Perlas", "Coiba", "Malpelo",
                                       "Rapa Nui", "Ducie", "Gambier", 
                                       "Moorea", "Aitutaki", "Niue",
                                        "Upolu", "Guam")) +
          scale_shape_manual(name = '', values=c(23,21,24,22,25,1), #16,17,18,15,4
                           labels=c("P. effusa", "P. meandrina","P. verrucosa", 
                                    "P. grandis", "SSH5_pver","Unknown")) +
          # scale_x_continuous(limits=c(-3.2,3.2), breaks=seq(from=-3, to=3, by=1), labels=seq(from=-3, to=3, by=1)) +
          # scale_y_continuous(limits=c(-3.2,3.2), breaks=seq(from=-3, to=3, by=1), labels=seq(from=-3, to=3, by=1)) +
          theme_bw() +
          theme(text = element_text(size=20, face="bold")) +
          #scale_x_continuous(limits = c(-2,5), breaks = seq(from=-2, to=3, by=1)) +
          #scale_y_continuous(limits = c(-2.5,2.52), breaks = seq(from=-2, to=2, by=1)) +
          scale_size(guide="none") +
          # guides(fill=F) +
          guides(fill=guide_legend(override.aes=list(color=mycols, size=4))) +
          guides(shape=guide_legend(override.aes=list(fill=c(mycols_clade,"white"), size=4))) +
          theme(legend.direction = 'horizontal',
              legend.position = 'bottom',
              legend.text = element_text(size = 14)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          coord_fixed()
        plot3
        
        plot3a <- ggplot() +
          stat_ellipse(data=subset(ford_POC, Score == 'sites'), geom="polygon", 
                       aes(x=CCA1,y=CCA2,
                           color=factor(vlookup(lookup_value = rownames(data_cts_vsd_env),
                                                dict = metaData_host,
                                                lookup_column = "Colony",
                                                result_column = "Clade")),
                           fill=factor(vlookup(lookup_value = rownames(data_cts_vsd_env),
                                                dict = metaData_host,
                                                lookup_column = "Colony",
                                                result_column = "Clade"))), 
                   type="t", level=0.68, alpha=0.5, show.legend=F,) +
          # scale_color_manual(name='', values = mycols, labels = c("Las Perlas", "Coiba", "Malpelo", "Rapa Nui",
          #                                                         "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
          #                                                         "Upolu", "Guam")) +
          # scale_shape_manual(name = '', values=c(16,17,15,4), labels=c("Site 01", "Site 02",
          #                                                              "Site 03", "Site 04")) +
          scale_fill_manual(name='', values=mycols_clade,
                            labels=c("SVD 1", "SVD 2","SVD 3",
                                    "SVD 4", "SVD 5")) +
          scale_color_manual(name='', values=mycols_clade,
                            labels=c("SVD 1", "SVD 2","SVD 3", 
                                    "SVD 4", "SVD 5")) +
          new_scale("fill") +
          new_scale("color") +
          geom_point(data = subset(ford_POC, Score == 'sites'),
                     mapping = aes(x = CCA1, y = CCA2,
                                   shape=factor(vlookup(lookup_value = rownames(data_cts_vsd_env),
                                                dict = metaData_host,
                                                lookup_column = "Colony",
                                                result_column = "Clade")),
                                   fill=factor(vlookup(lookup_value = rownames(data_cts_vsd_env),
                                                dict = metaData_host,
                                                lookup_column = "Colony",
                                                result_column = "Ile")), size=3)) +
          geom_segment(data = arrows_POC,
                       mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
                       arrow = arrow(length = unit(0.01, "npc"))) +
          geom_text(data = arrows_POC, # crudely push labels away arrow heads
                    mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1)) +
          scale_fill_manual(name='', values = mycols,
                            labels = c("Las Perlas", "Coiba", "Malpelo","Rapa Nui",
                                       "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
                                       "Upolu", "Guam")) +
          # scale_fill_manual(name='', values = mycols,
          #                   labels = c("Las Perlas", "Coiba", "Malpelo", "Rapa Nui",
          #                              "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
          #                               "Upolu", "Guam")) +
          scale_shape_manual(name = '', values=c(23,21,24,22,25,1), #16,17,18,15,4
                           labels=c("P. effusa", "P. meandrina","P. verrucosa", 
                                    "P. grandis", "SSH5_pver","Unknown")) +
          # scale_x_continuous(limits=c(-3.2,3.2), breaks=seq(from=-3, to=3, by=1), labels=seq(from=-3, to=3, by=1)) +
          # scale_y_continuous(limits=c(-3.2,3.2), breaks=seq(from=-3, to=3, by=1), labels=seq(from=-3, to=3, by=1)) +
          theme_bw() +
          theme(text = element_text(size=20, face="bold")) +
          scale_size(guide="none") +
          guides(fill=F) +
          guides(shape=guide_legend(override.aes=list(fill=c(mycols_clade,"white"), size=4))) +
          # guides(fill=guide_legend(override.aes=list(color=mycols[c(1:2,4:11)], size=2))) +
          # guides(fill=guide_legend(override.aes=list(color=mycols, size=2))) +
          # guides(shape=guide_legend(override.aes=list(fill=mycols2, size=2))) +
          theme(legend.direction = 'horizontal',
              legend.position = 'bottom',
              legend.text = element_text(size = 14)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          coord_fixed()
        plot3a
        
        
        plot3sym <- ggplot() +
          stat_ellipse(data=subset(ford_SYM, Score == 'sites'), geom="polygon", 
                       aes(x=CCA1,y=CCA2,
                           color=factor(vlookup(lookup_value = rownames(data_cts_sym_vsd_env),
                                                dict = metaData_sym,
                                                lookup_column = "Colony",
                                                result_column = "SymClade")),
                           fill=factor(vlookup(lookup_value = rownames(data_cts_sym_vsd_env),
                                                dict = metaData_sym,
                                                lookup_column = "Colony",
                                                result_column = "SymClade"))), 
                   type="t", level=0.68, alpha=0.5, show.legend=F,) +
          # scale_color_manual(name='', values = mycols, labels = c("Las Perlas", "Coiba", "Malpelo", "Rapa Nui",
          #                                                         "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
          #                                                         "Upolu", "Guam")) +
          # scale_shape_manual(name = '', values=c(16,17,15,4), labels=c("Site 01", "Site 02",
          #                                                              "Site 03", "Site 04")) +
          scale_fill_manual(name='', values=mycols_symclade,
                            labels=c("Group 1", "Group 2","Group 3",
                                    "Group 4", "Group 5", "Unassigned")) +
          scale_color_manual(name='', values=mycols_symclade,
                            labels=c("Group 1", "Group 2","Group 3", 
                                    "Group 4", "Group 5", "Unassigned")) +
          new_scale("fill") +
          new_scale("color") +
          geom_point(data = subset(ford_SYM, Score == 'sites'),
                     mapping = aes(x = CCA1, y = CCA2,
                                   shape=factor(vlookup(lookup_value = rownames(data_cts_sym_vsd_env),
                                                dict = metaData_sym,
                                                lookup_column = "Colony",
                                                result_column = "SymClade")),
                                   fill=factor(vlookup(lookup_value = rownames(data_cts_sym_vsd_env),
                                                dict = metaData_sym,
                                                lookup_column = "Colony",
                                                result_column = "Ile")), size=3)) +
          geom_segment(data = arrows_SYM,
                       mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
                       arrow = arrow(length = unit(0.01, "npc"))) +
          geom_text(data = arrows_SYM, # crudely push labels away arrow heads
                    mapping = aes(#label = vlookup(lookup_value = arrows_SYM$Label,
                                                  #dict = arrowlabels,
                                                  #lookup_column = "Label",
                                                  #result_column = "Number"),
                                  label = Label,
                                  x = CCA1 * 1.1, y = CCA2 * 1.1)) +
          scale_fill_manual(name='', values = mycols,
                            labels = c("Las Perlas", "Coiba", "Malpelo","Rapa Nui",
                                       "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
                                       "Upolu", "Guam")) +
          scale_shape_manual(name = '', values=c(23,21,24,22,25), #16,17,18,15,4
                           labels=c("C. goreaui", "C. latusorum (L2)","C. latusorum (L3)", 
                                    "C. pacificum (L4)", "C. pacificum (L5)", "Unassigned")) +
          scale_x_continuous(limits=c(-2,3), breaks=seq(from=-2, to=3, by=1), labels=seq(from=-2, to=3, by=1)) +
          # scale_y_continuous(limits=c(-3.2,1.5), breaks=seq(from=-3, to=1, by=1), labels=seq(from=-3, to=1, by=1)) +
          theme_bw() +
          theme(text = element_text(size=20, face="bold")) +
          scale_size(guide="none") +
          guides(fill=F) +
          guides(shape=guide_legend(override.aes=list(fill=mycols_symclade[c(1:5)], size=4))) +
          # guides(fill=guide_legend(override.aes=list(color=mycols, size=2))) +
          # guides(shape=guide_legend(override.aes=list(fill=mycols2, size=2))) +
          theme(legend.direction = 'horizontal',
              legend.position = 'bottom',
              legend.text = element_text(size = 14)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          coord_fixed()
        plot3sym
        
        library(cowplot)
        
        legend_ile <- get_legend(
          # create some space to the left of the legend
          plot3 + guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
            guides(shape=F) +
            theme(legend.position = "bottom")
        )
        
        plot3_both <-  plot_grid(
          plot3a,
          plot3sym,
          align = 'vh',
          labels = c("A", "B"),
          rel_widths = c(1,1),
          hjust = -1,
          nrow = 1)
        plot3_both
        
        pdf(file = paste0(outputPrefix, "_Both_CCA_TopEnv_HistCurr_SiteAvg_biplot.pdf"), w = 14, h = 8.5)
        plot_grid(plot3_both, legend_ile, nrow=2, rel_heights = c(1,0.1), labels=c("",""))
        dev.off()
        
        png(file = paste0(outputPrefix, "_Both_CCA_TopEnv_HistCurr_SiteAvg_biplot.png"), w = 14, h = 8.5, units="in", res=300)
        plot_grid(plot3_both, legend_ile, nrow=2, rel_heights = c(1,0.1), labels=c("",""))
        dev.off()
        
        
        
        
        plot3b <- ggplot() +
          stat_ellipse(data=subset(ford_POC, Score == 'sites'), geom="polygon", 
                       aes(x=CCA1,y=CCA2,color=factor(metaData_env_keep_west$Clade), fill=factor(metaData_env_keep_west$Clade)), 
                   type="t", level=0.68, alpha=0.5, show.legend=F,) +
          # scale_color_manual(name='', values = mycols, labels = c("Las Perlas", "Coiba", "Malpelo", "Rapa Nui",
          #                                                         "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
          #                                                         "Upolu", "Guam")) +
          # scale_shape_manual(name = '', values=c(16,17,15,4), labels=c("Site 01", "Site 02",
          #                                                              "Site 03", "Site 04")) +
          scale_fill_manual(name='', values=mycols2[2:3],
                            labels=c("SVD 2","SVD 3")) +
          scale_color_manual(name='', values=mycols2[2:3],
                            labels=c("SVD 2","SVD 3")) +
          new_scale("fill") +
          new_scale("color") +
          geom_point(data = subset(ford_POC, Score == 'sites'),
                     mapping = aes(x = CCA1, y = CCA2,
                                   shape=factor(metaData_env_keep_west$Clade),
                                   fill=factor(metaData_env_keep_west$Ile), size=3)) +
          geom_segment(data = arrows_POC,
                       mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
                       arrow = arrow(length = unit(0.01, "npc"))) +
          geom_text(data = arrows_POC, # crudely push labels away arrow heads
                    mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1)) +
          scale_fill_manual(name='', values = mycols[c(8:11)],
                            labels = c("Aitutaki", "Niue",
                                       "Upolu", "Guam")) +
          # scale_fill_manual(name='', values = mycols,
          #                   labels = c("Las Perlas", "Coiba", "Malpelo", "Rapa Nui",
          #                              "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
          #                               "Upolu", "Guam")) +
          scale_shape_manual(name = '', values=c(21,24), #16,17,18,15,4
                           labels=c("SVD 2","SVD 3")) +
          # scale_x_continuous(limits=c(-3.2,3.2), breaks=seq(from=-3, to=3, by=1), labels=seq(from=-3, to=3, by=1)) +
          # scale_y_continuous(limits=c(-3.2,3.2), breaks=seq(from=-3, to=3, by=1), labels=seq(from=-3, to=3, by=1)) +
          theme_bw() +
          theme(text = element_text(size=20, face="bold")) +
          scale_size(guide="none") +
          guides(fill=F) +
          guides(shape=guide_legend(override.aes=list(fill=mycols2[c(2:3)], size=2))) +
          # guides(fill=guide_legend(override.aes=list(color=mycols[c(8:11)], size=2))) +
          # guides(fill=guide_legend(override.aes=list(color=mycols, size=2))) +
          # guides(shape=guide_legend(override.aes=list(fill=mycols2, size=2))) +
          theme(legend.direction = 'horizontal',
              legend.position = 'bottom',
              legend.text = element_text(size = 14)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          coord_fixed()
        plot3b
        
        plot3bsym <- ggplot() +
          stat_ellipse(data=subset(ford_SYM, Score == 'sites'), geom="polygon", 
                       aes(x=CCA1,y=CCA2,color=factor(metaData_env_sym_keep_west$SymClade), fill=factor(metaData_env_sym_keep_west$SymClade)), 
                   type="t", level=0.68, alpha=0.5, show.legend=F,) +
          # scale_color_manual(name='', values = mycols, labels = c("Las Perlas", "Coiba", "Malpelo", "Rapa Nui",
          #                                                         "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
          #                                                         "Upolu", "Guam")) +
          # scale_shape_manual(name = '', values=c(16,17,15,4), labels=c("Site 01", "Site 02",
          #                                                              "Site 03", "Site 04")) +
          scale_fill_manual(name='', values=mycols3[c(3,5)],
                            labels=c("Group 3","Group 5")) +
          scale_color_manual(name='', values=mycols3[c(3,5)],
                            labels=c("Group 3","Group 5")) +
          new_scale("fill") +
          new_scale("color") +
          geom_point(data = subset(ford_SYM, Score == 'sites'),
                     mapping = aes(x = CCA1, y = CCA2,
                                   shape=factor(metaData_env_sym_keep_west$SymClade),
                                   fill=factor(metaData_env_sym_keep_west$Ile), size=3)) +
          geom_segment(data = arrows_SYM,
                       mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
                       arrow = arrow(length = unit(0.01, "npc"))) +
          geom_text(data = arrows_SYM, # crudely push labels away arrow heads
                    mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1)) +
          scale_fill_manual(name='', values = mycols[c(8:11)],
                            labels = c("Aitutaki", "Niue",
                                       "Upolu", "Guam")) +
          # scale_fill_manual(name='', values = mycols,
          #                   labels = c("Las Perlas", "Coiba", "Malpelo", "Rapa Nui",
          #                              "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
          #                               "Upolu", "Guam")) +
          scale_shape_manual(name = '', values=c(24,25), #16,17,18,15,4
                           labels=c("Group 3","Group 5")) +
          # scale_x_continuous(limits=c(-3.2,1.5), breaks=seq(from=-3, to=1, by=1), labels=seq(from=-3, to=1, by=1)) +
          # scale_y_continuous(limits=c(-3.2,1.5), breaks=seq(from=-3, to=1, by=1), labels=seq(from=-3, to=1, by=1)) +
          theme_bw() +
          theme(text = element_text(size=20, face="bold")) +
          scale_size(guide="none") +
          guides(fill=F) +
          guides(shape=guide_legend(override.aes=list(fill=mycols3[c(2,5)], size=2))) +
          # guides(fill=guide_legend(override.aes=list(color=mycols, size=2))) +
          # guides(shape=guide_legend(override.aes=list(fill=mycols2, size=2))) +
          theme(legend.direction = 'horizontal',
              legend.position = 'bottom',
              legend.text = element_text(size = 14)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          coord_fixed()
        plot3bsym
        
        
        library(cowplot)
        
        legend_ileb <- get_legend(
          # create some space to the left of the legend
          plot3b + guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
            guides(shape=F) +
            theme(legend.position = "bottom")
        )
        
        plot3b_both <-  plot_grid(
          plot3b,
          plot3bsym,
          align = 'vh',
          labels = c("A", "B"),
          rel_widths = c(1,1),
          hjust = -1,
          nrow = 1)
        plot3b_both
        
        pdf(file = paste0(outputPrefix, "_Both_CCA_TopEnv_HistCurr_SiteAvg_WestPac_biplot.pdf"), w = 14, h = 8.5)
        plot_grid(plot3b_both, legend_ileb, nrow=2, rel_heights = c(1,0.1), labels=c("",""))
        dev.off()
        
        png(file = paste0(outputPrefix, "_Both_CCA_TopEnv_HistCurr_SiteAvg_WestPac_biplot.png"), w = 14, h = 8.5, units="in", res=300)
        plot_grid(plot3b_both, legend_ileb, nrow=2, rel_heights = c(1,0.1), labels=c("",""))
        dev.off()


# Plot PCA by SVD Clade

# Choose colors
mycols2 <- c("tomato2","azure2","turquoise4", "palegreen", "darkorchid")
mycols3 <- c(mycols2,"grey60")

plot4 <- ggbiplot(pca.POC, choices = c(1,2), obs.scale = 1, var.scale = 1, groups = metaData_keep$Clade, 
                  ellipse = TRUE, ellipse.prob = 0.68, var.axes = FALSE) +
  geom_point(aes(shape=factor(metaData_keep$Site), color=factor(metaData_keep$Clade)), size=4) +
  #ggtitle("Pocillopora Principle Component Anlaysis") +
  scale_shape_manual(name = '', values=c(16,17,15,4)) +
  scale_color_manual(name='', values= mycols2) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme_bw() +
  theme(text = element_text(size=16, face="bold")) +
  theme(legend.direction = 'horizontal', 
        legend.position = 'top') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot4

plot5 <- ggbiplot(pca.POC, choices = c(2,3), obs.scale = 1, var.scale = 1, groups = metaData_keep$Clade, 
                  ellipse = TRUE, ellipse.prob = 0.68, var.axes = FALSE) +
  geom_point(aes(shape=factor(metaData_keep$Site), color=factor(metaData_keep$Clade)), size=4) +
  #ggtitle("Pocillopora Principle Component Anlaysis") +
  scale_shape_manual(name = '', values=c(16,17,15,4)) +
  scale_color_manual(name='', values= mycols2) +
  scale_x_continuous(limits=c(-80,130), breaks=c(-50,0,50,100), labels=c(-50,0,50,100)) +
  scale_y_continuous(limits=c(-80,90), breaks=c(-50,0,50), labels=c(-50,0,50)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme_bw() +
  theme(text = element_text(size=16, face="bold")) +
  #theme(legend.direction = 'horizontal', 
  #      legend.position = 'top') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot5

plot6 <-  plot_grid(
  plot4 + theme(legend.position="none"),
  plot5 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)
plot6


legend2 <- get_legend(
  # create some space to the left of the legend
  plot4 + guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    theme(legend.position = "right")
)

pdf(file = paste0(outputPrefix, "_SVDClade_PCA_biplot.pdf"), w = 12, h = 8.5)
plot_grid(plot6, legend2, ncol = 1, rel_heights = c(1, .1))
dev.off()

png(file = paste0(outputPrefix, "_SVDClade_PCA_biplot.png"), w = 12, h = 8.5, units="in", res=300)
plot_grid(plot6, legend2, ncol = 1, rel_heights = c(1, .1))
dev.off()


# Both hosts ile and clade
legend <- get_legend(
  # create some space to the left of the legend
  plot1 + guides(color = guide_legend(ncol = 3), shape = guide_legend(ncol = 2)) +
    theme(legend.position = "top", legend.direction ='horizontal')
)

legend2 <- get_legend(
  # create some space to the left of the legend
  plot4 + guides(color = guide_legend(ncol = 3), shape = guide_legend(ncol = 2)) +
    theme(legend.position = "top", legend.direction ='horizontal')
)

plot_6A_col1 <- plot_grid(plot1 + theme(legend.position="none"),
                          plot2 + theme(legend.position="none"),
                          legend,
                          align='vh', labels=c("(A)","",""),
                          ncol = 1, rel_heights = c(10,10,3))
plot_6A_col1

plot_6A_col2 <- plot_grid(plot4 + theme(legend.position="none"),
                          plot5 + theme(legend.position="none"),
                          legend2,
                          align='vh', labels=c("(B)","",""),
                          ncol = 1, rel_heights = c(10,10,3))
plot_6A_col2

plot6A <-  plot_grid(
  plot_6A_col1,
  plot_6A_col2,
  align = 'vh',
  labels = c("", ""),
  hjust = -1,
  nrow = 1
)
plot6A


pdf(file = paste0(outputPrefix, "_PCA_biplot_Both.pdf"), w = 12, h = 8.5)
plot6A
dev.off()

png(file = paste0(outputPrefix, "_PCA_biplot_Both.png"), w = 12, h = 8.5, units="in", res=300)
plot6A
dev.off()


plot4SymC <- ggbiplot(pca.SYM, choices = c(1,2), obs.scale = 1, var.scale = 1, groups = metaData_keep_sym$SymClade, 
                  ellipse = TRUE, ellipse.prob = 0.68, var.axes = FALSE) +
  geom_point(aes(shape=factor(metaData_keep_sym$Site), color=factor(metaData_keep_sym$SymClade)), size=4) +
  #ggtitle("Pocillopora Principle Component Anlaysis") +
  scale_shape_manual(name = '', values=c(16,17,15,4)) +
  scale_color_manual(name='', values= mycols3) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme_bw() +
  theme(text = element_text(size=16, face="bold")) +
  scale_x_continuous(limits=c(-60,210), breaks=c(-50,0,50,100,150,200), labels=c(-50,0,50,100,150,200)) +
  scale_y_continuous(limits=c(-110,110), breaks=c(-100,-50,0,50,100), labels=c(-100,-50,0,50,100)) +
  #theme(legend.direction = 'horizontal', 
  #      legend.position = 'none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot4SymC

plot5SymC <- ggbiplot(pca.SYM, choices = c(2,3), obs.scale = 1, var.scale = 1, groups = metaData_keep_sym$Clade, 
                  ellipse = TRUE, ellipse.prob = 0.68, var.axes = FALSE) +
  geom_point(aes(shape=factor(metaData_keep_sym$Site), color=factor(metaData_keep_sym$Clade)), size=4) +
  #ggtitle("Pocillopora Principle Component Anlaysis") +
  scale_shape_manual(name = '', values=c(16,17,15,4)) +
  scale_color_manual(name='', values= mycols2) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme_bw() +
  theme(text = element_text(size=16, face="bold")) +
  scale_x_continuous(limits=c(-60,210), breaks=c(-50,0,50,100,150,200), labels=c(-50,0,50,100,150,200)) +
  scale_y_continuous(limits=c(-110,110), breaks=c(-100,-50,0,50,100), labels=c(-100,-50,0,50,100)) +
  #theme(legend.direction = 'horizontal', 
  #      legend.position = 'none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot5SymC


pdf(file = paste0(outputPrefix, "_SymC_PCA_biplot.pdf"), w = 12, h = 8.5)
plot_grid(plot1SymC, plot4SymC, ncol = 2)
dev.off()

png(file = paste0(outputPrefix, "_SymC_PCA_biplot.png"), w = 12, h = 8.5, units="in", res=300)
plot_grid(plot1SymC, plot4SymC, ncol = 2)
dev.off()



## Figure for Abstract


plot_Ile_Abs <-  plot_grid(
  plot2 + theme(legend.position="none"),
  #NULL,
  plot1SymC,
  #align = 'vh',
  labels = NULL,
  rel_widths = c(1, 1.15),
  rel_heights = c(1,1),
  #rel_widths = c(1, 0.02, 1),
  hjust = -1,
  nrow = 1)
plot_Ile_Abs


#plot_Abs1 <- plot_grid(plot_Ile_Abs, legend, nrow = 1, rel_widths = c(1, 0.1), labels = c("A",""))
#plot_Abs1

plot_Clade_Abs <- plot_grid(
  plot5 + theme(legend.position="none"),
  #NULL,
  plot4SymC,
  #align = 'vh',
  labels = NULL,
  rel_widths = c(1, 1.15),
  rel_heights = c(1,1),
  #rel_widths = c(1, 0.02, 1),
  hjust = -1,
  nrow = 1)
plot_Clade_Abs

#plot_Abs2 <- plot_grid(plot_Clade_Abs, legend2, nrow = 1, rel_widths = c(1, 0.1), labels = c("B",""))
#plot_Abs2

plot_final <- plot_grid(
  plot_Ile_Abs,
  NULL,
  plot_Clade_Abs,
  #align='vh',
  labels =c ("(A)","","(B)"),
  rel_heights = c(1,0.2,1),
  hjust = -1,
  nrow = 3)
plot_final

library(magick)
Poc_image <- c("C:/Users/uax75/OneDrive/Documents/R/TaraCoral_2020/Coral_PCA/Pocillopora_Silhouette.png")

#plot_final <- plot_grid(plot_Ile_Abs, NULL, plot_Clade_Abs, ncol = 3, rel_widths = c(1,-0.3,1.2))

pdf(file = paste0(outputPrefix, "_Abstract_PCA_biplot.pdf"), w = 12, h = 10)
ggdraw(plot_final) +
  draw_image(Poc_image, x = 0.44, y = 0.7, hjust = 1, vjust = 1, width = 0.075, height = 0.1) +
  draw_image(Poc_image, x = 0.44, y = 0.15, hjust = 1, vjust = 1, width = 0.075, height = 0.1)
dev.off()

png(file = paste0(outputPrefix, "_Abstract_PCA_biplot.png"), w = 12, h = 10, units="in", res=300)
ggdraw(plot_final) +
  draw_image(Poc_image, x = 0.44, y = 0.7, hjust = 1, vjust = 1, width = 0.075, height = 0.1) +
  draw_image(Poc_image, x = 0.44, y = 0.15, hjust = 1, vjust = 1, width = 0.075, height = 0.1)
dev.off()



##### Zoom in on SVD 5 #####
plotSVD5 <- ggbiplot(pca.POC, choices = c(1,2), obs.scale = 1, var.scale = 1, groups = clade.poc, 
                  ellipse = TRUE, ellipse.prob = 0.68, var.axes = FALSE) +
  geom_point(aes(shape=factor(sites.poc), color=factor(clade.poc)), size=4) +
  #ggtitle("Pocillopora Principle Component Anlaysis") +
  scale_shape_manual(name = '', values=c(16,17,15,4)) +
  scale_color_manual(name='', values= c("transparent","transparent","transparent",
                                        "transparent","darkorchid","transparent")) +
  #scale_x_continuous(limits=c(-80,130), breaks=c(-50,0,50,100), labels=c(-50,0,50,100)) +
  #scale_y_continuous(limits=c(-80,90), breaks=c(-50,0,50), labels=c(-50,0,50)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme_bw() +
  theme(text = element_text(size=16, face="bold")) +
  #theme(legend.direction = 'horizontal', 
  #      legend.position = 'top') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotSVD5

library(ggnewscale)
plotSVD5_2 <- ggbiplot(pca.POC, choices = c(1,2), obs.scale = 1, var.scale = 1, groups = clade.poc, 
                     ellipse = TRUE, ellipse.prob = 0.68, var.axes = FALSE) +
  geom_point(aes(shape=factor(sites.poc), color=factor(ile.poc)), size=4) +
  scale_color_manual(name='', values= c("transparent","transparent","transparent",
                                        mycols[4],mycols[5],"transparent",
                                        mycols[7],"transparent","transparent","transparent",
                                        "transparent","transparent","transparent","transparent",
                                        "transparent","transparent")) +

  #ggtitle("Pocillopora Principle Component Anlaysis") +
  new_scale("color") +
  scale_shape_manual(name = '', values=c(16,17,15,4)) +
  scale_color_manual(name='', values= c("transparent","transparent","transparent",
                                        "transparent","darkorchid","transparent")) +
  #scale_x_continuous(limits=c(-80,130), breaks=c(-50,0,50,100), labels=c(-50,0,50,100)) +
  #scale_y_continuous(limits=c(-80,90), breaks=c(-50,0,50), labels=c(-50,0,50)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme_bw() +
  theme(text = element_text(size=16, face="bold")) +
  #theme(legend.direction = 'horizontal', 
  #      legend.position = 'top') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotSVD5_2

pdf(file = paste0(outputPrefix, "_PCA_biplot_SVD5.pdf"), w = 12, h = 8.5)
plotSVD5
dev.off()

png(file = paste0(outputPrefix, "_PCA_SVD5.png"), w = 12, h = 8.5, units="in", res=300)
plotSVD5
dev.off()


pdf(file = paste0(outputPrefix, "_PCA_biplot_SVD5_2.pdf"), w = 12, h = 8.5)
plotSVD5_2
dev.off()

png(file = paste0(outputPrefix, "_PCA_SVD5_2.png"), w = 12, h = 8.5, units="in", res=300)
plotSVD5_2
dev.off()

Host_cts <- cts_raw_host_keep[,colnames(cts_raw_sym_keep)]
Host_MetaData <- data.frame("Colony" = colnames(Host_cts),
                            "Clade" = vlookup(colnames(Host_cts), CladeTable, result_column = "PocilloSVD", lookup_column = "Colony"),
                            "SymClade" = vlookup(colnames(Host_cts), CladeTable, result_column = "PopMetaT", lookup_column = "Colony"),
                            "Ile" = substr(colnames(Host_cts), 1,3),
                            "HeatIndex" =  gsub("I04|I05|I06|I07|I08|I09","ambient",gsub("I01|I02|I03|I10|I15","Warm",substr(colnames(Host_cts), 1,3))))

dds_Host <- DESeqDataSetFromMatrix(countData = Host_cts, 
                                      colData = Host_MetaData, 
                                      design = ~ Ile)
vsd_mat_Host <- varianceStabilizingTransformation(dds_Host, blind=TRUE)
    cts_vsd_Host <- assay(vsd_mat_Host)
data_cts_vsd_Host <- data.frame(t(cts_vsd_Host))


host_ANOVAData <- cbind(Host_MetaData, data_cts_vsd_Host)

# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


# PerMANOVA - partitioning the euclidean distance matrix by Clade
  library(vegan)
  adonis(host_ANOVAData[,-c(1:5)] ~ Ile * Clade, data = host_ANOVAData, method='eu')
  library(pairwiseAdonis)
  pairwise.adonis(x=host_ANOVAData[,-c(1:5)],factors=host_ANOVAData$Clade,
                  sim.function='vegdist', sim.method='euclidian', p.adjust.m='holm')


# PerMANOVA - partitioning the euclidean distance matrix by SymClade
  adonis(host_ANOVAData[,-c(1:5)] ~ SymClade, data = host_ANOVAData, method='eu')
  library(pairwiseAdonis)
  pairwise.adonis(x=host_ANOVAData[,-c(1:5)],factors=host_ANOVAData$SymClade,
                          sim.function='vegdist', sim.method='euclidian', p.adjust.m='holm')

# PerMANOVA - partitioning the euclidean distance matrix by Ile
  adonis(host_ANOVAData[,-c(1:5)] ~ Ile, data = host_ANOVAData, method='eu')
  library(pairwiseAdonis)
  pairwise.adonis(x=host_ANOVAData[,-c(1:5)],factors=host_ANOVAData$Ile,
                          sim.function='vegdist', sim.method='euclidian', p.adjust.m='holm')

# PerMANOVA - partitioning the euclidean distance matrix by HeatIndex
  adonis(host_ANOVAData[,-c(1:5)] ~ HeatIndex, data = host_ANOVAData, method='eu')
  library(pairwiseAdonis)
  pairwise.adonis(x=host_ANOVAData[,-c(1:5)],factors=host_ANOVAData$HeatIndex,
                          sim.function='vegdist', sim.method='euclidian', p.adjust.m='holm')

Host_cts_all <- highcts_raw_host_keep
Host_MetaData_all <- data.frame("Colony" = colnames(Host_cts_all),
                            "Clade" = vlookup(colnames(Host_cts_all), CladeTable, result_column = "PocilloSVD", lookup_column = "Colony"),
                            "SymClade" = vlookup(colnames(Host_cts_all), CladeTable, result_column = "PopMetaT", lookup_column = "Colony"),
                            "Ile" = substr(colnames(Host_cts_all), 1,3))

Host_MetaData_east <- subset(Host_MetaData_all, Ile %in% c("I01","I02","I03") & SymClade %in% c("D1","Group5"))
Host_cts_east <- Host_cts_all[,Host_MetaData_east$Colony]

dds_Host_east <- DESeqDataSetFromMatrix(countData = Host_cts_east, 
                                      colData = Host_MetaData_east, 
                                      design = ~ Ile)
vsd_mat_Host_east <- varianceStabilizingTransformation(dds_Host_east, blind=TRUE)
    cts_vsd_Host_east <- assay(vsd_mat_Host_east)
data_cts_vsd_Host_east <- data.frame(t(cts_vsd_Host_east))

host_ANOVAData_east <- cbind(Host_MetaData_east, data_cts_vsd_Host_east)


Host_MetaData_west <- subset(Host_MetaData_all, Ile %in% c("I08","I09","I10","I15") & SymClade %in% c("Group3","Group5"))
Host_cts_west <- Host_cts_all[,Host_MetaData_west$Colony]

dds_Host_west <- DESeqDataSetFromMatrix(countData = Host_cts_west, 
                                      colData = Host_MetaData_west, 
                                      design = ~ Ile)
vsd_mat_Host_west <- varianceStabilizingTransformation(dds_Host_west, blind=TRUE)
    cts_vsd_Host_west <- assay(vsd_mat_Host_west)
data_cts_vsd_Host_west <- data.frame(t(cts_vsd_Host_west))

host_ANOVAData_west <- cbind(Host_MetaData_west, data_cts_vsd_Host_west)


library(MASS)

# metaMDS Ordination

host_mds_east <- metaMDS(data_cts_vsd_Host_east)
plot(host_mds_east$points, pch=as.integer(as.factor(Host_MetaData_east$SymClade)),col=as.integer(as.factor(Host_MetaData_east$Ile)))

ford_POC_east <- fortify(host_mds_east, axes= 1:2) # fortify the ordination

library(ggplot2)
ploteast <- ggplot() +
          stat_ellipse(data=subset(ford_POC_east, Score == 'sites'), geom="polygon", 
                       aes(x=NMDS1,y=NMDS2,color=factor(Host_MetaData_east$Ile), fill=factor(Host_MetaData_east$Ile)), 
                   type="t", level=0.68, alpha=0.5, show.legend=F,) +
          # scale_color_manual(name='', values = mycols, labels = c("Las Perlas", "Coiba", "Malpelo", "Rapa Nui",
          #                                                         "Ducie", "Gambier", "Moorea", "Aitutaki", "Niue",
          #                                                         "Upolu", "Guam")) +
          # scale_shape_manual(name = '', values=c(16,17,15,4), labels=c("Site 01", "Site 02",
          #                                                              "Site 03", "Site 04")) +
          scale_fill_manual(name='', values=mycols[1:3],
                            labels=c("Las Perlas", "Coiba","Malpelo")) +
          scale_color_manual(name='', values=mycols[1:3],
                            labels=c("Las Perlas", "Coiba","Malpelo")) +
          geom_point(data = subset(ford_POC_east, Score == 'sites'),
                     mapping = aes(x = NMDS1, y = NMDS2,
                                   shape=factor(Host_MetaData_east$SymClade), 
                                   color=factor(Host_MetaData_east$Ile),
                                   fill=factor(Host_MetaData_east$Ile), size=3)) +
          scale_shape_manual(name = '', values=c(17,16), #16,17,18,15,4
                           labels=c("Durusdinium", "Cladocopium Group 5")) +
          # scale_x_continuous(limits=c(-3.2,3.2), breaks=seq(from=-3, to=3, by=1), labels=seq(from=-3, to=3, by=1)) +
          # scale_y_continuous(limits=c(-3.2,3.2), breaks=seq(from=-3, to=3, by=1), labels=seq(from=-3, to=3, by=1)) +
          theme_bw() +
          theme(text = element_text(size=20, face="bold")) +
          scale_size(guide="none") +
          theme(legend.direction = 'horizontal',
              legend.position = 'bottom',
              legend.text = element_text(size = 14)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          coord_fixed()
ploteast

pdf(file = paste0(outputPrefix, "_MDS_SVD4_CvsD.pdf"), w = 12, h = 8.5)
ploteast
dev.off()

png(file = paste0(outputPrefix, "_MDS_SVD4_CvsD.png"), w = 12, h = 8.5, units="in", res=300)
ploteast
dev.off()



# PerMANOVA - partitioning the euclidean distance matrix by SymClade
  library(vegan)  
  adonis(host_ANOVAData_east[,-c(1:4)] ~ Ile, data = host_ANOVAData_east, method='eu')
  library(pairwiseAdonis)
  pairwise.adonis(x=host_ANOVAData_east[,-c(1:4)],factors=host_ANOVAData_east$Ile,
                  sim.function='vegdist', sim.method='euclidian', p.adjust.m='holm')

  SST <- vlookup(lookup_value = host_ANOVAData_west$Colony, dict = metaData_env,
                 lookup_column = "Colony", result_column = "SST")
  adonis(host_ANOVAData[,-c(1:4)] ~ SST, data = host_ANOVAData, method = 'eu')

  
### Genetic vs Expression Distances
  
  genData <- read.table("Poc_G111_Debug_RAxML_distances.txt", header = F)
    colnames(genData)[1] <- "Colony"
    rownames(genData) <- genData$Colony
    genData <- genData[,2:ncol(genData)]
    colnames(genData) <- rownames(genData)
    library(reshape)
    Pair_GenDst <- melt(as.matrix(genData))
    colnames(Pair_GenDst) <- c("ColonyA","ColonyB","PairGenDst")
    
  Host_cts_all <- highcts_raw_host_keep
  Host_MetaData_all <- data.frame("Colony" = colnames(Host_cts_all),
                               "Clade" = vlookup(colnames(Host_cts_all), CladeTable, result_column = "PocilloSVD", lookup_column = "Colony"),
                                "SymClade" = vlookup(colnames(Host_cts_all), CladeTable, result_column = "PopMetaT", lookup_column = "Colony"),
                                "Ile" = substr(colnames(Host_cts_all), 1,3))

  dds_Host_all <- DESeqDataSetFromMatrix(countData = Host_cts_all, 
                                          colData = Host_MetaData_all, 
                                          design = ~ Ile)
  vsd_mat_Host_all <- varianceStabilizingTransformation(dds_Host_all, blind=TRUE)
      cts_vsd_Host_all <- assay(vsd_mat_Host_all)
  data_cts_vsd_Host_all <- data.frame(t(cts_vsd_Host_all))
  
  host_ANOVAData_all <- cbind(Host_MetaData_all, data_cts_vsd_Host_all)
