# Selection of top environmental variables and constrained correspondence analysis (CCA)

To produce Figure 5 of the article, the following files are required: <br>
-the associations between host, symbiont, and island: Variables11Islands.txt in this directory <br>
-raw gene expression for the host: Pocillopora_MetaT_ReadCount.tab available in the folder for Figure_4 or at https://doi.org/10.5281/zenodo.6341761 <br>
-raw gene expression for the symbiont: CladocopiumC1_MetaT_ReadCount.tab available in the folder for Figure_4 or at https://doi.org/10.5281/zenodo.6341761 <br>
-the historical and in situ environmental data: TaraPacific_Nutrents_SST_timeseries_mean_products-20220317_11Islands.xlsx in this directory <br>
<br>
```r
## 1. INITIALIZATION ------------------------------------------------------------------------------------------------------------------------------------

  # Step 1 - Set Defaults

    # Set the prefix for each output file name
      outputPrefix <- "TaraPacific_Holobiont_CCA_GitHub_2022-12"


## 2. LOAD DATA & METADATA -------------------------------------------------------------------------------------------------------------------------------

  # Step 1 - Load Pocillopora Host and Photosymiont raw count data (reads mapped to the predicted coding sequences)
    cts_raw_host <- read.table("Pocillopora_MetaT_ReadCount.tab", header=TRUE, com='', sep='\t', row.names=1, check.names=FALSE)
    cts_raw_sym <- read.table("CladocopiumC1_MetaT_ReadCount.tab", header=TRUE, com='', sep='\t', row.names=1, check.names=FALSE)


  # Step 2 - Load table with lineage assignations
    CladeTable <- read.table("Variables11Islands.txt", header = T)
    
  
  # Step 3 - Load historical and in situ environmental data
    library(xlsx)
    envData <- read.xlsx("TaraPacific_Nutrents_SST_timeseries_mean_products-20220317_11Islands.xlsx", sheetName = "Sheet1")
        
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
```
