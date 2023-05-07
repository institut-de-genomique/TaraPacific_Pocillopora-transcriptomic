# Biogeography of Host and Photosymbiont Lineages
<br>
To produce Figure 1 of the article, the following files are required:
<br>
-the associations between host, symbiont, and island: ITS2-Kmeans-clustering_Pocillopora_JL.xlsx in this directory <br>
-raw gene expression for the host: Pocillopora_MetaT_ReadCount.tab available in this directory or at https://doi.org/10.5281/zenodo.6341761 <br>
-the historical and in situ environmental data: TaraPacific_SST_timeseries_mean_products_mai2021.xlsx and TaraPacific_EnvironmentContext_condensed-20210618_preliminary_Barbara.xlsx in this directory <br>
-the mean SST climatology data (1981 - 2010) from NOAA: X90.79.167.222.277.1.2.3.nc in this directory <br>
-the proportion of Durusdinium in each colony: Pocillopora_Figures_2021-10_SymbiontClade_PropData_D.csv in this directory <br>
<br>

```r

  # Set the prefix for each output file name
    outputPrefix <- "Pocillopora_Figures_2021-10"

############################
##### Load Count Data  #####
    
  # Step 1 - Load Pocillopora and Cladocopium raw count file with reads mapped to the predicted coding sequence
    cts_raw_poc <- read.table("Pocillopora_MetaT_ReadCount.tab", header=TRUE, com='', sep='\t', row.names=1, check.names=FALSE)
    colnames(cts_raw_poc) <- gsub("C0","C",colnames(cts_raw_poc))
    
  # Step 2 - Merge Host and Symbiont read count dataframes by colony
    cols_keep <- colnames(cts_raw_poc)
    cts_raw <- cts_raw_poc[,cols_keep]
  
  # Step 3 - Add coral host clade information from D. Forcioli's SVD tree
    Colony <- colnames(cts_raw)
    Ile <- substr(Colony, 1, 3)
    IleSite <- substr(Colony, 1, 6)
    Site <- substr(Colony, 4, 6)
    
    library(xlsx)
    CladeTable <- read.xlsx("ITS2-Kmeans-clustering_Pocillopora_JL.xlsx", sheetName = "RData", header = T)
    
    library(expss)
    Clade <- vlookup(Colony, CladeTable, result_column = "PocilloSVD", lookup_column = "Colony")
    IleClade <- paste(substr(unique(Colony), 1, 3), Clade, sep="")
    
    #3b - Add symbiont clade information from metaT-based assignations (G1 - G6)
      SymClade <- vlookup(Colony, CladeTable, result_column = "PopMetaT", lookup_column = "Colony")
      IleSymClade <- paste(substr(unique(Colony), 1, 3), SymClade, sep="")
      Holobiotype <- paste(Clade, SymClade, sep="_")
      IleHolobiotype <- paste(Ile, Holobiotype, sep='_')

  # Step 4 - Create table containing all pertinent metadata
    sampleTable <- data.frame("Colony" = Colony,
                       "Ile" = Ile, "Site" = Site, "IleSite" = IleSite, 
                       "Clade" = Clade, "SymClade" = SymClade, 
                       "Holobiotype" = Holobiotype, "IleHolobiotype" = IleHolobiotype)
    
  # Step 6 - Convert count table to DESeq friendly counts matrix
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = cts_raw,
                                  colData = sampleTable,
                                  design = ~ Ile)
                                # design = ~ CladeSymClade)
    
    
###############################################################  
##### Collection site map and Symbiodinaciae assignations #####
    
    
    library("pheatmap")
    library(RColorBrewer)
    library(rcartocolor)
    mycols <- brewer.pal(11, "Spectral")
    mycols[6] <- "black"
    mycarto_cols <- carto_pal(12, "Safe")
    mybrew_cols <- brewer.pal(12, "Paired")
    
    mycols2 <- c(mybrew_cols[5], mybrew_cols[1], mybrew_cols[3], mybrew_cols[2], mybrew_cols[4])
    mycols3 <- c(mycarto_cols[3], 
                 mycarto_cols[4], mycarto_cols[2], mycarto_cols[7], 
                 mycarto_cols[5],mycarto_cols[11])
    
    # Didier Colors: mycols2 <- c(mybrew_cols[10], mybrew_cols[4], mybrew_cols[6], mybrew_cols[8], mybrew_cols[2])
    # Didier Colors: mycols3 <- c(mybrew_cols[11], mybrew_cols[1], mybrew_cols[9],
    #              mybrew_cols[5], mybrew_cols[3], mybrew_cols[7])
    
    # mybrew_cols <- brewer.pal(12, "Paired")
    # mycols2 <- c("tomato2","azure2","turquoise4", "palegreen", "darkorchid","white")
    # mycols3 <- c("chartreuse","darkred","#F768A1", "#0066CC", "#CCCC33", "#FF9933", "#009900") #brewer.pal(7,"Set3")
    
  # Step 6 - Add select environmental data to sample info table
    library('xlsx')
    library('dplyr')
    
    histData <- read.xlsx("TaraPacific_SST_timeseries_mean_products_mai2021.xlsx",
                         # sheetName = "TaraPacific_SST_timeseries_mean")
                         sheetName = "RData")
    # rownames(histData) <- histData$Colony
    # currData <- read.xlsx("TaraPacific_EnvironmentContext_condensed-20210618_preliminary.xlsx",
    #                      sheetName = "IleSiteData")
    currData <- read.xlsx("TaraPacific_EnvironmentContext_condensed-20210618_preliminary_Barbara.xlsx",
                         sheetName = "RData")
    
    # rownames(currData) <- currData$Colony
    
    # histData_sitesumm <- aggregate(. ~ IleSite, data = histData[,c(2,6:ncol(histData))], FUN = mean)
    # currData_sitesumm <- aggregate(. ~ IleSite, data = currData[,c(2,6:ncol(currData))], FUN = mean)
     
    library(dplyr)
    histData_sitesumm <- histData %>% group_by(IleSite) %>% summarise_all(funs(mean)) 
    currData_sitesumm <- currData %>% group_by(IleSite) %>% summarise_all(funs(mean))
    
    # CommonColonies <- intersect(histData$Colony, currData$Colony)
    # histData_common <- histData[CommonColonies,]
    # currData_common <- currData[CommonColonies,]
    
    envData <- left_join(histData_sitesumm, currData_sitesumm, by="IleSite")
    colnames(envData)
    envData <- envData[,!(colnames(envData) %in% c("ISCOA","lat.x","lon.x","station","dt","lat.y","lon.y","Colony"))]
    
    
    sampleTable$SST <- vlookup(lookup_value = substr(sampleTable$Colony,1,6),
                               dict = envData,
                               lookup_column = "IleSite",
                               result_column = "SST_mean_DegC")
    

library(ggplot2)
library(cowplot)
    
  # Plot on map
      library(mapdata)
      library(RColorBrewer)
      library(rcartocolor)
      display_carto_all(colorblind_friendly = TRUE)
      
      library(xlsx)
      #sitedata <- read.xlsx("IslandSiteData.xlsx", sheetName="SiteData")
      # Create Visual Map with Data
      Lat = c(3.942909,7.49999,9.00973,-27.112722,-24.680725,-23.11813,-17.48992,
              -21.236736,-19.054445,-13.759029,13.444304)
      Lon360 = c((360-81.562814),(360-80.449997),(360-79.51939),(360-109.349686),
              (360-124.787888),(360-134.97124),(360-149.7687),
              (360-159.777664),(360-169.867233),(360-172.10463),144.793732)
      LIle = c("I03","I02","I01","I04",
               "I05","I06","I07","I08",
               "I09","I10","I15")
      LName = c("Malpelo","Coiba","Las Perlas","Rapa Nui",
                "Ducie","Gambier","Moorea","Aitutaki",
                "Niue","Upolu","Guam")
      labframe2 = data.frame(Lat,Lon360,LIle,LName)
      
      pacific <- map_data("world2")

      library(ggrepel)
      sampleTable$Ile <- factor(sampleTable$Ile, levels=c('I15','I10','I09','I08',
                                     'I07','I06','I05','I04',
                                     'I03','I02','I01'))
      ile_names <- list(
        'I15'="Guam",
        'I10'="Upolu",
        'I09'="Niue",
        'I08'="Aitutaki",
        'I07'="Moorea",
        'I06'="Gambier",
        'I05'="Ducie",
        'I04'="Rapa Nui",
        'I03'="Malpelo",
        'I02'="Coiba",
        'I01'="Las Perlas")
      
      ile_labeller <- function(variable,value){
        return(ile_names[value])}
      
      
#################################################
##### Load Mean SST Climatology 1981 - 2010 #####  
  library(ncdf4)
  library(reshape2)
  library(rgdal)
  library(rgeos)
      
  sstFile <- nc_open("X90.79.167.222.277.1.2.3.nc")
      
      ret <- list()
      names(sstFile$var)
      ret$lat <- ncvar_get(sstFile, "lat")
      ret$lon <- ncvar_get(sstFile, "lon")
      ret$sst <- ncvar_get(sstFile, "sst")
      nc_close(sstFile)
      
      
      # Prepare to Subset Data to Include only Tropical Pacific
      
      melt_sst <- function(L) {
        # dimnames(L$dhw0) <- list(long = L$lon, lat = L$lat)
        # ret <- melt(L$dhw0, value.name = "dhw0")
        dimnames(L$sst) <- list(long = L$lon, lat = L$lat)
        ret <- melt(L$sst, value.name = "sst")
        #dimnames(L$dhw8) <- list(long = L$lon, lat = L$lat)
        #ret <- melt(L$dhw8, value.name = "dhw8")
      }
      
      # Convert Data to Dataframe for Subsetting
      msst <-  melt_sst(ret)      
      

plot_map <- ggplot() +
  ggtitle("Annual Mean SST Climatology 1981 - 2010") +
  geom_raster(data = msst, aes(x = long, y = lat, fill = sst),interpolate = TRUE) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  # ggtitle("NOAA Coral Reef Watch 5km 12-wk DHW: 12-Sept to 03-Dec 2016") +
  geom_map(data=pacific, map=pacific,aes(map_id=region),  fill="grey90", color="grey") +
  # geom_raster(data = mdhw_sub, aes(x = long, y = lat, fill = dhw),interpolate = TRUE) +
  geom_point(data = labframe2, aes(x=Lon360, y=Lat), shape=1) +
  geom_text_repel(data = labframe2, aes(x=Lon360, y=Lat, label=LName),
                  min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_fill_fermenter(palette = "RdBu", breaks = c(5,10,15,20,24,26,27,28,29,30), na.value = "transparent") +
  # scale_fill_gradientn(limits=c(0,10), breaks = c(0,5,10),
  #                      labels=c("0","5","10+"), #colours=c("transparent","yellow","tomato4","black"),
  #                      colours = c(myPalette(100)), oob=squish, na.value = "transparent") +
  ylab(expression(bold(paste("Latitude (°N)")))) +
  xlab(expression(bold(paste("Longitude (°E)")))) +
  scale_x_continuous(limits = c(100, 300), breaks = seq(from = 100, to = 300, by = 50), expand = c(0,0)) +
  scale_y_continuous(limits = c(-50, 50), breaks = seq(from = -50, to = 500, by = 25), expand = c(0,0)) +
  # xlim(100, 300) +
  # ylim(-50,50) +
  theme_bw() +
  labs(fill="Annual \nMean SST (?C)") +
  # labs(fill="Cumulative DHW \n (?C - Weeks)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))
plot_map
  
library(forcats)

plot_SVD <- ggplot(subset(sampleTable, !Clade %in% c("Unassigned") & !SymClade %in% c("Unassigned")))+
  geom_bar(mapping = aes(x = fct_reorder(Colony, Clade), y = 1, fill = Clade), 
           stat = "identity", 
           width = 0.7) +
  facet_grid(.~Ile, scales = "free_x", labeller = ile_labeller) +
  scale_x_discrete(label = function(x) substr(x, 4,6)) +
  scale_fill_manual(values=mycols2) +
  guides(fill=guide_legend(ncol=1)) +
  #theme_void() +
  theme(panel.spacing.x = unit(1, "mm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.key.size = unit(0.4, "cm"))
plot_SVD

propData2 <- read.csv("Pocillopora_Figures_2021-10_SymbiontClade_PropData_D.csv", header = T)
propData2$Ile <- factor(propData2$Ile, levels=c('I15','I10','I09','I08',
                                     'I07','I06','I05','I04',
                                     'I03','I02','I01'))

plot_SymGG <- ggplot(subset(propData2, !Clade %in% c("Unassigned")))+
  geom_bar(mapping = aes(x = fct_reorder(Colony, Clade), y = Prop, fill = SymClade), 
           stat = "identity", position="fill",
           width = 0.7) +
  facet_grid(.~Ile, scales = "free_x", labeller = ile_labeller) +
  #scale_x_discrete(limits = levels(sampleTable$Ile)) +
  scale_fill_manual(values=mycols3, labels = c("Durusdinium D1",
                               "Cladocopium C1/C42 Lineage 1",
                               "Cladocopium C1/C42 Lineage 2",
                               "Cladocopium C1/C42 Lineage 3",
                               "Cladocopium C1/C42 Lineage 4",
                               "Cladocopium C1/C42 Lineage 5")) +
  guides(fill=guide_legend(ncol=1)) +
  #theme_void() +
  theme(panel.spacing.x = unit(1, "mm"),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.key.size = unit(0.4, "cm"))
plot_SymGG

plot_SST <- ggplot(subset(sampleTable, !Clade %in% c("Unassigned") & !SymClade %in% c("Unassigned")))+
  geom_bar(mapping = aes(x = fct_reorder(Colony, Clade), y = 1, fill = SST), 
           stat = "identity", 
           width = 0.7) +
  facet_grid(.~Ile, scales = "free_x", labeller = ile_labeller) +
  scale_x_discrete(label = function(x) substr(x, 4,6)) +
  scale_fill_distiller(palette = "RdBu", name = "SST (°C)",
                       breaks = c(5,10,15,20,24,26,27,28,29,30), na.value = "transparent") +
  # guides(fill=guide_legend(title="SST (?C)")) +
  #theme_void() +
  theme(panel.spacing.x = unit(1, "mm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        strip.text.x = element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.key.size = unit(0.4, "cm"))
plot_SST

legend <- plot_grid(get_legend(plot_map), get_legend(plot_SymGG), get_legend(plot_SVD),  ncol = 1, align='v',axis='tb') # get_legend(plot_SST)
plot_map_noleg <- plot_map + theme(legend.position = "none")
plot_SymGG_noleg <- plot_SymGG + theme(legend.position = "none",
                                       plot.margin = unit(c(0,0,0,0), "cm"))
plot_SVD_noleg <- plot_SVD + theme(legend.position = "none",
                                   plot.margin = unit(c(0,0,0,0), "cm"))
plot_SST_noleg <- plot_SST + theme(legend.position = "none",
                                   plot.margin = unit(c(0,0,0,0), "cm"))

plot <- plot_grid(plot_SymGG_noleg, NULL, plot_SVD_noleg, NULL, plot_SST_noleg,
                  align = "vh", ncol = 1, rel_heights = c(0.5,-0.18,0.5,-0.15,0.3)) #,-0.1,0.3))
plot

plot2 <- plot_grid(plot_map_noleg, plot,
                   align = 'vh', ncol = 1, axis = 'tb', 
                   rel_heights = c(3,2),
                   labels = c("A","B"))
plot2

plot3 <- plot_grid(plot2, legend, nrow = 1, rel_widths = c(10, 1.8))
plot3

pdf(file = paste0(outputPrefix, "_SymbiontClade_NewCols_SST_Map.pdf"), w = 12, h = 8.5)
plot3
dev.off()
    
png(file = paste0(outputPrefix, "_SymbiontClade_NewCols_SST_Map.png"), w = 12, h = 8.5, units = "in", res=300)
plot3
dev.off()
