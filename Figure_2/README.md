
# Tanglegram and PACo analysis 
This file contains the code to produce a tanglgram between the host and the symbiont (Figure 2a) and the code for the Procrustean Approach to Cophylogeny (PACo)  analysis (Figure 2b)  

## 1. Tanglegram (Figure 2a)

```r
library(dendextend)
library(DECIPHER)
library(ape)
library(ggtree)

#Genetic clades for the host and the symbiont based on the SNP.
Variable<-read.table("Variables11Islands.txt",sep="\t",h=T)

PocTree<-ReadDendrogram(file = "Poc_Debug_RAxML.meg_rooted.nwk",internalLabels = F,keepRoot = T)

PocTree2<- prune(PocTree, c(setdiff(labels(PocTree),labels(dendSNP2)),"I09S03C010POC"))

PocTree3 <- as.dendrogram(PocTree2) %>%
  hang.dendrogram(hang = -1) %>%
  set("labels_cex", c(.5)) %>%
  set("labels_cex", c(.5)) %>% 
  set("branches_lwd",c(0.5)) %>%
  set("labels_col",c("#fb9a99", "#a6cee3" ,"#b2df8a", "#1f78b4","#33a02c")[as.factor(unlist(lapply(labels(PocTree2),function(x){Variable[Variable$Samples==x,"PocilloGG"]})))])


CladoTree<-ReadDendrogram(file = "Cladocopium_SNP.dendrogram.nwk",internalLabels = F,keepRoot = T)

CladoTree2 <- as.dendrogram(CladoTree) %>%
  set("labels_cex", c(.5)) %>%
  set("labels_cex", c(.5)) %>% 
  set("branches_lwd",c(0.5)) %>%
  set("labels_col",c("#117733", "#cc6677" ,"#44aa99", "#332288","#6699cc")[as.factor(unlist(lapply(dendSNP$labels[dendSNP$order],function(x){Variable[Variable$Samples==x,"SymbioGG"]})))]) 


CladoTree3 <- CladoTree2 %>% prune(setdiff(labels(CladoTree2),labels(PocTree3)))

dendCombined<-dendlist(PocTree3, CladoTree3)

tanglegram(dendCombined)

pdf(file="tangledram.pdf",height=12,width=9)
dendCombined2<-dendCombined %>%
  untangle(method = "step2side") %>%
  tanglegram(common_subtrees_color_lines = F,
             highlight_distinct_edges = F,
             highlight_branches_lwd = F,
             common_subtrees_color_lines_default_single_leaf_color = "black",
             lwd = 1.5,
             sort = F)
dev.off()
```

## 2. PACo analysis  (Figure 2b)

```r
library(vegan)
library(ape)
library(paco)
library(ggplot2)
library(dplyr)
library(stringr)
library(xlsx)

#Establish PACo function: adjustment prior to procustes analysis
    PACo <- function(H.dist, P.dist, HP.bin)
      { HP.bin <- which(HP.bin > 0, arr.in=TRUE)
      H.PCo <- pcoa(H.dist, correction="cailliez") #Performs PCo of Host distances 
      P.PCo <- pcoa(P.dist, correction="cailliez") #Performs PCo of Parasite distances
      if (is.null(H.PCo$vectors.cor)==TRUE) H.PCo <- H.PCo$vectors else
        H.PCo <- H.PCo$vectors.cor      # returns corrected pcoord 
      if (is.null(P.PCo$vectors.cor)==TRUE) P.PCo <- P.PCo$vectors else
        P.PCo <- P.PCo$vectors.cor
      H.PCo <- H.PCo[HP.bin[,1],]  #adjust Host PCo vectors 
      P.PCo <- P.PCo[HP.bin[,2],]  #adjust Symbiont PCo vectors
      list (H.PCo = H.PCo, P.PCo = P.PCo)}

#Genetic clades for the host and the symbiont based on the SNP.
#Variable<-read.table("Variables11Islands.txt",sep="\t",h=T)
Variable<-read.table("C:/Users/uax75/OneDrive/Documents/R/TaraCoral_2020/Final_Results/Pocillopora/ManuscriptDocs/Scripts/Variables11Islands.txt",sep="\t",h=T)

#Remove hybrid host colony
library(tidyr)
Variable <- Variable %>% drop_na(PocilloGG)

# Create data matrix
HP <- data.frame(matrix(, nrow = length(unique(Variable$PocilloGG)), ncol = nrow(Variable)))
colnames(HP) <- Variable$Samples
rownames(HP) <- sort(unique(Variable$PocilloGG))

library(expss)
for (j in 1:ncol(HP)) {
  HP[,j] <- vlookup(lookup_value = colnames(HP)[j],
                        dict = Variable,
                        lookup_column = "Samples",
                        result_column = "PocilloGG")
  for (k in 1:nrow(HP)) {
    if (HP[k,j] == rownames(HP)[k]) {
      HP[k,j] <- 1
    } else {
      HP[k,j] <- 0
    }
  }
}

HP[] <- lapply(HP, function(x) type.convert(as.character(x)))

# Rename colonies so that they match with phylogeny
#colnames(HP) <- gsub("POC0","POC",gsub("C","C0",colnames(HP)))

#Phylogenetic trees:
#Read in Pocillopora host tree
TreeH <- read.tree("POC_SVDSub_outTree_midpointRooted_RaxMLBootstrap.nwk")

#Read in Cladocopium tree
TreeP <- read.tree("Cladocopium_SNP.dendrogram.newick")

host.D <- cophenetic(TreeH)
para.D <- cophenetic(TreeP)

#Sort host and photosymbiont taxa in distance matrices to match the HP matrix:
# Subset HP to include only colonies with full assignment for host and symbiont
SymbSamples <- TreeP$tip.label

HP <- HP[,c(intersect(SymbSamples,colnames(HP)))] 
#HP <- HP[-c(6),]
host.D <- host.D[rownames(HP), rownames(HP)]
para.D <- para.D[colnames(HP), colnames(HP)]

#Apply PACo Function
PACo.fit <- PACo(host.D, para.D, HP)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo) #Procrustes Ordination 
NLinks = sum(HP) #Number of H-P links; needed for further computations

#Goodness-of-fit-test 
# Note: Takes a long time to run
m2.obs <- HP.proc$ss #observed sum of squares
N.perm = 100000 #set number of permutations for testing
P.value = 0
set.seed(8765) ### use this option to obtain reproducible randomizations
    
    for (n in c(1:N.perm)){ 
      if (NLinks <= nrow(HP) | NLinks <= ncol(HP)) 	#control statement to avoid all symbionts being associated to a single host 
      {	flag2 <- TRUE 
      while (flag2 == TRUE)	{ 
        HP.perm <- t(apply(HP,1,sample))
        if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE
      }  
      } else { HP.perm <- t(apply(HP,1,sample))} #permutes each HP row independently
      PACo.perm <- PACo(host.D, para.D, HP.perm)
      m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss #randomized sum of squares
      if (m2.perm <= m2.obs)
      {P.value = P.value + 1} 
    }

P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", "P-value = ", P.value, " based on ", N.perm," permutations.")

#Contribution of individual links
HP.ones <- which(HP > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks) # empty matrix of jackknifed squared residuals
colnames (SQres.jackn) <- paste(rownames(HP.proc$X),rownames(HP.proc$Yrot), sep="-") # colnames identify the H-P link
t.critical = qt(0.975,NLinks-1) # needed to compute 95% confidence intervals.

    for(i in c(1:NLinks)){
      HP.ind <- HP
      HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
      PACo.ind <- PACo(host.D, para.D, HP.ind)
      Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo) 
      res.Proc.ind <- c(residuals(Proc.ind))
      res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
      SQres.jackn [i, ] <- res.Proc.ind	#Append residuals to matrix of jackknifed squared residuals
    }
    
SQres.jackn <- SQres.jackn**2 # Jackknifed residuals are squared
SQres <- (residuals (HP.proc))**2 # Vector of original square residuals

#jackknife calculations:
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres)) # apply jackknife function to matrix
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE) # mean jackknife estimate per link
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE) # standard deviation of estimates
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks) # upper 95% confidence interval

#Convert vectors into dataframes
phi.UCI.df <- as.data.frame(phi.UCI)
phi.mean.df <- as.data.frame(phi.mean)
phi.UCI.df <- cbind(UCI_links = rownames(phi.UCI.df), phi.UCI.df)
rownames(phi.UCI.df) <- 1:nrow(phi.UCI.df)
phi.mean.df <- cbind(mean_links = rownames(phi.mean.df), phi.mean.df)
rownames(phi.mean.df) <- 1:nrow(phi.mean.df)

#Merge dataframes
phi.UCI.mean <- merge(phi.UCI.df, phi.mean.df, by.x=c("UCI_links"), by.y=c("mean_links"))
    
#Add new column for Pocillopora species/haplotype
phi.UCI.mean$Pocillopora <- phi.UCI.mean$UCI_links
phi.UCI.mean$Pocillopora <- gsub("\\-.*","",phi.UCI.mean$Pocillopora)

#Merge with Cladocopium clades
PACo_clades <- phi.UCI.mean
PACo_clades$Clad_clades <- paste(vlookup(lookup_value = substr(phi.UCI.mean$UCI_links,6,18),
                             dict = Variable,
                             lookup_column = "Samples",
                             result_column = "SymbioGG"),
                           substr(phi.UCI.mean$UCI_links,6,18),
                           sep = '_')

#Find mean of phi.mean
median(phi.UCI.mean$phi.mean) # median = 3.472969e-05
Value_phi.mean <- median(phi.UCI.mean$phi.mean)
      
PACo_clades$Pocillopora <- factor(PACo_clades$Pocillopora,
                                  levels=c("SVD1","SVD2","SVD3","SVD4","SVD5"))

#Plot
library(RColorBrewer)
mybrew_cols <- brewer.pal(12, "Paired")
mycols_clade <- c(mybrew_cols[5], mybrew_cols[1], mybrew_cols[3],
                  mybrew_cols[2], mybrew_cols[4])
    # these correspond (in order) to: SVD1, SVD2, SVD3, SVD4, and SVD5
      
library(rcartocolor)
mycarto_cols <- carto_pal(12, "Safe")
mycols_symclade <- c(mycarto_cols[1], mycarto_cols[3], mycarto_cols[4],
                     mycarto_cols[2], mycarto_cols[7], mycarto_cols[5],
                     mycarto_cols[11], "white")
SymCladeColorTable <- data.frame("SymClade" = c("Group1","Group2","Group3","Group4","Group5"),
                                 "SymColor" = c(mycarto_cols[4],mycarto_cols[2],
                                                mycarto_cols[7], mycarto_cols[5],
                                                mycarto_cols[11]))
    # these correspond (in order) to: Cladocopium C1 (undefined), D1, Lineage 1, Lineage 2, Lineage 3, Lineage 4, Lineage 5, Unassigned
      
procrustplot <- ggplot(PACo_clades, aes(x=reorder(Clad_clades, phi.mean), y=phi.mean, fill=Pocillopora))+
  geom_bar(stat='identity')+
  geom_hline(yintercept=Value_phi.mean, linetype='dashed', col = 'red')+
  theme_minimal(base_size=15)+
  geom_errorbar(aes(ymin=phi.mean, ymax=phi.UCI), width=.1)+
  scale_fill_manual(name="Pocillopora lineage",
                    values = c("SVD1"= mybrew_cols[5],
                               "SVD2"= mybrew_cols[1],
                               "SVD3"= mybrew_cols[3],
                               "SVD4"= mybrew_cols[2],
                               "SVD5"= mybrew_cols[4]))+
  theme(axis.text=element_text(size=8))+
        # axis.text.y = element_text(colour = rev(vlookup(lookup_value = substr(reorder(PACo_clades$Clad_clades, PACo_clades$phi.mean),1,6),
        #                                             dict = SymCladeColorTable,
        #                                             lookup_column = "SymClade",
        #                                             result_column = "SymColor"))))+
  coord_flip()+
  labs(title="",
       x ="Pocillopora - Cladocopium links", y = "Squared residuals")
procrustplot

pdf(file = "ProcrusteanFit_Plot.pdf", h = 11, w = 8.5)
procrustplot
dev.off()
```
