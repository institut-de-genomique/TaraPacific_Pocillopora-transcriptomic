
# Tanglegram and PACo analysis 
This file contains the code to produce a tanglgram between the host and the symbiont (Figure 2a) and the code for the Procrustean Approach to Cophylogeny (PACo)  analysis (Figure 2b)  

1. Tanglegram (Figure 2a)

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


CladoTree<-ReadDendrogram(file = "Cladocopium_SNP.dendrogram.nex,internalLabels = F,keepRoot = T)

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

2. PACo analysis  (Figure 2b)
