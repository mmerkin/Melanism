library(GenotypePlot)
library(tidyverse)
library(vcfR)
library(cowplot)

vcftest <- read.vcfR("FILE_25_sorted.vcf.gz")

SPECIES_popmap <- read.csv("POPMAP.csv")
G_plot <- genotype_plot(vcf    = "FILE_25_sorted.vcf.gz", 
                         chr = "CHROMOSOME",
                         start  = START_POSITION,
                         end    = END_POSITION,
                         popmap = SPECIES_popmap,
                         cluster        = T,
                         plot_phased=F,
                         colour_scheme=c("#332288","#88CCEE","#AA4499"))


G_meta_order <- SPECIES_popmap[match(G_plot$dendro_labels, SPECIES_popmap$Ind),]
G_seg <- ggplot() + geom_tile(aes(y=1:length(G_meta_order$Ind), x=0.1, 
                                  fill=G_meta_order$pop), show.legend = F) + 
  scale_fill_manual(values = c("black", "grey")) +
  theme_void()

plot_grid(G_plot$dendrogram,
          G_seg, 
          G_plot$genotypes + guides(fill="none"), rel_widths = c(1,0.5,5), 
          axis = "tblr", align = "h", ncol = 3, nrow = 1)
