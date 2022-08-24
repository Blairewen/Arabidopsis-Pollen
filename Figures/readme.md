# Distinct chromatin signatures in Arabidopsis male gametophyte

## Requirements
```
homer v4.11.1
deeptools v3.3.1
R v4.0.3
```

## 1. venn diagram (using Spm bivalency as an example)
```bash
library(VennDiagram)
library(ggplot2)

Spm_bivalency <- draw.pairwise.venn(10761,10841,7843,
                                    category = c('Spm_H3K27me3','Spm_H3K4me3'),
                                    col = "transparent",
                                    fill = c("skyblue3", "hotpink3"),
                                    alpha = 0.5,
                                    label.col = c("darkblue", "white", "darkred"),
                                    cex = 2,
                                    fontfamily = "serif",
                                    fontface = "bold",
                                    cat.default.pos = "outer",
                                    cat.col = c("darkblue", "darkred"),
                                    cat.cex = 1.5,
                                    cat.fontfamily = "serif",
                                    cat.pos  = c(10,350),
                                    rotation.degree=180,
                                    inverted = F)

ggsave("Figures/Spm_bivalency_venn.png",plot = Spm_bivalency,width = 7,height = 7)
```
![Spm_bivalency_venn](https://github.com/Blairewen/Arabidopsis-Pollen/blob/main/Figures/Spm_bivalency_venn.png)

Using "shuffleBed(bedtools)" to take the random regions in the same chromsome with equal histone peaks' length 10,000 times, and using "mergePeaks" of homer to calculate the overlap between random  regions, the times that the overlap regions count more than 7843 devide 10,000 as the venn significant test p-value.

```bash
shuffleBed -i Spm_H3K27me3.bed -g tair10_genome_size.txt -chrom > Spm_H3K27me3_shuffle_1.bed

shuffleBed -i Spm_H3K4me3.bed -g tair10_genome_size.txt -chrom > Spm_H3K4me3_shuffle_1.bed

mergePeaks -d given \
           -venn Spm_bivalency_shuffle_1.venn \
           -matrix Spm_bivalency_shuffle_1 \
           -prefix Spm_bivalency_shuffle_1 \
           Spm_H3K27me3_shuffle_1.bed Spm_H3K4me3_shuffle_1.bed

# The peak count of "Spm_bivalency_shuffle_1_Spm_H3K27me3_shuffle_1.bed_Spm_H3K4me3_shuffle_1.bed" is the overlap of the 1st shuffle.
```

## 2. box plot and scatter diagram (using HBG 4 quantiles genes' FPKM as the example)
```bash
library(ggplot2)
library(ggpubr)
HBG <- read.table("Figures/HBG_SortedBy_PeakRPKM_H3K27me3_4quantile.data",sep = "\t",header = T,quote = "")
HBG_boxplot <- ggplot(HBG,aes(x = quantile, y = log10(RS_SN_median+1), fill = quantile)) + 
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = c('none'))+
  ggtitle("H3K27me3 quartile")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab("log10(FPKM)")+
  scale_fill_manual(values= c("#af2337", "#8B6508","#2967a0","#2f3c28"))+
  scale_color_manual(values= c("#000000", "#000000","#000000","#000000"))+
  stat_compare_means(comparisons = list(c("1_quantile","2_quantile"),c("2_quantile","3_quantile"),c("3_quantile","4_quantile")))
HBG_boxplot
ggsave(filename = "Figures/HBG_SortedBy_PeakRPKM_H3K27me3_4quantile_boxplot.png",plot = HBG_boxplot,width = 6,height = 9)

HBG_scatter <- ggplot(data=HBG, aes(x = log10(peak_H3K27me3_RPKM+1), y = log10(peak_H3K4me3_RPKM+1))) + 
  geom_point()+ 
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = c('none'))+
  xlab("log10(H3K27me3_RPKM)")+
  ylab("log10(H3K4me3_RPKM)")
HBG_scatter
ggsave("Figures/HBG_histone_scatter.png",plot = HBG_scatter,width = 4,height = 4)
```
![HBG_SortedBy_PeakRPKM_H3K27me3_4quantile_boxplot](https://github.com/Blairewen/Arabidopsis-Pollen/blob/main/Figures/HBG_SortedBy_PeakRPKM_H3K27me3_4quantile_boxplot.png)
![HBG_histone_scatter](https://github.com/Blairewen/Arabidopsis-Pollen/blob/main/Figures/HBG_histone_scatter.png)

## 3. GO enrichment analysis (using HBG as the example)
```bash
library(clusterProfiler)
library(tidyverse)
library(org.At.tair.db)
library(GOplot)
library(ggplot2)
gene_FPKM <- read.table("Figures/gene_FPKM_in_12tissues.data",header = T,sep = "\t",quote = "")
SN_exp_genes <- subset(gene_FPKM, Spm >=0.5)
HBG <- read.table("Figures/HBG_SortedBy_PeakRPKM_H3K27me3_4quantile.data",sep = "\t",header = T,quote = "")
HBG_gene <- unique(HBG$annogene_ID)

HBG_GO = enrichGO(HBG_gene, 
                  universe = SN_exp_genes$gene_ID, 
                  OrgDb = "org.At.tair.db", 
                  keyType="TAIR", 
                  ont="BP",
                  qvalueCutoff = 0.1)
write.table(HBG_GO@result,file = "Figures/HBG_GO.txt",sep = "\t", col.names = T, row.names = F,quote = F)

HBG_GO_result <- read.table("Figures/HBG_GO.txt",sep = "\t", header = T, quote = "")
HBG_GO_plot <- 
  ggplot(HBG_GO_result[1:10,],aes(-log10(qvalue),reorder(Description,-log10(qvalue)))) +
  geom_point()+
  geom_point(aes(size=Count,color=-1*log10(qvalue)))+
  scale_color_gradient(low="blue",high = "red")+
  labs(color=expression(-log[10](q-value)),size="Count",x="-log10(q-value)",y="GO Term")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("HBG")
HBG_GO_plot
ggsave("Figures/HBG_GO.png",plot = HBG_GO_plot,width = 9,height = 10)
```
![HBG_GO](https://github.com/Blairewen/Arabidopsis-Pollen/blob/main/Figures/HBG_GO.png)

## 4. Histone signal meta plot (using HBG as the example, deeptools: computeMatrix)
```bash
computeMatrix scale-regions -p 8 \
                            -R HBG.bed \
                            -S Spm_H3K27me3_reps_treat_pileup.bw \
                               MS_H3K27me3_reps_treat_pileup.bw \
                               Sdl_H3K27me3_reps_treat_pileup.bw \
                               Spm_H3K4me3_reps_treat_pileup.bw \
                               MS_H3K4me3_reps_treat_pileup.bw \
                               Sdl_H3K4me3_reps_treat_pileup.bw \
                            -b 500 -a 500 -bs 50 \
                            --skipZeros 
                            --samplesLabel H3K27me3_Spm H3K27me3_MS H3K27me3_Sdl \  
                                           H3K4me3_Spm H3K4me3_MS H3K4me3_Sdl \
                            --outFileNameMatrix HBG_metaplot.tab \
                            -o HBG_metaplot.matrix.mat.gz

plotHeatmap -m HBG_metaplot.matrix.mat.gz \
            --outFileNameMatrix HBG_metaplot_heatmap.data \
            -o HBG_metaplot_heatmap.pdf \
            --colorMap Reds \
            --xAxisLabel "" \
            --regionsLabel HBG \
            --samplesLabel H3K27me3_Spm H3K27me3_MS H3K27me3_Sdl \
                           H3K4me3_Spm H3K4me3_MS H3K4me3_Sdl \
            --plotTitle "4597 Hot Bivalent Gene (HBG)"
# The "HBG_metaplot_heatmap.pdf" is the Fig3A we showed in this paper.
```

## 5. Heatmap (using Spm reset genes as the example)
```bash
library(pheatmap)
library(ggplot2)
Spm_reset_gene <- read.table("Figures/Spm_reset_gene_heatmap.data",header = T,sep = "\t",quote = "")
row.names(Spm_reset_gene) <- Spm_reset_gene$gene
Spm_reset_gene <- Spm_reset_gene[,-1]
Spm_reset_gene <- Spm_reset_gene[rowSums(Spm_reset_gene)!=0,]
Spm_reset_gene_heatmap <- pheatmap(Spm_reset_gene,
                          cluster_cols = F,
                          show_rownames = F,
                          cluster_rows  = T,
                          scale = "row",
                          angle_col = "45",
                          fontsize_row = 5,
                          cellwidth = NA,
                          cellheight = NA,
                          border_color = "NA",
                          main = "Resetting genes in Spm")
ggsave("Figures/Spm_reset_gene_heatmap.png",plot = Spm_reset_gene_heatmap,width = 10,height = 5)
```
![Spm_reset_gene_heatmap](https://github.com/Blairewen/Arabidopsis-Pollen/blob/main/Figures/Spm_reset_gene_heatmap.png)

## 6. Peak density plot in chromosomes (using ATAC data as the example)
```bash
# Berfor run "RIdeogram", user should prepared some files:
# (1) Karyotype file, which including the start, end, centromere region of each chromosome;
# (2) Peak_centre.gff file, which including chromosome number, genome name, 
#     peak class, peak centre location, peak score, peak strand, score, peak source.

require(RIdeogram)
Veg_ATAC_density <- GFFex(input = "Figures/ATAC_centre.gff", 
                          karyotype = "Files/TAIR10_karyotype.txt", 
                          feature = "Veg_ATAC", 
                          window = 10000)
summary(Veg_ATAC_density$Value)
write.table(Veg_ATAC_density,"Figures/Veg_ATAC_density.txt",col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")

Spm_ATAC_density <- GFFex(input = "Figures/ATAC_centre.gff", 
                          karyotype = "Files/TAIR10_karyotype.txt", 
                          feature = "Spm_ATAC", 
                          window = 10000)
summary(Spm_ATAC_density$Value)
Spm_ATAC_density$Value[which(Spm_ATAC_density$Value == 7)] = 10    # To make sure the color gradient is consistent with Spm 
write.table(Spm_ATAC_density,"Figures/Spm_ATAC_density.txt",col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")

tair10_karyotype <- read.table("Files/TAIR10_karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
ideogram(karyotype = tair10_karyotype, 
         overlaid = Veg_ATAC_density, 
         label = Spm_ATAC_density, 
         label_type = "heatmap", 
         width = 120, 
         Lx = 80, Ly = 25, 
         colorset1 = c("#f7f7f7", "#7e331f"), 
         colorset2 = c("#f7f7f7", "#7e331f"),
         output = "Figures/ATAC_ChrDensity.svg") 
convertSVG("Figures/ATAC_ChrDensity.svg", device = "png",file = "Figures/ATAC_ChrDensity.png")

# In the result of "Figures/ATAC_ChrDensity.png", there are 5 chromosomes. 
# In each chromosome, the left one represents the Veg_ATAC_density, and the right one represents the Spm_ATAC_density.
```
![ATAC_ChrDensity](https://github.com/Blairewen/Arabidopsis-Pollen/blob/main/Figures/ATAC_ChrDensity.png)