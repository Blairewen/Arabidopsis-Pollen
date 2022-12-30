# Distinct chromatin signatures in Arabidopsis male gametophyte

## Introduction
The code in this repository contains the commands and pipeline to reproduce the data from the paper "Distinct chromatin signatures in Arabidopsis male gametophyte" (under review). The code of main figures in this paper can be found in `Figures` directory, and some necessary files can be found in `Files` directory  

## Raw Sequencing Data
The raw data (fastq files) are deposited into ArrayExpress under the accession number E-MTAB-10965, E-MTAB-10966, E-MTAB-10967, E-MTAB-12137.
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10965
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10966
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10967
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12137

## Requirements

The software packages needed for the anlaysis are:

```
python2
R v4.0.3
fastp v0.20.0
cutadapt v2.10
histat2 v2.1.0
samtools v1.9
picard v2.20.8
macs2 v2.1.4
homer v4.11.1
bedtools v2.29.2
bdg2bw_RemoveOverlap.sh
cuffdiff v2.2.1
STAR v2.7.9a
```

All of them can be found in the `Packages` directory in this repository.

## 1. ChIP-seq pipeline (using Spm_H3K4me3_rep1 as an example)
### 1.1 Data quality control
```bash
fastp -i Spm_H3K4me3_rep1_R1.fastq.gz \
      -I Spm_H3K4me3_rep1_R2.fastq.gz \
      -o Spm_H3K4me3_rep1_R1.fastp.fastq.gz \
      -O Spm_H3K4me3_rep1_R2.fastp.fastq.gz \
      -w 8 -l 25 --detect_adapter_for_pe \
      -j Spm_H3K4me3_rep1.fastp.json \
      -h Spm_H3K4me3_rep1.fastp.html
```

### 1.2 Genome mapping
```bash
hisat2 -p 8 \
       -x tair10_index \
       --no-temp-splicesite \
       --no-spliced-alignment \
       --summary-file Spm_H3K4me3_rep1.summary \
       -1 Spm_H3K4me3_rep1_R1.fastp.fastq.gz \
       -2 Spm_H3K4me3_rep1_R2.fastp.fastq.gz | \
       samtools view -ShuF 4 -q 30 -f 2 -@ 2 - | \
       samtools sort -@ 2 -o Spm_H3K4me3_rep1.sorted.bam -
```

### 1.3 Redundancy removal
```bash
java -Xmx5g \
     -XX:ParallelGCThreads=8 \
     -jar picard.jar MarkDuplicates \
     I= Spm_H3K4me3_rep1.sorted.bam \
     O= Spm_H3K4me3_rep1.sorted.picardMD.bam \
     M= Spm_H3K4me3_rep1.sorted.picardMD.txt \
     REMOVE_DUPLICATES=true 
```

### 1.4 Peak calling
```bash
macs2 callpeak -c Spm_H3K4me3_input_rep1.sorted.picardMD.bam \
               -t Spm_H3K4me3_rep1.sorted.picardMD.bam \
               -f BAMPE -g 1.053e8 --keep-dup all \
               -q 0.01 -n Spm_H3K4me3_rep1 \
               --outdir macs2/ -B --SPMR

# Note for H3K27me3 and H3K9me2, we using “--broad --nolambda” to replace “-q 0.01”
macs2 callpeak -c Spm_H3K27me3_input_rep1.sorted.picardMD.bam \
               -t Spm_H3K27me3_rep1.sorted.picardMD.bam \
               -f BAMPE -g 1.053e8 --keep-dup all \
               --broad --nolambda -n Spm_H3K27me3_rep1 \
               --outdir macs2/ -B --SPMR 


macs2 callpeak -c Spm_H3K9me2_input_rep1.sorted.picardMD.bam \
               -t Spm_H3K9me2_rep1.sorted.picardMD.bam \
               -f BAMPE -g 1.053e8 --keep-dup all \
               --broad --nolambda -n Spm_H3K9me2_rep1 \
               --outdir macs2/ -B --SPMR 
```

### 1.5 Find high confidence peaks (homer: mergePeaks)
```bash
# We do this by only keeping peaks that are presentin in all replicates
mergePeaks -d given \
           -venn mergePeaks_narrow_H3K4me3_Spm_given.venn \
           -prefix mergePeaks_narrow_H3K4me3_Spm_given \
           -matrix mergePeaks_narrow_H3K4me3_Spm_given \
           Spm_H3K4me3_rep1.narrowPeak Spm_H3K4me3_rep2.narrowPeak Spm_H3K4me3_rep 3.narrowPeak
           
# Note change the overlap peak file's name 
# "mergePeaks_narrow_H3K4me3_Spm_given_Spm_H3K4me3_rep1.narrowPeak_Spm_H3K4me3_rep2.narrowPeak_Spm_H3K4me3_rep 3.narrowPeak" 
# as "merged_narrow_Spm_H3K4me3.peak", which is the high confidence peak file we used below

# H3K9me2 often exists as very broad and continuous domains,
# therefore, we first stitch peaks in each replicates within 5000bp as one peaks, which using “mergePeaks -d 5000”,
# and then only keep the stitched peaks that are present in all replicates.
mergePeaks -d 5000 \
           -venn stitched_Spm_H3K9me2_rep1.venn \
           -prefix stitched_Spm_H3K9me2_rep1 \
           -matrix stitched_Spm_H3K9me2_rep1 \
           Spm_H3K9me2_rep1.broadPeak

mergePeaks -d 5000 \
           -venn stitched_Spm_H3K9me2_rep2.venn \
           -prefix stitched_Spm_H3K9me2_rep2 \
           -matrix stitched_Spm_H3K9me2_rep2 \
           Spm_H3K9me2_rep2.broadPeak

mergePeaks -d given \
           -venn mergePeaks_stitched_Spm_H3K9me2_given.venn \
           -prefix mergePeaks_stitched_Spm_H3K9me2_given \
           -matrix mergePeaks_stitched_Spm_H3K9me2_given \
           stitched_Spm_H3K9me2_rep1.peak stitched_Spm_H3K9me2_rep2.peak
```

### 1.6 Annotate peaks (homer: annotatePeaks.pl)
```bash
# Find the nearest gene of each peak
annotatePeaks.pl merged_narrow_Spm_H3K4me3.peak \
                 tair10 -gtf TAIR10.gtf \
                 -annStats merged_narrow_Spm_H3K4me3_tair10.annStats > \
                 merged_narrow_Spm_H3K4me3_tair10.annotate

# The result "merged_narrow_Spm_H3K4me3_tair10.annotate" contain the annotation information of each peaks, 
# then we remove the genes which annotated as "Intergenic",
# and remove the genes which located in mitochondrial chromosome and chloroplast chromosome,
# the rest of the genes as the histone marked genes we used below.
```

### 1.7 Prepare bigwig files for track visualization
```bash
# We using all replicates to call peak together, and then convert the bedgraph file to the bigwig file for visualization
macs2 callpeak -c Spm_H3K4me3_input_rep1.sorted.picardMD.bam \
                  Spm_H3K4me3_input_rep2.sorted.picardMD.bam \
                  Spm_H3K4me3_input_rep3.sorted.picardMD.bam \
               -t Spm_H3K4me3_rep1.sorted.picardMD.bam \
                  Spm_H3K4me3_rep2.sorted.picardMD.bam \
                  Spm_H3K4me3_rep3.sorted.picardMD.bam \
               -f BAMPE -g 1.053e8 --keep-dup all \
               -q 0.01 -n Spm_H3K4me3_reps \
               --outdir macs2/ -B --SPMR

bdg2bw_RemoveOverlap.sh Spm_H3K4me3_reps_treat_pileup.bdg tair10.chrom.sizes
```

## 2. ATAC-seq pipeline (using Spm_ATAC_rep1 as an example)
### 2.1 Data quality control
```bash
fastp -i Spm_ATAC_rep1_R1.fastq.gz \
      -I Spm_ATAC_rep1_R2.fastq.gz \
      -o Spm_ATAC_rep1_R1.fastp.fastq.gz \
      -O Spm_ATAC_rep1_R2.fastp.fastq.gz \
      -w 8 -l 25 --detect_adapter_for_pe \
      -j Spm_ATAC_rep1.fastp.json \
      -h Spm_ATAC_rep1.fastp.html
```

### 2.2 Genome mapping
```bash
hisat2 -p 8 \
       -x tair10_index \
       -X 2000 \
       --no-temp-splicesite \
       --no-spliced-alignment \
       --summary-file Spm_ATAC.summary \
       -1 Spm_ATAC_rep1_R1.fastp.fastq.gz \
       -2 Spm_ATAC_rep1_R2.fastp.fastq.gz | \
       samtools view -ShuF 4 -q 30 -f 2 -@ 2 - | \
       samtools sort -@ 2 -o Spm_ATAC_rep1.sorted.bam -
```

### 2.3 Redundancy removal
```bash
java -Xmx5g \
     -XX:ParallelGCThreads=8 \
     -jar picard.jar MarkDuplicates \
     I= Spm_ATAC_rep1.sorted.bam \
     O= Spm_ATAC_rep1.sorted.picardMD.bam \
     M= Spm_ATAC_rep1.sorted.picardMD.txt \
     REMOVE_DUPLICATES=true 
```

### 2.4 Convert the BAM format to BED format (bedtools: bamToBed)
```bash
bamToBed -i Spm_ATAC_rep1.sorted.picardMD.bam > Spm_ATAC_rep1.sorted.picardMD.bed
```

### 2.5 Peak calling
```bash
macs2 callpeak -c SRR4000479_gDNA_ATAC.sorted.picardMD.bed \
               -t Spm_ATAC_rep1.sorted.picardMD.bed \
               -f BED -g 1.053e8 --keep-dup all \
               -n Spm_ATAC_rep1 --outdir macs2/ \
               -B --SPMR --nomodel \
               --shift -100 --extsize 200
```

### 2.6 Filter for validated peaks (homer: mergePeaks)
```bash
# We do this by only keeping peaks that are presentin in all replicates
mergePeaks -d given \
           -venn mergePeaks_narrow_ATAC_Spm_given.venn \
           -prefix mergePeaks_narrow_ATAC_Spm_given \
           -matrix mergePeaks_narrow_ATAC_Spm_given \
           Spm_ATAC_rep1.narrowPeak Spm_ATAC_rep2.narrowPeak

# Note change the overlap peak file's name 
# "mergePeaks_narrow_ATAC_Spm_given_Spm_ATAC_rep1.narrowPeak_Spm_ATAC_rep2.narrowPeak" 
# as "merged_narrow_Spm_ATAC.peak", which is the high confidence peak file we used below
```

### 2.7 Annotate peaks (homer: annotatePeaks.pl)
```bash
annotatePeaks.pl merged_narrow_Spm_ATAC.peak \
                 tair10 -gtf TAIR10.gtf \
                 -annStats merged_narrow_Spm_ATAC_tair10.annStats > \
                 merged_narrow_Spm_ATAC_tair10.annotate
```

### 2.8 Motif enrichment analysis (homer: findMotifsGenome.pl)
```bash
findMotifsGenome.pl merged_narrow_Spm_ATAC.peak \
                    -p 8 tair10 Spm_ATAC_motif/ 
```

### 2.9 Prepare bigwig files for track visualization
```bash
# We using all replicates to call peak together, and then convert the bedgraph file to the bigwig file for visualization
macs2 callpeak -c SRR4000479_gDNA_ATAC.sorted.picardMD.bed \
               -t Spm_ATAC_rep1.sorted.picardMD.bed \
                  Spm_ATAC_rep2.sorted.picardMD.bed \
               -f BED -g 1.053e8 --keep-dup all \
               -n Spm_ATAC_reps --outdir macs2/ \
               -B --SPMR --nomodel \
               --shift -100 --extsize 200


# Because the chromosomes of vegetative cells are looser than those of sperm cells,
# so that the average level of chromosome accessibility of Veg is weakened at the same sequencing depth.
# Therefore, it is necessary to normalize the track of ATAC to avoid visual errors as much as possible,
# the normalization method is is to refer to the article (https://www.nature.com/articles/nature19360)
# "Broad histone H3K4me3 domains in mouse oocytes modulate maternal-to-zygotic transition",
# the normalization details as follows: (using Spm_ATAC_reps_treat_pileup.bdg as the example)
# Step 1.Calculate the RPKM of all the Spm_ATAC peaks (method see the section of "Region RPKM" below);
# Step 2.Pick the top 1000th peak’s RPKM 41.37 as normalized value;
# Step 3.The value in the bdg file divided by 41.37;
# Step 4.Prepare bigwig files for track visualization by normalized bdg files: Spm_ATAC_reps_treat_pileup_normalized.bdg

bdg2bw_RemoveOverlap.sh Spm_ATAC_reps_treat_pileup_normalized.bdg tair10.chrom.sizes
```

## 3. RNA-seq pipeline (using Spm_rep1 as an example)
### 3.1 Data quality control
```bash
fastp -i Spm_rep1_R1.fastq.gz \
      -I Spm_rep1_R2.fastq.gz \
      -o Spm_rep1_R1.fastp.fastq.gz \
      -O Spm_rep1_R2.fastp.fastq.gz \
      -w 8 -q 20 –u 20 \
      -j Spm_rep1.fastp.json \
      -h Spm_rep1.fastp.html

cutadapt -j 24 \
         -g ^ATTGCGCAATGNNNNNNNNGGG \
         -o Spm_rep1_R1.fastp.cutadapt.fastq.gz \
         -p Spm_rep1_R2.fastp.cutadapt.fastq.gz \
         Spm_rep1_R1.fastp.fastq.gz Spm_rep1_R2.fastp.fastq.gz
```

### 3.2 Genome mapping
```bash
hisat2 -p 8 \
       -x tair10_index \
       --summary-file Spm.summary \
       --dta-cufflinks \
       -1 Spm_rep1_R1.fastp.cutadapt.fastq.gz \
       -2 Spm_rep1_R2.fastp.cutadapt.fastq.gz | \
       samtools view -ShuF 4 -q 30 -f 2 -@ 2 - | \
       samtools sort -@ 2 -o Spm_rep1.sorted.bam -
```

### 3.3 Convert the BAM format to bigwig format for track visualization
```bash
genomeCoverageBed -ibam Spm_rep1.sorted.bam -bg -split > RS_Spm_rep1.bdg

bdg2bw_RemoveOverlap.sh RS_Spm_rep1.bdg tair10.chrom.sizes
```

### 3.4 Get the gene expression matrix
```bash
cuffdiff -o cuffdiff_results/ \
         -L Spm,Veg,MS,Sdl -p 16 \
         -b tair10_genome.fa -u \
         --library-type fr-unstranded \
         TAIR10.gtf \
         Spm_rep1.sorted.bam,Spm_rep2.sorted.bam,Spm_rep3.sorted.bam \
         Veg_rep1.sorted.bam,Veg_rep2.sorted.bam,Veg_rep3.sorted.bam \
         MS_rep1.sorted.bam,MS_rep2.sorted.bam,MS_rep3.sorted.bam \ 
         Sdl_rep1.sorted.bam,Sdl_rep2.sorted.bam

# The file "gene_exp.diff" contains the genes average FPKM of different cell types,
# and the file "genes.read_group_tracking" contains genes FPKM of each replicates.         
```

## 4. scRNA-seq pipeline
### 4.1 Build the STAR index of TAIR10
```bash
STAR-2.7.9a/source/STAR --runThreadN 20 \
                        --runMode genomeGenerate \
                        --genomeDir STAR_index/tair10/ \
                        --genomeFastaFiles tair10_genome.fa \
                        --sjdbGTFfile TAIR10.gtf
```

### 4.2 Genome mapping
```bash
# PlateA
STAR-2.7.9a/source/STAR --genomeDir STAR_index/tair10/ \
                        --readFilesCommand zcat \
                        --readFilesIn PlateA_R2.fastq.gz PlateA_R1.fastq.gz \
                        --soloCBstart 1 \
                        --soloCBlen 8 \
                        --soloUMIstart 9 \
                        --soloUMIlen 10 \
                        --soloType CB_UMI_Simple \
                        --soloCBwhitelist whitelist.csv \
                        --runThreadN 16 \
                        --soloBarcodeReadLength 0 \
                        --clip3pNbases 116 \
                        --outSAMattributes CB UB \
                        --outSAMtype BAM SortedByCoordinate

# PlateB 
STAR-2.7.9a/source/STAR --genomeDir STAR_index/tair10/ \
                        --readFilesCommand zcat \
                        --readFilesIn PlateB_R2.fastq.gz PlateB_R1.fastq.gz \
                        --soloCBstart 1 \
                        --soloCBlen 8 \
                        --soloUMIstart 9 \
                        --soloUMIlen 10 \
                        --soloType CB_UMI_Simple \
                        --soloCBwhitelist whitelist.csv \
                        --runThreadN 16 \
                        --soloBarcodeReadLength 0 \
                        --clip3pNbases 116 \
                        --outSAMattributes CB UB \
                        --outSAMtype BAM SortedByCoordinate
```

### 4.3 Analysis scRNA-seq mapping results
Firstly, saving the results of genome mapping by samples. Here we have 2 samples: PlateA and PlateB. In this paper, we used 384 well plate to perform the single-cell RNA-seq. In each plate, vegtative nuclei were on the left side of the plate, and sperm nuclei were on the right side.

```bash
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(dplyr)
library(readr)
library(stringr)
library(ggcorrplot)
library(tidyverse)
library(harmony)

set.seed(1234)
# 1.Creat seurat object
samples=list.files("scRNA-seq/")
dir <- file.path('./scRNA-seq',samples)
counts <- Read10X(data.dir = dir)
scRNA = CreateSeuratObject(counts, min.cells=1)

# 2.Quality control
scRNA[["log10_nCount_RNA"]] <- log10(scRNA$nCount_RNA)
scRNA[["log10_nFeature_RNA"]] <- log10(scRNA$nFeature_RNA)
scRNA <- subset(scRNA,subset = log10_nFeature_RNA >= 2.5 )
afterQC_vlnplot <- VlnPlot(scRNA,
                           features = c("log10_nFeature_RNA","log10_nCount_RNA"),
                           ncol = 3, 
                           pt.size = 0.1)
print(afterQC_vlnplot)
ggsave("scRNA-seq/after_scRNAQC_vlnplot.png",plot = afterQC_vlnplot)

# 3.NormalizeData
scRNA <- NormalizeData(scRNA,normalization.method = "LogNormalize")
all.genes = rownames(scRNA)
scRNA = ScaleData(scRNA,features = all.genes)

# 4.Linear dimension reduction
scRNA <- FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000)
scRNA <- RunPCA(scRNA,features = VariableFeatures(object = scRNA),verbose = FALSE)

# 5.Clustering
scRNA <- RunHarmony(scRNA, "orig.ident" , verbose = T)
scRNA <- FindNeighbors(scRNA,dims = 1:20)
scRNA <- FindClusters(scRNA,resolution = 0.4)
scRNA <- RunUMAP(scRNA, reduction = "harmony" , dims = 1:20)
harmony_plot <- DimPlot(scRNA , reduction = "umap")
harmony_plot
ggsave("scRNA-seq/UMAP.png",plot = harmony_plot)

# 6.Mark genes plot
genes_sample <- FeaturePlot(scRNA, ncol = 3,
                            features = c("AT4G30860","AT1G08600","AT5G55490","AT4G35700","AT3G62230","AT5G02110",
                                         "AT3G04690","AT5G28680","AT4G39110","AT3G12690","AT5G60740","AT4G24580",
                                         "log10_nCount_RNA"))
ggsave("scRNA-seq/MarkGenes_UMAP.pdf",plot = genes_sample,width = 25,height = 30)

saveRDS(scRNA,file = "scRNA-seq/scRNA.rds")
```
![QC](https://github.com/Blairewen/Arabidopsis-Pollen/blob/main/Figures/scRNA-seq/after_scRNAQC_vlnplot.png)
![umap](https://github.com/Blairewen/Arabidopsis-Pollen/blob/main/Figures/scRNA-seq/UMAP.png)

## 5. Region RPKM (calculating peaks' RPKM of Spm_ATAC in Spm as the example)
### 5.1 SampleReads: 
```bash
samtools merge -@ 16 \
               Spm_ATAC.sorted.picardMD.merged.bam \
               Spm_ATAC_rep1.sorted.picardMD.bam \
               Spm_ATAC_rep2.sorted.picardMD.bam

samtools view -c Spm_ATAC.sorted.picardMD.merged.bam
> 29833118
# The result of "29833118" is the "SampleReads" 
```

$$
SampleReadsDepth = \cfrac{SampleReads}{1,000,000}
$$

The SampleReadsDepth of Spm_ATAC is 29.83

### 5.2 RegionLength

$$
RegionLength = PeakEnd - PeakStart
$$

### 5.3 ReadsWithinRegion
```bash
# Counting the reads on each peaks
pos2bed.pl merged_narrow_Spm_ATAC.peak \
           -o merged_narrow_Spm_ATAC.bed

coverageBed -a merged_narrow_Spm_ATAC.bed \
            –b Spm_ATAC.sorted.picardMD.merged.bam > \
            Spm_ATAC_coverageBed.txt

# The value of 7th column in "Spm_ATAC_coverageBed.txt" is the reads count on the peaks,
# which is labeled as "ReadsWithinRegion" below.
```

### 5.4 The RPKM of peaks were calculated as:
$$
PeakRPKM= \cfrac{ReadsWithinRegion}{\cfrac{RegionLength}{1000} * SampleReadsDepth} 
$$

Note: Users also can calculate RPKM by using script "peak_RPKM.pl" we provide in the `Packages` directory:

```bash
peak_RPKM.pl -depth 29.83 \
             --col 7 \
             --input Spm_ATAC_coverageBed.txt \
             --output Spm_ATAC_coverageBed_RPKM.txt

# The RPKM of Spm_ATAC peaks in Spm were saved in the last column of output file "Spm_ATAC_coverageBed_RPKM.txt".
```

## 6. HBG and CBG identification
In this paper, we defined HBG and CBG as following steps:
### (1) Find the bivalent peaks in Spm  (homer: mergePeaks)
```bash
mergePeaks -d given \
           -venn mergePeaks_bivalent_Spm_given.venn \
           -prefix mergePeaks_bivalent_Spm_given \
           -matrix mergePeaks_bivalent_Spm_given \
           merged_narrow_Spm_H3K4me3.peak merged_broad_Spm_H3K27me3.peak

# Note change the overlap peak file's name 
# "mergePeaks_bivalent_Spm_given_merged_narrow_Spm_H3K4me3.peak_merged_broad_Spm_H3K27me3.peak" 
# as "Spm_bivalent.peak", which is the bivalent peak file we used below.
```

### (2) Annotate bivalent peaks (homer: annotatePeaks.pl)
```bash
annotatePeaks.pl Spm_bivalent.peak \
                 tair10 -gtf TAIR10.gtf \
                 -annStats Spm_bivalent_tair10.annStats > \
                 Spm_bivalent_tair10.annotate
```

### (3) Filter bivalent genes
Accroding annotate information in "Spm_bivalent_tair10.annotate", we removed the genes which annotated as "Intergenic", and removed the genes which located in mitochondrial chromosome and chloroplast chromosome. The rest genes were Spm bivalent genes(BG).

### (4) Define HBG and CBG combined with Spm RNA-seq
HBG(Hot Bivalent Genes): BG_FPKM >= 0.5

CBG(Cold Bivalent Genes): BG_FPKM < 0.5

## 7. H3K27me3 reset genes idenfitication
In this paper, we defined H3K27me3 reset genes as following steps:
### (1) Find high confidence H3K27me3 peaks of seedling(Sdl)  
```bash
mergePeaks -d given \
           -venn mergePeaks_broad_H3K27me3_Sdl_given.venn \
           -prefix mergePeaks_broad_H3K27me3_Sdl_given \
           -matrix mergePeaks_broad_H3K27me3_Sdl_given \
           Sdl_H3K27me3_rep1.broadPeak Sdl_H3K27me3_rep2.broadPeak Sdl_H3K27me3_rep3.broadPeak
           
# Note change the overlap peak file's name 
# "mergePeaks_broad_H3K27me3_Sdl_given_Sdl_H3K27me3_rep1.broadPeak_Sdl_H3K27me3_rep2.broadPeak_Sdl_H3K27me3_rep3.broadPeak" 
# as "merged_broad_Sdl_H3K27me3.peak", which is the high confidence peak file we used below

pos2bed.pl merged_broad_Sdl_H3K27me3.peak \
           -o merged_broad_Sdl_H3K27me3.bed
```

### (2) Calculate the H3K27me3_RPKM of all the Sdl H3K27me3 peaks (method see the section of "Region RPKM")
```bash
# Peaks' H3K27me3_RPKM in Sdl
coverageBed -a merged_broad_Sdl_H3K27me3.bed \
            –b Sdl_H3k27me3.sorted.picardMD.merged.bam > \
            merged_broad_Sdl_H3K27me3_to_Sdlbam_coverageBed.txt

peak_RPKM.pl -depth 51.20 \
             --col 7 \
             --input merged_broad_Sdl_H3K27me3_to_Sdlbam_coverageBed.txt \
             --output merged_broad_Sdl_H3K27me3_to_Sdlbam_coverageBed_RPKM.txt


# Peaks' H3K27me3_RPKM in Spm
coverageBed -a merged_broad_Sdl_H3K27me3.bed \
            –b Spm_H3k27me3.sorted.picardMD.merged.bam > \
            merged_broad_Sdl_H3K27me3_to_Spmbam_coverageBed.txt

peak_RPKM.pl -depth 36.75 \
             --col 7 \
             --input merged_broad_Sdl_H3K27me3_to_Spmbam_coverageBed.txt \
             --output merged_broad_Sdl_H3K27me3_to_Spmbam_coverageBed_RPKM.txt
```

### (3) Filter peaks with stable H3K27me3 modification in Sdl
Here, we picked up the Sdl_H3K27me3_peaks with H3K27me3_RPKM in Sdl higher than 20, and called them as Sdl_high_H3K27me3_peaks below.

### (4) Define the reset peaks
Find the Sdl_high_H3K27me3_peaks which Sdl_H3K27me3_RPKM higher 4 times than its Spm_H3K27me3_RPKM as the H3K27me3 reset peaks.

$$
FC = \cfrac{SdlH3K27me3RPKM}{SpmH3K27me3RPKM}
$$

### (5) Annotate the reset peaks
```bash
# Find the nearest gene of each peak
annotatePeaks.pl reset.peaks \
                 tair10 -gtf TAIR10.gtf \
                 -annStats reset_tair10.annStats > \
                 reset_tair10.annotate
# The annotated genes in "reset_tair10.annotate" are the reset genes we defined in Spm.
```
