# ====================================== #
# CO-ACCESSIBILITY LINKS FROM scATAC-seq #
# ====================================== #


# 0. Packages ------------------------------------------------------------------
.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")
# install.packages("Seurat")
# install.packages("Signac")


library(Seurat)
library(Signac)
library(SeuratWrappers)
library(cicero)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(tibble)
library(regioneR)
library(stringr)



# 1. Load single cell dataset --------------------------------------------------
so <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/PBMC/data/final_pbmc_atac.rds")
DefaultAssay(so) <- "ATAC"





# 2. Data format ---------------------------------------------------------------
## To Monocle3
cds <- SeuratWrappers::as.cell_data_set(x = so)
## To cicero CDS
cds.cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP.ATAC)

## Get the chromosome sizes from the Seurat object
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
genome <- genome[paste0("chr", 1:22)]
genome.df <- data.frame("chr" = names(genome), "length" = genome) # to dataframe





# 3. Run cicero ----------------------------------------------------------------
conns <- run_cicero(cds.cicero, genomic_coords = genome.df, sample_num = 100)
saveRDS(conns, "/homes/users/gfuentes/scratch/projects/spicey_paper/PBMC/data/PBMC_COACC_LINKS.rds")
