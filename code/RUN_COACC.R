# ====================================== #
# CO-ACCESSIBILITY LINKS FROM scATAC-seq #
# ====================================== #



# 0. Packages ------------------------------------------------------------------
.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")
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
so <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/SUBSET_CTRL_HPAP_2.rds")
DefaultAssay(so) <- "ATAC"





# 2. Data format ---------------------------------------------------------------
## To Monocle3
cds <- SeuratWrappers::as.cell_data_set(x = so)
## To cicero CDS
cds.cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$UMAP.HAR)
## Get the chromosome sizes from the Seurat object
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
genome <- genome[paste0("chr", 1:22)]
genome.df <- data.frame("chr" = names(genome), "length" = genome) # to dataframe





# 3. Run cicero ----------------------------------------------------------------
conns <- run_cicero(cds.cicero, genomic_coords = genome.df, sample_num = 100)
saveRDS(conns, "/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/CTRL_LINKS_2.rds")


ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(so[["ATAC"]]) <- links
saveRDS(so, "/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/HPAP_SUBSET_CLEAN_LINKS_TH0.8.rds")
