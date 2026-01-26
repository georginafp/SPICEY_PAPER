## Differential analysis ----

### General setup ----
.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")

# Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(regioneR)
library(purrr)
library(ggplot2)
library(Signac)
library(future)
library(future.apply)   # more memory-friendly than purrr::map with future
library(DESeq2)

# Increase maximum global size for future
options(future.globals.maxSize = 6 * 1024^3)  # 6 GB

# Parallel plan
plan(multisession, workers = 4)  # safer on cluster than multicore

# Helper function: collapse cell types
collapse_celltypes <- function(re) {
  re$cell_type <- as.character(re$cell_type)
  collapse_map <- c(
    "CD14 Mono" = "Mono",
    "CD16 Mono" = "Mono",
    "CD4 Naive" = "CD4 T",
    "CD4 TCM"   = "CD4 T",
    "CD4 TEM"   = "CD4 T",
    "Treg"      = "CD4 T",
    "CD8 Naive" = "CD8 T",
    "CD8 TEM_1" = "CD8 T",
    "CD8 TEM_2" = "CD8 T",
    "NK"        = "NK",
    "Naive B"        = "B",
    "Memory B"       = "B",
    "Intermediate B" = "B",
    "gdT"  = "other T",
    "MAIT" = "other T",
    "cDC" = "DC",
    "pDC" = "DC"
  )
  re$cell_type <- dplyr::recode(re$cell_type, !!!collapse_map, .default = NA_character_)
  re <- subset(re, cells = which(!is.na(re$cell_type)))
  Idents(re) <- re$cell_type
  return(re)
}

### -----------------------------
### ATAC differential analysis ----
re <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/PBMC/data/final_pbmc_atac.rds")
re$cell_type <- re$seurat_annotations

# Keep cell types with >100 cells
celltypes_to_keep <- names(which(table(re$cell_type) > 100))
Idents(re) <- re$cell_type
re <- subset(re, idents = celltypes_to_keep)
DefaultAssay(re) <- "ATAC"

# Collapse cell types
re <- collapse_celltypes(re)

# Function to run differential analysis per cell type
findmarkers_atac <- function(ct, seurat_obj) {
  n_cells <- sum(as.character(Idents(seurat_obj)) == ct, na.rm = TRUE)
  if (n_cells < 3) return(NULL)

  # Run FindMarkers (LR test)
  res <- FindMarkers(
    object = seurat_obj,
    test.use = "negbinom",
    ident.1 = ct,
    logfc.threshold = 0,
    min.pct = 0.01
  )

  # Convert to GRanges
  res <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "region") %>%
    separate(region, into = c("chr", "start", "end"), sep = "-", convert = TRUE) %>%
    mutate(cell_type = ct) %>%
    dplyr::select(chr, start, end, cell_type, p_val_adj, p_val, avg_log2FC) %>%
    regioneR::toGRanges()

  return(res)
}

# Valid cell types
valid_cell_types <- names(which(table(re$cell_type) >= 3))

# Run in parallel using future_lapply (memory-friendly)
gr_list_atac <- future_lapply(valid_cell_types, findmarkers_atac, seurat_obj = re)
gr_list_atac <- gr_list_atac[!sapply(gr_list_atac, is.null)]

# Save results
saveRDS(gr_list_atac, "/homes/users/gfuentes/scratch/projects/spicey_paper/PBMC/data/NB_DA_ATAC_PBMC.rds")


### -----------------------------
### RNA differential analysis ----
re <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/PBMC/data/final_pbmc_rna.rds")
re$cell_type <- re$seurat_annotations

# Keep cell types with >100 cells
celltypes_to_keep <- names(which(table(re$cell_type) > 100))
Idents(re) <- re$cell_type
re <- subset(re, idents = celltypes_to_keep)
DefaultAssay(re) <- "RNA"

# Collapse cell types
re <- collapse_celltypes(re)

# Function to run differential expression per cell type
findmarkers_rna <- function(ct, seurat_obj) {
  n_cells <- sum(as.character(Idents(seurat_obj)) == ct, na.rm = TRUE)
  if (n_cells < 3) return(NULL)

  res <- FindMarkers(
    object = seurat_obj,
    test.use = "negbinom",
    ident.1 = ct,
    logfc.threshold = 0,
    min.pct = 0.01
  )

  res <- res %>%
    as.data.frame() %>%
    mutate(cell_type = ct)

  return(res)
}

# Valid cell types
valid_cell_types <- names(which(table(re$cell_type) >= 3))

# Run in parallel
gr_list_rna <- future_lapply(valid_cell_types, findmarkers_rna, seurat_obj = re)
gr_list_rna <- gr_list_rna[!sapply(gr_list_rna, is.null)]
names(gr_list_rna) <- valid_cell_types

# Save results
saveRDS(gr_list_rna, "/homes/users/gfuentes/scratch/projects/spicey_paper/PBMC/data/NB_DA_RNA_PBMC.rds")
