## Differential analysis ----

### ATAC ----
.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")

library(Seurat)
library(dplyr)
library(tidyr)
library(regioneR)
library(purrr)
library(ggplot2)
library(Signac)
library(future)
library(DESeq2)

# Load dataset
re <- pbmc.atac
re$cell_type <- re$seurat_annotations

# Keep only cell types with >100 cells
celltypes_to_keep <- names(which(table(re$cell_type) > 100))
Idents(re) <- re$cell_type
re <- subset(re, idents = celltypes_to_keep)
DefaultAssay(re) <- "ATAC"

# -----------------------------
# Collapse cell types (17 -> 7)
# -----------------------------
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

# Apply mapping
re$cell_type <- dplyr::recode(re$cell_type, !!!collapse_map, .default = NA_character_)

# Remove cells that became NA
re <- subset(re, cells = which(!is.na(re$cell_type)))

# Set collapsed cell type as identity
Idents(re) <- re$cell_type

# -----------------------------
# Parallel setup
# -----------------------------
plan("multicore", workers = 4)

# -----------------------------
# Function to run differential analysis per cell type
# -----------------------------
findmarkers <- function(re, celltype) {

  # Skip if celltype not present
  if (!(celltype %in% as.character(Idents(re)))) {
    message(paste("Skipping:", celltype, "because it is not present"))
    return(NULL)
  }

  # Skip if fewer than 3 cells
  n_cells <- sum(as.character(Idents(re)) == celltype, na.rm = TRUE)
  if (n_cells < 3) {
    message(paste("Skipping:", celltype, "due to insufficient cells:", n_cells))
    return(NULL)
  }

  # Run FindMarkers
  t <- FindMarkers(
    re,
    ident.1 = celltype,
    logfc.threshold = 0,
    min.pct = 0.01
  )

  t <- t %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "region") %>%
    separate(region, into = c("chr", "start", "end"), sep = "-", convert = TRUE) %>%
    mutate(
      cell_type = celltype
    ) %>%
    dplyr::select(chr, start, end, cell_type, p_val_adj, p_val, avg_log2FC) %>%
    regioneR::toGRanges()

  return(t)
}

# -----------------------------
# Define valid cell types (≥3 cells)
# -----------------------------
valid_cell_types <- re$cell_type %>%
  table() %>%
  .[. >= 3] %>%
  names()

# -----------------------------
# Run differential analysis
# -----------------------------
gr_list <- purrr::map(valid_cell_types, ~ findmarkers(re, .x)) %>%
  purrr::compact()  # remove NULLs

# Save results
saveRDS(gr_list, "/homes/users/gfuentes/scratch/projects/spicey_paper/PBMC/data/DA_ATAC_PBMC.rds")


### RNA ----
.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")

library(Seurat)
library(dplyr)
library(purrr)
library(future)
library(ggplot2)

# Load dataset
re <- pbmc.rna
re$cell_type <- re$seurat_annotations

# Keep only cell types with >100 cells
celltypes_to_keep <- names(which(table(re$cell_type) > 100))
Idents(re) <- re$cell_type
re <- subset(re, idents = celltypes_to_keep)
DefaultAssay(re) <- "RNA"

# -----------------------------
# Collapse cell types (17 -> 7)
# -----------------------------
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

# Apply mapping
re$cell_type <- dplyr::recode(re$cell_type, !!!collapse_map, .default = NA_character_)

# Remove cells that became NA
re <- subset(re, cells = which(!is.na(re$cell_type)))

# Set collapsed cell type as identity
Idents(re) <- re$cell_type

# -----------------------------
# Parallel setup
# -----------------------------
plan("multicore", workers = 4)

# -----------------------------
# Function to run differential expression per cell type
# -----------------------------
findmarkers <- function(re, celltype) {

  # Skip if celltype not present
  if (!(celltype %in% as.character(Idents(re)))) {
    message(paste("Skipping:", celltype, "because it is not present"))
    return(NULL)
  }

  # Skip if fewer than 3 cells
  n_cells <- sum(as.character(Idents(re)) == celltype, na.rm = TRUE)
  if (n_cells < 3) {
    message(paste("Skipping:", celltype, "due to insufficient cells:", n_cells))
    return(NULL)
  }

  # Run FindMarkers
  t <- FindMarkers(
    re,
    ident.1 = celltype,
    logfc.threshold = 0,
    min.pct = 0.01
  )

  # Add cell type column
  t <- t %>%
    as.data.frame() %>%
    mutate(cell_type = celltype)

  return(t)
}

# -----------------------------
# Define valid cell types (≥3 cells)
# -----------------------------
valid_cell_types <- re$cell_type %>%
  table() %>%
  .[. >= 3] %>%
  names()

# -----------------------------
# Run differential expression
# -----------------------------
gr_list <- purrr::map(valid_cell_types, ~ findmarkers(re, .x)) %>%
  purrr::compact()  # remove NULLs
names(gr_list) <- valid_cell_types

# Save results
saveRDS(gr_list, "/homes/users/gfuentes/scratch/projects/spicey_paper/PBMC/data/DA_RNA_PBMC.rds")
