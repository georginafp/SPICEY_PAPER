# ================================== #
# Differential analysis - ATAC & RNA #
# ================================== #


# 0. Packages ------------------------------------------------------------------
.libPaths("/scratch/lab_lpasquali/shared_data/rstudio-singularity/packages/4.4.2")
library(Seurat)
library(dplyr)
library(tidyr)
library(regioneR)
library(purrr)
library(Signac)
library(future)



# 1. Load single cell dataset --------------------------------------------------
re <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/SUBSET_CTRL_HPAP_2.rds")




# 2. Helper function: clean and filter cell types ------------------------------
clean_and_filter <- function(obj, min_cells = 100) {
  obj$cell_type <- obj$integrated_clusters_names
  obj$cell_type <- gsub("_(A|B)", "", obj$cell_type)
  obj$cell_type <- gsub("(Act|Quies)_stellate", "Stellate", obj$cell_type)
  Idents(obj) <- obj$cell_type
  keep <- names(which(table(obj$cell_type) > min_cells))
  subset(obj, idents = keep)
}



# 3. Wraper FindMarkers function -----------------------------------------------
run_findmarkers <- function(obj, assay = c("ATAC","RNA"), min_cells = 3) {
  assay <- match.arg(assay)
  DefaultAssay(obj) <- assay
  valid_cell_types <- names(which(table(obj$cell_type) >= min_cells)) # valid cell types
  # Function to compute markers for a single cell type
  find_markers <- function(celltype) {
    n_cells <- sum(obj$cell_type == celltype)
    if(n_cells < min_cells) return(NULL)
    t <- FindMarkers(obj, ident.1 = celltype, logfc.threshold = 0, min.pct = 0.01) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "region") %>%
      mutate(cell_type = celltype)
    if(assay == "ATAC") {
      # For ATAC: split region string into chr, start, end
      t <- t %>%
        tidyr::separate(region, into = c("chr","start","end"), sep = "-", convert = TRUE) %>%
        dplyr::select(chr, start, end, cell_type, p_val_adj, p_val, avg_log2FC) %>%
        regioneR::toGRanges()
    }
    return(t)
  }
  future::plan("multicore", workers = 4)   # Run in parallel
  gr_list <- purrr::map(valid_cell_types, find_markers) %>% purrr::compact()
  names(gr_list) <- valid_cell_types
  return(gr_list)
}



# 4. Clean and filter cells ----------------------------------------------------
re <- clean_and_filter(re, min_cells = 100)




# 5. Run differential analysis -------------------------------------------------
# ATAC
da_atac <- run_findmarkers(re, assay = "ATAC")
saveRDS(da_atac, "/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/DA_ATAC_HPAP-CTRL.rds")

# RNA
da_rna <- run_findmarkers(re, assay = "RNA")
saveRDS(da_rna, "/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/DA_RNA_HPAP-CTRL.rds")
