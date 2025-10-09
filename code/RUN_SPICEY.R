# =============================== #
# Tissue specificity using SPICEY #
# =============================== #


# 0. Packages ------------------------------------------------------------------
library(GenomicRanges)
library(SPICEY)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


# 1. Load data -----------------------------------------------------------------
## 1.1. Differential analysis data derived from code/RUN_DIFF_ANALYSIS.R
atac <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/DA_ATAC_HPAP-CTRL.rds")
rna <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/DA_RNA_HPAP-CTRL.rds")

## 1.2. Co-accessible links derived from code/RUN_COACC.R
links <- readRDS("/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/CTRL_LINKS_2.rds")




# 2. Data formatting -----------------------------------------------------------
atac <- setNames(
  lapply(atac, function(x) {
    x <- as.data.frame(x)
    x$region_id <- paste0(x$seqnames, "-", x$start, "-", x$end)
    GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  }),
  vapply(atac, function(x) unique(as.character(mcols(x)$cell_type)), character(1))
)
rna <- lapply(rna, function(df) { rownames_to_column(df, var = "gene_id")})
links <- links %>%
  dplyr::filter(coaccess > 0.5)
peaks <- unlist(GRangesList(atac))
names(peaks) <- NULL




# 3. Annotate regions to target gene -------------------------------------------
## 3.1. Through proximity
annotation_near <- annotate_with_nearest(peaks = peaks,
                                         txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                         annot_dbi = org.Hs.eg.db,
                                         protein_coding_only = TRUE,
                                         verbose = TRUE,
                                         add_tss_annotation = FALSE,
                                         upstream = 2000,
                                         downstream = 2000)
## 3.2. Through co-accessible links
annotation_coacc <- annotate_with_coaccessibility(peaks = peaks,
                                                  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                                  links_df=links,
                                                  annot_dbi = org.Hs.eg.db,
                                                  protein_coding_only = TRUE,
                                                  verbose = TRUE,
                                                  add_tss_annotation = TRUE,
                                                  upstream = 2000,
                                                  downstream = 2000)

# 4. Run SPICEY ----------------------------------------------------------------
SPICEY_ANNOTATED_COACC <- SPICEY(atac=atac,
                                 rna=rna,
                                 annotation = annotation_coacc)
# add annotation of the region (promoter/distal)
SPICEY_ANNOTATED_COACC <- lapply(SPICEY_ANNOTATED_COACC, function(x) {
  if ("region_id" %in% colnames(x)) {
    x |>
      data.frame() |>
      left_join(annotation_near[, c("region_id", "annotation")], by = "region_id")
  } else {
    x  # leave unchanged if no region_id
  }
})

saveRDS(SPICEY_ANNOTATED_COACC, "/homes/users/gfuentes/scratch/projects/spicey_paper/HPAP_CTRL/data/HPAP_CTRL_SPICEY_COACC.rds")
