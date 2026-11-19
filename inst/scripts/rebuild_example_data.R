# Rebuild all example containers from the package source root.
# Usage: Rscript inst/scripts/rebuild_example_data.R [output_dir] [franzosa_rds] [cache_dir]
# Preparation requires SeuratObject to read the PBMC source object. The saved
# examples use Bioconductor containers and require no Seurat dependency to load.
# Sources, source licences, checksums, and transformations are documented in
# inst/scripts/README.md and inst/extdata/DATA_LICENSES.md.

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else file.path("inst", "extdata")
source_rds <- if (length(args) >= 2L) args[[2L]] else tempfile(fileext = ".rds")
cache_dir <- if (length(args) >= 3L) args[[3L]] else tempfile("tweedieverse_sources_")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

fetch_source <- function(filename, url, md5) {
  path <- file.path(cache_dir, filename)
  if (!file.exists(path)) {
    utils::download.file(url, path, mode = "wb")
  }
  if (unname(tools::md5sum(path)) != md5) {
    stop("Source checksum mismatch for ", filename)
  }
  path
}

variable_features <- function(counts, n = 20L) {
  n_samples <- ncol(counts)
  totals <- Matrix::rowSums(counts)
  squared_totals <- Matrix::rowSums(counts^2)
  variances <- (squared_totals - totals^2 / n_samples) / (n_samples - 1)
  order(-variances, rownames(counts))[seq_len(min(n, nrow(counts)))]
}

# HMP2 baseline microbiomes, distributed in MaAsLin2 under MIT.
maaslin_revision <- "39eda3be7d801f0eb594eb7b09b26fad7d5c5875"
maaslin_url <- paste0(
  "https://raw.githubusercontent.com/biobakery/Maaslin2/",
  maaslin_revision, "/inst/extdata/"
)
hmp_counts_file <- fetch_source(
  "HMP2_taxonomy.tsv", paste0(maaslin_url, "HMP2_taxonomy.tsv"),
  "78231446ddedd5c080cea53bbf37ac38"
)
hmp_metadata_file <- fetch_source(
  "HMP2_metadata.tsv", paste0(maaslin_url, "HMP2_metadata.tsv"),
  "30dce926858e6feda6fbde8a149fd2b4"
)
hmp_abundance <- utils::read.delim(hmp_counts_file, row.names = 1, check.names = FALSE)
hmp_metadata <- utils::read.delim(hmp_metadata_file, row.names = 1, check.names = FALSE)
hmp_metadata <- hmp_metadata[order(hmp_metadata$collection, rownames(hmp_metadata)), ]
hmp_metadata <- hmp_metadata[!duplicated(hmp_metadata$subject), ]
hmp_metadata <- hmp_metadata[
  hmp_metadata$collection < 3 &
    stats::complete.cases(hmp_metadata[, c("age", "antibiotics", "diagnosis")]),
  , drop = FALSE
]
hmp_groups <- c("nonIBD", "CD", "UC")
hmp_samples <- unlist(lapply(hmp_groups, function(group) {
  head(rownames(hmp_metadata)[hmp_metadata$diagnosis == group], 16L)
}), use.names = FALSE)
hmp_metadata <- hmp_metadata[hmp_samples, , drop = FALSE]
hmp_metadata$diagnosis <- factor(hmp_metadata$diagnosis, levels = hmp_groups)
hmp_counts <- t(as.matrix(hmp_abundance[hmp_samples, , drop = FALSE]))
hmp_counts <- hmp_counts[variable_features(hmp_counts), , drop = FALSE]
microbiome_tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = list(counts = hmp_counts), colData = hmp_metadata
)

# PBMC3k raw UMI counts and published cell labels from the Seurat tutorial object.
# Original 10x Genomics data are distributed under CC BY 4.0.
pbmc_file <- fetch_source(
  "pbmc3k_final.rds",
  "https://www.dropbox.com/s/63gnlw45jf7cje8/pbmc3k_final.rds?dl=1",
  "ca57284952b2a17553cee5b87a71e7a7"
)
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("SeuratObject is needed only to prepare the PBMC source data.")
}
pbmc <- readRDS(pbmc_file)
pbmc_counts <- methods::slot(methods::slot(pbmc, "assays")[["RNA"]], "counts")
pbmc_labels <- methods::slot(pbmc, "active.ident")
pbmc_groups <- ifelse(
  as.character(pbmc_labels) == "CD8 T", "CD8",
  ifelse(as.character(pbmc_labels) %in% c("Naive CD4 T", "Memory CD4 T"), "CD4", NA)
)
pbmc_cells <- unlist(lapply(c("CD8", "CD4"), function(group) {
  head(sort(names(pbmc_labels)[which(pbmc_groups == group)]), 30L)
}), use.names = FALSE)
pbmc_metadata <- methods::slot(pbmc, "meta.data")[pbmc_cells, , drop = FALSE]
pbmc_metadata$CellType <- factor(
  pbmc_groups[match(pbmc_cells, names(pbmc_labels))], levels = c("CD8", "CD4")
)
pbmc_metadata$library_size <- unname(Matrix::colSums(pbmc_counts[, pbmc_cells]))
pbmc_counts <- pbmc_counts[, pbmc_cells, drop = FALSE]
pbmc_counts <- pbmc_counts[variable_features(pbmc_counts), , drop = FALSE]
umi_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = pbmc_counts), colData = pbmc_metadata
)
rm(pbmc)

# Brain read counts from SC2P (GPL), originating from GEO GSE67835.
brain_file <- fetch_source(
  "brain_scRNAseq.rda",
  paste0(
    "https://raw.githubusercontent.com/haowulab/SC2P/",
    "bdbb526cb81d598da99b63f55cc64c826863b954/data/brain_scRNAseq.rda"
  ),
  "15562abf4745bc7aa8147f878b210fe9"
)
brain_source <- new.env()
load(brain_file, envir = brain_source)
brain_counts <- brain_source$Y
brain_metadata <- droplevels(brain_source$design)
colnames(brain_counts) <- rownames(brain_metadata)
brain_cells <- unlist(lapply(c("astrocytes", "oligodendrocytes"), function(group) {
  head(sort(rownames(brain_metadata)[brain_metadata$celltype == group]), 30L)
}), use.names = FALSE)
brain_counts <- brain_counts[, brain_cells, drop = FALSE]
brain_metadata <- droplevels(brain_metadata[brain_cells, , drop = FALSE])
brain_full <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = brain_counts), colData = brain_metadata
)
brain_full <- scran::computeSumFactors(
  brain_full, clusters = brain_metadata$celltype, sizes = c(10, 20)
)
brain_metadata$source_size_factor <- SingleCellExperiment::sizeFactors(brain_full)
brain_counts <- brain_counts[variable_features(brain_counts), , drop = FALSE]
nonumi_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = brain_counts), colData = brain_metadata
)

# HuMMANet v1 archive, DOI 10.5281/zenodo.22053943, CC BY 4.0.
# Verma S, Arora N, Mallick H, Ghosh TS (2026).
source_url <- paste0(
  "https://zenodo.org/records/22053943/files/",
  "franzosae_2019.rds?download=1"
)
if (!file.exists(source_rds)) {
  utils::download.file(source_url, source_rds, mode = "wb")
}
expected_md5 <- "9ab747a9c36003e9b76f4d9d58d67e77"
if (unname(tools::md5sum(source_rds)) != expected_md5) {
  stop("The source RDS does not match the pinned HuMMANet v1 archive file.")
}
source_data <- readRDS(source_rds)
study_metadata <- source_data$MetadataProfile
rownames(study_metadata) <- study_metadata$sample_id

# Convert samples-by-features tables to features-by-samples assays.
# Retain the archive's abundance values and feature names without rescaling.
layer_names <- c("taxaAbundanceProfile", "harmonizedMetaboliteProfile")
layers <- lapply(layer_names, function(layer_name) {
  table <- source_data[[layer_name]]
  rownames(table) <- table$sample_id
  table$sample_id <- NULL
  abundance <- t(as.matrix(table))
  abundance <- abundance[, rownames(study_metadata), drop = FALSE]

  if (layer_name == "taxaAbundanceProfile") {
    annotations <- data.frame(
      taxon_id = rownames(abundance), row.names = rownames(abundance)
    )
  } else {
    mapping <- source_data$metaboliteIdentifierMapping
    annotations <- mapping[
      match(rownames(abundance), mapping$Query_Name), , drop = FALSE
    ]
    annotations$feature_id <- rownames(abundance)
    rownames(annotations) <- rownames(abundance)
  }
  SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = abundance),
    rowData = annotations,
    colData = study_metadata
  )
})
names(layers) <- layer_names
FranzosaE_2019_data <- list(
  FranzosaE_2019 = MultiAssayExperiment::MultiAssayExperiment(
    experiments = layers, colData = study_metadata
  )
)
# Preserve the historical filename spelling for existing source references.
save(
  FranzosaE_2019_data,
  file = file.path(output_dir, "FransozaE_2019_MAE.RData"),
  compress = "xz", version = 2
)

# Fixed, compact vignette subset: eight gastrointestinal and eight healthy
# samples, four microbial species, and four metabolites. Selection is explicit
# and independent of model results. Layer annotations and metadata are retained.
samples <- c(
  "PRISM.7122", "PRISM.7147", "PRISM.7150", "PRISM.7153",
  "PRISM.7184", "PRISM.7238", "PRISM.7406", "PRISM.7408",
  "PRISM.7662", "PRISM.7847", "PRISM.7855", "PRISM.7858",
  "PRISM.7860", "PRISM.7861", "PRISM.7862", "PRISM.7870"
)
features <- list(
  taxaAbundanceProfile = c(
    "Bacteroides_vulgatus", "Ruminococcus_gnavus",
    "Eubacterium_rectale", "Faecalibacterium_prausnitzii"
  ),
  harmonizedMetaboliteProfile = c(
    "C54:3 TAG", "3-Carboxypropyl trimethylammonium", "C54:4 TAG", "C54:5 TAG"
  )
)
subset_layers <- lapply(layer_names, function(layer_name) {
  layers[[layer_name]][features[[layer_name]], samples, drop = FALSE]
})
names(subset_layers) <- layer_names
franzosa_example_mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = subset_layers,
  colData = study_metadata[samples, c("diseaseCat_clean", "Age"), drop = FALSE]
)
save(
  microbiome_tse, umi_sce, nonumi_sce, franzosa_example_mae,
  file = file.path(output_dir, "tweedieverse_data.RData"),
  compress = "xz", version = 2
)

example_objects <- list(
  microbiome_tse = microbiome_tse,
  umi_sce = umi_sce,
  nonumi_sce = nonumi_sce,
  franzosa_example_mae = franzosa_example_mae
)
for (object_name in names(example_objects)) {
  saveRDS(
    example_objects[[object_name]],
    file = file.path(output_dir, paste0(object_name, ".rds")),
    compress = "xz", version = 2
  )
}
