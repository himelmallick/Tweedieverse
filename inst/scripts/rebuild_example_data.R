# Rebuild all example containers from the package source root.
# Usage: Rscript inst/scripts/rebuild_example_data.R [output_dir] [franzosa_rds]
# [cache_dir]
# Preparation requires SeuratObject to read the PBMC source object. The saved
# examples use Bioconductor containers and require no Seurat dependency to load.
# Sources, source licences, checksums, and transformations are documented in
# inst/scripts/README.md and inst/extdata/DATA_LICENSES.md.

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else file.path("inst",
    "extdata")
source_rds <- if (length(args) >= 2L) args[[2L]] else tempfile(fileext = ".rds")
cache_dir <- if (length(args) >= 3L) args[[3L]] else tempfile(
    "tweedieverse_sources_")
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
hmp_abundance <- utils::read.delim(hmp_counts_file, row.names = 1, check.names =
    FALSE)
hmp_metadata <- utils::read.delim(hmp_metadata_file, row.names = 1,
    check.names = FALSE)
hmp_metadata <- hmp_metadata[order(hmp_metadata$collection, rownames(
    hmp_metadata)), ]
hmp_metadata <- hmp_metadata[!duplicated(hmp_metadata$subject), ]
hmp_metadata <- hmp_metadata[
    hmp_metadata$collection < 3 &
        stats::complete.cases(hmp_metadata[, c("age", "antibiotics",
            "diagnosis")]), ,
    drop = FALSE
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

# PBMC3k raw UMI counts and published cell labels from the Seurat tutorial
# object.
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
    ifelse(as.character(pbmc_labels) %in% c("Naive CD4 T", "Memory CD4 T"),
        "CD4", NA)
)
pbmc_cells <- unlist(lapply(c("CD8", "CD4"), function(group) {
    head(sort(names(pbmc_labels)[which(pbmc_groups == group)]), 30L)
}), use.names = FALSE)
pbmc_metadata <- methods::slot(pbmc, "meta.data")[pbmc_cells, , drop = FALSE]
pbmc_metadata$CellType <- factor(
    pbmc_groups[match(pbmc_cells, names(pbmc_labels))],
    levels = c("CD8", "CD4")
)
pbmc_metadata$library_size <- unname(Matrix::colSums(pbmc_counts[, pbmc_cells]))
pbmc_counts <- pbmc_counts[, pbmc_cells, drop = FALSE]
pbmc_counts <- pbmc_counts[variable_features(pbmc_counts), , drop = FALSE]
umi_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = pbmc_counts), colData = pbmc_metadata
)
rm(pbmc)

# Tabula Muris Smart-seq2 read counts, Figshare v8, CC BY 4.0.
brain_archive <- fetch_source(
    "FACS.zip", "https://ndownloader.figshare.com/files/10700143",
    "ef6e44ca2170155052591b5ab9d5d778"
)
brain_annotations <- fetch_source(
    "annotations_facs.csv",
    "https://ndownloader.figshare.com/files/13088129",
    "8ba7b81d91fe569cea07c384c2d34912"
)
brain_metadata <- utils::read.csv(brain_annotations, check.names = FALSE)
brain_groups <- c("astrocyte", "oligodendrocyte")
brain_metadata <- brain_metadata[
    brain_metadata$tissue == "Brain_Non-Myeloid" &
        brain_metadata$cell_ontology_class %in% brain_groups, , drop = FALSE
]
# Keep both cell types within one mouse; this is a cell-type illustration.
brain_group_counts <- table(
    brain_metadata$mouse.id, brain_metadata$cell_ontology_class
)
brain_mouse <- sort(rownames(brain_group_counts)[
    apply(brain_group_counts[, brain_groups], 1L, min) >= 30L
])[1L]
brain_metadata <- brain_metadata[
    brain_metadata$mouse.id == brain_mouse, , drop = FALSE
]
brain_cells <- unlist(lapply(brain_groups, function(group) {
    head(sort(brain_metadata$cell[
        brain_metadata$cell_ontology_class == group
    ]), 30L)
}), use.names = FALSE)
brain_metadata <- brain_metadata[
    match(brain_cells, brain_metadata$cell), , drop = FALSE
]
rownames(brain_metadata) <- brain_cells
brain_metadata$celltype <- factor(
    brain_metadata$cell_ontology_class, levels = brain_groups
)
brain_member <- "FACS/Brain_Non-Myeloid-counts.csv"
brain_directory <- file.path(cache_dir, "tabula_muris_v8")
utils::unzip(brain_archive, files = brain_member, exdir = brain_directory)
brain_csv <- file.path(brain_directory, brain_member)
brain_header <- names(utils::read.csv(
    brain_csv, nrows = 0L, check.names = FALSE
))
# Read only the selected cells, retaining every source gene for size factors.
brain_classes <- rep("NULL", length(brain_header))
brain_classes[1L] <- "character"
brain_classes[match(brain_cells, brain_header)] <- "numeric"
brain_table <- utils::read.csv(
    brain_csv, check.names = FALSE, colClasses = brain_classes
)
brain_counts <- as.matrix(brain_table[, -1L, drop = FALSE])
rownames(brain_counts) <- brain_table[[1L]]
brain_counts <- brain_counts[, brain_cells, drop = FALSE]
brain_full <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = brain_counts), colData = brain_metadata
)
brain_full <- scran::computeSumFactors(
    brain_full, clusters = brain_metadata$celltype, sizes = c(10, 20)
)
brain_metadata$source_size_factor <- SingleCellExperiment::sizeFactors(
    brain_full
)
brain_counts <- brain_counts[variable_features(brain_counts), , drop = FALSE]
nonumi_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = brain_counts), colData = brain_metadata
)
S4Vectors::metadata(nonumi_sce) <- list(
    dataset = "Tabula Muris FACS brain, Figshare v8",
    source_doi = "10.6084/m9.figshare.5829687.v8",
    licence = "CC BY 4.0",
    licence_url = "https://creativecommons.org/licenses/by/4.0/",
    protocol = "Smart-seq2 (non-UMI)",
    mouse = brain_mouse,
    source_md5 = c(
        FACS.zip = "ef6e44ca2170155052591b5ab9d5d778",
        annotations_facs.csv = "8ba7b81d91fe569cea07c384c2d34912"
    ),
    processing = paste(
        "First 30 sorted cells per type within the first eligible mouse;",
        "SCRAN factors use all source genes before selecting 20 by variance;",
        "source counts are unchanged."
    )
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
            match(rownames(abundance), mapping$Query_Name), ,
            drop = FALSE
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
        "C54:3 TAG", "3-Carboxypropyl trimethylammonium", "C54:4 TAG",
            "C54:5 TAG"
    )
)
subset_layers <- lapply(layer_names, function(layer_name) {
    layers[[layer_name]][features[[layer_name]], samples, drop = FALSE]
})
names(subset_layers) <- layer_names
franzosa_example_mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = subset_layers,
    colData = study_metadata[samples, c("diseaseCat_clean", "Age"), drop =
        FALSE]
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
