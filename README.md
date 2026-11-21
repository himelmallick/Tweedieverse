# Tweedieverse <img src="man/figures/logo.png" align="right" height="139" alt="Tweedieverse logo" />

**Differential analysis of omics data using Tweedie models.**

Tweedieverse fits one model per feature, supports multiple covariates and
optional random effects, and returns association tables with effect estimates
and adjusted p-values. It accepts microbiome, single-cell, bulk RNA-seq, and
other omics measurements through Bioconductor containers.

[Input](#input) · [Quick start](#quick-start) · [Output](#output) ·
[Model options](#model-options) · [Multi-omics](#multi-omics) ·
[Documentation](#documentation-and-example-data)

## Installation

Once Tweedieverse is available in your Bioconductor version, install it with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("Tweedieverse")
```

See the [Bioconductor installation guide](https://bioconductor.org/install/)
for R/Bioconductor version compatibility and setup.

Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("himelmallick/Tweedieverse")
```

## Input

Supply a container whose assay has **features in rows and samples in columns**.
Sample covariates belong in `colData`; sample identifiers are used to align the
assay and metadata. Select the assay with `assay_name` (default `"counts"`).

| Data / `domain` | Preferred container | Default normalization |
| --- | --- | --- |
| Microbiome / `"microbiome"` | `TreeSummarizedExperiment` | `TSS` |
| Single-cell / `"single_cell"` | `SingleCellExperiment` | `SCRAN` |
| Bulk RNA-seq / `"bulk_rnaseq"` | `SummarizedExperiment` | `TMM` |
| Other omics / `"custom"` | A supported SummarizedExperiment-like container | `MEDIAN` |

For multiple assays, use a `MultiAssayExperiment` (see [below](#multi-omics)).
To use your own count matrix and metadata, first wrap them in a supported
container. Plain feature tables, matrices, lists, and file paths are not accepted
as `input_features`. An optional `input_metadata` data frame can override
container metadata.

## Quick start

This example uses packaged PBMC3k UMI counts: 20 genes across 60 CD4/CD8 cells.
The library-size column was calculated from all source genes before subsetting.

```r
library(Tweedieverse)

umi_sce <- readRDS(system.file(
    "extdata", "umi_sce.rds", package = "Tweedieverse"
))

fit <- Tweedieverse(
    input_features = umi_sce,
    output = NULL,
    assay_name = "counts",
    domain = "single_cell",
    fixed_effects = "CellType",
    reference = "CellType,CD8",
    scale_factor = "library_size",
    tweedie_p = 1,
    return_all = TRUE
)

head(fit)
discoveries <- subset(fit, !is.na(qval) & qval <= 0.05)
```

The example fixes `p = 1` to illustrate the Poisson variance case. Leaving
`tweedie_p = NULL` estimates a compound Poisson-gamma power for each feature.
This is a cell-type illustration from one donor; population-level inference
requires biological replication and an appropriate study design.

## Output

For a single assay, `Tweedieverse()` returns a `data.frame`, ordered by q-value.
By default, it returns associations with `qval <= max_significance` (0.05).
Use `return_all = TRUE` to inspect the full available table, including with
`output = NULL`. Missing q-values are excluded from significance selection.

| Columns | Meaning |
| --- | --- |
| `feature`, `metadata`, `value` | Feature, covariate, and tested level or predictor |
| `coef`, `stderr` | Coefficient and standard error on the model's link scale |
| `pval`, `qval` | Test p-value and adjusted p-value; correction defaults to BH |
| `base.model`, `tweedie.index` | Model label and fitted or fixed variance power |
| `N`, `N.not.zero`, `percent.zero` | Aligned sample count, positive measurements, and rounded percentage of zeros |

In the quick start, `CellType` coefficients compare CD4 with CD8. With the
log link, a positive coefficient indicates a higher fitted mean in CD4 cells,
and `exp(coef)` is the mean ratio at equal library size. Continuous covariates
are standardized by default; use `standardize = FALSE` for their original units.

Set `output` to a directory to write `all_results.tsv`,
`significant_results.tsv`, logs, and model-index plots. Add
`plot_heatmap = TRUE` or `plot_scatter = TRUE` for optional association plots.
`output = NULL` creates no output directory or files.

## Model options

### Variance power

| `tweedie_p` | Model choice |
| --- | --- |
| `NULL` | Estimate a feature-specific compound Poisson-gamma power, `1 < p < 2`, for non-negative data |
| `0` | Gaussian case; MaAsLin2 linear modeling by default |
| `1` | Poisson variance case |
| Between `1` and `2` | Fixed compound Poisson-gamma power |
| `2` | Gamma variance case |
| `3` | Inverse Gaussian variance case |

The default link is log. At `p = 0`, the package uses an identity link and
disables log offsets. Negative measurements require `p = 0`; an unspecified
power resolves to zero when negative values remain after filtering. Powers
between zero and one are undefined. See the function help for other supported
powers and random-effect restrictions.

At `p = 0`, use `Maaslin2_run = FALSE` to select the Tweedie GLM implementation.
Optional `p0_transform` choices are `NONE`, `CLR`, `RCLR`, `LOG`, and
`ARC_SIGNED_SQRT`. These operate on an analysis copy of the feature matrix.

### Normalization and offsets

Normalization estimates sample size factors and enters them as
`offset(log(size_factor))`; it does not divide or replace the response matrix.
Choose a strategy with `normalization`, supply a `colData` column with
`scale_factor`, or disable offsets with `adjust_offset = FALSE`.
A supplied `scale_factor` overrides the normalization strategy.

For already scaled relative abundances, use `normalization = "NONE"`.
Filtering uses `abd_threshold`, `prev_threshold`, and `var_threshold`;
`entropy_threshold` controls metadata filtering.

Sparse assays are preserved through filtering and Tweedie fitting. MaAsLin2,
the edgeR TMM backend, selected transformations, and association plots can
require dense representations; see the tutorial for details.

### Optional analyses

- **Presence/absence:** set `run_presence_absence_model = TRUE` to combine
  abundance and prevalence evidence using a Cauchy combination test. By default,
  prevalence fitting requires at least three detected and three undetected
  samples; adjust this with `min_presence_absence_count`. Below either minimum,
  `pval_presence` is missing and the abundance p-value is used directly. Exact
  component p-values of 0 and 1 combine to 0. Component results are retained
  alongside the joint p-value and q-value; `coef` remains the abundance effect.
- **Median comparison:** for relative-abundance analyses, enable
  `median_comparison = TRUE` to test effects against a covariate-specific median.
  Set `median_comparison_seed` for reproducible Monte Carlo calculations. A
  supplied seed restores the previous RNG state afterward.
- **Study design and execution:** use `fixed_effects`, `random_effects`, and
  `reference` to specify the design, `method_args` for engine-specific options,
  and `BPPARAM` for BiocParallel execution. The default backend is serial.

## Multi-omics

A `MultiAssayExperiment` returns a named list of association tables, each with
an `omics` column. Layers are fitted independently and q-values are adjusted
within each layer. Sample alignment uses `sampleMap`.

```r
mae <- readRDS(system.file(
    "extdata", "franzosa_example_mae.rds", package = "Tweedieverse"
))

multi_fit <- Tweedieverse(
    input_features = mae,
    output = NULL,
    assay_name = "abundance",
    fixed_effects = c("diseaseCat_clean", "Age"),
    domain = c(
        taxaAbundanceProfile = "microbiome",
        harmonizedMetaboliteProfile = "custom"
    ),
    normalization = c(
        taxaAbundanceProfile = "NONE",
        harmonizedMetaboliteProfile = "MEDIAN"
    ),
    return_all = TRUE
)

names(multi_fit)
head(multi_fit$taxaAbundanceProfile)
```

Most arguments accept a shared value or named values keyed by layer. Use a named
list when a layer needs `NULL`, such as a mix of estimated and fixed powers.

## Documentation and example data

The [tutorial](https://himelmallick.github.io/Tweedieverse/articles/Tweedieverse.html) covers microbiome abundances, UMI and
non-UMI single-cell data, pseudobulk aggregation with muscat, bulk RNA-seq with
airway, multi-omics analysis, result interpretation, and visualization.

```r
vignette("Tweedieverse", package = "Tweedieverse")
?Tweedieverse
```

Packaged examples run without downloads. See the
[data provenance and licences](https://github.com/himelmallick/Tweedieverse/blob/master/inst/extdata/DATA_LICENSES.md) and
[preparation instructions](https://github.com/himelmallick/Tweedieverse/blob/master/inst/scripts/README.md) for their sources and
reconstruction steps.

## Citation and support

Mallick H et al. (2022). [Differential expression of single-cell RNA-seq data
using Tweedie models](https://doi.org/10.1002/sim.9430).
*Statistics in Medicine*, 41(18), 3492–3510.

Use `citation("Tweedieverse")` for the package citation. Please report problems
through [GitHub issues](https://github.com/himelmallick/Tweedieverse/issues).
