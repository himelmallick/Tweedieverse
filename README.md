Tweedieverse: Differential Analysis of Omics Data Based on the Tweedie Distribution
================
Himel Mallick, Ali Rahnavard
<img src="docs/logo.jpg" align="right" width="365px"/>

- [Introduction](#introduction)
- [Installation](#installation)
- [What Tweedieverse Fits](#what-tweedieverse-fits)
- [Domain-Specific Normalization and Offsets](#domain-specific-normalization-and-offsets)
- [Tweedie Variance Power](#tweedie-variance-power)
- [MaAsLin2 Linear Model Option](#maaslin2-linear-model-option)
- [Median Comparison](#median-comparison)
- [Presence-Absence Model and CCT Ranking](#presence-absence-model-and-cct-ranking)
- [Method-Specific Arguments](#method-specific-arguments)
- [Basic Usage](#basic-usage)
- [Input](#input)
- [Output](#output)
- [Getting Started with Tweedieverse](#getting-started-with-tweedieverse)
- [Citation](#citation)
- [Issues](#issues)

Introduction
------------

Tweedieverse is an R package for differential analysis of omics data using generalized linear models built around the [Tweedie distribution](https://en.wikipedia.org/wiki/Tweedie_distribution).

The package is designed for high-throughput data types where the outcome can be counts, non-negative continuous measurements, sparse measurements with many zeros, or approximately continuous abundance-like values. This includes microbiome taxonomic or functional profiles, bulk and single-cell RNA-seq measurements, metabolomics peak intensities, and other omics tables where many features are tested against the same metadata.

Tweedieverse fits one model per feature, extracts coefficient-level inference for each metadata variable, adjusts p-values across tested associations, and returns a ranked results table. It supports multiple covariates, optional random effects, configurable Tweedie variance powers, optional MaAsLin2 linear modeling for the normal case, optional median-comparison adjustment, and optional presence-absence logistic regression with Cauchy combination test ranking.

Installation
------------

Install the development version from GitHub:

```r
install.packages("devtools")
devtools::install_github("himelmallick/Tweedieverse")
library(Tweedieverse)
```

MaAsLin2 is an imported dependency because it is the default engine for the `tweedie_p = 0` linear-model path. Standard package installation will install required dependencies.

What Tweedieverse Fits
----------------------

The main function is `Tweedieverse()`. At a high level, it:

1. reads and aligns the feature and metadata tables,
2. filters features by abundance, prevalence, and variance,
3. builds a per-feature model formula from the selected metadata,
4. fits one model per feature,
5. returns coefficient estimates, standard errors, p-values, q-values, prevalence summaries, and model labels.

The default base model is `CPLM`, using compound Poisson Tweedie modeling for abundance-style outcomes. The default `tweedie_p = NULL` lets the model estimate the Tweedie index in the usual compound Poisson range between 1 and 2 when all filtered feature values are non-negative.

Domain-Specific Normalization and Offsets
-----------------------------------------

Tweedieverse normalization is **offset-only**. The input feature table is not normalized, divided, log-transformed, or otherwise replaced by the normalization step. Instead, `domain` and `normalization` are used to estimate a sample-level size factor, and that size factor enters the model as:

```r
offset(log(size_factor))
```

The response remains on the original input scale:

```r
raw_feature_value ~ metadata + offset(log(size_factor))
```

Domain defaults and supported strategies are:

| Domain | Default | Supported normalization strategies |
| --- | --- | --- |
| `microbiome` | `TSS` | `TSS`, `GMPR`, `CSS`, `MEDIAN`, `NONE` |
| `single_cell` | `SCRAN` | `SCRAN`, `MEDIAN`, `NONE` |
| `bulk_rnaseq` | `TMM` | `TMM`, `RLE` / `DESEQ2`, `CPM`, `MEDIAN`, `NONE` |
| `custom` | `TSS` | all supported strategies |

Examples:

```r
# Microbiome default: TSS-derived size factors as model offsets.
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  domain = "microbiome",
  fixed_effects = "diagnosis"
)

# Microbiome GMPR size factors. The feature table is still untouched.
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  domain = "microbiome",
  normalization = "GMPR",
  fixed_effects = "diagnosis"
)

# Bulk RNA-seq median-ratio size factors.
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  domain = "bulk_rnaseq",
  normalization = "RLE",
  fixed_effects = "condition"
)

# Use a user-supplied metadata column as the size factor.
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  scale_factor = "library_size",
  fixed_effects = "group"
)
```

`scale_factor` always overrides `domain` and `normalization`. Set `adjust_offset = FALSE` to fit without a model offset.

Tweedie Variance Power
----------------------

The `tweedie_p` argument controls the Tweedie variance power:

| `tweedie_p` | Interpretation | Tweedieverse behavior |
| --- | --- | --- |
| `NULL` | Estimated compound Poisson Tweedie index, usually `1 < p < 2` | Default CPLM path |
| `0` | Normal variance case | MaAsLin2 linear model by default, or Tweedie GLM when `Maaslin2_run = FALSE` |
| `1` | Poisson | Fixed-power Tweedie GLM for fixed effects; Poisson GLMM for random effects |
| `(1, 2)` | Compound Poisson with non-negative mass at zero | Estimated by default when `tweedie_p = NULL`; fixed values can be passed for fixed-effect GLMs |
| `2` | Gamma | Fixed-power Tweedie GLM for fixed effects; Gamma GLMM for random effects |
| `3` | Inverse Gaussian | Fixed-power Tweedie GLM for fixed effects; inverse Gaussian GLMM for random effects |
| `> 2` | Stable positive-support variance powers | Fixed-power Tweedie GLM for fixed effects |

The index is resolved after the feature table has been oriented, aligned with metadata, and filtered. `tweedie_p = NA` is treated the same as `tweedie_p = NULL`.

Values between 0 and 1 are undefined and are rejected. Negative supplied Tweedie powers are allowed for non-negative data with a warning, but they require strictly positive fitted means. Negative feature values are supported only at `p = 0`. If `tweedie_p` is left unspecified and the filtered data contain negative values, Tweedieverse resolves the index to `p = 0`, the Gaussian case.

When the resolved index is `p = 0`, Tweedieverse automatically uses `link = "identity"` and disables the log-offset because log, sqrt, and inverse links are not appropriate for responses that can span the negative half line.

For random-effects models, fixed powers are supported at `p = 0`, `p = 1`, `1 < p < 2`, `p = 2`, and `p = 3`. Fixed random-effect powers inside `(1, 2)` are pinned in `glmmTMB`; boundary values use the corresponding Gaussian, Poisson, Gamma, or inverse Gaussian family.

For `p = 0`, Tweedieverse can optionally transform a **copy** of the filtered feature table for the Gaussian abundance model:

```r
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = "demo_output/p0_log",
  fixed_effects = "diagnosis",
  tweedie_p = 0,
  p0_transform = "LOG",
  p0_transform_pseudocount = 1
)
```

Supported `p0_transform` values are:

| `p0_transform` | Behavior |
| --- | --- |
| `NONE` | Use the filtered feature table as-is |
| `CLR` | Centered log-ratio transform with `p0_transform_pseudocount` |
| `RCLR` | Robust CLR using positive values only; zeros remain zero |
| `LOG` | `log(x + p0_transform_pseudocount)` |
| `ARC_SIGNED_SQRT` | `sign(x) * sqrt(abs(x))` |

The original input object is not modified. When `output` is provided and a non-`NONE` transform is used, the analysis copy is written to `transformed_features.tsv`.

MaAsLin2 Linear Model Option
----------------------------

For `tweedie_p = 0`, MaAsLin2 is the default:

```r
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "diagnosis",
  tweedie_p = 0,
  Maaslin2_run = TRUE
)
```

This route uses MaAsLin2's linear-model workflow and returns MaAsLin2-labeled results with `tweedie.index = 0`.

MaAsLin2 can also be used together with Tweedieverse median comparison by explicitly enabling `median_comparison`:

```r
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "diagnosis",
  tweedie_p = 0,
  Maaslin2_run = TRUE,
  median_comparison = TRUE,
  median_subtraction = TRUE
)
```

To force the Tweedie GLM implementation for the same normal variance power, set `Maaslin2_run = FALSE`:

```r
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "diagnosis",
  tweedie_p = 0,
  Maaslin2_run = FALSE,
  link = "identity"
)
```

Median Comparison
-----------------

`median_comparison` is optional and defaults to `FALSE`.

When enabled, Tweedieverse compares fitted coefficients against a metadata-specific median effect rather than the default null of zero. This can be useful for relative-abundance or compositional settings where a global shift in coefficients may be expected.

```r
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "diagnosis",
  median_comparison = TRUE,
  median_subtraction = TRUE
)
```

When `median_comparison = TRUE`, the returned `coef`, `pval`, and `qval` columns are replaced by the median-comparison adjusted values.

Presence-Absence Model and CCT Ranking
--------------------------------------

Tweedieverse can also fit a presence-absence logistic regression model in addition to the abundance model:

```r
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "diagnosis",
  run_presence_absence_model = TRUE
)
```

The presence-absence model tests whether each feature is detected (`feature > 0`) as a function of the metadata. When `run_presence_absence_model = TRUE`, Tweedieverse runs this logistic model alongside the selected abundance model, including the MaAsLin2 abundance path used by default for `tweedie_p = 0`. Tweedieverse keeps the individual abundance-model and presence-absence-model p-values, then combines them with a Cauchy combination test (CCT) to produce a joint ranking.

The combined output includes columns such as:

- `pval_abundance` and `qval_abundance`
- `pval_presence` and `qval_presence`
- `base.model_abundance`
- `base.model_presence`
- combined `pval`, `qval`, and `base.model = "CCT"`

For example, this runs MaAsLin2 for the `p = 0` abundance model, then also runs the presence-absence logistic regression model and reports the CCT joint ranking:

```r
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "diagnosis",
  tweedie_p = 0,
  Maaslin2_run = TRUE,
  run_presence_absence_model = TRUE
)
```

Method-Specific Arguments
-------------------------

Use `method_args` to pass options to method-specific engines. The main current use is forwarding MaAsLin2 options when `tweedie_p = 0` and `Maaslin2_run = TRUE`.

```r
fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "diagnosis",
  tweedie_p = 0,
  method_args = list(
    Maaslin2 = list(
      normalization = "NONE",
      transform = "NONE"
    )
  )
)
```

The alias `method.args` is also supported for compatibility, but only one of `method_args` or `method.args` should be supplied in a single call.

Basic Usage
-----------

```r
library(Tweedieverse)

fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "group",
  cores = 1
)
```

A small reproducible example:

```r
set.seed(123)

features <- as.data.frame(matrix(rpois(60, lambda = 5), nrow = 20, ncol = 3))
colnames(features) <- paste0("feature", seq_len(ncol(features)))
rownames(features) <- paste0("sample", seq_len(nrow(features)))

metadata <- data.frame(
  group = rep(c("A", "B"), each = 10),
  row.names = rownames(features)
)

fit <- Tweedieverse(
  input_features = features,
  input_metadata = metadata,
  output = NULL,
  fixed_effects = "group",
  tweedie_p = 1,
  median_comparison = FALSE,
  max_significance = 1,
  cores = 1
)

head(fit)
```

Input
-----

Tweedieverse requires a feature table and a metadata table.

The feature table should contain samples and omics features such as taxa, genes, transcripts, pathways, metabolites, or peaks. The metadata table should contain the sample-level variables to test. Sample identifiers must overlap between the two inputs; Tweedieverse will align the tables before fitting models.

`input_features` can be:

- a data frame,
- a tab-delimited file path,
- a domain-appropriate Bioconductor container,
- a `MultiAssayExperiment` containing multiple omics layers.

The domain-specific Bioconductor container expectations are:

| Domain | Preferred container |
| --- | --- |
| `microbiome` | `TreeSummarizedExperiment` |
| `single_cell` | `SingleCellExperiment` |
| `bulk_rnaseq` | `SummarizedExperiment` |
| `custom` | any supported SummarizedExperiment-like container |

When a Bioconductor object is supplied, `assay_name` selects the assay to use, and `colData` is used as metadata unless `input_metadata` is provided separately.

### Multi-omics input with `MultiAssayExperiment`

When `input_features` is a `MultiAssayExperiment`, Tweedieverse iterates over each experiment/omics layer and returns a named list of omics-specific result tables. If `output` is provided, each omics layer is written to its own subdirectory.

Most Tweedieverse arguments can be supplied either as a single value used for every omics layer or as a named list/vector keyed by experiment name. Use a named list when one layer needs `NULL`, because atomic vectors cannot reliably store per-layer `NULL` values.

```r
fit <- Tweedieverse(
  input_features = mae,
  output = "tweedieverse_multiomics",
  fixed_effects = "group",
  domain = c(
    microbiome = "microbiome",
    rnaseq = "bulk_rnaseq"
  ),
  normalization = c(
    microbiome = "TSS",
    rnaseq = "TMM"
  ),
  tweedie_p = list(
    microbiome = NULL,
    rnaseq = 1.5
  ),
  run_presence_absence_model = c(
    microbiome = TRUE,
    rnaseq = FALSE
  ),
  cores = 1
)

names(fit)
fit$microbiome
fit$rnaseq
```

Output
------

For single-omics input, `Tweedieverse()` returns a data frame ordered by increasing q-value. For `MultiAssayExperiment` input, it returns a named list of omics-specific data frames, each with an `omics` column. The main columns are:

- `feature`: the tested feature,
- `metadata`: the metadata variable,
- `value`: the coefficient level or contrast,
- `coef`: coefficient estimate,
- `stderr`: standard error,
- `pval`: p-value,
- `qval`: multiple-testing adjusted p-value,
- `base.model`: model used for the reported ranking,
- `tweedie.index`: estimated or fixed Tweedie variance power,
- `N`: number of samples,
- `N.not.zero`: number of nonzero samples,
- `percent.zero`: percent zero values.

When `run_presence_absence_model = TRUE`, additional abundance-specific and presence-specific columns are returned so the joint CCT ranking can be inspected alongside the individual model results.

Getting Started with Tweedieverse
---------------------------------

Check out the [Tweedie Labs](https://github.com/himelmallick/TweedieLabs/) repository for walkthrough tutorials on applying Tweedieverse to different omics data types.

For full function options, see:

```r
?Tweedieverse
```

Citation
--------

To cite **Tweedieverse** in publications, please use:

Mallick, H, Chatterjee, S, Chowdhury, S, Chatterjee, S, Rahnavard, A, Hicks, SC. [Differential expression of single-cell RNA-seq data using Tweedie models](https://onlinelibrary.wiley.com/doi/10.1002/sim.9430). Statistics in Medicine. 2022; 41(18): 3492-3510. doi:10.1002/sim.9430

To cite the **Tweedieverse** software, please use:

Mallick H et al. (2021). [Tweedieverse - A Unified Statistical Framework for Differential Analysis of Multi-omics Data](https://github.com/himelmallick/Tweedieverse). R package, <https://github.com/himelmallick/Tweedieverse>.

Issues
------

We are happy to troubleshoot issues with the package. Please contact the maintainer by email or [open an issue](https://github.com/himelmallick/Tweedieverse/issues) in the GitHub repository.
