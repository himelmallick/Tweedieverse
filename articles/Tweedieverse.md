# Tweedieverse

## Overview

`Tweedieverse` provides a unified differential analysis workflow for
omics measurements whose empirical variance structure is well
represented by the Tweedie family. The framework is intended for
heterogeneous high-throughput data, including microbiome abundances,
single-cell RNA-seq measurements, bulk RNA-seq summaries, metabolite
intensities, and other continuous or count-like assays.

The package uses Bioconductor containers as the primary interface. A
single omics assay is represented as a `SummarizedExperiment`,
`TreeSummarizedExperiment`, or `SingleCellExperiment`, and matched
multi-omics assays are represented as a `MultiAssayExperiment`. This
keeps feature assays, sample metadata, and experiment-level annotations
together and makes the workflow compatible with the broader Bioconductor
ecosystem.

This vignette demonstrates six Bioconductor-oriented workflows:

- microbiome relative abundances stored in a `TreeSummarizedExperiment`,
- single-cell RNA-seq UMI counts stored in a `SingleCellExperiment`,
- single-cell RNA-seq non-UMI/read counts stored in a
  `SingleCellExperiment`,
- pseudobulk RNA-seq counts generated with
  [`muscat::aggregateData()`](https://rdrr.io/pkg/muscat/man/aggregateData.html),
- bulk RNA-seq counts from the Bioconductor `airway` dataset stored in a
  `SummarizedExperiment`,
- matched microbiome and metabolomics assays stored in a
  `MultiAssayExperiment`.

## Motivation

In omics data, biological processes or technical factors may induce
distinct statistical properties (such as overdispersion and
zero-inflation) across experimental platforms (or library preparation
protocols). Most published differential analysis methods are
omics-specific and assume a single probabilistic model for all features.
`DESeq2` and `edgeR` use the negative binomial distribution for
sequencing counts. `MAST` uses a two-part hurdle model for single-cell
expression. `MaAsLin2` supports several microbiome association models,
including linear and compound Poisson models, while `MaAsLin3` models
prevalence and abundance separately. `limma` fits linear models to
expression data, including transformed RNA-seq measurements.

These choices ignore the fact that data distributions are
feature-specific and vary considerably both within and between platforms
or technologies. A typical multi-omics study therefore runs a different
tool for each assay, each with its own assumptions, preprocessing, and
output format, and the resulting effect sizes and significance calls are
difficult to compare.

`Tweedieverse` provides a single home for differential analysis across
omics data types. It is built on Tweedie generalized linear models,
whose flexibility enables feature-specific modeling of omics data by
adjusting a single configuration parameter per feature: the Tweedie
index `p`, which spans the Gaussian, Poisson, compound Poisson-gamma,
gamma, and inverse Gaussian distributions. `Tweedieverse` estimates `p`
separately for each feature, so every feature is modeled with a
mean-variance relationship learned from its own data. This yields a
unified framework that adapts to biological and technological
variability in multi-omics profiles. Model adequacy and study design
determine the validity and power of the resulting tests; these examples
demonstrate the workflow rather than benchmark comparative performance.

With `Tweedieverse`, users can:

- run differential analysis on microbiome, single-cell, pseudobulk, bulk
  RNA-seq, and metabolomics data with one function,
  [`Tweedieverse()`](https://himelmallick.github.io/Tweedieverse/reference/Tweedieverse.md);
- pass data directly as a `SummarizedExperiment`,
  `SingleCellExperiment`, `TreeSummarizedExperiment`, or
  `MultiAssayExperiment`, with sample metadata and feature annotations
  kept intact;
- compare results across assays on a common scale, since every assay is
  analyzed with the same model family, testing procedure, and output
  format.

The sections below walk through each of these workflows.

## Relationship to existing tools

Bioconductor provides established differential-analysis workflows:
[DESeq2](https://bioconductor.org/packages/DESeq2) and
[edgeR](https://bioconductor.org/packages/edgeR) use negative-binomial
models for sequencing counts;
[limma](https://bioconductor.org/packages/limma) uses linear models with
empirical Bayes inference, including precision-weighted RNA-seq analysis
with voom. [MAST](https://bioconductor.org/packages/MAST) uses a
two-part hurdle model for single-cell detection and positive expression.
[MaAsLin2](https://bioconductor.org/packages/Maaslin2) offers
multivariable meta-omic association analysis with several model choices,
including compound Poisson models, while
[MaAsLin3](https://bioconductor.org/packages/release/bioc/html/maaslin3.html)
combines prevalence and abundance analyses.

Tweedie modeling is also available through general-purpose CRAN tools.
[cplm](https://CRAN.R-project.org/package=cplm) fits compound Poisson
regression models, and
[glmmTMB](https://CRAN.R-project.org/package=glmmTMB) supports Tweedie
generalized linear mixed models.
[statmod](https://CRAN.R-project.org/package=statmod) supplies Tweedie
GLM families, while
[tweedie](https://CRAN.R-project.org/package=tweedie) provides
distribution and likelihood calculations. Tweedieverse uses cplm,
glmmTMB, and statmod as modeling engines.

Tweedieverse adds a common omics workflow around these engines: direct
Bioconductor-container input, sample alignment, feature filtering,
domain-specific size-factor offsets, feature-wise fitting, optional
prevalence tests, multiple-testing correction, and consistently
structured association tables. It supports both fixed and estimated
variance powers and analyzes `MultiAssayExperiment` layers through the
same interface. Each layer is fitted independently; effect sizes retain
their assay- and link-specific interpretation.

## Installation

After Bioconductor acceptance, the release version can be installed
with:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("Tweedieverse")
```

The development version can be installed from GitHub:

``` r

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("himelmallick/Tweedieverse")
```

``` r

library(Tweedieverse)
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
library(SingleCellExperiment)
library(TreeSummarizedExperiment)
#> Loading required package: Biostrings
#> Loading required package: XVector
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
library(MultiAssayExperiment)
library(muscat)
library(airway)
```

## Loading Packaged Example Data

The example containers used in this vignette are distributed with the
package under `inst/extdata`. This allows the vignette to run during
Bioconductor checks without network access.

The microbiome and UMI examples are compact subsets of the real HMP2 and
PBMC3k datasets used in the original [TweedieLabs
tutorials](https://github.com/himelmallick/TweedieLabs). The non-UMI
example uses CC BY 4.0 Tabula Muris Smart-seq2 brain data. The Franzosa
multi-omics example comes from the [HuMMANet v1
archive](https://doi.org/10.5281/zenodo.22053943). The preparation code,
source details, and per-object licences are included in
scripts/rebuild_example_data.R, scripts/README.md, and
extdata/DATA_LICENSES.md inside the installed package.

``` r

microbiome_tse <- readRDS(system.file(
    "extdata", "microbiome_tse.rds",
    package = "Tweedieverse"
))
umi_sce <- readRDS(system.file(
    "extdata", "umi_sce.rds",
    package = "Tweedieverse"
))
nonumi_sce <- readRDS(system.file(
    "extdata", "nonumi_sce.rds",
    package = "Tweedieverse"
))
franzosa_example_mae <- readRDS(system.file(
    "extdata", "franzosa_example_mae.rds",
    package = "Tweedieverse"
))
```

## Microbiome Relative Abundance Example

For microbiome data,
[`Tweedieverse()`](https://himelmallick.github.io/Tweedieverse/reference/Tweedieverse.md)
accepts a `TreeSummarizedExperiment`. The assay contains taxa/features
in rows and samples in columns; `colData` contains sample metadata.

This subset contains the 20 most variable taxa across 48 baseline HMP2
samples, with 16 samples each from the nonIBD, CD, and UC groups. The
source tables are distributed with MaAsLin2. Baseline selection and
metadata preparation follow the [original microbiome
tutorial](https://github.com/himelmallick/TweedieLabs) .

The example table contains relative abundances, so this analysis
disables additional offset normalization by setting
`normalization = "NONE"`. For raw microbiome counts, omit that argument
to use the microbiome default `TSS` size-factor offset.

``` r

microbiome_tse
#> class: TreeSummarizedExperiment 
#> dim: 20 48 
#> metadata(0):
#> assays(1): counts
#> rownames(20): Prevotella copri Bacteroides vulgatus ... Veillonella
#>   unclassified Ruminococcus torques
#> rowData names(0):
#> colnames(48): CSM67UH7 CSM6J2H9_P ... MSM5LLIQ_P MSM79HBZ
#> colData names(11): site age ... dysbiosisCD collection
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowLinks: NULL
#> rowTree: NULL
#> colLinks: NULL
#> colTree: NULL

microbiome_fit <- Tweedieverse(
    input_features = microbiome_tse,
    output = NULL,
    assay_name = "counts",
    domain = "microbiome",
    normalization = "NONE",
    fixed_effects = c("diagnosis", "antibiotics", "age"),
    reference = "diagnosis,nonIBD",
    return_all = TRUE,
    tweedie_p = NULL,
    median_comparison = FALSE
)
#> The specified assay has been extracted
#> 2026-11-21 21:00:47.846237 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:47.853296 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:47.860777 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:47.861206 INFO::Total samples in data: 48
#> 2026-11-21 21:00:47.86158 INFO::Min samples required with min abundance for a feature not to be filtered: 4.800000
#> 2026-11-21 21:00:47.862215 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:47.862681 INFO::Filtered feature names:
#> 2026-11-21 21:00:47.863389 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:47.86382 INFO::Filtered feature names:
#> 2026-11-21 21:00:47.86454 INFO::Offset size factors computed with normalization: NONE
#> 2026-11-21 21:00:47.865748 INFO::Excluded metadata with entropy less or equal to 0: collection
#> 2026-11-21 21:00:47.866473 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-21 21:00:47.878552 INFO::Formula for fixed effects: expr ~  diagnosis + antibiotics + age
#> 2026-11-21 21:00:47.879287 INFO::Running selected analysis method: CPLM
#> 2026-11-21 21:00:47.900184 INFO::Fitting model to feature number 1, Prevotella copri
#> 2026-11-21 21:00:49.02098 INFO::Fitting model to feature number 2, Bacteroides vulgatus
#> 2026-11-21 21:00:49.202194 INFO::Fitting model to feature number 3, Bacteroides dorei
#> 2026-11-21 21:00:49.364125 INFO::Fitting model to feature number 4, Faecalibacterium prausnitzii
#> 2026-11-21 21:00:49.5204 INFO::Fitting model to feature number 5, Akkermansia muciniphila
#> 2026-11-21 21:00:49.664987 INFO::Fitting model to feature number 6, Bacteroides uniformis
#> 2026-11-21 21:00:49.807489 INFO::Fitting model to feature number 7, Bacteroides fragilis
#> 2026-11-21 21:00:49.936823 INFO::Fitting model to feature number 8, Bacteroides stercoris
#> 2026-11-21 21:00:50.084427 INFO::Fitting model to feature number 9, Bacteroides cellulosilyticus
#> 2026-11-21 21:00:50.218735 INFO::Fitting model to feature number 10, Eubacterium rectale
#> 2026-11-21 21:00:50.364851 INFO::Fitting model to feature number 11, Escherichia coli
#> Warning: Model fit failed for feature 'Escherichia coli': inner loop 1; cannot
#> correct step size
#> 2026-11-21 21:00:50.378738 INFO::Fitting model to feature number 12, Bacteroides ovatus
#> 2026-11-21 21:00:50.502841 INFO::Fitting model to feature number 13, Clostridium bolteae
#> 2026-11-21 21:00:50.668323 INFO::Fitting model to feature number 14, Bacteroides thetaiotaomicron
#> 2026-11-21 21:00:50.808132 INFO::Fitting model to feature number 15, Alistipes putredinis
#> 2026-11-21 21:00:50.93421 INFO::Fitting model to feature number 16, Subdoligranulum unclassified
#> 2026-11-21 21:00:51.12592 INFO::Fitting model to feature number 17, Ruminococcus gnavus
#> 2026-11-21 21:00:51.269919 INFO::Fitting model to feature number 18, Bacteroides massiliensis
#> 2026-11-21 21:00:51.406082 INFO::Fitting model to feature number 19, Veillonella unclassified
#> Warning: Model fit failed for feature 'Veillonella unclassified': NA/NaN/Inf in
#> 'x'
#> 2026-11-21 21:00:51.421452 INFO::Fitting model to feature number 20, Ruminococcus torques
#> 2026-11-21 21:00:51.569103 INFO::Counting prevalence for each feature

head(microbiome_fit)
#>                        feature  metadata value       coef    stderr
#> 1          Clostridium bolteae diagnosis    CD  3.1337650 0.8303239
#> 2     Bacteroides massiliensis diagnosis    CD -3.8494500 1.1246612
#> 3        Bacteroides uniformis       age   age  0.5623203 0.1993187
#> 4 Bacteroides cellulosilyticus       age   age  1.3947464 0.4531786
#> 5 Subdoligranulum unclassified       age   age  0.3964299 0.1323723
#> 6          Ruminococcus gnavus       age   age -1.2215268 0.4317158
#>           pval       qval base.model tweedie.index  N N.not.zero percent.zero
#> 1 0.0004870498 0.03506759       CPLM         1.780 48         34           29
#> 2 0.0013718752 0.04938751       CPLM         1.784 48         21           56
#> 3 0.0072123785 0.07418446       CPLM         1.759 48         42           12
#> 4 0.0036238469 0.07418446       CPLM         1.799 48         19           60
#> 5 0.0045416657 0.07418446       CPLM         1.912 48         48            0
#> 6 0.0070577406 0.07418446       CPLM         1.799 48         32           33
```

### Reading the result table

Each row describes one feature–covariate coefficient. Factor variables
can produce several rows per feature, one per non-reference level. These
examples use `median_comparison = FALSE` and the default log link.

- `feature` , `metadata` , `value` : Feature identifier, tested
  covariate, and factor level or continuous predictor name. Factor
  effects are relative to the reference level.
- `coef` , `stderr` : Estimated coefficient and its standard error on
  the link scale. With a log link, `exp(coef)` is the fitted mean ratio
  at equal offsets. Continuous covariates are standardized by default,
  so their coefficients describe a one-standard-deviation increase
  unless `standardize = FALSE` .
- `pval` , `qval` : Unadjusted test p-value and multiple-testing
  adjusted p-value (Benjamini–Hochberg by default across the reported
  coefficient tests within an assay).
- `base.model` , `tweedie.index` : Model engine/family label and
  estimated or fixed variance power.
- `N` , `N.not.zero` , `percent.zero` : Aligned sample count, number
  with positive measurements, and rounded percentage of observed
  measurements equal to zero. Missing model covariates can reduce a
  particular fit’s effective sample count.

`return_all = TRUE` keeps the full result table in memory. The default
returns only rows with `qval <= max_significance` (default 0.05). An
empty default table means no returned associations meet that cutoff,
rather than no models were attempted. Failed fits produce feature-named
warnings and are absent from the result table; inspect those warnings
before interpreting the results. The number of hypotheses here covers
only the selected features; these q-values are not results from testing
the entire source dataset.

``` r

diagnosis_results <- subset(microbiome_fit, metadata == "diagnosis")
diagnosis_results$mean_ratio <- exp(diagnosis_results$coef)
head(diagnosis_results[, c(
    "feature", "value", "coef", "mean_ratio", "pval", "qval"
)], 2)
#>                    feature value      coef  mean_ratio         pval       qval
#> 1      Clostridium bolteae    CD  3.133765 22.96026215 0.0004870498 0.03506759
#> 2 Bacteroides massiliensis    CD -3.849450  0.02129144 0.0013718752 0.04938751
```

The first displayed diagnosis association is **Clostridium bolteae** in
**CD** relative to nonIBD: its coefficient is 3.13 and its adjusted
p-value is 0.0351. The second is **Bacteroides massiliensis**, with
coefficient -3.85 and adjusted p-value 0.0494. A positive coefficient
indicates higher fitted relative abundance and a negative coefficient
lower fitted relative abundance, adjusting for age and antibiotic use.
These are conditional associations in this subset, not evidence of
causality.

### Visualizing associations

For an exploratory overview, the heatmap below includes associations
meeting `qval <= 0.25`, an explicitly relaxed 25% FDR threshold. It
shows up to five features and all retained contrasts for those features.
The package heatmap uses signed `-log(qval)`, clipped to \[-20, 20\]:
red denotes positive abundance coefficients and blue negative
coefficients. Larger magnitudes indicate stronger evidence. Empty
feature–contrast combinations are displayed as zero; they do not
establish absence of an effect. The plot’s significance label refers to
the chosen 0.25 cutoff.

``` r

heatmap_results <- subset(microbiome_fit, !is.na(qval) & qval <= 0.25)
Tweedieverse_heatmap(
    output_results = heatmap_results,
    cell_value = "qval",
    first_n = 5,
    title = "HMP2 associations (q <= 0.25)",
    data_label = "Taxon",
    metadata_label = "Covariate / contrast"
)
```

![](Tweedieverse_files/figure-html/microbiome-heatmap-1.png)

## Single-Cell RNA-seq UMI Count Example

For single-cell RNA-seq data,
[`Tweedieverse()`](https://himelmallick.github.io/Tweedieverse/reference/Tweedieverse.md)
accepts a `SingleCellExperiment`. The UMI example uses integer molecule
counts.

This subset contains 20 genes and 60 PBMC3k cells, with 30 CD4 and 30
CD8 cells. Cell labels and raw UMI counts come from the processed object
used in the [original UMI
tutorial](https://github.com/himelmallick/TweedieLabs) . Genes are
selected by count variance across the retained cells.

UMI counts are count-like. Here we pin the variance power to `p = 1`,
the Poisson variance case, and use the single-cell domain. The library
sizes in colData were calculated using all source genes before selecting
the 20 genes. Passing `scale_factor = "library_size"` uses those
full-library sizes as offsets.

``` r

umi_sce
#> class: SingleCellExperiment 
#> dim: 20 60 
#> metadata(0):
#> assays(1): counts
#> rownames(20): MALAT1 TMSB4X ... RPS27A RPL10A
#> rowData names(0):
#> colnames(60): AAACGCTGACCAGT AAACGCTGGTTCTT ... AACCTACTGTGTTG
#>   AACCTTACGAGACG
#> colData names(8): orig.ident nCount_RNA ... CellType library_size
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

umi_fit <- Tweedieverse(
    input_features = umi_sce,
    output = NULL,
    assay_name = "counts",
    domain = "single_cell",
    scale_factor = "library_size",
    fixed_effects = "CellType",
    reference = "CellType,CD8",
    tweedie_p = 1,
    median_comparison = FALSE
)
#> The specified assay has been extracted
#> 2026-11-21 21:00:52.109301 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:52.110266 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:52.133683 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:52.134252 INFO::Total samples in data: 60
#> 2026-11-21 21:00:52.13463 INFO::Min samples required with min abundance for a feature not to be filtered: 6.000000
#> 2026-11-21 21:00:52.182647 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:52.183498 INFO::Filtered feature names:
#> 2026-11-21 21:00:52.201443 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:52.202062 INFO::Filtered feature names:
#> 2026-11-21 21:00:52.204972 INFO::Offset size factors computed with normalization: USER
#> 2026-11-21 21:00:52.207433 INFO::Excluded metadata with entropy less or equal to 0: orig.ident
#> 2026-11-21 21:00:52.207948 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-21 21:00:52.209833 INFO::Formula for fixed effects: expr ~  CellType
#> 2026-11-21 21:00:52.210289 INFO::Running selected analysis method: CPLM
#> 2026-11-21 21:00:52.2239 INFO::Fitting model to feature number 1, MALAT1
#> 2026-11-21 21:00:52.227526 INFO::Fitting model to feature number 2, TMSB4X
#> 2026-11-21 21:00:52.230321 INFO::Fitting model to feature number 3, B2M
#> 2026-11-21 21:00:52.233559 INFO::Fitting model to feature number 4, RPL10
#> 2026-11-21 21:00:52.236387 INFO::Fitting model to feature number 5, RPS6
#> 2026-11-21 21:00:52.239197 INFO::Fitting model to feature number 6, RPL13
#> 2026-11-21 21:00:52.241916 INFO::Fitting model to feature number 7, RPL3
#> 2026-11-21 21:00:52.24468 INFO::Fitting model to feature number 8, RPL19
#> 2026-11-21 21:00:52.247566 INFO::Fitting model to feature number 9, RPLP1
#> 2026-11-21 21:00:52.250336 INFO::Fitting model to feature number 10, RPS2
#> 2026-11-21 21:00:52.253042 INFO::Fitting model to feature number 11, RPS18
#> 2026-11-21 21:00:52.255808 INFO::Fitting model to feature number 12, RPL13A
#> 2026-11-21 21:00:52.258479 INFO::Fitting model to feature number 13, RPS12
#> 2026-11-21 21:00:52.261225 INFO::Fitting model to feature number 14, RPL14
#> 2026-11-21 21:00:52.264064 INFO::Fitting model to feature number 15, RPL11
#> 2026-11-21 21:00:52.266869 INFO::Fitting model to feature number 16, RPL32
#> 2026-11-21 21:00:52.269627 INFO::Fitting model to feature number 17, MPC2
#> 2026-11-21 21:00:52.272709 INFO::Fitting model to feature number 18, RPS4X
#> 2026-11-21 21:00:52.275466 INFO::Fitting model to feature number 19, RPS27A
#> 2026-11-21 21:00:52.278246 INFO::Fitting model to feature number 20, RPL10A
#> 2026-11-21 21:00:52.286453 INFO::Counting prevalence for each feature

head(umi_fit)
#>   feature metadata value       coef     stderr         pval        qval
#> 1    RPS6 CellType   CD4  0.4006848 0.09630737 0.0001063576 0.001920124
#> 2   RPL32 CellType   CD4  0.3784633 0.09501197 0.0001920124 0.001920124
#> 3     B2M CellType   CD4 -0.3007403 0.08640942 0.0009575612 0.003830245
#> 4   RPLP1 CellType   CD4  0.4495490 0.12908761 0.0009513840 0.003830245
#> 5  RPL13A CellType   CD4  0.2438439 0.06731234 0.0006148881 0.003830245
#> 6   RPS18 CellType   CD4  0.3306069 0.10218562 0.0020085431 0.006695144
#>    base.model tweedie.index  N N.not.zero percent.zero
#> 1 Tweedie GLM             1 60         60            0
#> 2 Tweedie GLM             1 60         60            0
#> 3 Tweedie GLM             1 60         60            0
#> 4 Tweedie GLM             1 60         60            0
#> 5 Tweedie GLM             1 60         60            0
#> 6 Tweedie GLM             1 60         60            0
```

The `CellType` rows compare CD4 with the CD8 reference group. For
example, **RPS6** has coefficient 0.4 and q-value 0.00192. At equal
library size, its fitted CD4/CD8 mean ratio is 1.49. This is a cell-type
illustration from one donor. Cells are not independent biological
replicates for population-level treatment inference; use replicated
samples and an appropriate design for that question.

## Single-Cell RNA-seq Non-UMI / Read Count Example

The non-UMI example uses read-count-like measurements. These can have a
more flexible mean-variance relationship than UMI counts, so we leave
`tweedie_p = NULL` and let Tweedieverse estimate the index after
filtering.

The source is the original Tabula Muris Smart-seq2 FACS dataset,
distributed under CC BY 4.0 in [Figshare version
8](https://doi.org/10.6084/m9.figshare.5829687.v8). This example retains
20 genes and 60 cells from mouse 3_10_M, balanced between astrocytes and
oligodendrocytes. The SCRAN factors in colData were computed on all
source genes for these cells before selecting the 20 most variable
genes. The source counts are unchanged. This is a cell-type illustration
within one mouse; it does not establish population-level differential
expression across biological replicates. See the [original
study](https://doi.org/10.1038/s41586-018-0590-4).

``` r

nonumi_sce
#> class: SingleCellExperiment 
#> dim: 20 60 
#> metadata(8): dataset source_doi ... source_md5 processing
#> assays(1): counts
#> rownames(20): Rn45s Plp1 ... Gsn Mfge8
#> rowData names(0):
#> colnames(60): A13.MAA000560.3_10_M.1.1 A15.MAA000560.3_10_M.1.1 ...
#>   B19.MAA000550.3_10_M.1.1 B19.MAA000564.3_10_M.1.1
#> colData names(26): Neurog3>0_raw Neurog3>0_scaled ... celltype
#>   source_size_factor
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

nonumi_fit <- Tweedieverse(
    input_features = nonumi_sce,
    output = NULL,
    assay_name = "counts",
    domain = "single_cell",
    scale_factor = "source_size_factor",
    fixed_effects = "celltype",
    reference = "celltype,astrocyte",
    standardize = FALSE,
    tweedie_p = NULL,
    median_comparison = FALSE
)
#> The specified assay has been extracted
#> 2026-11-21 21:00:52.492244 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:52.493211 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:52.494334 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:52.494723 INFO::Total samples in data: 60
#> 2026-11-21 21:00:52.495121 INFO::Min samples required with min abundance for a feature not to be filtered: 6.000000
#> 2026-11-21 21:00:52.495527 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:52.495955 INFO::Filtered feature names:
#> 2026-11-21 21:00:52.49655 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:52.497008 INFO::Filtered feature names:
#> 2026-11-21 21:00:52.497417 INFO::Offset size factors computed with normalization: USER
#> 2026-11-21 21:00:52.500109 INFO::Excluded metadata with entropy less or equal to 0: Neurog3.0_raw,Neurog3.0_scaled,free_annotation,mouse.id,mouse.sex,subsetB,subsetB_cluster.ids,subsetC,subsetC_cluster.ids,subsetD,subsetD_cluster.ids,subsetE,subsetE_cluster.ids,tissue
#> 2026-11-21 21:00:52.500611 INFO::Bypass z-score application to metadata
#> 2026-11-21 21:00:52.501054 INFO::Formula for fixed effects: expr ~  celltype
#> 2026-11-21 21:00:52.501464 INFO::Running selected analysis method: CPLM
#> 2026-11-21 21:00:52.51472 INFO::Fitting model to feature number 1, Rn45s
#> 2026-11-21 21:00:52.694451 INFO::Fitting model to feature number 2, Plp1
#> 2026-11-21 21:00:52.837879 INFO::Fitting model to feature number 3, Cnp
#> 2026-11-21 21:00:52.995595 INFO::Fitting model to feature number 4, Malat1
#> 2026-11-21 21:00:53.199742 INFO::Fitting model to feature number 5, Trf
#> 2026-11-21 21:00:53.341814 INFO::Fitting model to feature number 6, Mbp
#> 2026-11-21 21:00:53.439664 INFO::Fitting model to feature number 7, Ptgds
#> 2026-11-21 21:00:53.589407 INFO::Fitting model to feature number 8, Apoe
#> 2026-11-21 21:00:53.725608 INFO::Fitting model to feature number 9, Mag
#> 2026-11-21 21:00:53.84647 INFO::Fitting model to feature number 10, Mal
#> 2026-11-21 21:00:53.966659 INFO::Fitting model to feature number 11, Clu
#> 2026-11-21 21:00:54.113257 INFO::Fitting model to feature number 12, Fth1
#> 2026-11-21 21:00:54.2657 INFO::Fitting model to feature number 13, Tubb4a
#> 2026-11-21 21:00:54.376868 INFO::Fitting model to feature number 14, Aldoc
#> 2026-11-21 21:00:54.517546 INFO::Fitting model to feature number 15, Aplp1
#> 2026-11-21 21:00:54.652581 INFO::Fitting model to feature number 16, Cst3
#> 2026-11-21 21:00:54.848462 INFO::Fitting model to feature number 17, Cd81
#> 2026-11-21 21:00:54.995076 INFO::Fitting model to feature number 18, Mobp
#> 2026-11-21 21:00:55.099393 INFO::Fitting model to feature number 19, Gsn
#> 2026-11-21 21:00:55.262994 INFO::Fitting model to feature number 20, Mfge8
#> 2026-11-21 21:00:55.417593 INFO::Counting prevalence for each feature

head(nonumi_fit)
#>   feature metadata           value      coef    stderr         pval
#> 1     Trf celltype oligodendrocyte  9.970411 0.4825781 2.001878e-28
#> 2     Cnp celltype oligodendrocyte  4.140469 0.2190457 1.806553e-26
#> 3    Cst3 celltype oligodendrocyte -4.010578 0.2415055 1.040934e-23
#> 4     Mbp celltype oligodendrocyte  6.254016 0.3965322 1.215910e-22
#> 5   Aplp1 celltype oligodendrocyte  2.225031 0.1469565 8.236202e-22
#> 6     Mal celltype oligodendrocyte 10.409580 0.7187001 6.354322e-21
#>           qval base.model tweedie.index  N N.not.zero percent.zero
#> 1 4.003755e-27       CPLM         1.693 60         40           33
#> 2 1.806553e-25       CPLM         1.506 60         52           13
#> 3 6.939561e-23       CPLM         1.741 60         58            3
#> 4 6.079551e-22       CPLM         1.486 60         38           37
#> 5 3.294481e-21       CPLM         1.364 60         57            5
#> 6 2.118107e-20       CPLM         1.586 60         33           45
```

Here the contrast is oligodendrocytes relative to astrocytes. The first
returned gene, **Trf**, has coefficient 9.97, q-value 4^{-27}, and
fitted variance power 1.693. The coefficient’s sign describes the
direction after accounting for the precomputed SCRAN factor. The fitted
index describes the variance model; it is not an effect size or
significance measure.

## Pseudobulk Single-Cell RNA-seq Example

Pseudobulk analysis aggregates single-cell counts to sample-level counts
before modeling. This is useful when the biological replicate, rather
than the cell, is the experimental unit. The `muscat` package provides a
Bioconductor-native workflow for this aggregation. Its
[`aggregateData()`](https://rdrr.io/pkg/muscat/man/aggregateData.html)
function summarizes a `SingleCellExperiment` by variables in `colData`,
commonly `cluster_id` and `sample_id`.

This example follows the workflow used by the `muscat` vignette:
multi-subject single-cell data are first represented as a prepared
`SingleCellExperiment`, cells are aggregated into cluster-by-sample
pseudobulk counts, and the resulting sample-level count matrix is
analyzed with
[`Tweedieverse()`](https://himelmallick.github.io/Tweedieverse/reference/Tweedieverse.md).
The package-provided `example_sce` object contains the required
`cluster_id`, `sample_id`, and `group_id` fields and is a lightweight
version of the Kang et al. interferon stimulation example used
throughout `muscat`.

``` r

data("example_sce", package = "muscat")

pseudobulk_input <- example_sce[seq_len(100), ]

pseudobulk_sce <- muscat::aggregateData(
    pseudobulk_input,
    assay = "counts",
    fun = "sum",
    by = c("cluster_id", "sample_id"),
    verbose = FALSE
)

SummarizedExperiment::assayNames(pseudobulk_sce)
#> [1] "B cells"           "CD14+ Monocytes"   "CD4 T cells"      
#> [4] "CD8 T cells"       "FCGR3A+ Monocytes"

pseudobulk_assay <- "CD4 T cells"
pseudobulk_counts <- SummarizedExperiment::assay(
    pseudobulk_sce,
    pseudobulk_assay
)
top_pseudobulk_genes <- order(
    rowSums(pseudobulk_counts),
    decreasing = TRUE
)[seq_len(20)]

pseudobulk_metadata <- unique(as.data.frame(
    SummarizedExperiment::colData(pseudobulk_input)[, c("sample_id",
        "group_id")]
))
rownames(pseudobulk_metadata) <- pseudobulk_metadata$sample_id
pseudobulk_metadata <- pseudobulk_metadata[
    colnames(pseudobulk_counts), ,
    drop = FALSE
]
pseudobulk_metadata$condition <- pseudobulk_metadata$group_id

pseudobulk_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = pseudobulk_counts[top_pseudobulk_genes, , drop =
        FALSE]),
    colData = pseudobulk_metadata
)

pseudobulk_fit <- Tweedieverse(
    input_features = pseudobulk_se,
    output = NULL,
    assay_name = "counts",
    domain = "bulk_rnaseq",
    fixed_effects = c("condition"),
    reference = "condition,ctrl",
    tweedie_p = 1,
    median_comparison = FALSE,
    max_significance = 1
)
#> The specified assay has been extracted
#> 2026-11-21 21:00:55.961173 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:55.962157 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:55.963079 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:55.963481 INFO::Total samples in data: 4
#> 2026-11-21 21:00:55.963864 INFO::Min samples required with min abundance for a feature not to be filtered: 0.400000
#> 2026-11-21 21:00:55.964258 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:55.964661 INFO::Filtered feature names:
#> 2026-11-21 21:00:55.96529 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:55.965698 INFO::Filtered feature names:
#> calcNormFactors has been renamed to normLibSizes
#> 2026-11-21 21:00:55.969991 INFO::Offset size factors computed with normalization: TMM
#> 2026-11-21 21:00:55.970693 INFO::Excluded metadata with entropy less or equal to 0:
#> 2026-11-21 21:00:55.971222 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-21 21:00:55.97329 INFO::Formula for fixed effects: expr ~  condition
#> 2026-11-21 21:00:55.973743 INFO::Running selected analysis method: CPLM
#> 2026-11-21 21:00:55.986892 INFO::Fitting model to feature number 1, RPL11
#> 2026-11-21 21:00:55.989587 INFO::Fitting model to feature number 2, RPS27
#> 2026-11-21 21:00:55.991977 INFO::Fitting model to feature number 3, RPS8
#> 2026-11-21 21:00:55.994381 INFO::Fitting model to feature number 4, ISG15
#> 2026-11-21 21:00:55.996878 INFO::Fitting model to feature number 5, RPL5
#> 2026-11-21 21:00:55.999347 INFO::Fitting model to feature number 6, IFI6
#> 2026-11-21 21:00:56.00185 INFO::Fitting model to feature number 7, SH3BGRL3
#> 2026-11-21 21:00:56.004249 INFO::Fitting model to feature number 8, YBX1
#> 2026-11-21 21:00:56.006697 INFO::Fitting model to feature number 9, S100A11
#> 2026-11-21 21:00:56.009053 INFO::Fitting model to feature number 10, ENO1
#> 2026-11-21 21:00:56.011501 INFO::Fitting model to feature number 11, CD53
#> 2026-11-21 21:00:56.013952 INFO::Fitting model to feature number 12, S100A10
#> 2026-11-21 21:00:56.016298 INFO::Fitting model to feature number 13, S100A6
#> 2026-11-21 21:00:56.018769 INFO::Fitting model to feature number 14, TMEM59
#> 2026-11-21 21:00:56.021235 INFO::Fitting model to feature number 15, CD2
#> 2026-11-21 21:00:56.023629 INFO::Fitting model to feature number 16, NDUFS5
#> 2026-11-21 21:00:56.025998 INFO::Fitting model to feature number 17, CD52
#> 2026-11-21 21:00:56.047689 INFO::Fitting model to feature number 18, LAPTM5
#> 2026-11-21 21:00:56.050522 INFO::Fitting model to feature number 19, SYF2
#> 2026-11-21 21:00:56.053 INFO::Fitting model to feature number 20, RPL22
#> 2026-11-21 21:00:56.061074 INFO::Counting prevalence for each feature

head(pseudobulk_fit)
#>    feature  metadata value       coef     stderr        pval      qval
#> 1     IFI6 condition  stim  3.8839623 0.27615612 0.005017428 0.1003486
#> 2    ISG15 condition  stim  3.8805689 0.57303845 0.021117751 0.2111775
#> 3     RPL5 condition  stim -0.4307140 0.12650534 0.076498419 0.3059937
#> 4     ENO1 condition  stim  0.4636025 0.13026340 0.070681909 0.3059937
#> 5    RPL22 condition  stim -0.2870616 0.07274647 0.058629532 0.3059937
#> 6 SH3BGRL3 condition  stim -0.1499796 0.07721882 0.191591787 0.4414440
#>    base.model tweedie.index N N.not.zero percent.zero
#> 1 Tweedie GLM             1 4          4            0
#> 2 Tweedie GLM             1 4          4            0
#> 3 Tweedie GLM             1 4          4            0
#> 4 Tweedie GLM             1 4          4            0
#> 5 Tweedie GLM             1 4          4            0
#> 6 Tweedie GLM             1 4          4            0
```

The `condition` coefficient compares stimulated samples with control
samples within CD4 T cells. For **IFI6**, the coefficient is 3.88 and
q-value is 0.1. The rows represent sample-level pseudobulks rather than
individual cells, so `N` counts pseudobulk samples.
`max_significance = 1` displays all available tests, including weak
evidence; it does not declare them all significant. This compact example
omits donor blocking; a paired stimulation study should include donor in
the model.

## Bulk RNA-seq Count Example

For bulk RNA-seq data,
[`Tweedieverse()`](https://himelmallick.github.io/Tweedieverse/reference/Tweedieverse.md)
accepts a `SummarizedExperiment`. This example uses the Bioconductor
`airway` dataset, which contains RNA-seq read counts from airway smooth
muscle cells treated with dexamethasone.

The full `airway` object contains tens of thousands of genes. To keep
the vignette fast, the example keeps a small set of high-count genes
while preserving the original `airway` sample metadata. Bulk RNA-seq
uses `TMM` as the default normalization strategy in `Tweedieverse`; this
computes library-size factors for the model offset without modifying the
count assay.

Set `run_presence_absence_model = TRUE` to include a test for whether
the feature is detected. By default, the prevalence fit requires at
least three samples with nonzero counts and three with zero counts for
each gene. Use `min_presence_absence_count` to change both minimums. In
bulk and pseudobulk data, genes with only one or two zeros can still
enter the abundance analysis. For retained genes below either minimum,
`pval_presence` is `NA` and `pval_abundance` becomes the joint `pval`.
When the two reported p-values are exactly 0 and 1, the joint `pval` is
0. Tweedieverse applies the selected multiple-testing correction to the
joint p-values to compute `qval`.

Use `return_all = TRUE` to retrieve the full association table with
`output = NULL`, including associations above `max_significance`. The
default returns associations passing the threshold. Missing q-values
after median comparison are retained in the full table and excluded from
significant results.

``` r

data("airway", package = "airway")

airway_counts <- SummarizedExperiment::assay(airway, "counts")
top_airway_genes <- order(rowSums(airway_counts), decreasing = TRUE)[seq_len(
    20)]

bulk_rnaseq_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = airway_counts[top_airway_genes, , drop = FALSE]),
    colData = SummarizedExperiment::colData(airway)
)

bulk_rnaseq_se
#> class: SummarizedExperiment 
#> dim: 20 8 
#> metadata(0):
#> assays(1): counts
#> rownames(20): ENSG00000115414 ENSG00000011465 ... ENSG00000137801
#>   ENSG00000142173
#> rowData names(0):
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample

bulk_rnaseq_fit <- Tweedieverse(
    input_features = bulk_rnaseq_se,
    output = NULL,
    assay_name = "counts",
    domain = "bulk_rnaseq",
    fixed_effects = c("dex"),
    reference = "dex,untrt",
    tweedie_p = 1,
    median_comparison = FALSE,
    run_presence_absence_model = TRUE,
    max_significance = 0.05,
    return_all = TRUE
)
#> The specified assay has been extracted
#> 2026-11-21 21:00:57.311591 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:57.312578 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:57.313604 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:57.314022 INFO::Total samples in data: 8
#> 2026-11-21 21:00:57.314389 INFO::Min samples required with min abundance for a feature not to be filtered: 0.800000
#> 2026-11-21 21:00:57.314808 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:57.315223 INFO::Filtered feature names:
#> 2026-11-21 21:00:57.315815 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:57.316239 INFO::Filtered feature names:
#> calcNormFactors has been renamed to normLibSizes
#> 2026-11-21 21:00:57.319399 INFO::Offset size factors computed with normalization: TMM
#> 2026-11-21 21:00:57.320291 INFO::Excluded metadata with entropy less or equal to 0: albut
#> 2026-11-21 21:00:57.320781 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-21 21:00:57.32288 INFO::Formula for fixed effects: expr ~  dex
#> 2026-11-21 21:00:57.323333 INFO::Running selected analysis method: CPLM
#> 2026-11-21 21:00:57.336726 INFO::Fitting model to feature number 1, ENSG00000115414
#> 2026-11-21 21:00:57.339471 INFO::Fitting model to feature number 2, ENSG00000011465
#> 2026-11-21 21:00:57.341968 INFO::Fitting model to feature number 3, ENSG00000198804
#> 2026-11-21 21:00:57.344415 INFO::Fitting model to feature number 4, ENSG00000156508
#> 2026-11-21 21:00:57.346855 INFO::Fitting model to feature number 5, ENSG00000164692
#> 2026-11-21 21:00:57.349362 INFO::Fitting model to feature number 6, ENSG00000116260
#> 2026-11-21 21:00:57.351952 INFO::Fitting model to feature number 7, ENSG00000087086
#> 2026-11-21 21:00:57.35449 INFO::Fitting model to feature number 8, ENSG00000115461
#> 2026-11-21 21:00:57.357 INFO::Fitting model to feature number 9, ENSG00000077942
#> 2026-11-21 21:00:57.359417 INFO::Fitting model to feature number 10, ENSG00000166923
#> 2026-11-21 21:00:57.361977 INFO::Fitting model to feature number 11, ENSG00000087245
#> 2026-11-21 21:00:57.36446 INFO::Fitting model to feature number 12, ENSG00000091986
#> 2026-11-21 21:00:57.366953 INFO::Fitting model to feature number 13, ENSG00000103888
#> 2026-11-21 21:00:57.369529 INFO::Fitting model to feature number 14, ENSG00000210082
#> 2026-11-21 21:00:57.372106 INFO::Fitting model to feature number 15, ENSG00000124942
#> 2026-11-21 21:00:57.37464 INFO::Fitting model to feature number 16, ENSG00000108821
#> 2026-11-21 21:00:57.377207 INFO::Fitting model to feature number 17, ENSG00000196205
#> 2026-11-21 21:00:57.379719 INFO::Fitting model to feature number 18, ENSG00000150991
#> 2026-11-21 21:00:57.382287 INFO::Fitting model to feature number 19, ENSG00000137801
#> 2026-11-21 21:00:57.384941 INFO::Fitting model to feature number 20, ENSG00000142173
#> 2026-11-21 21:00:57.393151 INFO::Counting prevalence for each feature

head(bulk_rnaseq_fit)
#>           feature metadata value       coef     stderr         pval        qval
#> 1 ENSG00000108821      dex   trt -0.8287049 0.08821054 8.260095e-05 0.001652019
#> 2 ENSG00000137801      dex   trt  1.3972630 0.23161755 9.374555e-04 0.009374555
#> 3 ENSG00000124942      dex   trt  0.2768458 0.07273457 8.902531e-03 0.059350204
#> 4 ENSG00000164692      dex   trt -0.2559281 0.08206528 2.062187e-02 0.095415880
#> 5 ENSG00000150991      dex   trt  0.2134059 0.07101706 2.385397e-02 0.095415880
#> 6 ENSG00000087086      dex   trt -0.3890967 0.14496888 3.634185e-02 0.121139513
#>   coef_abundance stderr_abundance pval_abundance qval_abundance coef_presence
#> 1     -0.8287049       0.08821054   8.260095e-05    0.001652019            NA
#> 2      1.3972630       0.23161755   9.374555e-04    0.009374555            NA
#> 3      0.2768458       0.07273457   8.902531e-03    0.059350204            NA
#> 4     -0.2559281       0.08206528   2.062187e-02    0.095415880            NA
#> 5      0.2134059       0.07101706   2.385397e-02    0.095415880            NA
#> 6     -0.3890967       0.14496888   3.634185e-02    0.121139513            NA
#>   stderr_presence pval_presence qval_presence base.model base.model_abundance
#> 1              NA            NA            NA        CCT          Tweedie GLM
#> 2              NA            NA            NA        CCT          Tweedie GLM
#> 3              NA            NA            NA        CCT          Tweedie GLM
#> 4              NA            NA            NA        CCT          Tweedie GLM
#> 5              NA            NA            NA        CCT          Tweedie GLM
#> 6              NA            NA            NA        CCT          Tweedie GLM
#>   base.model_presence tweedie.index
#> 1                  NA             1
#> 2                  NA             1
#> 3                  NA             1
#> 4                  NA             1
#> 5                  NA             1
#> 6                  NA             1
```

The `dex` contrast compares dexamethasone-treated with untreated
samples. For **ENSG00000108821**, `coef_abundance` is -0.83 and the
joint q-value is 0.00165. A positive abundance coefficient indicates
higher offset-adjusted mean counts in treated samples.

With the presence model enabled, `coef` and `stderr` retain
abundance-model meaning; the joint CCT p-value is not a test of a
combined effect size. `coef_presence` is a logistic-model log-odds
coefficient, while `pval_abundance`/`qval_abundance` and
`pval_presence`/`qval_presence` describe component tests. Inspect
`pval_presence` to see which genes lacked enough zeros or detections to
fit the prevalence model. For those genes, the joint p-value falls back
to the abundance p-value. `return_all = TRUE` includes rows above 0.05;
use `subset(bulk_rnaseq_fit, !is.na(qval) & qval <= 0.05)` for discovery
selection. The full airway experiment is paired by cell line; include
`cell` alongside `dex` when analyzing that study design.

## Multi-omics Container

For multi-omics data,
[`Tweedieverse()`](https://himelmallick.github.io/Tweedieverse/reference/Tweedieverse.md)
accepts a `MultiAssayExperiment`. The packaged Franzosa et al. example
is stored as an individual `.rds` file in `inst/extdata`. It contains
matched microbiome species abundances and metabolite profiles. Each
omics layer is analyzed independently, and per-omics parameters can be
supplied as named vectors or lists.

``` r

franzosa_example_mae
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] taxaAbundanceProfile: SummarizedExperiment with 4 rows and 16 columns
#>  [2] harmonizedMetaboliteProfile: SummarizedExperiment with 4 rows and 16 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files

multi_fit <- Tweedieverse(
    input_features = franzosa_example_mae,
    output = NULL,
    assay_name = "abundance",
    domain = c(
        taxaAbundanceProfile = "microbiome",
        harmonizedMetaboliteProfile = "custom"
    ),
    normalization = c(
        taxaAbundanceProfile = "NONE",
        harmonizedMetaboliteProfile = "MEDIAN"
    ),
    fixed_effects = c("diseaseCat_clean", "Age"),
    tweedie_p = list(
        taxaAbundanceProfile = NULL,
        harmonizedMetaboliteProfile = NULL
    ),
    median_comparison = FALSE,
    max_significance = 1
)
#> The specified assay has been extracted
#> 2026-11-21 21:00:58.443863 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:58.445118 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:58.446095 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:58.446523 INFO::Total samples in data: 16
#> 2026-11-21 21:00:58.446929 INFO::Min samples required with min abundance for a feature not to be filtered: 1.600000
#> 2026-11-21 21:00:58.447334 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:58.447727 INFO::Filtered feature names:
#> 2026-11-21 21:00:58.448238 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:58.44864 INFO::Filtered feature names:
#> 2026-11-21 21:00:58.449045 INFO::Offset size factors computed with normalization: NONE
#> 2026-11-21 21:00:58.449684 INFO::Excluded metadata with entropy less or equal to 0:
#> 2026-11-21 21:00:58.450188 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-21 21:00:58.452952 INFO::Formula for fixed effects: expr ~  diseaseCat_clean + Age
#> 2026-11-21 21:00:58.453412 INFO::Running selected analysis method: CPLM
#> 2026-11-21 21:00:58.467813 INFO::Fitting model to feature number 1, Bacteroides_vulgatus
#> 2026-11-21 21:00:58.579274 INFO::Fitting model to feature number 2, Ruminococcus_gnavus
#> 2026-11-21 21:00:58.678505 INFO::Fitting model to feature number 3, Eubacterium_rectale
#> 2026-11-21 21:00:58.857907 INFO::Fitting model to feature number 4, Faecalibacterium_prausnitzii
#> 2026-11-21 21:00:58.965734 INFO::Counting prevalence for each feature
#> The specified assay has been extracted
#> 2026-11-21 21:00:58.97619 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:58.977138 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:58.97804 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:58.97845 INFO::Total samples in data: 16
#> 2026-11-21 21:00:58.978836 INFO::Min samples required with min abundance for a feature not to be filtered: 1.600000
#> 2026-11-21 21:00:58.979231 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:58.979631 INFO::Filtered feature names:
#> 2026-11-21 21:00:58.980081 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:58.980479 INFO::Filtered feature names:
#> 2026-11-21 21:00:58.98144 INFO::Offset size factors computed with normalization: MEDIAN
#> 2026-11-21 21:00:58.982114 INFO::Excluded metadata with entropy less or equal to 0:
#> 2026-11-21 21:00:58.982583 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-21 21:00:58.985115 INFO::Formula for fixed effects: expr ~  diseaseCat_clean + Age
#> 2026-11-21 21:00:58.985567 INFO::Running selected analysis method: CPLM
#> 2026-11-21 21:00:58.998823 INFO::Fitting model to feature number 1, C54:3 TAG
#> 2026-11-21 21:00:59.125011 INFO::Fitting model to feature number 2, 3-Carboxypropyl trimethylammonium
#> 2026-11-21 21:00:59.278584 INFO::Fitting model to feature number 3, C54:4 TAG
#> 2026-11-21 21:00:59.32457 INFO::Fitting model to feature number 4, C54:5 TAG
#> 2026-11-21 21:00:59.514212 INFO::Counting prevalence for each feature

names(multi_fit)
#> [1] "taxaAbundanceProfile"        "harmonizedMetaboliteProfile"
head(multi_fit$taxaAbundanceProfile)
#>                  omics                      feature         metadata   value
#> 1 taxaAbundanceProfile          Eubacterium_rectale              Age     Age
#> 2 taxaAbundanceProfile          Eubacterium_rectale diseaseCat_clean Healthy
#> 3 taxaAbundanceProfile          Ruminococcus_gnavus              Age     Age
#> 4 taxaAbundanceProfile Faecalibacterium_prausnitzii diseaseCat_clean Healthy
#> 5 taxaAbundanceProfile         Bacteroides_vulgatus              Age     Age
#> 6 taxaAbundanceProfile Faecalibacterium_prausnitzii              Age     Age
#>         coef    stderr        pval        qval base.model tweedie.index  N
#> 1 -2.7824200 0.6723543 0.001166404 0.009331229       CPLM         1.959 16
#> 2 -4.0777339 1.2742125 0.006964917 0.027859668       CPLM         1.959 16
#> 3  1.6620788 0.8757181 0.080128960 0.160257919       CPLM         1.990 16
#> 4  2.2109010 1.1622964 0.079534477 0.160257919       CPLM         1.869 16
#> 5  0.5819870 0.5659511 0.322546636 0.446136909       CPLM         1.951 16
#> 6  0.5975592 0.5963310 0.334602682 0.446136909       CPLM         1.869 16
#>   N.not.zero percent.zero
#> 1         16            0
#> 2         16            0
#> 3         16            0
#> 4         16            0
#> 5         16            0
#> 6         16            0
head(multi_fit$harmonizedMetaboliteProfile)
#>                         omics                           feature
#> 1 harmonizedMetaboliteProfile                         C54:3 TAG
#> 2 harmonizedMetaboliteProfile 3-Carboxypropyl trimethylammonium
#> 3 harmonizedMetaboliteProfile                         C54:5 TAG
#> 4 harmonizedMetaboliteProfile                         C54:3 TAG
#> 5 harmonizedMetaboliteProfile                         C54:5 TAG
#> 6 harmonizedMetaboliteProfile 3-Carboxypropyl trimethylammonium
#>           metadata   value       coef    stderr       pval      qval base.model
#> 1              Age     Age  0.4576428 0.2540086 0.09482193 0.2528585       CPLM
#> 2 diseaseCat_clean Healthy -3.6533834 1.9329772 0.08125656 0.2528585       CPLM
#> 3              Age     Age -0.3291460 0.1750951 0.08273022 0.2528585       CPLM
#> 4 diseaseCat_clean Healthy  0.5723166 0.4908363 0.26456010 0.5291202       CPLM
#> 5 diseaseCat_clean Healthy -0.2784898 0.3378094 0.42458092 0.6793295       CPLM
#> 6              Age     Age -0.4177465 1.0004811 0.68309312 0.6830931       CPLM
#>   tweedie.index  N N.not.zero percent.zero
#> 1          1.99 16         16            0
#> 2          1.99 16         16            0
#> 3          1.99 16         16            0
#> 4          1.99 16         16            0
#> 5          1.99 16         16            0
#> 6          1.99 16         16            0
```

Each list element is a separate association table, and its `omics`
column identifies the layer. `diseaseCat_clean` rows compare the named
`value` with the reference disease category; `Age` rows describe a
one-standard-deviation increase because continuous predictors are
standardized by default. Within either layer, a positive log-link
coefficient means a higher fitted mean, conditional on the other
covariate and offset. Taxon abundances and metabolite intensities have
different measurement scales, so coefficients should be interpreted
within layers. q-values are adjusted separately per layer, not across
the entire multi-omics list. As above, `max_significance = 1` displays
tests without requiring a discovery threshold.

Median comparison uses Monte Carlo simulations when
`median_comparison = TRUE`. Set `median_comparison_seed = 123` to
reproduce those simulations. The seed is scoped with `withr`, restoring
the preceding random-number state after the median-comparison
calculation. When the seed is `NULL`, results can vary between runs.

## Session Info

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] airway_1.32.0                   muscat_1.26.0                  
#>  [3] MultiAssayExperiment_1.38.0     TreeSummarizedExperiment_2.20.0
#>  [5] Biostrings_2.80.2               XVector_0.52.0                 
#>  [7] SingleCellExperiment_1.34.0     SummarizedExperiment_1.42.0    
#>  [9] Biobase_2.72.0                  GenomicRanges_1.64.0           
#> [11] Seqinfo_1.2.0                   IRanges_2.46.0                 
#> [13] S4Vectors_0.50.3                BiocGenerics_0.58.1            
#> [15] generics_0.1.4                  MatrixGenerics_1.24.0          
#> [17] matrixStats_1.5.0               Tweedieverse_0.99.3            
#> [19] BiocStyle_2.40.0               
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.6.2            bitops_1.1-0             tibble_3.3.1            
#>   [4] lifecycle_1.0.5          Rdpack_2.6.6             edgeR_4.10.7            
#>   [7] doParallel_1.0.17        lattice_0.23-1           MASS_7.3-66             
#>  [10] backports_1.5.1          magrittr_2.0.5           limma_3.68.5            
#>  [13] sass_0.5.0               rmarkdown_2.32           jquerylib_0.1.4         
#>  [16] yaml_2.3.12              otel_0.3.0               DBI_1.3.0               
#>  [19] minqa_1.2.8              RColorBrewer_1.1-3       abind_1.4-8             
#>  [22] EnvStats_3.1.0           cplm_0.7-12.1            glmmTMB_1.1.15.2        
#>  [25] purrr_1.2.2              yulab.utils_0.2.6        rappdirs_0.3.4          
#>  [28] sandwich_3.1-3           circlize_0.4.18          ggrepel_0.9.8           
#>  [31] pbkrtest_0.5.5           irlba_2.4.1              tidytree_0.4.9          
#>  [34] pheatmap_1.0.13          tweedie_3.2.2            pkgdown_2.2.1           
#>  [37] codetools_0.2-20         DelayedArray_0.38.2      scuttle_1.22.0          
#>  [40] tidyselect_1.2.1         shape_1.4.6.1            farver_2.1.2            
#>  [43] lme4_2.1-0               ScaledMatrix_1.20.0      viridis_0.6.5           
#>  [46] jsonlite_2.0.0           GetoptLong_1.1.1         BiocNeighbors_2.6.0     
#>  [49] scater_1.40.2            iterators_1.0.14         systemfonts_1.3.3       
#>  [52] foreach_1.5.2            tools_4.6.2              progress_1.2.3          
#>  [55] treeio_1.36.1            ragg_1.5.2               Rcpp_1.1.2              
#>  [58] blme_1.0-8               glue_1.8.1               gridExtra_2.3.1         
#>  [61] SparseArray_1.12.3       xfun_0.62                mgcv_1.9-4              
#>  [64] dplyr_1.2.1              withr_3.0.3              numDeriv_2016.8-1.1     
#>  [67] BiocManager_1.30.27      fastmap_1.2.0            boot_1.3-32             
#>  [70] caTools_1.18.4           digest_0.6.39            rsvd_2.0.0              
#>  [73] R6_2.6.1                 textshaping_1.1.0        colorspace_2.1-4        
#>  [76] gtools_3.9.5             RhpcBLASctl_0.23-42      tidyr_1.3.2             
#>  [79] variancePartition_1.42.0 corpcor_1.6.10           prettyunits_1.2.0       
#>  [82] S4Arrays_1.12.1          pkgconfig_2.0.3          gtable_0.3.6            
#>  [85] ComplexHeatmap_2.28.0    S7_1.0.0                 remaCor_0.0.20          
#>  [88] htmltools_0.5.9          bookdown_0.48            TMB_1.9.25              
#>  [91] clue_0.3-68              scales_1.4.0             logging_0.10-111        
#>  [94] png_0.1-9                fANCOVA_0.6-1            reformulas_0.4.4        
#>  [97] knitr_1.52               reshape2_1.4.5           rjson_0.2.23            
#> [100] coda_0.19-4.1            nlme_3.1-171             nloptr_2.2.1            
#> [103] cachem_1.1.0             zoo_1.9-1                GlobalOptions_0.1.4     
#> [106] stringr_1.6.0            KernSmooth_2.23-27       parallel_4.6.2          
#> [109] vipor_0.4.7              desc_1.4.3               pillar_1.11.1           
#> [112] grid_4.6.2               vctrs_0.7.3              gplots_3.3.0            
#> [115] BiocSingular_1.28.1      beachmat_2.28.0          cluster_2.1.8.3         
#> [118] beeswarm_0.4.0           evaluate_1.0.5           RTMB_2.0                
#> [121] mvtnorm_1.4-2            cli_3.6.6                locfit_1.5-9.12         
#> [124] compiler_4.6.2           rlang_1.3.0              crayon_1.5.3            
#> [127] plyr_1.8.9               fs_2.1.0                 ggbeeswarm_0.7.3        
#> [130] stringi_1.8.9            viridisLite_0.4.3        BiocParallel_1.46.0     
#> [133] lmerTest_3.2-1           lazyeval_0.2.3           aod_1.3.3               
#> [136] Matrix_1.7-6             hms_1.1.4                ggplot2_4.0.4           
#> [139] statmod_1.5.2            rbibutils_2.4.1          broom_1.0.13            
#> [142] bslib_0.12.0             biglm_0.9-3              ape_5.8-1
```
