# Differential analysis of multi-omics data using Tweedie GLMs

Fit a per-feature Tweedie generalized linear model to omics features.

## Usage

``` r
Tweedieverse(
    input_features,
    input_metadata = NULL,
    output = NULL,
    assay_name = "counts",
    abd_threshold = 0,
    prev_threshold = 0.1,
    var_threshold = 0,
    entropy_threshold = 0,
    base_model = "CPLM",
    link = "log",
    tweedie_p = NULL,
    p0_transform = "NONE",
    p0_transform_pseudocount = 1,
    fixed_effects = NULL,
    random_effects = NULL,
    domain = "microbiome",
    normalization = NULL,
    adjust_offset = TRUE,
    scale_factor = NULL,
    max_significance = 0.05,
    correction = "BH",
    Maaslin2_run = TRUE,
    median_comparison = FALSE,
    median_subtraction = FALSE,
    run_presence_absence_model = FALSE,
    method_args = NULL,
    method.args = NULL,
    standardize = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    optimizer = "nlminb",
    na.action = na.exclude,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    heatmap_first_n = 50,
    reference = NULL,
    min_presence_absence_count = 3L,
    return_all = FALSE,
    median_comparison_seed = NULL
)
```

## Arguments

- input_features:

  A domain-appropriate Bioconductor container, such as
  `TreeSummarizedExperiment` for `domain = "microbiome"`,
  `SingleCellExperiment` for `domain = "single_cell"`, and
  `SummarizedExperiment` for `domain = "bulk_rnaseq"`. The `assays` slot
  contains the expression or abundance matrix and defaults to an assay
  named `"counts"`. A `MultiAssayExperiment` object is also accepted. In
  that case, `Tweedieverse()` iterates over each experiment/omics layer
  and returns a named list of omics-specific results. This matrix should
  have one row for each feature and one sample for each column. The
  `colData` slot should contain a data frame with one row per sample and
  columns that contain metadata for each sample. Additional information
  about the experiment can be contained in the `metadata` slot as a
  list.

- input_metadata:

  Optional R data frame of metadata. Samples are expected to have
  matching sample names with `input_features`. If `NULL`, metadata are
  taken from the Bioconductor container's `colData`.

- output:

  The output folder to write results.

- assay_name:

  If the input is provided as one of the accepted Bioconductor objects,
  this argument selects the name of the assay slot in the input object
  that contains the omics measurements. For `MultiAssayExperiment`
  input, most arguments can be supplied as a single value used for every
  omics layer, or as a named list/vector keyed by experiment name.

- abd_threshold:

  If prevalence-abundance filtering is desired, only features that are
  present (or detected) in at least `prev_threshold` percent of samples
  at `abd_threshold` minimum abundance (read count or proportion) are
  retained. Default value for `abd_threshold` is `0.0`. To disable
  prevalence-abundance filtering, set `abd_threshold = -Inf`.

- prev_threshold:

  If prevalence-abundance filtering is desired, only features that are
  present (or detected) in at least `prev_threshold` percent of samples
  at `abd_threshold` minimum abundance (read count or proportion) are
  retained. Default value for `prev_threshold` is `0.1`.

- var_threshold:

  If variance filtering is desired, only features that have variances
  greater than `var_threshold` are retained. This step is done after the
  prevalence-abundance filtering. Default value for `var_threshold` is
  `0.0` (i.e. no variance filtering).

- entropy_threshold:

  Metadata columns are retained when their entropy, computed by
  [`tweedie_entropy()`](https://himelmallick.github.io/Tweedieverse/reference/tweedie_entropy.md),
  exceeds this threshold. Default is 0, which removes metadata columns
  with zero entropy, including constant columns.

- base_model:

  The per-feature base model. Only "CPLM" is supported.

- link:

  A specification of the GLM link function. Default is "log". Must be
  one of "log", "identity", "sqrt", or "inverse".

- tweedie_p:

  Numeric Tweedie variance power. Default is `NULL`, which estimates the
  compound-Poisson Tweedie index in the usual `1 < p < 2` range when all
  filtered feature values are non-negative. If `tweedie_p` is `NULL` or
  `NA` and filtered features contain negative values, it is resolved to
  `0`, the Gaussian case. Values between 0 and 1 are undefined and
  rejected. Other supplied values pin the model to that fixed index.
  Negative supplied powers are allowed with a warning for non-negative
  data, but negative feature values are supported only at `p = 0`.
  Random-effect models support fixed `p = 0`, `p = 1`, `1 < p < 2`,
  `p = 2`, and `p = 3`.

- p0_transform:

  Transformation to apply to a copy of the filtered feature table when
  the resolved Tweedie index is `0`. The original input data are not
  modified. Must be one of `"NONE"`, `"CLR"`, `"RCLR"`, `"LOG"`, or
  `"ARC_SIGNED_SQRT"`. The transformed copy is used for the `p = 0`
  abundance model and written to `transformed_features.tsv` when
  `output` is provided. Default is `"NONE"`.

- p0_transform_pseudocount:

  Numeric pseudocount used for `p0_transform = "CLR"` and
  `p0_transform = "LOG"`. Default is 1.

- fixed_effects:

  Metadata variable(s) describing the fixed effects coefficients.

- random_effects:

  Metadata variable(s) describing the random effects part of the model.

- domain:

  Domain used to choose the default size-factor normalization. Must be
  one of `"microbiome"`, `"single_cell"`, `"bulk_rnaseq"`, or
  `"custom"`. The default is `"microbiome"`.

- normalization:

  Size-factor strategy used only to compute the model offset. The input
  feature table is never normalized or transformed. If `NULL`, the
  domain default is used: `"TSS"` for microbiome, `"SCRAN"` for
  single-cell, `"TMM"` for bulk RNA-seq, and `"MEDIAN"` for custom/other
  omics. Supported strategies are `"TSS"`, `"GMPR"`, `"CSS"`, `"SCRAN"`,
  `"TMM"`, `"RLE"` / `"DESEQ2"`, `"CPM"`, `"MEDIAN"`, and `"NONE"`.

- adjust_offset:

  If TRUE (default), an offset term will be included as the logarithm of
  the size factor estimated by `normalization`, or by `scale_factor`
  when supplied.

- scale_factor:

  Name of a numerical metadata variable containing user-supplied sample
  size factors. When supplied, this overrides `domain` and
  `normalization` for offset calculation.

- max_significance:

  The q-value threshold for significance. Default is 0.05.

- correction:

  The correction method for computing the q-value (see
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for options,
  default is 'BH').

- Maaslin2_run:

  Logical. For `tweedie_p = 0`, run the MaAsLin2 linear-model path
  instead of the Tweedie GLM path. Set to FALSE to fit `tweedie_p = 0`
  with the Tweedie GLM. Default is TRUE. When the resolved Tweedie index
  is `0`, non-identity links are changed to `"identity"` and
  `adjust_offset` is disabled.

- median_comparison:

  If TRUE, coefficients will be tested against a null value
  corresponding to the median coefficient for a covariate in the
  `metadata`. Default is FALSE. Should only be used for relative
  abundance data.

- median_subtraction:

  If TRUE, coefficients minus median will be used for compositionality
  adjustment.

- run_presence_absence_model:

  If TRUE, also fit a DAssemble-style presence-absence logistic
  regression model alongside the selected abundance model, including the
  MaAsLin2 path used by default for `tweedie_p = 0`. Retains individual
  abundance and presence-absence model results, and ranks features by a
  Cauchy combination test of the two p-values. Default is FALSE.

- method_args:

  Optional named list of method-specific arguments. For `tweedie_p = 0`,
  use `method_args = list(Maaslin2 = list(...))` to pass arguments to
  [`Maaslin2`](https://rdrr.io/pkg/Maaslin2/man/Maaslin2.html);
  user-supplied values override Tweedieverse defaults.

- method.args:

  Alias for `method_args`.

- standardize:

  Should continuous metadata be standardized? Default is TRUE. Bypassed
  for categorical variables.

- BPPARAM:

  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object controlling Bioconductor-native parallel evaluation. Default is
  [`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

- optimizer:

  The optimization routine to be used for estimating the parameters of
  the Tweedie model. Possible choices are `"nlminb"` (the default, see
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html)), `"bobyqa"`
  ([`bobyqa`](https://rdrr.io/pkg/minqa/man/bobyqa.html)), and
  `"L-BFGS-B"` ([`optim`](https://rdrr.io/r/stats/optim.html)). Ignored
  for random effects modeling which uses an alternative Template Model
  Builder (TMB) approach
  ([`glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)).

- na.action:

  How to handle missing values? See
  [`na.action`](https://rdrr.io/r/stats/na.action.html). Default is
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html).

- plot_heatmap:

  Logical. If TRUE (default is FALSE), generate a heatmap of the (top
  `heatmap_first_n`) significant associations.

- plot_scatter:

  Logical. If TRUE (default is FALSE), generate scatter/box plots of
  individual associations.

- heatmap_first_n:

  In heatmap, plot top N features with significant associations (default
  is 50).

- reference:

  The factor to use as a reference for a variable with more than two
  levels provided as a string of 'variable,reference' semi-colon
  delimited for multiple variables (default is NULL).

- min_presence_absence_count:

  Minimum number of detected samples and minimum number of undetected
  samples required for the presence-absence fit. Default is 3 for each.
  For retained features below either minimum, the presence-absence
  p-value is missing and the joint p-value uses the abundance p-value
  directly.

- return_all:

  Logical. If TRUE, return the full result table before the significance
  threshold is applied, including when output is NULL. Default is FALSE,
  which returns associations passing max_significance.

- median_comparison_seed:

  Optional non-negative integer seed for the median-comparison
  simulations. Default is NULL, which uses the current random-number
  stream. A supplied seed is scoped with withr and restores the
  preceding random-number state after the median-comparison calculation.

## Value

For single-omics input, a data frame containing coefficient estimates,
p-values, and q-values (multiplicity-adjusted p-values) is returned. For
`MultiAssayExperiment` input, a named list of omics-specific result data
frames is returned. With `return_all = TRUE`, each table contains the
results before filtering by `max_significance`.

## Author

Himel Mallick, <him4004@med.cornell.edu>

## Examples

``` r

set.seed(123)
features <- as.data.frame(matrix(rpois(24, lambda = 5), nrow = 12, ncol = 2))
colnames(features) <- c("feature1", "feature2")
rownames(features) <- paste0("sample", seq_len(nrow(features)))
metadata <- data.frame(
    group = rep(c("A", "B"), each = 6),
    row.names = rownames(features)
)
experiment <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(features))),
    colData = metadata
)
fit <- Tweedieverse(
    input_features = experiment,
    output = NULL,
    fixed_effects = "group",
    median_comparison = FALSE,
    BPPARAM = BiocParallel::SerialParam()
)
#> The specified assay has been extracted
#> 2026-11-21 21:00:28.855442 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:28.862613 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:28.869947 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:28.870357 INFO::Total samples in data: 12
#> 2026-11-21 21:00:28.870721 INFO::Min samples required with min abundance for a feature not to be filtered: 1.200000
#> 2026-11-21 21:00:28.871334 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:28.871782 INFO::Filtered feature names:
#> 2026-11-21 21:00:28.872376 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:28.872777 INFO::Filtered feature names:
#> 2026-11-21 21:00:28.873875 INFO::Offset size factors computed with normalization: TSS
#> 2026-11-21 21:00:28.874632 INFO::Excluded metadata with entropy less or equal to 0:
#> 2026-11-21 21:00:28.875246 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-21 21:00:28.883012 INFO::Formula for fixed effects: expr ~  group
#> 2026-11-21 21:00:28.883679 INFO::Running selected analysis method: CPLM
#> 2026-11-21 21:00:28.904745 INFO::Fitting model to feature number 1, feature1
#> 2026-11-21 21:00:29.700823 INFO::Fitting model to feature number 2, feature2
#> 2026-11-21 21:00:29.762852 INFO::Counting prevalence for each feature
head(fit)
#>  [1] feature       metadata      value         coef          stderr       
#>  [6] pval          qval          base.model    tweedie.index N            
#> [11] N.not.zero    percent.zero 
#> <0 rows> (or 0-length row.names)

maaslin2_median_fit <- Tweedieverse(
    input_features = experiment,
    output = NULL,
    fixed_effects = "group",
    tweedie_p = 0,
    Maaslin2_run = TRUE,
    median_comparison = TRUE,
    median_subtraction = TRUE,
    BPPARAM = BiocParallel::SerialParam()
)
#> The specified assay has been extracted
#> 2026-11-21 21:00:29.776884 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:29.777818 INFO::Aligning assay columns with metadata rows
#> 2026-11-21 21:00:29.7787 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:29.779147 INFO::Total samples in data: 12
#> 2026-11-21 21:00:29.77951 INFO::Min samples required with min abundance for a feature not to be filtered: 1.200000
#> 2026-11-21 21:00:29.779913 INFO::Total filtered features with prevalence-abundance filtering: 0
#> 2026-11-21 21:00:29.780333 INFO::Filtered feature names:
#> 2026-11-21 21:00:29.780755 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:29.78119 INFO::Filtered feature names:
#> tweedie_p = 0 uses the Gaussian variance case; setting link = 'identity' because the requested link is undefined for negative responses.
#> tweedie_p = 0 with identity link does not use a log(scale_factor) offset; setting adjust_offset = FALSE.
#> 2026-11-21 21:00:29.782203 INFO::Offset adjustment disabled; feature data remain on the input scale.
#> 2026-11-21 21:00:29.782828 INFO::Excluded metadata with entropy less or equal to 0:
#> 2026-11-21 21:00:29.783314 INFO::Applying z-score to standardize continuous metadata
#> 2026-11-21 21:00:29.785256 INFO::Formula for fixed effects: expr ~  group
#> 2026-11-21 21:00:29.785877 INFO::Running selected analysis method: CPLM
#> [1] "Creating output feature tables folder"
#> [1] "Creating output fits folder"
#> 2026-11-21 21:00:30.374977 INFO::Writing function arguments to log file
#> 2026-11-21 21:00:30.378526 INFO::Verifying options selected are valid
#> 2026-11-21 21:00:30.405643 INFO::Determining format of input files
#> 2026-11-21 21:00:30.406245 INFO::Input format is data samples as rows and metadata samples as rows
#> 2026-11-21 21:00:30.407582 INFO::Formula for fixed effects: expr ~  group
#> 2026-11-21 21:00:30.408157 INFO::Factor detected for categorial metadata 'group'. Provide a reference argument or manually set factor ordering to change reference level.
#> 2026-11-21 21:00:30.40862 INFO::Filter data based on min abundance and min prevalence
#> 2026-11-21 21:00:30.409148 INFO::Total samples in data: 12
#> 2026-11-21 21:00:30.409608 INFO::Min samples required with min abundance for a feature not to be filtered: 1.200000
#> 2026-11-21 21:00:30.410319 INFO::Total filtered features: 0
#> 2026-11-21 21:00:30.410856 INFO::Filtered feature names from abundance and prevalence filtering:
#> 2026-11-21 21:00:30.411489 INFO::Total filtered features with variance filtering: 0
#> 2026-11-21 21:00:30.412011 INFO::Filtered feature names from variance filtering:
#> 2026-11-21 21:00:30.412461 INFO::Running selected normalization method: NONE
#> 2026-11-21 21:00:30.413018 INFO::Bypass z-score application to metadata
#> 2026-11-21 21:00:30.4135 INFO::Running selected transform method: NONE
#> 2026-11-21 21:00:30.414052 INFO::Running selected analysis method: LM
#> 2026-11-21 21:00:30.417465 INFO::Fitting model to feature number 1, feature1
#> 2026-11-21 21:00:30.420569 INFO::Fitting model to feature number 2, feature2
#> 2026-11-21 21:00:30.424859 INFO::Counting total values for each feature
#> 2026-11-21 21:00:30.425981 INFO::Writing filtered data to file /tmp/RtmpBS4ToU/maaslin2_55023642/features/filtered_data.tsv
#> 2026-11-21 21:00:30.427272 INFO::Writing filtered, normalized data to file /tmp/RtmpBS4ToU/maaslin2_55023642/features/filtered_data_norm.tsv
#> 2026-11-21 21:00:30.428197 INFO::Writing filtered, normalized, transformed data to file /tmp/RtmpBS4ToU/maaslin2_55023642/features/filtered_data_norm_transformed.tsv
#> 2026-11-21 21:00:30.429101 INFO::Writing residuals to file /tmp/RtmpBS4ToU/maaslin2_55023642/fits/residuals.rds
#> 2026-11-21 21:00:30.429712 INFO::Writing fitted values to file /tmp/RtmpBS4ToU/maaslin2_55023642/fits/fitted.rds
#> 2026-11-21 21:00:30.430305 INFO::Writing all results to file (ordered by increasing q-values): /tmp/RtmpBS4ToU/maaslin2_55023642/all_results.tsv
#> 2026-11-21 21:00:30.431352 INFO::Writing the significant results (those which are less than or equal to the threshold of 1.000000 ) to file (ordered by increasing q-values): /tmp/RtmpBS4ToU/maaslin2_55023642/significant_results.tsv
#> 2026-11-21 21:00:30.434696 INFO::Counting prevalence for each feature

# For Bioconductor-container workflow examples, see:
# vignette("Tweedieverse", package = "Tweedieverse")
```
