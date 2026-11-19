test_that("Tweedieverse rejects the legacy feature and metadata list input", {
  withr::local_seed(123)
  features <- matrix(stats::runif(800, -2, 5), ncol = 40)
  features[features < 0] <- 0
  rownames(features) <- paste0("feature", seq_len(nrow(features)))
  colnames(features) <- paste0("subject", seq_len(ncol(features)))
  metadata <- data.frame(
    grouping = sample(0:1, ncol(features), replace = TRUE),
    row.names = colnames(features)
  )
  expect_error(
    Tweedieverse(list(feature = features, group = metadata)),
    "input_features must be a supported Bioconductor container",
    fixed = TRUE
  )
})

make_test_se <- function(features, metadata) {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(features))),
    colData = metadata
  )
}

test_that("size-factor normalization does not modify feature data", {
  features <- data.frame(
    feature1 = c(10, 20, 30),
    feature2 = c(5, 0, 10),
    feature3 = c(1, 3, 9),
    row.names = paste0("sample", 1:3)
  )
  features_before <- features
  metadata <- data.frame(group = c("A", "B", "A"), row.names = rownames(features))

  sf <- compute_tweedieverse_size_factor(
    features = features,
    metadata = metadata,
    normalization = "TSS"
  )

  expect_equal(features, features_before)
  expect_equal(length(sf), nrow(features))
  expect_true(all(is.finite(sf)))
  expect_true(all(sf > 0))
})

test_that("domain defaults and normalization choices resolve as expected", {
  expect_equal(resolve_tweedieverse_normalization("microbiome", NULL), "TSS")
  expect_equal(resolve_tweedieverse_normalization("single_cell", NULL), "SCRAN")
  expect_equal(resolve_tweedieverse_normalization("bulk_rnaseq", NULL), "TMM")
  expect_equal(resolve_tweedieverse_normalization("custom", NULL), "MEDIAN")
  expect_equal(resolve_tweedieverse_normalization("bulk_rnaseq", "DESEQ2"), "RLE")
  expect_error(
    resolve_tweedieverse_normalization("single_cell", "TMM"),
    "not supported"
  )
})

test_that("Bioconductor container validation is domain specific", {
  map <- domain_bioc_container_map()
  expect_true("TreeSummarizedExperiment" %in% map$microbiome)
  expect_true("SingleCellExperiment" %in% map$single_cell)
  expect_true("SummarizedExperiment" %in% map$bulk_rnaseq)

  counts <- matrix(rpois(60, lambda = 5), nrow = 3, ncol = 20)
  rownames(counts) <- paste0("feature", seq_len(nrow(counts)))
  colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = data.frame(
      group = rep(c("A", "B"), each = 10),
      row.names = colnames(counts)
    )
  )

  expect_silent(validate_domain_bioc_container(se, "bulk_rnaseq"))
  expect_error(
    validate_domain_bioc_container(se, "single_cell"),
    "SingleCellExperiment"
  )
})

test_that("per-omics arguments resolve by experiment name", {
  omics_names <- c("microbiome", "rnaseq")

  expect_null(resolve_multiomics_arg(
    list(microbiome = NULL, rnaseq = 1.5),
    "microbiome",
    omics_names,
    "tweedie_p"
  ))
  expect_equal(resolve_multiomics_arg(
    c(microbiome = "TSS", rnaseq = "NONE"),
    "rnaseq",
    omics_names,
    "normalization"
  ), "NONE")
  expect_equal(resolve_multiomics_arg(
    list("microbiome", "bulk_rnaseq"),
    "rnaseq",
    omics_names,
    "domain"
  ), "bulk_rnaseq")
  expect_equal(resolve_multiomics_arg(
    c("age", "group"),
    "microbiome",
    omics_names,
    "fixed_effects"
  ), c("age", "group"))
})

test_that("MultiAssayExperiment recovers simulated effects and applies offsets", {
  skip_if_not_installed("MultiAssayExperiment")
  withr::local_seed(123)
  withr::local_options(lifecycle_verbosity = "quiet")
  samples <- paste0("sample", seq_len(80))
  group_b <- rep(c(FALSE, TRUE), each = 40)
  metadata <- data.frame(
    group = ifelse(group_b, "B", "A"),
    micro_size = ifelse(group_b, 2, 1),
    rna_size = ifelse(group_b, 3, 1),
    row.names = samples
  )
  micro_effects <- c(taxon1 = log(2), taxon2 = -log(2), taxon3 = 0)
  rna_effects <- c(gene1 = log(1.5), gene2 = -log(1.5), gene3 = 0, gene4 = log(3))
  micro_means <- 100 * exp(outer(micro_effects, as.numeric(group_b)))
  rna_means <- 100 * exp(outer(rna_effects, as.numeric(group_b)))
  micro_means <- sweep(micro_means, 2, metadata$micro_size, "*")
  rna_means <- sweep(rna_means, 2, metadata$rna_size, "*")
  microbiome_counts <- matrix(
    stats::rpois(length(micro_means), lambda = micro_means), nrow = 3,
    dimnames = list(names(micro_effects), samples)
  )
  rnaseq_counts <- matrix(
    stats::rpois(length(rna_means), lambda = rna_means), nrow = 4,
    dimnames = list(names(rna_effects), samples)
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(
      microbiome = SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = microbiome_counts), colData = metadata
      ),
      rnaseq = SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = rnaseq_counts), colData = metadata
      )
    ),
    colData = metadata
  )
  run_fit <- function(adjust_offset) Tweedieverse(
    input_features = mae,
    output = NULL,
    fixed_effects = "group",
    domain = c(microbiome = "microbiome", rnaseq = "bulk_rnaseq"),
    normalization = "NONE",
    scale_factor = c(microbiome = "micro_size", rnaseq = "rna_size"),
    tweedie_p = 1,
    adjust_offset = adjust_offset,
    standardize = FALSE,
    median_comparison = FALSE,
    return_all = TRUE
  )
  fit <- run_fit(TRUE)
  unadjusted <- run_fit(FALSE)
  expect_s3_class(fit, "TweedieverseMultiAssayResult")
  expect_equal(names(fit), c("microbiome", "rnaseq"))
  expected_effects <- list(microbiome = micro_effects, rnaseq = rna_effects)
  expected_depth_effects <- c(microbiome = log(2), rnaseq = log(3))
  for (omics in names(expected_effects)) {
    results <- fit[[omics]]
    effects <- expected_effects[[omics]]
    coefficients <- results$coef[match(names(effects), results$feature)]
    raw_coefficients <- unadjusted[[omics]]$coef[
      match(names(effects), unadjusted[[omics]]$feature)
    ]
    expect_equal(nrow(results), length(effects))
    expect_true(all(results$omics == omics))
    expect_true(all(results$metadata == "group" & results$value == "B"))
    expect_lt(max(abs(coefficients - effects)), 0.08)
    expect_equal(
      raw_coefficients - coefficients,
      rep(unname(expected_depth_effects[omics]), length(effects)),
      tolerance = 1e-6
    )
  }
})

test_that("scale_factor overrides domain normalization for offsets", {
  features <- data.frame(
    feature1 = c(10, 20, 30),
    feature2 = c(5, 0, 10),
    row.names = paste0("sample", 1:3)
  )
  metadata <- data.frame(
    group = c("A", "B", "A"),
    user_size = c(2, 4, 8),
    row.names = rownames(features)
  )

  expect_equal(
    resolve_tweedieverse_normalization("microbiome", "GMPR", "user_size"),
    "USER"
  )
  expect_equal(
    compute_tweedieverse_size_factor(features, metadata, "USER", "user_size"),
    metadata$user_size
  )
})

test_that("Tweedieverse keeps input data unchanged when normalization is selected", {
  set.seed(123)
  features <- as.data.frame(matrix(rpois(60, lambda = 5), nrow = 20, ncol = 3))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  features_before <- features
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    domain = "bulk_rnaseq",
    normalization = "CPM",
    tweedie_p = 1,
    median_comparison = FALSE,
    max_significance = 1
  )

  expect_equal(features, features_before)
  expect_true(nrow(fit) > 0)
})

test_that("BPPARAM controls BiocParallel execution", {
  set.seed(123)
  features <- as.data.frame(matrix(rpois(40, lambda = 5), nrow = 20, ncol = 2))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    domain = "bulk_rnaseq",
    normalization = "NONE",
    tweedie_p = 1,
    median_comparison = FALSE,
    max_significance = 1,
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_true(nrow(fit) > 0)
  expect_error(
    Tweedieverse(
      input_features = make_test_se(features, metadata),
      output = NULL,
      fixed_effects = "group",
      BPPARAM = list()
    ),
    "BiocParallelParam"
  )
})

test_that("vector fixed effects are parsed before metadata validation", {
  set.seed(123)
  features <- as.data.frame(matrix(rpois(80, lambda = 5), nrow = 20, ncol = 4))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    antibiotic = rep(c("no", "yes"), times = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = c("group", "antibiotic"),
    normalization = "NONE",
    tweedie_p = 1,
    median_comparison = FALSE,
    max_significance = 1
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(c("group", "antibiotic") %in% unique(fit$metadata)))
})

test_that("comma-separated fixed effects are rejected", {
  expect_error(
    parse_effect_names("group,antibiotic"),
    "character vectors"
  )
})

test_that("vector fixed effects are preserved", {
  expect_equal(
    parse_effect_names(c("group", "antibiotic")),
    c("group", "antibiotic")
  )
})

test_that("p equals 0 transformations return a copy and keep input data unchanged", {
  features <- data.frame(
    feature1 = c(1, 4, 9),
    feature2 = c(2, 0, 3),
    feature3 = c(4, 8, 0),
    row.names = paste0("sample", 1:3)
  )
  features_before <- features

  logged <- p0_transform_features(features, transform = "LOG", pseudocount = 1)
  clr <- p0_transform_features(features, transform = "CLR", pseudocount = 1)
  rclr <- p0_transform_features(features, transform = "RCLR")
  signed <- p0_transform_features(
    data.frame(feature1 = c(-4, 0, 9), row.names = paste0("sample", 1:3)),
    transform = "ARC_SIGNED_SQRT"
  )

  expect_equal(features, features_before)
  expect_equal(logged$feature1, log(features$feature1 + 1))
  expect_true(all(abs(rowMeans(clr)) < 1e-12))
  expect_equal(rclr$feature2[2], 0)
  expect_equal(signed$feature1, c(-2, 0, 3))
})

test_that("p equals 0 transformation is used for the analysis copy and written to output", {
  set.seed(123)
  features <- data.frame(
    feature1 = c(stats::rnorm(10, 1, 0.2), stats::rnorm(10, 3, 0.2)),
    feature2 = c(stats::rnorm(10, 2, 0.2), stats::rnorm(10, 4, 0.2))
  )
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  features_before <- features
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )
  output <- tempfile("tweedieverse_p0_transform_")

  suppressWarnings(
    fit <- Tweedieverse(
      input_features = make_test_se(features, metadata),
      output = output,
      fixed_effects = "group",
      tweedie_p = 0,
      Maaslin2_run = FALSE,
      p0_transform = "arc signed sqrt",
      median_comparison = FALSE,
      max_significance = 1
    )
  )

  transformed_file <- file.path(output, "transformed_features.tsv")
  transformed <- utils::read.delim(transformed_file, row.names = 1, check.names = FALSE)

  expect_equal(features, features_before)
  expect_true(file.exists(transformed_file))
  expect_equal(transformed$feature1, sign(features$feature1) * sqrt(abs(features$feature1)))
  expect_true(nrow(fit) > 0)
})

test_that("fixed Tweedie p values are supported for fixed-effect GLMs", {
  set.seed(123)
  features <- as.data.frame(matrix(rpois(60, lambda = 5), nrow = 20, ncol = 3))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    tweedie_p = 1,
    median_comparison = FALSE,
    max_significance = 1
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(fit$tweedie.index == 1))
  expect_true(all(fit$base.model == "Tweedie GLM"))
})

test_that("undefined Tweedie p values between 0 and 1 are rejected", {
  set.seed(123)
  features <- data.frame(
    feature1 = c(stats::rnorm(10, 1, 0.2), stats::rnorm(10, 3, 0.2)),
    feature2 = c(stats::rnorm(10, 2, 0.2), stats::rnorm(10, 4, 0.2))
  )
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  expect_error(
    Tweedieverse(
      input_features = make_test_se(features, metadata),
      output = NULL,
      fixed_effects = "group",
      tweedie_p = 0.5,
      median_comparison = FALSE
    ),
    "between 0 and 1"
  )
})

test_that("NA Tweedie p is treated as unspecified", {
  set.seed(123)
  features <- as.data.frame(matrix(rpois(60, lambda = 5), nrow = 20, ncol = 3))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    tweedie_p = NA,
    median_comparison = FALSE,
    max_significance = 1
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(fit$base.model == "CPLM"))
})

test_that("unspecified p with negative features resolves to p equals 0", {
  set.seed(123)
  features <- data.frame(
    feature1 = c(stats::rnorm(10, -1, 0.2), stats::rnorm(10, 1, 0.2)),
    feature2 = c(stats::rnorm(10, -2, 0.2), stats::rnorm(10, 2, 0.2))
  )
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  expect_message(
    fit <- Tweedieverse(
      input_features = make_test_se(features, metadata),
      output = NULL,
      fixed_effects = "group",
      tweedie_p = NULL,
      Maaslin2_run = FALSE,
      median_comparison = FALSE,
      max_significance = 1
    ),
    "Setting tweedie_p = 0"
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(fit$tweedie.index == 0))
  expect_true(all(fit$base.model == "Tweedie GLM"))
})

test_that("negative feature values are rejected for positive Tweedie powers", {
  features <- as.data.frame(matrix(rpois(40, lambda = 5), nrow = 20, ncol = 2))
  features[1, 1] <- -1
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  expect_error(
    Tweedieverse(
      input_features = make_test_se(features, metadata),
      output = NULL,
      fixed_effects = "group",
      tweedie_p = 1,
      median_comparison = FALSE
    ),
    "Negative feature values"
  )
})

test_that("negative supplied Tweedie p warns but runs for non-negative data", {
  set.seed(123)
  features <- data.frame(
    feature1 = c(stats::rnorm(10, 1, 0.2), stats::rnorm(10, 3, 0.2)),
    feature2 = c(stats::rnorm(10, 2, 0.2), stats::rnorm(10, 4, 0.2))
  )
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  expect_warning(
    fit <- Tweedieverse(
      input_features = make_test_se(features, metadata),
      output = NULL,
      fixed_effects = "group",
      tweedie_p = -1,
      link = "identity",
      median_comparison = FALSE,
      max_significance = 1
    ),
    "Negative Tweedie variance powers"
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(fit$tweedie.index == -1))
  expect_true(all(fit$base.model == "Tweedie GLM"))
})

test_that("presence-absence model adds individual and CCT results", {
  set.seed(123)
  features <- as.data.frame(matrix(1L + rpois(80, lambda = 5), nrow = 20, ncol = 4))
  features[1:3, 1] <- 0
  features[11:20, 2] <- features[11:20, 2] + 25L
  features[1:2, 3] <- 0
  features[3:20, 4] <- 0
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  run_fit <- function(...) Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    tweedie_p = 1,
    prev_threshold = 0,
    median_comparison = FALSE,
    run_presence_absence_model = TRUE,
    max_significance = 1,
    ...
  )
  fit <- run_fit()

  expect_true(all(c(
    "pval_abundance",
    "qval_abundance",
    "pval_presence",
    "qval_presence",
    "base.model_abundance",
    "base.model_presence"
  ) %in% colnames(fit)))
  expect_true(all(fit$base.model == "CCT"))
  expect_equal(fit$qval, sort(fit$qval))
  all_present <- fit[fit$feature == "feature2", ]
  expect_equal(nrow(all_present), 1L)
  expect_true(is.na(all_present$pval_presence))
  expect_equal(all_present$pval, all_present$pval_abundance)
  expect_lt(all_present$qval, 0.05)
  expect_false(is.na(fit$pval_presence[fit$feature == "feature1"]))
  rare <- fit[fit$feature %in% c("feature3", "feature4"), ]
  expect_equal(nrow(rare), 2L)
  expect_true(all(is.na(rare$pval_presence)))
  expect_equal(rare$pval, rare$pval_abundance)

  relaxed <- run_fit(min_presence_absence_count = 2L)
  expect_equal(sum(!is.na(relaxed$pval_presence[relaxed$feature %in% rare$feature])), 2L)
  expect_error(run_fit(min_presence_absence_count = 1.5), "must be a positive whole number")
})

test_that("p equals 0 uses the MaAsLin2 linear-model path", {
  skip_if_not_installed("Maaslin2")

  set.seed(123)
  features <- as.data.frame(matrix(rpois(60, lambda = 5), nrow = 20, ncol = 3))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    tweedie_p = 0,
    max_significance = 1
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(fit$tweedie.index == 0))
  expect_true(all(fit$base.model == "Maaslin2"))
})

test_that("presence-absence model can combine with the MaAsLin2 p equals 0 path", {
  skip_if_not_installed("Maaslin2")

  set.seed(123)
  features <- as.data.frame(matrix(rpois(80, lambda = 5), nrow = 20, ncol = 4))
  features[seq_len(10), 1] <- 0
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    tweedie_p = 0,
    run_presence_absence_model = TRUE,
    max_significance = 1
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(c(
    "pval_abundance",
    "qval_abundance",
    "pval_presence",
    "qval_presence",
    "base.model_abundance",
    "base.model_presence"
  ) %in% colnames(fit)))
  expect_true(all(fit$base.model == "CCT"))
  expect_true(all(fit$base.model_abundance == "Maaslin2"))
})

test_that("p equals 0 can run through the Tweedie GLM when Maaslin2_run is FALSE", {
  set.seed(123)
  features <- data.frame(
    feature1 = c(stats::rnorm(10, 1, 0.2), stats::rnorm(10, 3, 0.2)),
    feature2 = c(stats::rnorm(10, 2, 0.2), stats::rnorm(10, 4, 0.2))
  )
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    tweedie_p = 0,
    Maaslin2_run = FALSE,
    link = "identity",
    max_significance = 1
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(fit$tweedie.index == 0))
  expect_true(all(fit$base.model == "Tweedie GLM"))
})

test_that("MaAsLin2 method_args are accepted for the p equals 0 path", {
  skip_if_not_installed("Maaslin2")

  set.seed(123)
  features <- as.data.frame(matrix(rpois(60, lambda = 5), nrow = 20, ncol = 3))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    tweedie_p = 0,
    median_comparison = FALSE,
    max_significance = 1,
    method_args = list(
      Maaslin2 = list(
        normalization = "NONE",
        transform = "NONE"
      )
    )
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(fit$base.model == "Maaslin2"))
})

test_that("method.args is an alias but cannot be mixed with method_args", {
  skip_if_not_installed("Maaslin2")

  set.seed(123)
  features <- as.data.frame(matrix(rpois(60, lambda = 5), nrow = 20, ncol = 3))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = rownames(features)
  )

  fit <- Tweedieverse(
    input_features = make_test_se(features, metadata),
    output = NULL,
    fixed_effects = "group",
    tweedie_p = 0,
    median_comparison = FALSE,
    max_significance = 1,
    method.args = list(
      Maaslin2 = list(
        normalization = "NONE",
        transform = "NONE"
      )
    )
  )

  expect_true(nrow(fit) > 0)
  expect_error(
    Tweedieverse(
      input_features = make_test_se(features, metadata),
      output = NULL,
      fixed_effects = "group",
      tweedie_p = 0,
      method_args = list(Maaslin2 = list(transform = "NONE")),
      method.args = list(Maaslin2 = list(transform = "LOG"))
    ),
    "only one of method_args or method.args"
  )
})
