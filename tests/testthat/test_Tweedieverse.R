# Generate a matrix with feature and group data (from SDA package)
featureInfo <- matrix(runif(800, -2, 5), ncol = 40)
featureInfo[featureInfo<0] <- 0
feat.names <- paste("feature", 1:20, sep = '')
rownames(featureInfo) <- feat.names
sub.names <- paste('subject', 1:40, sep = '')
colnames(featureInfo) <- sub.names
groupInfo <- data.frame(grouping=matrix(sample(0:1, 40, replace = TRUE),
                                        ncol = 1))
rownames(groupInfo) <- colnames(featureInfo)

SEdata <- list(feature = featureInfo, group = groupInfo)

expect_error(Tweedieverse(SEdata))

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

test_that("MultiAssayExperiment input returns omics-specific results", {
  skip_if_not_installed("MultiAssayExperiment")

  set.seed(123)
  samples <- paste0("sample", seq_len(20))
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    row.names = samples
  )

  microbiome_counts <- matrix(rpois(60, lambda = 5) + 1, nrow = 3, ncol = 20)
  rnaseq_counts <- matrix(rpois(80, lambda = 8) + 1, nrow = 4, ncol = 20)
  rownames(microbiome_counts) <- paste0("taxon", seq_len(nrow(microbiome_counts)))
  rownames(rnaseq_counts) <- paste0("gene", seq_len(nrow(rnaseq_counts)))
  colnames(microbiome_counts) <- samples
  colnames(rnaseq_counts) <- samples

  microbiome <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = microbiome_counts),
    colData = metadata
  )
  rnaseq <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = rnaseq_counts),
    colData = metadata
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(microbiome = microbiome, rnaseq = rnaseq),
    colData = metadata
  )

  fit <- Tweedieverse(
    input_features = mae,
    output = NULL,
    fixed_effects = "group",
    domain = c(microbiome = "microbiome", rnaseq = "bulk_rnaseq"),
    normalization = c(microbiome = "NONE", rnaseq = "NONE"),
    tweedie_p = c(microbiome = 1, rnaseq = 1),
    median_comparison = FALSE,
    max_significance = 1,
    cores = 1
  )

  expect_s3_class(fit, "TweedieverseMultiAssayResult")
  expect_equal(names(fit), c("microbiome", "rnaseq"))
  expect_true(all(vapply(fit, nrow, integer(1)) > 0))
  expect_true(all(fit$microbiome$omics == "microbiome"))
  expect_true(all(fit$rnaseq$omics == "rnaseq"))
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
    max_significance = 1,
    cores = 1
  )

  expect_equal(features, features_before)
  expect_true(nrow(fit) > 0)
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
      max_significance = 1,
      cores = 1
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
    max_significance = 1,
    cores = 1
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
      median_comparison = FALSE,
      cores = 1
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
    max_significance = 1,
    cores = 1
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
      max_significance = 1,
      cores = 1
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
      median_comparison = FALSE,
      cores = 1
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
      max_significance = 1,
      cores = 1
    ),
    "Negative Tweedie variance powers"
  )

  expect_true(nrow(fit) > 0)
  expect_true(all(fit$tweedie.index == -1))
  expect_true(all(fit$base.model == "Tweedie GLM"))
})

test_that("presence-absence model adds individual and CCT results", {
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
    tweedie_p = 1,
    median_comparison = FALSE,
    run_presence_absence_model = TRUE,
    max_significance = 1,
    cores = 1
  )

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
    max_significance = 1,
    cores = 1
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
    max_significance = 1,
    cores = 1
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
    max_significance = 1,
    cores = 1
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
    ),
    cores = 1
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
    ),
    cores = 1
  )

  expect_true(nrow(fit) > 0)
  expect_error(
    Tweedieverse(
      input_features = make_test_se(features, metadata),
      output = NULL,
      fixed_effects = "group",
      tweedie_p = 0,
      method_args = list(Maaslin2 = list(transform = "NONE")),
      method.args = list(Maaslin2 = list(transform = "LOG")),
      cores = 1
    ),
    "only one of method_args or method.args"
  )
})
