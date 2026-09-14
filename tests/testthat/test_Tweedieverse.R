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
    input_features = features,
    input_metadata = metadata,
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
      input_features = features,
      input_metadata = metadata,
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
    input_features = features,
    input_metadata = metadata,
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
      input_features = features,
      input_metadata = metadata,
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
      input_features = features,
      input_metadata = metadata,
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
      input_features = features,
      input_metadata = metadata,
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
    input_features = features,
    input_metadata = metadata,
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
    input_features = features,
    input_metadata = metadata,
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
    input_features = features,
    input_metadata = metadata,
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
    input_features = features,
    input_metadata = metadata,
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
    input_features = features,
    input_metadata = metadata,
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
    input_features = features,
    input_metadata = metadata,
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
      input_features = features,
      input_metadata = metadata,
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
