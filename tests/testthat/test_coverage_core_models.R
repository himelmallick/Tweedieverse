make_core_features <- function(n = 18, p = 3, lambda = 6) {
  set.seed(321)
  features <- as.data.frame(matrix(
    stats::rpois(n * p, lambda = lambda) + 1,
    nrow = n,
    ncol = p
  ))
  colnames(features) <- paste0("feature", seq_len(ncol(features)))
  rownames(features) <- paste0("sample", seq_len(nrow(features)))
  features
}

make_core_metadata <- function(samples) {
  data.frame(
    group = rep(c("A", "B"), length.out = length(samples)),
    batch = rep(c("X", "Y", "Z"), length.out = length(samples)),
    age = seq_along(samples),
    site = rep(c("S1", "S2", "S3"), length.out = length(samples)),
    row.names = samples
  )
}

make_core_se <- function(features, metadata) {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(features))),
    colData = metadata
  )
}

make_row_oriented_se <- function(features) {
  feature_col_metadata <- data.frame(
    placeholder = seq_len(ncol(features)),
    row.names = colnames(features)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(features)),
    colData = feature_col_metadata
  )
}

test_that("Tweedieverse rejects unsupported input and metadata file paths", {
  features <- make_core_features()
  metadata <- make_core_metadata(rownames(features))
  se <- make_core_se(features, metadata)

  expect_error(
    Tweedieverse(input_features = data.frame(feature = 1)),
    "supported Bioconductor"
  )
  expect_error(
    Tweedieverse(input_features = se, input_metadata = tempfile()),
    "input_metadata file paths are not supported"
  )
})

test_that("Tweedieverse validates high-level options before fitting", {
  features <- make_core_features()
  metadata <- make_core_metadata(rownames(features))
  se <- make_core_se(features, metadata)

  expect_error(
    Tweedieverse(se, output = NULL, domain = "proteomics"),
    "Option not valid"
  )
  expect_error(
    Tweedieverse(se, output = NULL, link = "cloglog"),
    "Option not valid"
  )
  expect_error(
    Tweedieverse(se, output = NULL, base_model = "lm"),
    "Option not valid"
  )
  expect_error(
    Tweedieverse(se, output = NULL, p0_transform = "sqrt"),
    "Option not valid"
  )
  expect_error(
    Tweedieverse(se, output = NULL, p0_transform_pseudocount = -1),
    "non-negative finite numeric"
  )
  expect_error(
    Tweedieverse(se, output = NULL, correction = "fdr"),
    "Option not valid"
  )
  expect_error(
    Tweedieverse(se, output = NULL, optimizer = "optim"),
    "Option not valid"
  )
  expect_error(
    Tweedieverse(se, output = NULL, tweedie_p = c(1, 2)),
    "single finite numeric"
  )
  expect_error(
    Tweedieverse(se, output = NULL, run_presence_absence_model = NA),
    "TRUE or FALSE"
  )
  expect_error(
    Tweedieverse(se, output = NULL, Maaslin2_run = NA),
    "TRUE or FALSE"
  )
  expect_error(
    Tweedieverse(se, output = NULL, adjust_offset = NA),
    "TRUE or FALSE"
  )
  expect_error(
    Tweedieverse(se, output = NULL, method_args = "bad"),
    "named list"
  )
  expect_error(
    Tweedieverse(se, output = NULL, prev_threshold = 1.1),
    "outside \\[0, 1\\]"
  )
})

test_that("Tweedieverse handles all supported feature and metadata orientations", {
  features <- make_core_features(n = 16, p = 2)
  metadata <- make_core_metadata(rownames(features))

  row_row_fit <- Tweedieverse(
    input_features = make_row_oriented_se(features),
    input_metadata = metadata,
    output = NULL,
    fixed_effects = c("group"),
    domain = "custom",
    normalization = "NONE",
    tweedie_p = 1,
    max_significance = 1
  )
  column_column_fit <- Tweedieverse(
    input_features = make_core_se(features, metadata),
    input_metadata = as.data.frame(t(metadata)),
    output = NULL,
    fixed_effects = c("group"),
    domain = "custom",
    normalization = "NONE",
    tweedie_p = 1,
    max_significance = 1
  )
  row_column_fit <- Tweedieverse(
    input_features = make_row_oriented_se(features),
    input_metadata = as.data.frame(t(metadata)),
    output = NULL,
    fixed_effects = c("group"),
    domain = "custom",
    normalization = "NONE",
    tweedie_p = 1,
    max_significance = 1
  )

  expect_true(nrow(row_row_fit) > 0)
  expect_true(nrow(column_column_fit) > 0)
  expect_true(nrow(row_column_fit) > 0)
})

test_that("Tweedieverse handles output setup, references, and optional plots", {
  features <- make_core_features(n = 24, p = 1)
  metadata <- make_core_metadata(rownames(features))
  se <- make_core_se(features, metadata)
  output_file <- tempfile()
  writeLines("not a directory", output_file)

  expect_error(
    Tweedieverse(se, output = output_file, fixed_effects = "group"),
    "not a directory"
  )
  expect_error(
    Tweedieverse(
      se,
      output = NULL,
      fixed_effects = "site",
      domain = "custom",
      normalization = "NONE",
      tweedie_p = 1
    ),
    "Please provide the reference"
  )

  output_dir <- tempfile("tweedieverse_core_output_")
  dir.create(output_dir)
  writeLines("old log", file.path(output_dir, "Tweedieverse.log"))
  grDevices::pdf(file.path(output_dir, "coverage-device.pdf"))
  withr::defer(grDevices::dev.off())
  suppressWarnings(
    fit <- Tweedieverse(
      se,
      output = output_dir,
      fixed_effects = c("group", "site"),
      domain = "custom",
      normalization = "NONE",
      tweedie_p = 1,
      max_significance = 1,
      reference = "site,S2",
      plot_heatmap = TRUE,
      plot_scatter = TRUE,
      heatmap_first_n = 5
    )
  )

  expect_true(nrow(fit) > 0)
  expect_true(file.exists(file.path(output_dir, "all_results.tsv")))
  expect_true(file.exists(file.path(output_dir, "significant_results.tsv")))
  expect_true(dir.exists(file.path(output_dir, "figures")))
})

test_that("Tweedieverse drops missing effects and respects standardize false", {
  features <- make_core_features(n = 18, p = 2)
  metadata <- make_core_metadata(rownames(features))
  metadata$explicit <- factor(rep(c("low", "high"), length.out = nrow(metadata)))

  fit <- suppressWarnings(Tweedieverse(
    input_features = make_core_se(features, metadata),
    output = NULL,
    fixed_effects = "explicit",
    random_effects = "missing_random",
    domain = "custom",
    normalization = "NONE",
    adjust_offset = FALSE,
    standardize = FALSE,
    tweedie_p = 1,
    max_significance = 1
  ))

  expect_true(nrow(fit) > 0)
  expect_equal(unique(fit$metadata), "explicit")
})

test_that("fit.CPLM covers direct fixed-effect and error paths", {
  features <- make_core_features(n = 18, p = 2)
  metadata <- make_core_metadata(rownames(features))
  formula <- stats::as.formula("expr ~ group + age")

  cpglm_fit <- fit.CPLM(
    features = features,
    metadata = metadata,
    formula = formula,
    tweedie_p = NULL
  )
  glm_fit <- fit.CPLM(
    features = features,
    metadata = metadata,
    formula = formula,
    tweedie_p = 1
  )

  expect_true(nrow(cpglm_fit$results) > 0)
  expect_true(all(cpglm_fit$results$base.model == "CPLM"))
  expect_true(nrow(glm_fit$results) > 0)
  expect_true(all(glm_fit$results$base.model == "Tweedie GLM"))

  sqrt_fit <- fit.CPLM(
    features = features,
    metadata = metadata,
    link = "sqrt",
    tweedie_p = 1,
    formula = formula
  )
  inverse_fit <- fit.CPLM(
    features = features,
    metadata = metadata,
    link = "inverse",
    tweedie_p = 1,
    formula = formula
  )
  expect_true(nrow(sqrt_fit$results) > 0)
  expect_true(nrow(inverse_fit$results) > 0)

  bad_link <- suppressWarnings(
    fit.CPLM(
      features = features,
      metadata = metadata,
      link = "bad",
      tweedie_p = 1,
      formula = formula
    )
  )
  expect_true(all(is.na(bad_link$results$pval)))
})

test_that("fit.CPLM covers simple random-effect Tweedie GLMM paths", {
  skip_if_not_installed("glmmTMB")

  set.seed(456)
  features <- data.frame(
    feature1 = stats::rpois(20, lambda = 5) + 1,
    row.names = paste0("sample", seq_len(20))
  )
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 10),
    batch = rep(letters[1:5], each = 4),
    row.names = rownames(features)
  )

  fixed_p_fit <- fit.CPLM(
    features = features,
    metadata = metadata,
    formula = stats::as.formula("expr ~ group"),
    random_effects_formula = stats::as.formula("expr ~ (1 | batch)"),
    tweedie_p = 1
  )
  estimated_p_fit <- suppressWarnings(fit.CPLM(
    features = features,
    metadata = metadata,
    formula = stats::as.formula("expr ~ group"),
    random_effects_formula = stats::as.formula("expr ~ (1 | batch)"),
    tweedie_p = NULL
  ))
  fixed_between_one_and_two <- suppressWarnings(fit.CPLM(
    features = features,
    metadata = metadata,
    formula = stats::as.formula("expr ~ group"),
    random_effects_formula = stats::as.formula("expr ~ (1 | batch)"),
    tweedie_p = 1.5
  ))

  expect_true(nrow(fixed_p_fit$results) > 0)
  expect_true(all(fixed_p_fit$results$base.model == "Tweedie GLMM"))
  expect_true(nrow(estimated_p_fit$results) > 0)
  expect_true(nrow(fixed_between_one_and_two$results) > 0)
})

test_that("fit.CPLM reports failed feature fits and random-effect p restrictions", {
  features <- make_core_features(n = 18, p = 2)
  metadata <- make_core_metadata(rownames(features))
  bad_formula <- stats::as.formula("expr ~ missing_covariate")

  failed <- suppressWarnings(fit.CPLM(
    features = features,
    metadata = metadata,
    formula = bad_formula,
    tweedie_p = 1
  ))

  expect_true(all(is.na(failed$results$pval)))
  expect_error(
    fit.CPLM(
      features = features,
      metadata = metadata,
      formula = stats::as.formula("expr ~ group"),
      random_effects_formula = stats::as.formula("expr ~ (1 | batch)"),
      tweedie_p = 2.5
    ),
    "With random_effects"
  )
})
