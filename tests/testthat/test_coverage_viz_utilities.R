test_that("heatmap helpers handle data frames, files, top-n filtering, and output files", {
  results <- data.frame(
    feature = c("tax1", "tax2", "tax1", "tax2", "tax3", "tax3"),
    metadata = c("group", "group", "time", "time", "group", "time"),
    value = c("A", "B", "early", "late", "A", "late"),
    pval = c(0.01, 0.04, 0.02, 0.03, 0.20, 0.25),
    qval = c(0.02, 0.08, 0.03, 0.04, 0.30, 0.35),
    coef = c(1.2, -0.7, 0.9, -0.4, 0.2, -0.1)
  )
  out_dir <- tempfile("coverage_heatmap_")
  dir.create(out_dir)
  grDevices::pdf(file.path(out_dir, "heatmap-device.pdf"))
  withr::defer(grDevices::dev.off())
  results_file <- file.path(out_dir, "results.tsv")
  utils::write.table(results, results_file, sep = "\t", quote = FALSE, row.names = FALSE)

  p_q <- Tweedieverse_heatmap(results, first_n = 2, write_to = out_dir)
  p_p <- Tweedieverse_heatmap(results_file, cell_value = "pval", first_n = NA)
  p_c <- Tweedieverse_heatmap(results, cell_value = "coef", first_n = 3)

  expect_s3_class(p_q, "pheatmap")
  expect_s3_class(p_p, "pheatmap")
  expect_s3_class(p_c, "pheatmap")
  expect_true(file.exists(file.path(out_dir, "gg_heatmap.RDS")))

  pdf_file <- file.path(out_dir, "heatmap.pdf")
  Tweedieverse:::save_heatmap(
    results_file = results_file,
    heatmap_file = pdf_file,
    figures_folder = out_dir,
    first_n = 2
  )

  expect_true(file.exists(pdf_file))
  expect_true(file.exists(file.path(out_dir, "heatmap.jpg")))
})

test_that("heatmap helper returns NULL for insufficient associations", {
  one_row <- data.frame(
    feature = "tax1",
    metadata = "group",
    value = "A",
    pval = 0.01,
    qval = 0.02,
    coef = 1
  )
  one_feature <- data.frame(
    feature = c("tax1", "tax1"),
    metadata = c("group", "time"),
    value = c("A", "early"),
    pval = c(0.01, 0.02),
    qval = c(0.02, 0.03),
    coef = c(1, 0.5)
  )
  one_metadata <- data.frame(
    feature = c("tax1", "tax2"),
    metadata = c("group", "group"),
    value = c("A", "A"),
    pval = c(0.01, 0.02),
    qval = c(0.02, 0.03),
    coef = c(1, 0.5)
  )

  expect_message(expect_null(Tweedieverse_heatmap(one_row)), "no associations")
  expect_message(expect_null(Tweedieverse_heatmap(one_feature)), "not enough features")
  expect_message(expect_null(Tweedieverse_heatmap(one_metadata)), "not enough metadata")
})

test_that("association and Tweedie-index plots write expected artifacts", {
  out_dir <- tempfile("coverage_assoc_")
  dir.create(out_dir)
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir)

  samples <- paste0("sample", seq_len(8))
  features <- data.frame(
    feat1 = seq(2, 9),
    feat2 = c(4, 5, 6, 7, 8, 9, 10, 11),
    row.names = samples
  )
  metadata <- data.frame(
    age = seq(20, 27),
    group = factor(rep(c("A", "B"), each = 4)),
    row.names = samples
  )
  output <- data.frame(
    feature = c("feat1", "feat2"),
    metadata = c("age", "group"),
    value = c("age", "B"),
    coef = c(0.5, -0.4),
    stderr = c(0.1, 0.2),
    pval = c(0.01, 0.03),
    qval = c(0.02, 0.04),
    tweedie.index = c(1.2, 1.6)
  )

  features_file <- file.path(out_dir, "features.tsv")
  metadata_file <- file.path(out_dir, "metadata.tsv")
  output_file <- file.path(out_dir, "output.tsv")
  utils::write.table(features, features_file, sep = "\t", quote = FALSE)
  utils::write.table(metadata, metadata_file, sep = "\t", quote = FALSE)
  utils::write.table(output, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  suppressWarnings(
    Tweedieverse:::association_plots(
      metadata = metadata_file,
      features = features_file,
      output_results = output_file,
      write_to = out_dir,
      figures_folder = fig_dir,
      max_jpgs = 2
    )
  )
  Tweedieverse:::tweedie_index_plot(output_file, figures_folder = fig_dir)

  expect_true(file.exists(file.path(out_dir, "age.pdf")))
  expect_true(file.exists(file.path(out_dir, "group.pdf")))
  expect_true(file.exists(file.path(fig_dir, "age_1.jpg")))
  expect_true(file.exists(file.path(fig_dir, "group_1.jpg")))
  expect_true(file.exists(file.path(fig_dir, "age_gg_associations.RDS")))
  expect_true(file.exists(file.path(fig_dir, "group_gg_associations.RDS")))
  expect_true(file.exists(file.path(fig_dir, "tweedie_index_plot.pdf")))
  expect_true(file.exists(file.path(fig_dir, "tweedie_index_plot.jpg")))
  expect_s3_class(readRDS(file.path(fig_dir, "gg_tweedie_index_plot.RDS")), "ggplot")
})

test_that("plot helpers handle empty inputs and direct drawing branches", {
  expect_invisible(Tweedieverse:::draw_tweedieverse_plot(NULL))

  out_dir <- tempfile("coverage_empty_plot_")
  dir.create(out_dir)
  empty_output <- data.frame(
    feature = character(),
    metadata = character(),
    value = character(),
    coef = numeric(),
    stderr = numeric(),
    pval = numeric(),
    qval = numeric(),
    tweedie.index = numeric()
  )
  features <- data.frame(feat1 = c(1, 2), row.names = c("s1", "s2"))
  metadata <- data.frame(group = c("A", "B"), row.names = c("s1", "s2"))

  expect_message(
    expect_null(Tweedieverse:::association_plots(metadata, features, empty_output)),
    "no associations"
  )
  expect_message(
    expect_null(Tweedieverse:::tweedie_index_plot(empty_output, figures_folder = out_dir)),
    "no associations"
  )
})

test_that("information-gain and method-argument helpers cover numeric and discrete branches", {
  df <- data.frame(
    feature_num = c(1, 2, 3, 4, NA, 6),
    feature_group = c("low", "low", "high", "high", "low", "high"),
    target = c("A", "A", "B", "B", "A", "B")
  )

  expect_gt(Tweedieverse:::IG_numeric(df, "feature_num", "target", bins = 2), 0)
  expect_gt(Tweedieverse:::IG_discrete(df, "feature_group", "target"), 0)
  expect_equal(Tweedieverse:::merge_method_args(list(a = 1), list(b = 2)), list(a = 1, b = 2))
  expect_error(Tweedieverse:::merge_method_args(list(a = 1), "bad"), "list")
  expect_equal(Tweedieverse:::extract_method_args(list(Maaslin2 = list(a = 1)), "Maaslin2"), list(a = 1))
  expect_equal(Tweedieverse:::extract_method_args(list(maaslin2 = list(a = 2)), "Maaslin2"), list(a = 2))
  expect_equal(Tweedieverse:::extract_method_args(NULL, "Maaslin2"), list())
  expect_error(Tweedieverse:::extract_method_args("bad", "Maaslin2"), "named list")
  expect_error(Tweedieverse:::extract_method_args(list(Maaslin2 = "bad"), "Maaslin2"), "must be a list")
})

test_that("assay and MultiAssay helpers cover success and failure branches", {
  skip_if_not_installed("MultiAssayExperiment")

  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("gene", seq_len(3))
  colnames(counts) <- paste0("sample", seq_len(4))
  metadata <- data.frame(group = rep(c("A", "B"), each = 2), row.names = colnames(counts))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = metadata
  )

  expect_message(assay_df <- Tweedieverse:::extractAssay(se, "counts"), "extracted")
  expect_equal(dim(assay_df), c(3, 4))
  expect_message(expect_null(Tweedieverse:::extractAssay(se, "missing")), "not found")
  expect_true(Tweedieverse:::is_multiassay_experiment(
    MultiAssayExperiment::MultiAssayExperiment(experiments = list(rna = se), colData = metadata)
  ))
  expect_false(Tweedieverse:::is_multiassay_experiment(se))
  expect_identical(Tweedieverse:::coerce_multiassay_experiment_input(se), se)
  expect_error(Tweedieverse:::coerce_multiassay_experiment_input(list()), "not supported")

  metadata_file <- tempfile(fileext = ".tsv")
  utils::write.table(metadata, metadata_file, sep = "\t", quote = FALSE)
  expect_equal(Tweedieverse:::extract_multiassay_metadata(se, metadata), metadata)
  expect_equal(Tweedieverse:::extract_multiassay_metadata(se, metadata_file), metadata)
})

test_that("size-factor helpers cover additional normalization and validation branches", {
  x <- matrix(
    c(10, 20, 30, 5, 15, 25, 2, 4, 8, 1, 3, 9),
    nrow = 4,
    byrow = TRUE
  )
  rownames(x) <- paste0("sample", seq_len(nrow(x)))
  colnames(x) <- paste0("feature", seq_len(ncol(x)))
  features <- as.data.frame(x)
  metadata <- data.frame(size = c(1, 2, 4, 8), row.names = rownames(features))

  expect_equal(Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "NONE"), rep(1, 4))
  expect_true(all(Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "MEDIAN") > 0))
  expect_true(all(Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "RLE") > 0))
  expect_true(all(Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "GMPR") > 0))
  expect_true(all(Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "CSS") > 0))
  expect_equal(
    Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "USER", "size"),
    metadata$size
  )
  expect_error(Tweedieverse:::validate_tweedieverse_size_factor(c(1, NA)), "finite")
  expect_error(Tweedieverse:::validate_tweedieverse_size_factor(c(1, 0)), "strictly positive")
  expect_error(Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "USER", "missing"), "scale_factor")
  expect_error(Tweedieverse:::compute_tweedieverse_size_factor(replace(features, 1, -1), metadata, "TSS"), "non-negative")
  expect_error(Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "BAD"), "Unsupported")

  skip_if_not_installed("edgeR")
  expect_true(all(Tweedieverse:::compute_tweedieverse_size_factor(features, metadata, "TMM") > 0))
})

test_that("median comparison and CCT helpers cover edge cases", {
  set.seed(1)
  df <- data.frame(
    taxon = paste0("tax", 1:4),
    metadata = c("group", "group", "time", "time"),
    effect_size = c(0.4, 0.2, 0.1, NA),
    pval = c(0.01, 0.2, 0.99, 0.99),
    stderr = c(0.1, 0.1, 0, 0.2),
    qval = c(0.02, 0.25, 1, 1)
  )

  med <- median_comparison_tweedie(df, p_cutoff = 0.95, subtract_median = TRUE, n_sims = 50)
  no_use <- median_comparison_tweedie(df, p_cutoff = 0, n_sims = 10)

  expect_true(all(c("pval_median", "coef_median") %in% colnames(med)))
  expect_equal(no_use$coef_median, df$effect_size)
  expect_equal(Tweedieverse:::tweedieverse_cct(c(0, 0.5)), 0)
  expect_equal(Tweedieverse:::tweedieverse_cct(c(1, NA)), 1)
  expect_error(Tweedieverse:::tweedieverse_cct(c(0, 1)), "exact 0")
  expect_error(Tweedieverse:::tweedieverse_cct(c(-0.1, 0.5)), "between 0 and 1")
  expect_error(Tweedieverse:::tweedieverse_cct(c(0.2, 0.3), weights = 1), "same length")
  expect_error(Tweedieverse:::tweedieverse_cct(c(0.2, 0.3), weights = c(1, -1)), "non-negative")
  expect_true(Tweedieverse:::tweedieverse_cct(c(1e-20, 0.2)) < 1e-10)
  expect_equal(Tweedieverse:::tweedieverse_cct_rows(matrix(c(0.1, 0.2), ncol = 1)), c(0.1, 0.2))
})

test_that("presence-absence helper functions augment, fit, and combine results", {
  metadata <- data.frame(
    group = factor(rep(c("A", "B"), each = 5)),
    offset = seq(1, 2, length.out = 10)
  )
  data <- data.frame(expr = rep(c(0, 1), 5), metadata)
  formula <- stats::as.formula("expr ~ group + offset(log(offset))")

  augmented <- Tweedieverse:::augment_presence_data(formula, data)
  expect_equal(nrow(augmented$data), 30)
  expect_equal(length(augmented$weights), 30)
  expect_error(Tweedieverse:::augment_presence_data(formula, data, weights = 1), "weights")
  expect_error(Tweedieverse:::fit_augmented_presence_model(formula, data, offset = 1), "offset")
  fit <- Tweedieverse:::fit_augmented_presence_model(formula, data, offset = log(metadata$offset))
  expect_s3_class(fit, "glm")

  features <- data.frame(
    feat1 = rep(c(0, 1), 5),
    feat2 = rep(1, 10)
  )
  presence <- Tweedieverse:::fit_presence_absence_model(
    features = features,
    metadata = metadata,
    formula = formula
  )
  abundance <- data.frame(
    feature = c("feat1", "feat2"),
    metadata = c("group", "group"),
    value = c("B", "B"),
    coef = c(0.5, 0.2),
    stderr = c(0.1, 0.2),
    pval = c(0.01, 0.2),
    qval = c(0.02, 0.2),
    base.model = c("Tweedie GLM", "Tweedie GLM"),
    tweedie.index = c(1, 1)
  )
  combined <- Tweedieverse:::combine_abundance_presence_results(abundance, presence)

  expect_true(all(c("pval_presence", "base.model_presence", "pval_abundance") %in% colnames(combined)))
  expect_true(all(combined$base.model == "CCT"))
})
