test_that("return_all exposes the full table independently of significance and output", {
  features <- data.frame(
    gene1 = c(4, 6, 5, 7, 10, 12, 11, 13),
    gene2 = c(8, 6, 7, 9, 14, 12, 13, 15),
    row.names = paste0("sample", 1:8)
  )
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 4),
    row.names = rownames(features)
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(features))), colData = metadata
  )
  run_fit <- function(return_all = FALSE, output = NULL) Tweedieverse(
    se, fixed_effects = "group", normalization = "NONE",
    tweedie_p = 1, max_significance = 0,
    return_all = return_all, output = output
  )
  expect_equal(nrow(run_fit()), 0L)
  full <- run_fit(TRUE)
  expect_s3_class(full, "data.frame")
  expect_equal(nrow(full), 2L)
  expect_true(all(full$qval > 0))
  expect_equal(full$qval, sort(full$qval))
  output_dir <- tempfile("tweedieverse_return_all_")
  on.exit({
    unlink(output_dir, recursive = TRUE)
  }, add = TRUE)
  expect_equal(run_fit(TRUE, output_dir), full)
  saved <- utils::read.delim(file.path(output_dir, "all_results.tsv"))
  significant <- utils::read.delim(file.path(output_dir, "significant_results.tsv"))
  expect_equal(saved$feature, full$feature)
  expect_equal(saved$pval, full$pval)
  expect_equal(nrow(significant), 0L)
  expect_error(run_fit(NA), "return_all must be TRUE or FALSE")
  expect_error(run_fit(c(TRUE, FALSE)), "return_all must be TRUE or FALSE")
})

test_that("MultiAssayExperiment retains empty result tables and other layers", {
  counts <- rbind(
    gene1 = c(4, 6, 5, 7, 10, 12, 11, 13),
    gene2 = c(8, 6, 7, 9, 14, 12, 13, 15)
  )
  colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
  metadata <- S4Vectors::DataFrame(
    group = rep(c("A", "B"), each = 4), row.names = colnames(counts)
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = metadata
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(empty = se, retained = se), colData = metadata
  )
  run_fit <- function(max_significance, return_all = FALSE) Tweedieverse(
    input_features = mae, output = NULL, fixed_effects = "group",
    domain = "bulk_rnaseq", normalization = "NONE", tweedie_p = 1,
    median_comparison = FALSE, max_significance = max_significance,
    return_all = return_all
  )

  full <- run_fit(0, return_all = TRUE)
  expect_true(all(full$empty$qval > 0))
  mixed <- run_fit(c(empty = 0, retained = 1))
  expect_s3_class(mixed, "TweedieverseMultiAssayResult")
  expect_named(mixed, c("empty", "retained"))
  expect_equal(nrow(mixed$empty), 0L)
  expect_identical(mixed$empty, full$empty[FALSE, , drop = FALSE])
  expect_identical(mixed$retained, full$retained)
  expect_identical(mixed$empty$omics, character())
  expect_identical(mixed$retained$omics, rep("retained", 2L))

  all_empty <- run_fit(0)
  expect_named(all_empty, names(full))
  for (omics in names(full)) {
    expect_identical(all_empty[[omics]], full[[omics]][FALSE, , drop = FALSE])
  }
})
