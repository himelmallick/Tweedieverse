median_seed_example <- function() data.frame(
  taxon = paste0("gene", 1:4),
  metadata = rep("group", 4),
  effect_size = c(0.5, 0.2, -0.4, 0.1),
  pval = c(0.01, 0.03, 0.02, 0.04),
  stderr = c(0.2, 0.1, 0.25, 0.15),
  qval = c(0.04, 0.04, 0.04, 0.04)
)

test_that("median simulation seeds reproduce results and restore RNG state", {
  withr::local_seed(615)
  previous_state <- .Random.seed
  first <- median_comparison_tweedie(median_seed_example(), n_sims = 200, seed = 17)
  expect_identical(.Random.seed, previous_state)
  stats::runif(3)
  second <- median_comparison_tweedie(median_seed_example(), n_sims = 200, seed = 17)
  expect_identical(first, second)
  third <- median_comparison_tweedie(median_seed_example(), n_sims = 200, seed = 18)
  expect_false(identical(first$pval_median, third$pval_median))
  expect_error(
    median_comparison_tweedie(median_seed_example(), seed = 1.5),
    "seed must be NULL or a single non-negative integer"
  )
})

test_that("a scoped median seed restores the absence of an RNG state", {
  withr::local_preserve_seed()
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  median_comparison_tweedie(median_seed_example(), n_sims = 100, seed = 17)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("model failures are attempted once and warn with the feature and reason", {
  attempts <- 0L
  local_mocked_bindings(
    glm = function(...) {
      attempts <<- attempts + 1L
      stop("deliberate fitting failure")
    },
    .package = "stats"
  )
  features <- data.frame(gene1 = c(2, 3, 4, 5))
  metadata <- data.frame(group = c("A", "A", "B", "B"))
  expect_warning(
    fit <- fit.CPLM(
      features, metadata, tweedie_p = 1,
      formula = expr ~ group, BPPARAM = BiocParallel::SerialParam()
    ),
    "Model fit failed for feature 'gene1': deliberate fitting failure",
    fixed = TRUE
  )
  expect_equal(attempts, 1L)
  expect_true(all(is.na(fit$results$pval)))
})

test_that("median seeds reach simulations through single and multi-assay entry points", {
  counts <- rbind(
    gene1 = c(2, 4, 3, 5, 8, 10, 9, 11),
    gene2 = c(5, 7, 6, 8, 11, 13, 12, 14)
  )
  colnames(counts) <- paste0("sample", 1:8)
  metadata <- data.frame(
    group = rep(c("A", "B"), each = 4), row.names = colnames(counts)
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = metadata
  )
  original_median <- median_comparison_tweedie
  observed_seed <- NULL
  local_mocked_bindings(
    median_comparison_tweedie = function(..., seed = NULL) {
      observed_seed <<- seed
      original_median(..., seed = seed)
    },
    .package = "Tweedieverse"
  )
  run_fit <- function(input, seed) Tweedieverse(
    input, fixed_effects = "group", normalization = "NONE", tweedie_p = 1,
    median_comparison = TRUE, median_comparison_seed = seed, return_all = TRUE
  )
  first <- run_fit(se, 17)
  expect_identical(observed_seed, 17)
  expect_identical(first, run_fit(se, 17))
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rna = se), colData = metadata
  )
  multi <- run_fit(mae, c(rna = 17))$rna
  expect_identical(observed_seed, 17)
  expect_equal(multi[, names(first)], first)
})
