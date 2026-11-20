test_that(
    "container alignment uses sample names even when feature names overlap", {
    counts <- rbind(
        sample1 = c(2, 4, 3, 5, 8, 10, 9, 11),
        gene2 = c(5, 7, 6, 8, 11, 13, 12, 14)
    )
    colnames(counts) <- paste0("sample", 1:8)
    metadata <- data.frame(
        group = rep(c("A", "B"), each = 4),
        row.names = colnames(counts)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts), colData = metadata
    )
    run_fit <- function(metadata) {
        Tweedieverse(
            se,
            input_metadata = metadata, fixed_effects = "group",
            normalization = "NONE", tweedie_p = 1, max_significance = 1
        )
    }
    expected <- stats::glm(
        as.numeric(counts[1, ]) ~ group,
        data = metadata,
        family = statmod::tweedie(var.power = 1, link.power = 0)
    )
    actual <- run_fit(metadata[8:1, , drop = FALSE])
    expect_equal(
        actual$coef[actual$feature == "sample1"],
        unname(stats::coef(expected)["groupB"])
    )
    expect_equal(actual$N, c(8L, 8L))
    expect_equal(actual, run_fit(metadata))

    rownames(metadata) <- paste0("unknown", 1:8)
    expect_error(run_fit(metadata),
        "Assay column names must match sample row names")
    expect_error(
        Tweedieverse(
            se,
            fixed_effects = "constant",
            input_metadata = data.frame(
                constant = rep(1, ncol(counts)), row.names = colnames(counts)
            )
        ),
        "No fixed effects remain"
    )
    expect_error(
        Tweedieverse(se, assay_name = "missing"),
        "Assay 'missing' was not found. Available assays: counts",
        fixed = TRUE
    )
})

test_that("sparse assays remain sparse through normalization and fitting", {
    set.seed(512)
    counts <- matrix(stats::rpois(120, lambda = 4), nrow = 6)
    counts[matrix(stats::runif(120) < 0.45, nrow = 6)] <- 0
    rownames(counts) <- paste0("gene", 1:6)
    colnames(counts) <- paste0("sample", 1:20)
    metadata <- data.frame(
        group = rep(c("A", "B"), each = 10), row.names = colnames(counts)
    )
    sparse_counts <- Matrix::Matrix(counts, sparse = TRUE)
    dense_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts), colData = metadata
    )
    sparse_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = sparse_counts), colData = metadata
    )
    run_fit <- function(se) {
        Tweedieverse(
            se,
            fixed_effects = "group", normalization = "TSS",
            tweedie_p = 1, run_presence_absence_model = TRUE,
            max_significance = 1
        )
    }
    expected <- run_fit(dense_se)
    original_fit <- Tweedieverse:::fit.Tweedieverse
    local_mocked_bindings(
        fit.Tweedieverse = function(features, ...) {
            expect_s4_class(features, "sparseMatrix")
            original_fit(features = features, ...)
        },
        .package = "Tweedieverse"
    )
    expect_equal(run_fit(sparse_se), expected)
    expect_identical(SummarizedExperiment::assay(sparse_se), sparse_counts)
})

test_that("dense and sparse assays give equivalent size factors", {
    set.seed(513)
    counts <- matrix(stats::rpois(1200, lambda = 5), nrow = 20)
    counts[matrix(stats::runif(1200) < 0.3, nrow = 20)] <- 0
    sparse_counts <- Matrix::Matrix(counts, sparse = TRUE)
    metadata <- data.frame(group = rep(c("A", "B"), each = 10))
    for (method in c("NONE", "TSS", "CPM", "MEDIAN", "RLE", "GMPR", "CSS",
        "TMM", "SCRAN")) {
        expected <- compute_tweedieverse_size_factor(counts, metadata, method)
        actual <- compute_tweedieverse_size_factor(sparse_counts, metadata,
            method)
        expect_equal(actual, expected, tolerance = 1e-6, info = method)
    }
})
