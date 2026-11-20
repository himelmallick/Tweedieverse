logging_summary_se <- function() {
    counts <- rbind(
        gene1 = rep(c(4, 6, 5, 7, 10), 4),
        gene2 = rep(c(8, 6, 7, 9, 14), 4)
    )
    colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        colData = data.frame(
            group = rep(c("A", "B"), each = 10),
            row.names = colnames(counts)
        )
    )
}

test_that("output logging restores handlers and levels on success and error", {
    withr::local_options(lifecycle_verbosity = "quiet")
    logger <- logging::getLogger()
    original_handlers <- logger$handlers
    original_level <- logger$level
    output_root <- tempfile("tweedieverse_logging_")
    dir.create(output_root)
    withr::defer({
        logger$handlers <- original_handlers
        logger$level <- original_level
        unlink(output_root, recursive = TRUE)
    })
    sentinel <- file.path(output_root, "caller.log")
    logging::addHandler(
        "caller.file",
        action = logging::writeToFile, file = sentinel, level = "DEBUG"
    )
    logging::setLevel("WARN", logger)
    previous_handlers <- logger$handlers
    previous_level <- logger$level
    run_fit <- function(output, ...) {
        Tweedieverse(
            logging_summary_se(),
            output = output,
            fixed_effects = "group", normalization = "NONE",
            return_all = TRUE, ...
        )
    }
    first_dir <- file.path(output_root, "first")
    run_fit(first_dir, tweedie_p = 1)
    expect_identical(logger$handlers, previous_handlers)
    expect_identical(logger$level, previous_level)
    first_log <- readLines(file.path(first_dir, "Tweedieverse.log"))
    expect_true(any(grepl("Counting prevalence", first_log)))
    unlink(first_dir, recursive = TRUE)
    run_fit(file.path(output_root, "second"), tweedie_p = 1)
    expect_identical(logger$handlers, previous_handlers)
    expect_identical(logger$level, previous_level)
    expect_false(dir.exists(first_dir))
    expect_error(
        run_fit(file.path(output_root, "error"), tweedie_p = 1, optimizer =
            "bad"),
        "nlminb.*bobyqa.*L-BFGS-B"
    )
    expect_identical(logger$handlers, previous_handlers)
    expect_identical(logger$level, previous_level)

    maaslin_dir <- file.path(output_root, "maaslin")
    suppressMessages(run_fit(maaslin_dir, tweedie_p = 0))
    expect_identical(logger$handlers, previous_handlers)
    expect_identical(logger$level, previous_level)
    expect_true(any(grepl(
        "Counting prevalence", readLines(file.path(maaslin_dir,
            "Tweedieverse.log"))
    )))
    suppressMessages(run_fit(NULL, tweedie_p = 0))
    expect_identical(logger$handlers, previous_handlers)
    expect_identical(logger$level, previous_level)
    logging::logwarn("caller logging remains active")
    expect_true(any(grepl("caller logging remains active", readLines(
        sentinel))))
})

test_that(
    "feature summaries match by name with missing values and multiple coefficients", {
    features <- cbind(
        zero_heavy = c(0, 0, NA, 0, 1, 2, 0, 0, 3, 0, 0, 4),
        dense = 1:12,
        signed = c(-2, 0, 1, 3, NA, -1, 0, 2, 4, 0, 5, 6)
    )
    rownames(features) <- paste0("sample", seq_len(nrow(features)))
    metadata <- data.frame(
        group = rep(c("A", "B"), each = 6),
        age = rep(c(1, 3, 2), 4),
        row.names = rownames(features)
    )
    run_fit <- function(counts) {
        Tweedieverse(
            SummarizedExperiment::SummarizedExperiment(
                assays = list(counts = counts), colData = metadata
            ),
            fixed_effects = c("group", "age"), tweedie_p = 0,
            Maaslin2_run = FALSE, adjust_offset = FALSE,
            link = "identity", normalization = "NONE",
            return_all = TRUE, standardize = FALSE
        )
    }
    dense <- run_fit(t(features))
    sparse <- run_fit(Matrix::t(Matrix::Matrix(features, sparse = TRUE)))
    expect_equal(dense, sparse)
    expect_equal(nrow(dense), 6L)
    for (feature in colnames(features)) {
        values <- features[, feature]
        rows <- dense[dense$feature == feature, ]
        expect_equal(rows$N, rep(length(values), 2))
        expect_equal(rows$N.not.zero, rep(sum(values > 0, na.rm = TRUE), 2))
        expect_equal(
            rows$percent.zero, rep(round(mean(values == 0, na.rm = TRUE),
                2) * 100, 2)
        )
    }
})
