test_that(
    "direct entropy calculation preserves category and missing-value conventions", {
    expect_equal(tweedie_entropy(c("A", "A", "B", "B")), 1)
    expect_equal(tweedie_entropy(c(1, 1, 1, 2)), 0.8112781244591328)
    expect_equal(tweedie_entropy(rep("A", 4)), 0)
    expect_equal(
        tweedie_entropy(factor(c("A", "B", "A", "B"), levels = c("A", "B",
            "C"))),
        1
    )
    expect_equal(tweedie_entropy(c("A", "A", "B", NA)), 1)
    expect_equal(tweedie_entropy(c("A", NA)), 0.5)
    expect_equal(tweedie_entropy(rep(NA_real_, 4)), 0)
    expect_equal(tweedie_entropy(character()), 0)
    expect_equal(tweedie_entropy(c(0.1, 0.1 + .Machine$double.eps / 4)), 0)
})

test_that(
    "metadata entropy thresholds retain the same eligible model variables", {
    counts <- rbind(
        gene1 = c(2, 4, 3, 5, 8, 10, 9, 11),
        gene2 = c(5, 7, 6, 8, 11, 13, 12, 14)
    )
    colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
    metadata <- data.frame(
        group = rep(c("A", "B"), each = 4),
        age = seq_len(ncol(counts)),
        constant = 1,
        missing = NA_real_,
        row.names = colnames(counts)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts), colData = metadata
    )
    run_fit <- function(threshold) {
        Tweedieverse(
            se,
            fixed_effects = c("group", "age"),
            entropy_threshold = threshold,
            normalization = "NONE", tweedie_p = 1, return_all = TRUE
        )
    }
    expect_setequal(unique(run_fit(0)$metadata), c("group", "age"))
    expect_equal(unique(run_fit(1)$metadata), "age")
    expect_error(run_fit(4), "No fixed effects remain")
})
