review_test_se <- function(features, metadata) {
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = t(as.matrix(features))),
        colData = metadata
    )
}

test_that("significant results exclude missing median-comparison q-values", {
    features <- data.frame(
        gene1 = c(4, 6, 5, 7, 10, 12, 11, 13),
        gene2 = c(8, 6, 7, 9, 14, 12, 13, 15),
        row.names = paste0("sample", 1:8)
    )
    metadata <- data.frame(
        group = rep(c("A", "B"), each = 4),
        row.names = rownames(features)
    )
    all_missing <- FALSE
    local_mocked_bindings(
        median_comparison_tweedie = function(df, ...) {
            df$coef_median <- df$effect_size
            df$pval_median <- rep(NA_real_, nrow(df))
            if (!all_missing) {
                df$pval_median[1L] <- 0.01
            }
            df
        },
        .package = "Tweedieverse"
    )
    run_fit <- function(return_all = FALSE) {
        Tweedieverse(
            review_test_se(features, metadata),
            fixed_effects = "group",
            normalization = "NONE",
            tweedie_p = 1,
            median_comparison = TRUE,
            return_all = return_all
        )
    }
    fit <- run_fit()
    expect_equal(nrow(fit), 1L)
    expect_false(anyNA(fit$feature))
    expect_false(anyNA(fit$qval))
    full <- run_fit(TRUE)
    expect_equal(nrow(full), 2L)
    expect_equal(sum(is.na(full$qval)), 1L)
    expect_false(anyNA(full$feature))
    all_missing <- TRUE
    expect_equal(nrow(run_fit()), 0L)
    expect_equal(nrow(run_fit(TRUE)), 2L)
})

test_that("MultiAssayExperiment maps assay samples to primary metadata", {
    previous_options <- options(lifecycle_verbosity = "quiet")
    on.exit(options(previous_options), add = TRUE)
    primary <- paste0("person", 1:20)
    metadata <- data.frame(
        group = rep(c("A", "B"), each = 10),
        size = rep(c(1, 3), each = 10),
        row.names = primary
    )
    primary_indices <- c(20:1, 1L)
    assay_samples <- paste0("aliquot", seq_along(primary_indices))
    assay_metadata <- metadata[primary_indices, , drop = FALSE]
    rownames(assay_metadata) <- assay_samples
    features <- data.frame(
        gene1 = (10 + seq_along(primary_indices) %% 3) *
            assay_metadata$size * ifelse(assay_metadata$group == "B", 2, 1),
        gene2 = (15 + seq_along(primary_indices) %% 4) * assay_metadata$size,
        row.names = assay_samples
    )
    experiment <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = t(as.matrix(features)))
    )
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(rna = experiment),
        colData = metadata,
        sampleMap = S4Vectors::DataFrame(
            assay = factor(rep("rna", length(assay_samples))),
            primary = primary[primary_indices],
            colname = assay_samples
        )
    )
    args <- list(
        fixed_effects = "group", domain = "bulk_rnaseq",
        normalization = "NONE", scale_factor = "size",
        tweedie_p = 1, max_significance = 1, standardize = FALSE
    )
    expected <- do.call(Tweedieverse, c(
        list(input_features = review_test_se(features, assay_metadata)), args
    ))
    actual <- do.call(Tweedieverse, c(list(input_features = mae), args))$rna
    expect_equal(actual[, names(expected)], expected)
    expect_equal(actual$N, rep(length(assay_samples), 2))
    full <- do.call(Tweedieverse, c(
        list(input_features = mae),
        utils::modifyList(args, list(return_all = TRUE, max_significance = 0))
    ))$rna
    expect_equal(full, actual)

    override <- metadata[20:1, , drop = FALSE]
    overridden <- do.call(Tweedieverse, c(
        list(input_features = mae, input_metadata = override), args
    ))$rna
    expect_equal(overridden, actual)

    subset_metadata <- metadata[20:3, , drop = FALSE]
    keep <- primary[primary_indices] %in% rownames(subset_metadata)
    expected_subset <- do.call(Tweedieverse, c(
        list(input_features = review_test_se(
            features[keep, , drop = FALSE], assay_metadata[keep, , drop = FALSE]
        )), args
    ))
    actual_subset <- do.call(Tweedieverse, c(
        list(input_features = mae, input_metadata = subset_metadata), args
    ))$rna
    expect_equal(actual_subset[, names(expected_subset)], expected_subset)
})

test_that("abundance offsets enter fixed and estimated Tweedie fits", {
    set.seed(410)
    n <- 60L
    metadata <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)))
    metadata$offset <- exp(0.8 * (metadata$group == "B") + seq(-0.2, 0.2,
        length.out = n))
    data <- data.frame(
        expr = stats::rgamma(
            n,
            shape = 5,
            scale = exp(2 + 0.5 * (metadata$group == "B")) * metadata$offset / 5
        ),
        metadata
    )
    features <- data.frame(gene = data$expr)
    for (power in list(NULL, 1, 1.5, 2, 3)) {
        actual <- fit.Tweedieverse(
            features, metadata,
            formula = expr ~ group, tweedie_p = power
        )$results
        if (is.null(power)) {
            expected <- cplm::cpglm(
                expr ~ group + offset(log(offset)),
                data = data, link = "log"
            )
        } else {
            expected <- stats::glm(
                expr ~ group + offset(log(offset)),
                data = data,
                family = statmod::tweedie(var.power = power, link.power = 0)
            )
        }
        expect_equal(actual$coef, unname(stats::coef(expected)["groupB"]),
            tolerance = 1e-6)
    }
})

test_that("mixed abundance models preserve the log size-factor offset", {
    set.seed(411)
    n <- 60L
    metadata <- data.frame(
        group = factor(rep(rep(c("A", "B"), each = 3), 10)),
        subject = factor(rep(1:10, each = 6))
    )
    metadata$offset <- exp(0.8 * (metadata$group == "B") + rep(c(-0.2, 0, 0.2),
        20))
    subject_effect <- stats::rnorm(10, sd = 0.4)[metadata$subject]
    data <- data.frame(
        expr = stats::rpois(
            n, exp(3 + 0.5 * (metadata$group == "B") +
                subject_effect) * metadata$offset
        ),
        metadata
    )
    actual <- fit.Tweedieverse(
        data.frame(gene = data$expr), metadata,
        formula = expr ~ group, random_effects_formula = expr ~ (1 | subject),
        tweedie_p = 1
    )$results
    expected <- glmmTMB::glmmTMB(
        expr ~ group + offset(log(offset)) + (1 | subject),
        data = data, family = stats::poisson(link = "log")
    )
    expect_equal(
        actual$coef, unname(glmmTMB::fixef(expected)$cond["groupB"]),
        tolerance = 1e-6
    )
})

test_that("adjust_offset changes coefficients through the public entry point", {
    previous_options <- options(lifecycle_verbosity = "quiet")
    on.exit(options(previous_options), add = TRUE)
    features <- data.frame(
        gene = rep(c(9, 10, 11, 12), 5) * rep(c(1, 3), each = 10),
        row.names = paste0("sample", 1:20)
    )
    metadata <- data.frame(
        group = rep(c("A", "B"), each = 10),
        size = rep(c(1, 3), each = 10),
        row.names = rownames(features)
    )
    run_fit <- function(adjust_offset) {
        Tweedieverse(
            review_test_se(features, metadata),
            fixed_effects = "group",
            scale_factor = "size", adjust_offset = adjust_offset,
            tweedie_p = 1, max_significance = 1, standardize = FALSE
        )
    }
    adjusted <- run_fit(TRUE)
    unadjusted <- run_fit(FALSE)
    expect_equal(unadjusted$coef - adjusted$coef, log(3), tolerance = 1e-6)
})
