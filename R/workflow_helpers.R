# Internal stages of the single-assay analysis workflow.
# Options and intermediate data are passed explicitly between stages.

validate_tweedieverse_power <- function(tweedie_p) {
    if (!is.null(tweedie_p)) {
        if (length(tweedie_p) != 1L) {
            stop("tweedie_p must be NULL, NA, or a single ",
                "finite numeric value.")

        }
        if (is.na(tweedie_p)) {
            tweedie_p <- NULL
        } else if (!is.numeric(tweedie_p)) {
            stop("tweedie_p must be NULL, NA, or a single ",
                "finite numeric value.")

        } else if (!is.finite(tweedie_p)) {
            stop("tweedie_p must be NULL, NA, or a single ",
                "finite numeric value.")

        }
    }
    if (!is.null(tweedie_p)) {
        if (tweedie_p > 0 && tweedie_p < 1) {
            stop("Tweedie variance powers between 0 and 1 ",
                "are undefined and are not supported.")

        }
    }
    tweedie_p
}

validate_tweedieverse_transform <- function(options) {
    options$p0_transform <- gsub("[ -]", "_", toupper(options$p0_transform))
    if (options$p0_transform %in% c(
        "ARCSIGNEDSQRT", "ARC_SIGNED_SQUARE_ROOT",
        "SIGNED_SQRT"
    )) {
        options$p0_transform <- "ARC_SIGNED_SQRT"
    }
    if (!options$p0_transform %in% options$p0_transform_choices) {
        option_not_valid_error(
            "Please select a p0_transform from the list of available options",
            toString(options$p0_transform_choices)
        )
    }
    if (length(options$p0_transform_pseudocount) != 1L || !is.numeric(
        options$p0_transform_pseudocount) ||
        !is.finite(options$p0_transform_pseudocount) ||
            options$p0_transform_pseudocount <
        0) {
        stop("p0_transform_pseudocount must be a ",
            "single non-negative finite numeric ",
            "value.")
    }
    options$p0_transform
}

validate_tweedieverse_flags <- function(options) {
    if (!is.logical(options$run_presence_absence_model) || length(
        options$run_presence_absence_model) !=
        1L || is.na(options$run_presence_absence_model)) {
        stop("run_presence_absence_model must be TRUE or FALSE.")
    }
    if (!is.logical(options$return_all) || length(options$return_all) != 1L ||
        is.na(options$return_all)) {
        stop("return_all must be TRUE or FALSE.")
    }
    if (!is.logical(options$Maaslin2_run) || length(options$Maaslin2_run) !=
        1L || is.na(options$Maaslin2_run)) {
        stop("Maaslin2_run must be TRUE or FALSE.")
    }
    if (!is.logical(options$adjust_offset) || length(options$adjust_offset) !=
        1L || is.na(options$adjust_offset)) {
        stop("adjust_offset must be TRUE or FALSE.")
    }
}

validate_tweedieverse_thresholds <- function(options) {
    if (!is.numeric(options$min_presence_absence_count) || length(
        options$min_presence_absence_count) !=
        1L || !is.finite(options$min_presence_absence_count) ||
            options$min_presence_absence_count <
        1 || options$min_presence_absence_count != floor(
            options$min_presence_absence_count)) {
        stop("min_presence_absence_count must be a positive whole number.")
    }
    prop_options <- c(options$prev_threshold, options$max_significance)
    if (any(prop_options < 0) || any(prop_options > 1)) {
        stop("One of the following is outside [0, 1]: ",
            "prev_threshold, max_significance")

    }
}

resolve_tweedieverse_method_args <- function(options) {
    if (!is.null(options$method_args) && !is.null(options$method.args)) {
        stop("Please provide only one of method_args or method.args.")
    }
    if (is.null(options$method_args)) {
        options$method_args <- options$method.args
    }
    if (!is.null(options$method_args) && !is.list(options$method_args)) {
        stop("method_args must be NULL or a named list.")
    }
    options$method_args
}

validate_tweedieverse_options <- function(options) {
    model_choices <- c("CPLM")
    link_choices <- c("log", "identity", "sqrt", "inverse")
    p0_transform_choices <- c("NONE", "CLR", "RCLR", "LOG", "ARC_SIGNED_SQRT")
    domain_choices <- c("microbiome", "single_cell", "bulk_rnaseq", "custom")
    correction_choices <- c(
        "BH", "holm", "hochberg", "hommel", "bonferroni",
        "BY"
    )
    optimizer_choices <- c("nlminb", "bobyqa", "L-BFGS-B")
    options$domain <- match.arg(tolower(options$domain), domain_choices)
    options$link <- match.arg(options$link, link_choices)
    options$base_model <- match.arg(options$base_model, model_choices)
    options$correction <- match.arg(options$correction, correction_choices)
    options$optimizer <- match.arg(options$optimizer, optimizer_choices)
    options$p0_transform_choices <- p0_transform_choices
    options$p0_transform <- validate_tweedieverse_transform(options)
    options$p0_transform_choices <- NULL
    options$normalization <- resolve_tweedieverse_normalization(
        options$domain,
        options$normalization, options$scale_factor
    )
    options["tweedie_p"] <- list(validate_tweedieverse_power(options$tweedie_p))
    validate_tweedieverse_flags(options)
    validate_tweedieverse_thresholds(options)
    options["method_args"] <- list(resolve_tweedieverse_method_args(options))
    options$no_output <- is.null(options$output)
    options
}

extract_tweedieverse_input <- function(options) {
    valid_classes <- unique(unlist(domain_bioc_container_map()))
    is_supported_bioc_class <- inherits(options$input_features, valid_classes)
    if (is_supported_bioc_class) {
        validate_domain_bioc_container(options$input_features, options$domain)
        data <- extractAssay(options$input_features, options$assay_name)
        if (is.null(options$input_metadata)) {
            metadata <- data.frame(SummarizedExperiment::colData(
                options$input_features))
        } else {
            metadata <- options$input_metadata
        }
    } else {
        {
            condition_message1 <- paste0(
                "Input data of class <%s> not supported. ",
                "Please use a domain-appropriate ", "Bioconductor container.")
            stop(sprintf(condition_message1, class(options$input_features)[1]))
        }
    }
    list(data = data, metadata = metadata)
}

prepare_tweedieverse_output <- function(output) {
    if (is.null(output)) {
        return()
    }
    if (file.exists(output) && !dir.exists(output)) {
        stop(sprintf("Output path exists but is not a directory: %s",
            output))
    }
    if (!dir.exists(output)) {
        message(sprintf("Creating output folder: %s", output))
        dir.create(output, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(output)) {
        stop(sprintf("Unable to create output folder: %s", output))
    }
    figures_folder <- file.path(output, "figures")
    if (!dir.exists(figures_folder)) {
        message(sprintf("Creating output figures folder: %s",
            figures_folder))
        dir.create(figures_folder, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(figures_folder)) {
        stop(sprintf("Unable to create output figures folder: %s",
            figures_folder))
    }
    figures_folder
}

start_tweedieverse_logging <- function(output, .local_envir = parent.frame()) {
    if (is.null(output)) {
        return(invisible())
    }
    log_file <- file.path(output, "Tweedieverse.log")
    if (file.exists(log_file)) {
        warning(sprintf("Deleting existing log file: %s", log_file), call. =
            FALSE)
        unlink(log_file)
    }
    logger <- logging::getLogger()
    previous_handlers <- logger$handlers
    previous_level <- logger$level
    withr::defer(
        {
            logger$handlers <- previous_handlers
            logger$level <- previous_level
        },
        envir = .local_envir
    )
    logging::basicConfig(level = "FINEST")
    logging::addHandler(logging::writeToFile, file = log_file, level = "DEBUG")
    logging::setLevel(20, logging::getHandler("basic.stdout"))
}

log_tweedieverse_options <- function(options) {
    logging::loginfo("Writing function arguments to log file")
    logging::logdebug("Function arguments")
    logging::logdebug("Output folder: %s", options$output)
    logging::logdebug("Abundance threshold: %f", options$abd_threshold)
    logging::logdebug("Prevalence threshold: %f", options$prev_threshold)
    logging::logdebug("Variance threshold: %f", options$var_threshold)
    logging::logdebug("Base model: %s", options$base_model)
    logging::logdebug("Link function: %s", options$link)
    logging::logdebug("Tweedie variance power: %s", ifelse(is.null(
        options$tweedie_p),
        "NULL", options$tweedie_p
    ))
    logging::logdebug("p = 0 transform: %s", options$p0_transform)
    logging::logdebug("p = 0 transform pseudocount: %f",
        options$p0_transform_pseudocount)
    logging::logdebug("Fixed effects: %s", options$fixed_effects)
    logging::logdebug("Random effects: %s", options$random_effects)
    logging::logdebug("Domain: %s", options$domain)
    logging::logdebug("Normalization: %s", ifelse(is.null(
        options$normalization),
        "NULL", options$normalization
    ))
    logging::logdebug("Offset adjustment: %s", options$adjust_offset)
    logging::logdebug("Scale factor: %s", options$scale_factor)
    logging::logdebug("Max significance: %f", options$max_significance)
    logging::logdebug("Correction method: %s", options$correction)
    logging::logdebug("Run MaAsLin2 for tweedie_p = 0: %s",
        options$Maaslin2_run)
    logging::logdebug("Run presence-absence model: %s",
        options$run_presence_absence_model)
    logging::logdebug("Minimum presence-absence count: %s",
        options$min_presence_absence_count)
    logging::logdebug("Method-specific arguments provided: %s", !is.null(
        options$method_args) ||
        !is.null(options$method.args))
    logging::logdebug("Standardize: %s", options$standardize)
    logging::logdebug("BiocParallel backend: %s", class(options$BPPARAM)[1])
    logging::logdebug("Optimization routine: %s", options$optimizer)
}

align_tweedieverse_samples <- function(data, metadata) {
    logging::loginfo("Aligning assay columns with metadata rows")
    data <- Matrix::t(data)
    if (is.null(colnames(data))) {
        colnames(data) <- as.character(seq_len(ncol(data)))
    }
    metadata <- as.data.frame(metadata)
    if (is.null(rownames(data)) || is.null(rownames(metadata)) || !any(rownames(
        data) %in%
        rownames(metadata))) {
        stop("Assay column names must match sample ",
            "row names in colData or input_metadata.")

    }
    extra_feature_samples <- setdiff(rownames(data), rownames(metadata))
    if (length(extra_feature_samples) > 0) {
        logging::logdebug(
            paste0(
                "The following samples were found to ",
                "have features but no metadata. They ",
                "will be removed. %s"
            ),
            paste(extra_feature_samples, collapse = ",")
        )
    }
    extra_metadata_samples <- setdiff(rownames(metadata), rownames(data))
    if (length(extra_metadata_samples) > 0) {
        logging::logdebug(
            paste0(
                "The following samples were found to ",
                "have metadata but no features. They ",
                "will be removed. %s"
            ),
            paste(extra_metadata_samples, collapse = ",")
        )
    }
    intersect_samples <- intersect(rownames(data), rownames(metadata))
    logging::logdebug(
        "A total of %s samples were found in both the data and metadata",
        length(intersect_samples)
    )
    logging::logdebug("Reordering data/metadata to use same sample ordering")
    data <- data[intersect_samples, , drop = FALSE]
    metadata <- metadata[intersect_samples, , drop = FALSE]
    list(data = data, metadata = metadata)
}

prepare_tweedieverse_references <- function(metadata, fixed_effects,
    reference) {
    if (is.null(reference)) {
        reference <- ","
    }
    split_reference <- unlist(strsplit(reference, "[,;]"))
    for (i in fixed_effects) {
        if (is.numeric(metadata[, i])) {
            next
        }
        if (is.factor(metadata[, i]) && !(i %in% split_reference)) {
            logging::loginfo(
                paste0(
                    "Factor detected for categorial metadata ",
                    "'%s'. Provide a reference argument or ",
                    "manually set factor ordering to change ",
                    "reference level."
                ),
                i
            )
            next
        }
        metadata[, i] <- as.factor(metadata[, i])
        mlevels <- levels(metadata[, i])
        ref <- split_reference[match(i, split_reference) + 1]
        if (length(mlevels) < 2) {
            stop("Provided categorical metadata has fewer ",
                "than 2 unique, non-NA values.")

        }
        if (!is.na(ref)) {
            metadata[, i] <- stats::relevel(metadata[, i], ref = ref)
        } else if (length(mlevels) > 2) {
            {
                condition_message1 <- paste0(
                    "Please provide the reference for the ",
                    "variable '%s' which includes more than ", "2 levels: %s.")
                stop(sprintf(condition_message1, i, paste(as.character(mlevels),

                    collapse = ", ")))
            }
        }
    }
    metadata
}

filter_tweedieverse_prevalence <- function(data, abd_threshold,
    prev_threshold) {
    logging::loginfo("Filter data based on min abundance and min prevalence")
    total_samples <- nrow(data)
    logging::loginfo("Total samples in data: %d", total_samples)
    min_samples <- total_samples * prev_threshold
    logging::loginfo(
        paste0(
            "Min samples required with min abundance ",
            "for a feature not to be filtered: %f"
        ),
        min_samples
    )
    data_zeros <- data
    data_zeros[is.na(data_zeros)] <- 0
    filtered_data <- data[, Matrix::colSums(
        data_zeros > abd_threshold) > min_samples,
        drop = FALSE
    ]
    total_filtered_features <- ncol(data) - ncol(filtered_data)
    logging::loginfo(
        "Total filtered features with prevalence-abundance filtering: %d",
        total_filtered_features
    )
    filtered_feature_names <- setdiff(colnames(data), colnames(filtered_data))
    logging::loginfo("Filtered feature names: %s", toString(
        filtered_feature_names))
    filtered_data
}

filter_tweedieverse_variance <- function(filtered_data, var_threshold) {
    sds <- vapply(seq_len(ncol(filtered_data)), function(i) {
        stats::sd(as.numeric(filtered_data[, i]), na.rm = TRUE)
    }, numeric(1))
    final_features <- filtered_data[, which(sds > var_threshold), drop = FALSE]
    total_filtered_features_var <- ncol(filtered_data) - ncol(final_features)
    logging::loginfo("Total filtered features with variance filtering: %d",
        total_filtered_features_var)
    filtered_feature_names_var <- setdiff(colnames(filtered_data), colnames(
        final_features))
    logging::loginfo("Filtered feature names: %s", toString(
        filtered_feature_names_var))
    final_features
}

resolve_tweedieverse_feature_power <- function(final_features, tweedie_p) {
    has_negative_data <- any(final_features < 0, na.rm = TRUE)
    if (is.null(tweedie_p) && has_negative_data) {
        message("Negative feature values detected after filtering; ",
            "p = 0 is the only Tweedie index whose ",
            "support includes negative values. ", "Setting tweedie_p = 0.")
        tweedie_p <- 0
    } else if (!is.null(tweedie_p) && tweedie_p != 0 && has_negative_data) {
        stop("Negative feature values are supported only when tweedie_p = 0. ",
            "Set tweedie_p = 0, leave tweedie_p ",
                "unspecified, or transform the negatives ",
            "away.")
    }
    if (!is.null(tweedie_p) && tweedie_p < 0) {
        warning("Negative Tweedie variance powers are rarely used. ",
            "The model will run at the supplied ",
            "tweedie_p, but requires strictly ", "positive fitted means.",
            call. = FALSE)
    }
    tweedie_p
}

resolve_tweedieverse_model_options <- function(options, final_features) {
    options["tweedie_p"] <- list(resolve_tweedieverse_feature_power(
        final_features,
        options$tweedie_p
    ))
    {
        if (!is.null(options$tweedie_p) && options$tweedie_p == 0 &&
            options$link !=
            "identity") {
            message("tweedie_p = 0 uses the Gaussian variance case; ",
                "setting link = 'identity' because the ",
                "requested link is undefined for ", "negative responses.")
            options$link <- "identity"
        }
        if (!is.null(options$tweedie_p) && options$tweedie_p == 0 &&
            options$adjust_offset) {
            message("tweedie_p = 0 with identity link does ",
                "not use a log(scale_factor) offset; ",
                "setting adjust_offset = FALSE.")
            options$adjust_offset <- FALSE
        }
        options$link <- options$link
        options$adjust_offset <- options$adjust_offset
        options
    }
}

compute_tweedieverse_offset <- function(data, metadata, options) {
    offset <- NULL
    if (options$adjust_offset) {
        offset <- compute_tweedieverse_size_factor(
            features = data, metadata = metadata,
            normalization = options$normalization, scale_factor =
                options$scale_factor
        )
        logging::loginfo(
            "Offset size factors computed with normalization: %s",
            options$normalization
        )
    } else {
        logging::loginfo(paste0(
            "Offset adjustment disabled; feature ",
            "data remain on the input scale."
        ))
    }
    offset
}

filter_tweedieverse_metadata <- function(metadata, options) {
    if (!is.null(options$scale_factor) && options$scale_factor %in% colnames(
        metadata)) {
        metadata <- dplyr::select(metadata, -dplyr::all_of(
            options$scale_factor))
    }
    temp_filtered_metadata <- metadata[, apply(metadata, 2, tweedie_entropy) >
        options$entropy_threshold, drop = FALSE]
    excluded_metadata <- setdiff(colnames(metadata), colnames(
        temp_filtered_metadata))
    logging::loginfo(
        "Excluded metadata with entropy less or equal to %s: %s",
        options$entropy_threshold, paste(excluded_metadata, collapse = ",")
    )
    filtered_metadata <- temp_filtered_metadata
    filtered_metadata
}

resolve_tweedieverse_fixed_effects <- function(filtered_metadata,
    fixed_effects) {
    if (is.null(fixed_effects)) {
        fixed_effects <- colnames(filtered_metadata)
    } else {
        to_remove <- setdiff(fixed_effects, colnames(filtered_metadata))
        if (length(to_remove) > 0) {
            logging::logwarn(
                paste0(
                    "Feature name not found in metadata so ",
                    "not applied to formula as fixed effect: ",
                    "%s"
                ),
                paste(to_remove, collapse = " , ")
            )
        }
        fixed_effects <- setdiff(fixed_effects, to_remove)
        if (length(fixed_effects) == 0) {
            logging::logerror("No fixed effects included in formula.")
            stop("No fixed effects remain after matching ",
                "the requested effects to metadata.")

        }
    }
    fixed_effects
}

build_tweedieverse_random_formula <- function(random_effects) {
    random_effects_formula_text <- paste("expr ~ (1 | ", paste(random_effects,
        ")",
        sep = "", collapse = " + (1 | "
    ), sep = "")
    logging::loginfo("Formula for random effects: %s",
        random_effects_formula_text)
    random_effects_formula <- tryCatch(as.formula(random_effects_formula_text),
        error = function(e) {
            stop(sprintf("Invalid formula for random effects: %s",
                random_effects_formula_text))
        }
    )
    random_effects_formula
}

resolve_tweedieverse_random_effects <- function(
    filtered_metadata, fixed_effects,
    random_effects
) {
    random_effects_formula <- NULL
    if (!is.null(random_effects)) {
        fixed_effects <- setdiff(fixed_effects, random_effects)
        to_remove <- setdiff(random_effects, colnames(filtered_metadata))
        if (length(to_remove) > 0) {
            logging::logwarn(
                paste0(
                    "Feature name not found in metadata so ",
                    "not applied to formula as random ",
                    "effect: %s"
                ),
                paste(to_remove, collapse = " , ")
            )
        }
        random_effects <- setdiff(random_effects, to_remove)
        if (length(random_effects) > 0) {
            random_effects_formula <- build_tweedieverse_random_formula(
                random_effects)
        }
    }
    list(fixed = fixed_effects, random = random_effects, formula =
        random_effects_formula)
}

build_tweedieverse_fixed_formula <- function(fixed_effects) {
    formula_text <- paste("expr ~ ", paste(fixed_effects, collapse = " + "))
    logging::loginfo("Formula for fixed effects: %s", formula_text)
    formula <- tryCatch(as.formula(formula_text), error = function(e) {
        stop(sprintf(
            "Invalid formula. Please provide a different formula: %s",
            formula_text
        ))
    })
    formula
}

standardize_tweedieverse_metadata <- function(
    filtered_metadata, standardize,
    adjust_offset, offset
) {
    if (standardize) {
        logging::loginfo("Applying z-score to standardize continuous metadata")
        filtered_metadata <- filtered_metadata %>% dplyr::mutate_if(
            is.numeric,
            scale
        )
    } else {
        logging::loginfo("Bypass z-score application to metadata")
    }
    final_metadata <- as.data.frame(filtered_metadata)
    if (adjust_offset) {
        final_metadata$offset <- offset
    }
    final_metadata
}

prepare_tweedieverse_design <- function(metadata, offset, options) {
    filtered_metadata <- filter_tweedieverse_metadata(metadata, options)
    fixed <- resolve_tweedieverse_fixed_effects(filtered_metadata,
        options$fixed_effects)
    effects <- resolve_tweedieverse_random_effects(
        filtered_metadata, fixed,
        options$random_effects
    )
    filtered_metadata <- filtered_metadata[, union(effects$fixed,
        effects$random),
        drop = FALSE
    ]
    list(
        metadata = standardize_tweedieverse_metadata(
            filtered_metadata, options$standardize,
            options$adjust_offset, offset
        ), formula = build_tweedieverse_fixed_formula(effects$fixed),
        random_effects_formula = effects$formula
    )
}

prepare_tweedieverse_response <- function(final_features, options) {
    analysis_features <- final_features
    transformed_features <- NULL
    if (!is.null(options$tweedie_p) && options$tweedie_p == 0) {
        analysis_features <- p0_transform_features(
            features = final_features,
            transform = options$p0_transform, pseudocount =
                options$p0_transform_pseudocount
        )
        if (options$p0_transform != "NONE") {
            transformed_features <- analysis_features
            logging::loginfo(
                "Using p = 0 transformed feature copy with transform: %s",
                options$p0_transform
            )
            if (!options$no_output) {
                transformed_features_file <- file.path(options$output,
                    "transformed_features.tsv")
                logging::loginfo(
                    "Writing p = 0 transformed feature copy to file: %s",
                    transformed_features_file
                )
                write.table(transformed_features,
                    file = transformed_features_file,
                    sep = "\t", quote = FALSE, row.names = TRUE
                )
            }
        }
    } else if (options$p0_transform != "NONE") {
        message("p0_transform is only applied when the ",
            "resolved tweedie_p is 0; ignoring ",
            "p0_transform.")
    }
    analysis_features
}

fit_tweedieverse_abundance <- function(analysis_features, design, options) {
    logging::loginfo("Running selected analysis method: %s", options$base_model)
    fit_data <- fit.Tweedieverse(
        features = analysis_features, metadata = design$metadata,
        base_model = options$base_model, link = options$link, tweedie_p =
            options$tweedie_p,
        Maaslin2_run = options$Maaslin2_run, method_args = options$method_args,
        formula = design$formula, random_effects_formula =
            design$random_effects_formula,
        adjust_offset = options$adjust_offset, correction = options$correction,
        BPPARAM = options$BPPARAM, optimizer = options$optimizer, na.action =
            options$na.action
    )
    fit_data
}

summarize_tweedieverse_features <- function(results, final_features) {
    logging::loginfo("Counting prevalence for each feature")
    detected <- Matrix::colSums(final_features > 0, na.rm = TRUE)
    observed <- nrow(final_features) - Matrix::colSums(is.na(final_features))
    zeros <- observed - Matrix::colSums(final_features != 0, na.rm = TRUE)
    feature_index <- match(results$feature, colnames(final_features))
    results$N <- rep.int(nrow(final_features), nrow(results))
    results$N.not.zero <- unname(detected[feature_index])
    results$percent.zero <- unname(round(zeros / observed, 2)[feature_index] *
        100)
    results
}

adjust_tweedieverse_median <- function(ordered_results, options) {
    mc_input <- ordered_results
    names(mc_input)[names(mc_input) == "feature"] <- "taxon"
    names(mc_input)[names(mc_input) == "coef"] <- "effect_size"
    mc_out <- median_comparison_tweedie(mc_input,
        p_cutoff = 0.95, subtract_median = options$median_subtraction,
        n_sims = 10000, median_threshold = 0, seed =
            options$median_comparison_seed
    )
    ordered_results$coef <- mc_out$coef_median
    ordered_results$pval <- mc_out$pval_median
    ordered_results$qval <- p.adjust(mc_out$pval_median, method =
        options$correction)
    ordered_results <- ordered_results[order(ordered_results$qval), ]
    rownames(ordered_results) <- NULL
    ordered_results
}

adjust_tweedieverse_presence <- function(
    ordered_results, final_features, design,
    options
) {
    presence_results <- fit_presence_absence_model(
        features = final_features,
        metadata = design$metadata, formula = design$formula,
            random_effects_formula = design$random_effects_formula,
        correction = options$correction, BPPARAM = options$BPPARAM,
            min_presence_absence_count = options$min_presence_absence_count
    )
    ordered_results <- combine_abundance_presence_results(
        abundance_results = ordered_results,
        presence_results = presence_results, correction = options$correction
    )
    ordered_results <- ordered_results[order(ordered_results$qval), ]
    rownames(ordered_results) <- NULL
    ordered_results
}

prepare_tweedieverse_results <- function(results, final_features, design,
    options) {
    ordered_results <- summarize_tweedieverse_features(results, final_features)
    ordered_results <- ordered_results[order(ordered_results$qval), , drop =
        FALSE]
    ordered_results <- ordered_results[!is.na(ordered_results$qval), , drop =
        FALSE]
    if (options$median_comparison) {
        ordered_results <- adjust_tweedieverse_median(ordered_results, options)
    }
    if (options$run_presence_absence_model) {
        ordered_results <- adjust_tweedieverse_presence(
            ordered_results, final_features,
            design, options
        )
    }
    ordered_results <- dplyr::select(ordered_results, c(
        "feature", "metadata",
        "value", "coef", "stderr", "pval", "qval"
    ), everything())
    ordered_results
}

write_tweedieverse_results <- function(
    ordered_results, significant_results,
    output, max_significance
) {
    results_file <- file.path(output, "all_results.tsv")
    logging::loginfo(
        "Writing all results to file (ordered by increasing q-values): %s",
        results_file
    )
    write.table(ordered_results,
        file = results_file, sep = "\t", quote = FALSE,
        row.names = FALSE
    )
    significant_results_file <- file.path(output, "significant_results.tsv")
    logging::loginfo(
        paste0(
            "Writing the significant results (those ",
            "which are less than or equal to the ",
            "threshold of %f) to file (ordered by ",
            "increasing q-values): %s"
        ),
        max_significance, significant_results_file
    )
    write.table(significant_results,
        file = significant_results_file, sep = "\t",
        quote = FALSE, row.names = FALSE
    )
}

plot_tweedieverse_index <- function(ordered_results, figures_folder, output) {
    logging::loginfo("Writing Tweedie inxed plot to file: %s", output)
    tryCatch(
        {
            tweedie_index_plot(ordered_results, figures_folder)
        },
        error = function(err) {
            logging::logerror(paste0(
                "Unable to do make a Tweedie inxed plot ",
                "of results!!!"
            ))
            logging::logerror(err)
        }
    )
}

plot_tweedieverse_heatmap <- function(
    significant_results, figures_folder, output,
    options
) {
    significant_results_file <- file.path(options$output,
        "significant_results.tsv")
    if (options$plot_heatmap && nrow(significant_results) > 0 && length(unique(
        significant_results$metadata)) >
        1) {
        heatmap_file <- file.path(options$output, "Tweedieverse_Heatmap.pdf")
        logging::loginfo(
            "Writing heatmap of significant results to file: %s",
            heatmap_file
        )
        tryCatch(
            {
                save_heatmap(significant_results_file, heatmap_file,
                    figures_folder,
                    first_n = options$heatmap_first_n
                )
            },
            error = function(err) {
                logging::logerror("Unable to do make a hetamp of results!!!")
                logging::logerror(err)
            }
        )
    }
}

plot_tweedieverse_associations <- function(
    metadata, final_features, significant_results,
    offset, figures_folder, output, options
) {
    significant_results_file <- file.path(options$output,
        "significant_results.tsv")
    if (options$plot_scatter) {
        logging::loginfo(
            paste0(
                "Writing association plots (one for each ",
                "significant association) to output ",
                "folder: %s"
            ),
            options$output
        )
        plot_features <- as.data.frame(as.matrix(final_features[, unique(
            significant_results$feature),
            drop = FALSE
        ]))
        if (!is.null(offset)) {
            plot_features <- plot_features / offset
        }
        association_plots(
            metadata, plot_features, significant_results_file,
            options$output, figures_folder
        )
    }
}

output_tweedieverse_results <- function(
    ordered_results, final_features, metadata,
    offset, figures_folder, options
) {
    significant_results <- ordered_results[which(ordered_results$qval <=
        options$max_significance), ,
        drop = FALSE
    ]
    if (!is.null(options$output)) {
        write_tweedieverse_results(
            ordered_results, significant_results, options$output,
            options$max_significance
        )
        plot_tweedieverse_index(ordered_results, figures_folder, options$output)
        plot_tweedieverse_heatmap(
            significant_results, figures_folder, options$output,
            options
        )
        plot_tweedieverse_associations(
            metadata, final_features, significant_results,
            offset, figures_folder, options$output, options
        )
    }
    if (options$return_all) {
        ordered_results
    } else {
        significant_results
    }
}
