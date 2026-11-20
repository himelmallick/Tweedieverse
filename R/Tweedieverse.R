#' Differential analysis of multi-omics data using Tweedie GLMs
#'
#' Fit a per-feature Tweedie generalized linear model to omics features.

#' @param input_features A domain-appropriate Bioconductor container, such as
#' \code{TreeSummarizedExperiment} for \code{domain = "microbiome"},
#' \code{SingleCellExperiment} for \code{domain = "single_cell"}, and
#' \code{SummarizedExperiment} for \code{domain = "bulk_rnaseq"}. The
#' \code{assays} slot contains the expression
#' or abundance matrix and defaults to an assay named \code{"counts"}.
#' A \code{MultiAssayExperiment} object is also accepted. In that case,
#' \code{Tweedieverse()} iterates over each experiment/omics layer and returns a
#' named
#' list of omics-specific results.
#' This matrix should have one row for each feature and one sample for each
#' column.
#' The \code{colData} slot should contain a data frame with one row per
#' sample and columns that contain metadata for each sample.
#' Additional information about the experiment can be contained in the
#' \code{metadata} slot as a list.
#' @param input_metadata Optional R data frame of metadata. Samples are expected
#' to have
#' matching sample names with \code{input_features}. If \code{NULL}, metadata
#' are taken from
#' the Bioconductor container's \code{colData}.
#' @param output The output folder to write results.
#' @param assay_name If the input is provided as one of the accepted
#' Bioconductor objects,
#' this argument selects the name of the assay slot in the input object that
#' contains the omics measurements.
#' For \code{MultiAssayExperiment} input, most arguments can be supplied as a
#' single value
#' used for every omics layer, or as a named list/vector keyed by experiment
#' name.
#' @param abd_threshold If prevalence-abundance filtering is desired, only
#' features that are present (or detected)
#' in at least \code{prev_threshold} percent of samples at \code{abd_threshold}
#' minimum abundance (read count or proportion)
#' are retained. Default value for \code{abd_threshold} is \code{0.0}.
#' To disable prevalence-abundance filtering, set \code{abd_threshold = -Inf}.
#' @param prev_threshold If prevalence-abundance filtering is desired, only
#' features that are present (or detected)
#' in at least \code{prev_threshold} percent of samples at \code{abd_threshold}
#' minimum abundance (read count or proportion)
#' are retained. Default value for \code{prev_threshold} is \code{0.1}.
#' @param var_threshold If variance filtering is desired, only features that
#' have variances greater than
#' \code{var_threshold} are retained. This step is done after the
#' prevalence-abundance filtering.
#' Default value for \code{var_threshold} is \code{0.0} (i.e. no variance
#' filtering).
#' @param entropy_threshold Metadata columns are retained when their entropy,
#' computed by
#' \code{tweedie_entropy()}, exceeds this threshold. Default is 0, which removes
#' metadata columns with zero entropy, including constant columns.
#' @param base_model The per-feature base model. Only "CPLM" is supported.
#' @param link A specification of the GLM link function. Default is "log". Must
#' be one of "log", "identity", "sqrt", or "inverse".
#' @param tweedie_p Numeric Tweedie variance power. Default is \code{NULL},
#' which estimates the compound-Poisson
#' Tweedie index in the usual \code{1 < p < 2} range when all filtered feature
#' values are non-negative.
#' If \code{tweedie_p} is \code{NULL} or \code{NA} and filtered features contain
#' negative values, it is resolved
#' to \code{0}, the Gaussian case. Values between 0 and 1 are undefined and
#' rejected. Other supplied values
#' pin the model to that fixed index. Negative supplied powers are allowed with
#' a warning for non-negative data,
#' but negative feature values are supported only at \code{p = 0}.
#' Random-effect models support fixed \code{p = 0}, \code{p = 1}, \code{1 < p <
#' 2}, \code{p = 2}, and \code{p = 3}.
#' @param p0_transform Transformation to apply to a copy of the filtered feature
#' table when the
#' resolved Tweedie index is \code{0}. The original input data are not modified.
#' Must be one of
#' \code{"NONE"}, \code{"CLR"}, \code{"RCLR"}, \code{"LOG"}, or
#' \code{"ARC_SIGNED_SQRT"}.
#' The transformed copy is used for the \code{p = 0} abundance model and written
#' to
#' \code{transformed_features.tsv} when \code{output} is provided. Default is
#' \code{"NONE"}.
#' @param p0_transform_pseudocount Numeric pseudocount used for
#' \code{p0_transform = "CLR"} and
#' \code{p0_transform = "LOG"}. Default is 1.
#' @param fixed_effects Metadata variable(s) describing the fixed effects
#' coefficients.
#' @param random_effects Metadata variable(s) describing the random effects part
#' of the model.
#' @param domain Domain used to choose the default size-factor normalization.
#' Must be one of
#' \code{"microbiome"}, \code{"single_cell"}, \code{"bulk_rnaseq"}, or
#' \code{"custom"}.
#' The default is \code{"microbiome"}.
#' @param normalization Size-factor strategy used only to compute the model
#' offset. The input feature
#' table is never normalized or transformed. If \code{NULL}, the domain default
#' is used:
#' \code{"TSS"} for microbiome, \code{"SCRAN"} for single-cell, \code{"TMM"}
#' for bulk RNA-seq, and \code{"MEDIAN"} for custom/other omics.
#' Supported strategies are \code{"TSS"}, \code{"GMPR"}, \code{"CSS"},
#' \code{"SCRAN"}, \code{"TMM"}, \code{"RLE"} / \code{"DESEQ2"}, \code{"CPM"},
#' \code{"MEDIAN"}, and \code{"NONE"}.
#' @param adjust_offset If TRUE (default), an offset term will be included as
#' the logarithm of the
#' size factor estimated by \code{normalization}, or by \code{scale_factor} when
#' supplied.
#' @param scale_factor Name of a numerical metadata variable containing
#' user-supplied sample size
#' factors. When supplied, this overrides \code{domain} and \code{normalization}
#' for offset calculation.
#' @param max_significance The q-value threshold for significance. Default is
#' 0.05.
#' @param correction The correction method for computing the q-value (see
#' \code{\link[stats]{p.adjust}} for options, default is 'BH').
#' @param Maaslin2_run Logical. For \code{tweedie_p = 0}, run the MaAsLin2
#' linear-model path instead of the
#' Tweedie GLM path. Set to FALSE to fit \code{tweedie_p = 0} with the Tweedie
#' GLM. Default is TRUE.
#' When the resolved Tweedie index is \code{0}, non-identity links are changed
#' to \code{"identity"} and
#' \code{adjust_offset} is disabled.
#' @param median_comparison If TRUE, coefficients will be tested against a null
#' value corresponding to the median coefficient for a covariate in the
#' \code{metadata}. Default is FALSE. Should only be used for relative abundance
#' data.
#' @param median_subtraction If TRUE, coefficients minus median will be used for
#' compositionality adjustment.
#' @param median_comparison_seed Optional non-negative integer seed for the
#' median-comparison
#' simulations. Default is NULL, which uses the current random-number stream.
#' A supplied seed is scoped with withr and restores the preceding random-number
#' state
#' after the median-comparison calculation.
#' @param run_presence_absence_model If TRUE, also fit a DAssemble-style
#' presence-absence logistic regression model
#' alongside the selected abundance model, including the MaAsLin2 path used by
#' default for \code{tweedie_p = 0}.
#' Retains individual abundance and presence-absence model results, and ranks
#' features by a Cauchy combination test
#' of the two p-values. Default is FALSE.
#' @param min_presence_absence_count Minimum number of detected samples and
#' minimum number of
#' undetected samples required for the presence-absence fit. Default is 3 for
#' each.
#' For retained features below either minimum, the presence-absence p-value is
#' missing
#' and the joint p-value uses the abundance p-value directly.
#' @param method_args Optional named list of method-specific arguments. For
#' \code{tweedie_p = 0}, use
#' \code{method_args = list(Maaslin2 = list(...))} to pass arguments to
#' \code{\link[Maaslin2]{Maaslin2}};
#' user-supplied values override Tweedieverse defaults.
#' @param method.args Alias for \code{method_args}.
#' @param standardize Should continuous metadata be standardized? Default is
#' TRUE. Bypassed for categorical variables.
#' @param BPPARAM \code{\link[BiocParallel]{BiocParallelParam}} object
#' controlling Bioconductor-native parallel evaluation. Default is
#' \code{\link[BiocParallel]{SerialParam}}.
#' @param optimizer The optimization routine to be used for estimating the
#' parameters of the Tweedie model.
#' Possible choices are \code{"nlminb"} (the default, see
#' \code{\link[stats]{nlminb}}),
#' \code{"bobyqa"} (\code{\link[minqa]{bobyqa}}), and \code{"L-BFGS-B"}
#' (\code{\link[stats]{optim}}).
#' Ignored for random effects modeling which uses an alternative Template Model
#' Builder (TMB) approach (\code{\link[glmmTMB]{glmmTMB}}).
#' @param na.action How to handle missing values? See \code{\link{na.action}}.
#' Default is \code{\link{na.exclude}}.
#' @param plot_heatmap Logical. If TRUE (default is FALSE), generate a heatmap
#' of the (top \code{heatmap_first_n}) significant associations.
#' @param plot_scatter Logical. If TRUE (default is FALSE), generate scatter/box
#' plots of individual associations.
#' @param heatmap_first_n In heatmap, plot top N features with significant
#' associations (default is 50).
#' @param reference The factor to use as a reference for a variable with more
#' than two levels provided as a string of 'variable,reference' semi-colon
#' delimited for multiple variables (default is NULL).
#' @param return_all Logical. If TRUE, return the full result table before the
#' significance threshold is applied, including when output is NULL.
#' Default is FALSE, which returns associations passing max_significance.
#'
#' @importFrom grDevices colorRampPalette dev.off jpeg pdf
#' @importFrom stats coef fitted as.formula na.exclude offset p.adjust plogis
#'     relevel sd update
#' @importFrom utils capture.output read.table type.convert write.table
#' @importFrom SummarizedExperiment colData
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom dplyr %>% everything
#' @return For single-omics input, a data frame containing coefficient
#' estimates, p-values,
#' and q-values (multiplicity-adjusted p-values) is returned. For
#' \code{MultiAssayExperiment}
#' input, a named list of omics-specific result data frames is returned.
#' With \code{return_all = TRUE}, each table contains the results before
#' filtering by \code{max_significance}.
#'
#' @author Himel Mallick, \email{him4004@@med.cornell.edu}
#'
#' @examples
#'
#' set.seed(123)
#' features <- as.data.frame(matrix(rpois(24, lambda = 5), nrow = 12, ncol = 2))
#' colnames(features) <- c("feature1", "feature2")
#' rownames(features) <- paste0("sample", seq_len(nrow(features)))
#' metadata <- data.frame(
#'   group = rep(c("A", "B"), each = 6),
#'   row.names = rownames(features)
#' )
#' experiment <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = t(as.matrix(features))),
#'   colData = metadata
#' )
#' fit <- Tweedieverse(
#'   input_features = experiment,
#'   output = NULL,
#'   fixed_effects = "group",
#'   median_comparison = FALSE,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#' head(fit)
#'
#' maaslin2_median_fit <- Tweedieverse(
#'   input_features = experiment,
#'   output = NULL,
#'   fixed_effects = "group",
#'   tweedie_p = 0,
#'   Maaslin2_run = TRUE,
#'   median_comparison = TRUE,
#'   median_subtraction = TRUE,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#'
#' # For Bioconductor-container workflow examples, see:
#' # vignette("Tweedieverse", package = "Tweedieverse")
#' @keywords microbiome metagenomics multiomics scRNASeq tweedie singlecell
#' @export
Tweedieverse <- function(input_features,
    input_metadata = NULL,
    output = NULL,
    assay_name = "counts",
    abd_threshold = 0.0,
    prev_threshold = 0.1,
    var_threshold = 0.0,
    entropy_threshold = 0.0,
    base_model = "CPLM",
    link = "log",
    tweedie_p = NULL,
    p0_transform = "NONE",
    p0_transform_pseudocount = 1,
    fixed_effects = NULL,
    random_effects = NULL,
    domain = "microbiome",
    normalization = NULL,
    adjust_offset = TRUE,
    scale_factor = NULL,
    max_significance = 0.05,
    correction = "BH",
    Maaslin2_run = TRUE,
    median_comparison = FALSE,
    median_subtraction = FALSE,
    run_presence_absence_model = FALSE,
    method_args = NULL,
    method.args = NULL,
    standardize = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    optimizer = "nlminb",
    na.action = na.exclude,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    heatmap_first_n = 50,
    reference = NULL,
    min_presence_absence_count = 3L,
    return_all = FALSE,
    median_comparison_seed = NULL) {
    valid_input_classes <- c(
        "MultiAssayExperiment",
        unique(unlist(domain_bioc_container_map()))
    )
    if (!inherits(input_features, valid_input_classes)) {
        {
            condition_message1 <- paste0("input_features must be a supported ",
                "Bioconductor container: %s. Plain data ",
                    "frames, matrices, lists, and file paths ",
                "are not supported.")
            stop(sprintf(condition_message1, paste(valid_input_classes,
                collapse = ", ")),
                call. = FALSE)
        }
    }
    if (is.character(input_metadata)) {
        stop(
            "input_metadata file paths are not supported. ",
            "Use colData(input_features), colData on a MultiAssayExperiment, ",
            "or provide input_metadata as a data.frame.",
            call. = FALSE
        )
    }

    if (inherits(input_features, "MultiAssayExperiment")) {
        return(run_multiassay_tweedieverse(
            input_features = input_features,
            input_metadata = input_metadata,
            output = output,
            assay_name = assay_name,
            abd_threshold = abd_threshold,
            prev_threshold = prev_threshold,
            var_threshold = var_threshold,
            entropy_threshold = entropy_threshold,
            base_model = base_model,
            link = link,
            tweedie_p = tweedie_p,
            p0_transform = p0_transform,
            p0_transform_pseudocount = p0_transform_pseudocount,
            fixed_effects = fixed_effects,
            random_effects = random_effects,
            domain = domain,
            normalization = normalization,
            adjust_offset = adjust_offset,
            scale_factor = scale_factor,
            max_significance = max_significance,
            correction = correction,
            Maaslin2_run = Maaslin2_run,
            median_comparison = median_comparison,
            median_subtraction = median_subtraction,
            run_presence_absence_model = run_presence_absence_model,
            method_args = method_args,
            method.args = method.args,
            standardize = standardize,
            BPPARAM = BPPARAM,
            optimizer = optimizer,
            na.action = na.action,
            plot_heatmap = plot_heatmap,
            plot_scatter = plot_scatter,
            heatmap_first_n = heatmap_first_n,
            reference = reference,
            min_presence_absence_count = min_presence_absence_count,
            return_all = return_all,
            median_comparison_seed = median_comparison_seed
        ))
    }

    .Tweedieverse_single(
        input_features = input_features,
        input_metadata = input_metadata,
        output = output,
        assay_name = assay_name,
        abd_threshold = abd_threshold,
        prev_threshold = prev_threshold,
        var_threshold = var_threshold,
        entropy_threshold = entropy_threshold,
        base_model = base_model,
        link = link,
        tweedie_p = tweedie_p,
        p0_transform = p0_transform,
        p0_transform_pseudocount = p0_transform_pseudocount,
        fixed_effects = fixed_effects,
        random_effects = random_effects,
        domain = domain,
        normalization = normalization,
        adjust_offset = adjust_offset,
        scale_factor = scale_factor,
        max_significance = max_significance,
        correction = correction,
        Maaslin2_run = Maaslin2_run,
        median_comparison = median_comparison,
        median_subtraction = median_subtraction,
        run_presence_absence_model = run_presence_absence_model,
        method_args = method_args,
        method.args = method.args,
        standardize = standardize,
        BPPARAM = BPPARAM,
        optimizer = optimizer,
        na.action = na.action,
        plot_heatmap = plot_heatmap,
        plot_scatter = plot_scatter,
        heatmap_first_n = heatmap_first_n,
        reference = reference,
        min_presence_absence_count = min_presence_absence_count,
        return_all = return_all,
        median_comparison_seed = median_comparison_seed
    )
}

.Tweedieverse_single <- function(...) {
    options <- utils::modifyList(
        list(method.args = NULL), list(...),
        keep.null = TRUE
    )
    options <- validate_tweedieverse_options(options)
    input <- extract_tweedieverse_input(options)
    figures_folder <- prepare_tweedieverse_output(options$output)
    start_tweedieverse_logging(options$output)
    log_tweedieverse_options(options)
    aligned <- align_tweedieverse_samples(input$data, input$metadata)
    options["fixed_effects"] <- list(parse_effect_names(options$fixed_effects))
    options["random_effects"] <- list(parse_effect_names(
        options$random_effects))
    aligned$metadata <- prepare_tweedieverse_references(
        aligned$metadata, options$fixed_effects, options$reference
    )
    filtered <- filter_tweedieverse_prevalence(
        aligned$data, options$abd_threshold, options$prev_threshold
    )
    features <- filter_tweedieverse_variance(filtered, options$var_threshold)
    options <- resolve_tweedieverse_model_options(options, features)
    offset <- compute_tweedieverse_offset(aligned$data, aligned$metadata,
        options)
    design <- prepare_tweedieverse_design(aligned$metadata, offset, options)
    analysis_features <- prepare_tweedieverse_response(features, options)
    fit <- fit_tweedieverse_abundance(analysis_features, design, options)
    results <- prepare_tweedieverse_results(fit$results, features, design,
        options)
    output_tweedieverse_results(
        results, features, aligned$metadata, offset, figures_folder, options
    )
}


option_not_valid_error <- function(message, valid_options) {
    logging::logerror("%s: %s", message, toString(valid_options))
    stop("Option not valid", call. = FALSE)
}

## Quiets concerns of R CMD check
utils::globalVariables(c(
    "base.model",
    "base.model_abundance",
    "base.model_presence",
    "bin",
    "coef",
    "coef_abundance",
    "coef_presence",
    "data",
    "feature",
    "metadata",
    "name",
    "pval",
    "pval_abundance",
    "pval_presence",
    "qval",
    "qval_abundance",
    "qval_presence",
    "stderr",
    "stderr_abundance",
    "stderr_presence",
    "tweedie.index",
    "value",
    "xnames"
))
