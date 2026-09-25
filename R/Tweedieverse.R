#' Differential analysis of multi-omics data using Tweedie GLMs
#'
#' Fit a per-feature Tweedie generalized linear model to omics features.

#' @param input_features A domain-appropriate Bioconductor container, such as
#' \code{TreeSummarizedExperiment} for \code{domain = "microbiome"},
#' \code{SingleCellExperiment} for \code{domain = "single_cell"}, and
#' \code{SummarizedExperiment} for \code{domain = "bulk_rnaseq"}. The \code{assays} slot contains the expression
#' or abundance matrix and defaults to an assay named \code{"counts"}.
#' A \code{MultiAssayExperiment} object is also accepted. In that case,
#' \code{Tweedieverse()} iterates over each experiment/omics layer and returns a named
#' list of omics-specific results.
#' This matrix should have one row for each feature and one sample for each column.  
#' The \code{colData} slot should contain a data frame with one row per 
#' sample and columns that contain metadata for each sample. 
#' Additional information about the experiment can be contained in the
#' \code{metadata} slot as a list.
#' @param input_metadata Optional R data frame of metadata. Samples are expected to have
#' matching sample names with \code{input_features}. If \code{NULL}, metadata are taken from
#' the Bioconductor container's \code{colData}.
#' @param output The output folder to write results.
#' @param assay_name If the input is provided as one of the accepted Bioconductor objects,
#' this argument selects the name of the assay slot in the input object that contains the omics measurements.
#' For \code{MultiAssayExperiment} input, most arguments can be supplied as a single value
#' used for every omics layer, or as a named list/vector keyed by experiment name.
#' @param abd_threshold If prevalence-abundance filtering is desired, only features that are present (or detected)
#' in at least \code{prev_threshold} percent of samples at \code{abd_threshold} minimum abundance (read count or proportion)
#' are retained. Default value for \code{abd_threshold} is \code{0.0}.
#' To disable prevalence-abundance filtering, set \code{abd_threshold = -Inf}.
#' @param prev_threshold If prevalence-abundance filtering is desired, only features that are present (or detected)
#' in at least \code{prev_threshold} percent of samples at \code{abd_threshold} minimum abundance (read count or proportion)
#' are retained. Default value for \code{prev_threshold} is \code{0.1}.
#' @param var_threshold If variance filtering is desired, only features that have variances greater than
#' \code{var_threshold} are retained. This step is done after the prevalence-abundance filtering.
#' Default value for \code{var_threshold} is \code{0.0} (i.e. no variance filtering).
#' @param entropy_threshold If entropy-based filtering is desired for metadata, only features that have entropy greater than
#' \code{entropy_threshold} are retained. Default value for \code{entropy_threshold} is \code{0.0} (i.e. no entropy filtering).
#' @param base_model The per-feature base model. Only "CPLM" is supported.
#' @param link A specification of the GLM link function. Default is "log". Must be one of "log", "identity", "sqrt", or "inverse".
#' @param tweedie_p Numeric Tweedie variance power. Default is \code{NULL}, which estimates the compound-Poisson
#' Tweedie index in the usual \code{1 < p < 2} range when all filtered feature values are non-negative.
#' If \code{tweedie_p} is \code{NULL} or \code{NA} and filtered features contain negative values, it is resolved
#' to \code{0}, the Gaussian case. Values between 0 and 1 are undefined and rejected. Other supplied values
#' pin the model to that fixed index. Negative supplied powers are allowed with a warning for non-negative data,
#' but negative feature values are supported only at \code{p = 0}.
#' Random-effect models support fixed \code{p = 0}, \code{p = 1}, \code{1 < p < 2}, \code{p = 2}, and \code{p = 3}.
#' @param p0_transform Transformation to apply to a copy of the filtered feature table when the
#' resolved Tweedie index is \code{0}. The original input data are not modified. Must be one of
#' \code{"NONE"}, \code{"CLR"}, \code{"RCLR"}, \code{"LOG"}, or \code{"ARC_SIGNED_SQRT"}.
#' The transformed copy is used for the \code{p = 0} abundance model and written to
#' \code{transformed_features.tsv} when \code{output} is provided. Default is \code{"NONE"}.
#' @param p0_transform_pseudocount Numeric pseudocount used for \code{p0_transform = "CLR"} and
#' \code{p0_transform = "LOG"}. Default is 1.
#' @param fixed_effects Metadata variable(s) describing the fixed effects coefficients.
#' @param random_effects Metadata variable(s) describing the random effects part of the model.
#' @param domain Domain used to choose the default size-factor normalization. Must be one of
#' \code{"microbiome"}, \code{"single_cell"}, \code{"bulk_rnaseq"}, or \code{"custom"}.
#' The default is \code{"microbiome"}.
#' @param normalization Size-factor strategy used only to compute the model offset. The input feature
#' table is never normalized or transformed. If \code{NULL}, the domain default is used:
#' \code{"TSS"} for microbiome, \code{"SCRAN"} for single-cell, \code{"TMM"}
#' for bulk RNA-seq, and \code{"MEDIAN"} for custom/other omics.
#' Supported strategies are \code{"TSS"}, \code{"GMPR"}, \code{"CSS"},
#' \code{"SCRAN"}, \code{"TMM"}, \code{"RLE"} / \code{"DESEQ2"}, \code{"CPM"},
#' \code{"MEDIAN"}, and \code{"NONE"}.
#' @param adjust_offset If TRUE (default), an offset term will be included as the logarithm of the
#' size factor estimated by \code{normalization}, or by \code{scale_factor} when supplied.
#' @param scale_factor Name of a numerical metadata variable containing user-supplied sample size
#' factors. When supplied, this overrides \code{domain} and \code{normalization} for offset calculation.
#' @param max_significance The q-value threshold for significance. Default is 0.05.
#' @param correction The correction method for computing the q-value (see \code{\link[stats]{p.adjust}} for options, default is 'BH').
#' @param Maaslin2_run Logical. For \code{tweedie_p = 0}, run the MaAsLin2 linear-model path instead of the
#' Tweedie GLM path. Set to FALSE to fit \code{tweedie_p = 0} with the Tweedie GLM. Default is TRUE.
#' When the resolved Tweedie index is \code{0}, non-identity links are changed to \code{"identity"} and
#' \code{adjust_offset} is disabled.
#' @param median_comparison If TRUE, coefficients will be tested against a null value corresponding to the median coefficient for a covariate in the \code{metadata}. Default is FALSE. Should only be used for relative abundance data.
#' @param median_subtraction If TRUE, coefficients minus median will be used for compositionality adjustment. 
#' @param run_presence_absence_model If TRUE, also fit a DAssemble-style presence-absence logistic regression model
#' alongside the selected abundance model, including the MaAsLin2 path used by default for \code{tweedie_p = 0}.
#' Retains individual abundance and presence-absence model results, and ranks features by a Cauchy combination test
#' of the two p-values. Default is FALSE.
#' @param method_args Optional named list of method-specific arguments. For \code{tweedie_p = 0}, use
#' \code{method_args = list(Maaslin2 = list(...))} to pass arguments to \code{\link[Maaslin2]{Maaslin2}};
#' user-supplied values override Tweedieverse defaults.
#' @param method.args Alias for \code{method_args}.
#' @param standardize Should continuous metadata be standardized? Default is TRUE. Bypassed for categorical variables.
#' @param cores An integer that indicates the number of R processes to run in parallel. Default is 1.
#' @param optimizer The optimization routine to be used for estimating the parameters of the Tweedie model.
#' Possible choices are \code{"nlminb"} (the default, see \code{\link[stats]{nlminb}}),
#' \code{"bobyqa"} (\code{\link[minqa]{bobyqa}}), and \code{"L-BFGS-B"} (\code{\link[stats]{optim}}).
#' Ignored for random effects modeling which uses an alternative Template Model Builder (TMB) approach (\code{\link[glmmTMB]{glmmTMB}}).
#' @param na.action How to handle missing values? See \code{\link{na.action}}. Default is \code{\link{na.exclude}}.
#' @param plot_heatmap Logical. If TRUE (default is FALSE), generate a heatmap of the (top \code{heatmap_first_n}) significant associations.
#' @param plot_scatter Logical. If TRUE (default is FALSE), generate scatter/box plots of individual associations.
#' @param heatmap_first_n In heatmap, plot top N features with significant associations (default is 50).
#' @param reference The factor to use as a reference for a variable with more than two levels provided as a string of 'variable,reference' semi-colon delimited for multiple variables (default is NULL).
#'
#' @importFrom grDevices colorRampPalette dev.off jpeg pdf
#' @importFrom stats coef fitted as.formula na.exclude p.adjust plogis relevel sd update
#' @importFrom utils capture.output read.table type.convert write.table
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr %>% everything
#' @importFrom parallel clusterExport
#' @return For single-omics input, a data frame containing coefficient estimates, p-values,
#' and q-values (multiplicity-adjusted p-values) is returned. For \code{MultiAssayExperiment}
#' input, a named list of omics-specific result data frames is returned.
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
#'   cores = 1
#' )
#' head(fit)
#' 
#' \dontrun{
#' maaslin2_median_fit <- Tweedieverse(
#'   input_features = experiment,
#'   output = NULL,
#'   fixed_effects = "group",
#'   tweedie_p = 0,
#'   Maaslin2_run = TRUE,
#'   median_comparison = TRUE,
#'   median_subtraction = TRUE,
#'   cores = 1
#' )
#' }
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
                         cores = 1,
                         optimizer = "nlminb",
                         na.action = na.exclude,
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         heatmap_first_n = 50,
                         reference = NULL) {
  valid_input_classes <- c(
    "MultiAssayExperiment",
    unique(unlist(domain_bioc_container_map()))
  )
  if (!inherits(input_features, valid_input_classes)) {
    stop(
      paste(
        "input_features must be a supported Bioconductor container:",
        paste(valid_input_classes, collapse = ", "),
        ". Plain data frames, matrices, lists, and file paths are not supported."
      ),
      call. = FALSE
    )
  }
  if (is.character(input_metadata)) {
    stop(
      paste(
        "input_metadata file paths are not supported.",
        "Use colData(input_features), colData on a MultiAssayExperiment,",
        "or provide input_metadata as a data.frame."
      ),
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
      cores = cores,
      optimizer = optimizer,
      na.action = na.action,
      plot_heatmap = plot_heatmap,
      plot_scatter = plot_scatter,
      heatmap_first_n = heatmap_first_n,
      reference = reference
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
    cores = cores,
    optimizer = optimizer,
    na.action = na.action,
    plot_heatmap = plot_heatmap,
    plot_scatter = plot_scatter,
    heatmap_first_n = heatmap_first_n,
    reference = reference
  )
}

.Tweedieverse_single <- function(input_features,
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
                                 cores = 1,
                                 optimizer = "nlminb",
                                 na.action = na.exclude,
                                 plot_heatmap = FALSE,
                                 plot_scatter = FALSE,
                                 heatmap_first_n = 50,
                                 reference = NULL) {
  
  

  #################################
  # Specify all available options #
  #################################
  
  no_output <- is.null(output)
  
  model_choices <- c("CPLM")
  link_choices <- c("log", "identity", "sqrt", "inverse")
  p0_transform_choices <- c("NONE", "CLR", "RCLR", "LOG", "ARC_SIGNED_SQRT")
  domain_choices <- c("microbiome", "single_cell", "bulk_rnaseq", "custom")
  correction_choices <-
    c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY")
  optimizer_choices <- c("nlminb", "bobyqa", "L-BFGS-B")
  
  #######################################################################
  #=====================================================================#
  # Read in the data and metadata, create output folder, initialize log #
  #=====================================================================#
  #######################################################################
  
  #########################################
  # Support multiple Bioconductor classes #
  #########################################

  domain <- tolower(domain)
  if (!domain %in% domain_choices) {
    option_not_valid_error(
      "Please select a domain from the list of available options",
      toString(domain_choices)
    )
  }

  valid_classes <- unique(unlist(domain_bioc_container_map()))
  
  ##############################################################
  # Extract features and metadata based on user-provided input #
  ##############################################################
  
  is_supported_bioc_class <- inherits(input_features, valid_classes)
  if (is_supported_bioc_class) {
    validate_domain_bioc_container(input_features, domain)
    data <- extractAssay(input_features, assay_name)
    if (is.null(input_metadata)) {
      metadata <- data.frame(SummarizedExperiment::colData(input_features))
    } else {
      metadata <- input_metadata
    }
  } else {
    stop(
      sprintf(
        paste(
          "Input data of class <%s> not supported.",
          "Please use a domain-appropriate Bioconductor container."
        ),
        class(input_features)[1]
      )
    )
  }
    
  # create an output folder and figures folder if it does not exist
  if (!no_output) {
    if (file.exists(output) && !dir.exists(output)) {
      stop(sprintf("Output path exists but is not a directory: %s", output))
    }
    if (!dir.exists(output)) {
      message(sprintf("Creating output folder: %s", output))
      dir.create(output, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(output)) {
      stop(sprintf("Unable to create output folder: %s", output))
    }
    
    #if (plot_heatmap || plot_scatter) {
    figures_folder <- file.path(output, "figures")
    if (!dir.exists(figures_folder)) {
      message(sprintf("Creating output figures folder: %s", figures_folder))
      dir.create(figures_folder, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(figures_folder)) {
      stop(sprintf("Unable to create output figures folder: %s", figures_folder))
    }
    #}
    
    # Create log file (write info to stdout and debug level to log file)
    # Set level to finest so all log levels are reviewed
    log_file <- file.path(output, "Tweedieverse.log")
    # Remove log file if already exists (to avoid append)
    if (file.exists(log_file)) {
      print(paste("Warning: Deleting existing log file:", log_file))
      unlink(log_file)
    }
    logging::basicConfig(level = 'FINEST')
    logging::addHandler(logging::writeToFile,
                        file = log_file, level = "DEBUG")
    logging::setLevel(20, logging::getHandler('basic.stdout'))
  } else {
    # no_output mode: no folder, no figures, no log file on disk
    figures_folder <- NULL
  }
  
  #####################
  # Log the arguments #
  #####################
  
  logging::loginfo("Writing function arguments to log file")
  logging::logdebug("Function arguments")
  if (is.character(input_features)) {
    logging::logdebug("Input data file: %s", input_features)
  }
  if (is.character(input_metadata)) {
    logging::logdebug("Input metadata file: %s", input_metadata)
  }
  logging::logdebug("Output folder: %s", output)
  logging::logdebug("Abundance threshold: %f", abd_threshold)
  logging::logdebug("Prevalence threshold: %f", prev_threshold)
  logging::logdebug("Variance threshold: %f", var_threshold)
  logging::logdebug("Base model: %s", base_model)
  logging::logdebug("Link function: %s", link)
  logging::logdebug("Tweedie variance power: %s", ifelse(is.null(tweedie_p), "NULL", tweedie_p))
  logging::logdebug("p = 0 transform: %s", p0_transform)
  logging::logdebug("p = 0 transform pseudocount: %f", p0_transform_pseudocount)
  logging::logdebug("Fixed effects: %s", fixed_effects)
  logging::logdebug("Random effects: %s", random_effects)
  logging::logdebug("Domain: %s", domain)
  logging::logdebug("Normalization: %s", ifelse(is.null(normalization), "NULL", normalization))
  logging::logdebug("Offset adjustment: %s", adjust_offset)
  logging::logdebug("Scale factor: %s", scale_factor)
  logging::logdebug("Max significance: %f", max_significance)
  logging::logdebug("Correction method: %s", correction)
  logging::logdebug("Run MaAsLin2 for tweedie_p = 0: %s", Maaslin2_run)
  logging::logdebug("Run presence-absence model: %s", run_presence_absence_model)
  logging::logdebug("Method-specific arguments provided: %s",
                    !is.null(method_args) || !is.null(method.args))
  logging::logdebug("Standardize: %s", standardize)
  logging::logdebug("Cores: %d", cores)
  logging::logdebug("Optimization routine: %s", optimizer)
  
  
  #######################################
  # Check if valid options are selected #
  #######################################
  
  # Check if the selected link is valid
  if (!link %in% link_choices) {
    option_not_valid_error("Please select a link from the list of available options",
                           toString(link_choices))
  }
  
  # Check if the selected base_model is valid
  if (!base_model %in% model_choices) {
    option_not_valid_error(
      paste(
        "Please select an analysis method",
        "from the list of available options"
      ),
      toString(model_choices)
    )
  }

  p0_transform <- gsub("[ -]", "_", toupper(p0_transform))
  if (p0_transform %in% c("ARCSIGNEDSQRT", "ARC_SIGNED_SQUARE_ROOT", "SIGNED_SQRT")) {
    p0_transform <- "ARC_SIGNED_SQRT"
  }
  if (!p0_transform %in% p0_transform_choices) {
    option_not_valid_error(
      "Please select a p0_transform from the list of available options",
      toString(p0_transform_choices)
    )
  }
  if (length(p0_transform_pseudocount) != 1L ||
      !is.numeric(p0_transform_pseudocount) ||
      !is.finite(p0_transform_pseudocount) ||
      p0_transform_pseudocount < 0) {
    stop("p0_transform_pseudocount must be a single non-negative finite numeric value.")
  }

  normalization <- resolve_tweedieverse_normalization(
    domain = domain,
    normalization = normalization,
    scale_factor = scale_factor
  )
  
  # Check if the selected correction is valid
  if (!correction %in% correction_choices) {
    option_not_valid_error(
      paste(
        "Please select a correction method",
        "from the list of available options"
      ),
      toString(correction_choices)
    )
  }
  
  # Check if the selected optimizer is valid
  if (!optimizer %in% optimizer_choices) {
    option_not_valid_error(
      paste(
        "Please select an optimizer method",
        "from the list of available options"
      ),
      toString(correction_choices)
    )
  }

  if (!is.null(tweedie_p)) {
    if (length(tweedie_p) != 1L) {
      stop("tweedie_p must be NULL, NA, or a single finite numeric value.")
    }
    if (is.na(tweedie_p)) {
      tweedie_p <- NULL
    } else if (!is.numeric(tweedie_p)) {
      stop("tweedie_p must be NULL, NA, or a single finite numeric value.")
    } else if (!is.finite(tweedie_p)) {
      stop("tweedie_p must be NULL, NA, or a single finite numeric value.")
    }
  }

  if (!is.null(tweedie_p)) {
    if (tweedie_p > 0 && tweedie_p < 1) {
      stop("Tweedie variance powers between 0 and 1 are undefined and are not supported.")
    }
  }

  if (!is.logical(run_presence_absence_model) ||
      length(run_presence_absence_model) != 1L ||
      is.na(run_presence_absence_model)) {
    stop("run_presence_absence_model must be TRUE or FALSE.")
  }

  if (!is.logical(Maaslin2_run) ||
      length(Maaslin2_run) != 1L ||
      is.na(Maaslin2_run)) {
    stop("Maaslin2_run must be TRUE or FALSE.")
  }

  if (!is.logical(adjust_offset) ||
      length(adjust_offset) != 1L ||
      is.na(adjust_offset)) {
    stop("adjust_offset must be TRUE or FALSE.")
  }

  if (!is.null(method_args) && !is.null(method.args)) {
    stop("Please provide only one of method_args or method.args.")
  }
  if (is.null(method_args)) {
    method_args <- method.args
  }
  if (!is.null(method_args) && !is.list(method_args)) {
    stop("method_args must be NULL or a named list.")
  }
  
  ############################################################
  # Check if the selected numerical options are within range #
  ############################################################
  
  prop_options <- c(prev_threshold, max_significance)
  if (any(prop_options < 0) || any(prop_options > 1)) {
    stop(
      paste(
        "One of the following is outside [0, 1]:",
        "prev_threshold, max_significance"
      )
    )
  }
  
  ###############################################################
  # Determine orientation of data in input and reorder to match #
  ###############################################################
  
  logging::loginfo("Determining format of input files")
  samples_row_row <- intersect(rownames(data), rownames(metadata))
  if (length(samples_row_row) > 0) {
    # this is the expected formatting so do not modify data frames
    logging::loginfo(paste(
      "Input format is data samples",
      "as rows and metadata samples as rows"
    ))
  } else {
    samples_column_row <- intersect(colnames(data), rownames(metadata))
    if (length(samples_column_row) > 0) {
      logging::loginfo(paste(
        "Input format is data samples",
        "as columns and metadata samples as rows"
      ))
      # transpose data frame so samples are rows
      data <- utils::type.convert(as.data.frame(t(data)), as.is = TRUE)
      logging::logdebug("linked data so samples are rows")
    } else {
      samples_column_column <-
        intersect(colnames(data), colnames(metadata))
      if (length(samples_column_column) > 0) {
        logging::loginfo(
          paste(
            "Input format is data samples",
            "as columns and metadata samples as columns"
          )
        )
        data <- utils::type.convert(as.data.frame(t(data)), as.is = TRUE)
        metadata <- utils::type.convert(as.data.frame(t(metadata)), as.is = TRUE)
        logging::logdebug("linked data and metadata so samples are rows")
      } else {
        samples_row_column <-
          intersect(rownames(data), colnames(metadata))
        if (length(samples_row_column) > 0) {
          logging::loginfo(
            paste(
              "Input format is data samples",
              "as rows and metadata samples as columns"
            )
          )
          metadata <- utils::type.convert(as.data.frame(t(metadata)), as.is = TRUE)
          logging::logdebug("linked metadata so samples are rows")
        } else {
          logging::logerror(
            paste(
              "Unable to find samples in data and",
              "metadata files.",
              "Rows/columns do not match."
            )
          )
          logging::logdebug("input_features rows: %s",
                            paste(rownames(data), collapse = ","))
          logging::logdebug("input_features columns: %s",
                            paste(colnames(data), collapse = ","))
          logging::logdebug("Metadata rows: %s",
                            paste(rownames(metadata), collapse = ","))
          logging::logdebug("Metadata columns: %s",
                            paste(colnames(data), collapse = ","))
          stop()
        }
      }
    }
  }
  
  # Replace unexpected characters in feature names
  # colnames(data) <- make.names(colnames(data))
  
  # Check for samples without metadata
  extra_feature_samples <-
    setdiff(rownames(data), rownames(metadata))
  if (length(extra_feature_samples) > 0)
    logging::logdebug(
      paste(
        "The following samples were found",
        "to have features but no metadata.",
        "They will be removed. %s"
      ),
      paste(extra_feature_samples, collapse = ",")
    )
  
  # Check for metadata samples without features
  extra_metadata_samples <-
    setdiff(rownames(metadata), rownames(data))
  if (length(extra_metadata_samples) > 0)
    logging::logdebug(
      paste(
        "The following samples were found",
        "to have metadata but no features.",
        "They will be removed. %s"
      ),
      paste(extra_metadata_samples, collapse = ",")
    )
  
  # Get a set of the samples with both metadata and features
  intersect_samples <- intersect(rownames(data), rownames(metadata))
  logging::logdebug(
    "A total of %s samples were found in both the data and metadata",
    length(intersect_samples)
  )
  
  # Now order both data and metadata with the same sample ordering
  logging::logdebug("Reordering data/metadata to use same sample ordering")
  data <- data[intersect_samples, , drop = FALSE]
  metadata <- metadata[intersect_samples, , drop = FALSE]

  fixed_effects <- parse_effect_names(fixed_effects)
  random_effects <- parse_effect_names(random_effects)

  ########################################################################
  # Assign reference values to categorical metadata (fixed effects only) #
  ########################################################################
  
  if (is.null(reference)) {
    reference <- ","
  }
  split_reference <- unlist(strsplit(reference, "[,;]"))
  
  # for each fixed effect, check that a reference level has been set if necessary: number of levels > 2 and metadata isn't already an ordered factor
  for (i in fixed_effects) {
    # don't check for or require reference levels for numeric metadata
    if (is.numeric(metadata[,i])) {
      next
    }
    # respect ordering if a factor is explicitly passed in with no reference set
    if (is.factor(metadata[,i]) && !(i %in% split_reference)) {
      logging::loginfo(paste("Factor detected for categorial metadata '", 
                             i, "'. Provide a reference argument or manually set factor ordering to change reference level.", sep=""))
      next
    }
    
    # set metadata as a factor (ordered alphabetically)
    metadata[,i] <- as.factor(metadata[,i])
    mlevels <- levels(metadata[,i])
    
    # get reference level for variable being considered, returns NA if not found
    ref <- split_reference[match(i, split_reference)+1]
    
    # if metadata has 2 levels, allow but don't require setting reference level, otherwise require it
    if ((length(mlevels) == 2)) {
      if(!is.na(ref)) {
        metadata[, i] <- stats::relevel(metadata[, i], ref = ref)
      }
    } else if (length(mlevels) > 2) {
      if (!is.na(ref)) {
        metadata[, i] <- stats::relevel(metadata[, i], ref = ref)
      } else {
        stop(paste("Please provide the reference for the variable '",
                   i, "' which includes more than 2 levels: ",
                   paste(as.character(mlevels), collapse=", "), ".", sep=""))   
      } 
    } else {
      stop("Provided categorical metadata has fewer than 2 unique, non-NA values.")
    }
  }
  
  
  #########################################################
  # Non-specific filtering based on user-provided options #
  #########################################################
  
  unfiltered_data <- data
  unfiltered_metadata <- metadata
  
  # require at least total samples * min prevalence values
  # for each feature to be greater than min abundance
  logging::loginfo("Filter data based on min abundance and min prevalence")
  total_samples <- nrow(unfiltered_data)
  logging::loginfo("Total samples in data: %d", total_samples)
  min_samples <- total_samples * prev_threshold
  logging::loginfo(
    paste(
      "Min samples required with min abundance",
      "for a feature not to be filtered: %f"
    ),
    min_samples
  )
  
  # Filter by abundance using zero as value for NAs
  data_zeros <- unfiltered_data
  data_zeros[is.na(data_zeros)] <- 0
  filtered_data <-
    unfiltered_data[,
                    colSums(data_zeros > abd_threshold) > min_samples,
                    drop = FALSE]
  total_filtered_features <-
    ncol(unfiltered_data) - ncol(filtered_data)
  logging::loginfo(
    "Total filtered features with prevalence-abundance filtering: %d",
    total_filtered_features
  )
  filtered_feature_names <-
    setdiff(names(unfiltered_data), names(filtered_data))
  logging::loginfo("Filtered feature names: %s",
                   toString(filtered_feature_names))
  
  
  
  #################################
  # Filter data based on variance #
  #################################
  
  sds <- apply(filtered_data, 2, na.rm = TRUE, sd)
  final_features <-
    filtered_data[, which(sds > var_threshold), drop = FALSE]
  total_filtered_features_var <-
    ncol(filtered_data) - ncol(final_features)
  logging::loginfo("Total filtered features with variance filtering: %d",
                   total_filtered_features_var)
  filtered_feature_names_var <-
    setdiff(names(filtered_data), names(final_features))
  logging::loginfo("Filtered feature names: %s",
                   toString(filtered_feature_names_var))

  ############################################################
  # Resolve Tweedie index after filtering the feature matrix #
  ############################################################

  has_negative_data <- any(as.matrix(final_features) < 0, na.rm = TRUE)
  if (is.null(tweedie_p) && has_negative_data) {
    message(
      paste(
        "Negative feature values detected after filtering;",
        "p = 0 is the only Tweedie index whose support includes negative values.",
        "Setting tweedie_p = 0."
      )
    )
    tweedie_p <- 0
  } else if (!is.null(tweedie_p) && tweedie_p != 0 && has_negative_data) {
    stop(
      paste(
        "Negative feature values are supported only when tweedie_p = 0.",
        "Set tweedie_p = 0, leave tweedie_p unspecified, or transform the negatives away."
      )
    )
  }

  if (!is.null(tweedie_p) && tweedie_p < 0) {
    warning(
      paste(
        "Negative Tweedie variance powers are rarely used.",
        "The model will run at the supplied tweedie_p, but requires strictly positive fitted means."
      ),
      call. = FALSE
    )
  }

  if (!is.null(tweedie_p) && tweedie_p == 0 && link != "identity") {
    message(
      paste(
        "tweedie_p = 0 uses the Gaussian variance case;",
        "setting link = 'identity' because the requested link is undefined for negative responses."
      )
    )
    link <- "identity"
  }

  if (!is.null(tweedie_p) && tweedie_p == 0 && adjust_offset) {
    message(
      paste(
        "tweedie_p = 0 with identity link does not use a log(scale_factor) offset;",
        "setting adjust_offset = FALSE."
      )
    )
    adjust_offset <- FALSE
  }
  
  
  #############################################################
  # Compute size factors for the model offset without touching #
  # the feature table used as the response. ###################
  #############################################################

  offset <- NULL
  if (adjust_offset) {
    offset <- compute_tweedieverse_size_factor(
      features = unfiltered_data,
      metadata = unfiltered_metadata,
      normalization = normalization,
      scale_factor = scale_factor
    )
    logging::loginfo("Offset size factors computed with normalization: %s", normalization)
  } else {
    logging::loginfo("Offset adjustment disabled; feature data remain on the input scale.")
  }

  if (!is.null(scale_factor) && scale_factor %in% colnames(unfiltered_metadata)) {
    unfiltered_metadata <- dplyr::select(unfiltered_metadata, -scale_factor)
  }
  
  
  ####################################
  # Filter metadata based on entropy #
  ####################################
  
  # Reduce metadata to only include those pass entropy threshold
  temp_filtered_metadata <- unfiltered_metadata[, apply(unfiltered_metadata, 2, entropy) > entropy_threshold, drop = FALSE]
  excluded_metadata <- setdiff(colnames(unfiltered_metadata), colnames(temp_filtered_metadata))
  logging::loginfo(
    paste(
      "Excluded metadata with",
      "entropy less or equal to %s: %s"
    ),
    entropy_threshold, paste(excluded_metadata, collapse = ",")
  )
  filtered_metadata <- temp_filtered_metadata
  
  
  ###############################################
  # Compute the formula based on the user input #
  ###############################################

  #####################
  # Determine formula #
  #####################
  
  random_effects_formula <- NULL
  # Use all metadata if no fixed effects are provided
  if (is.null(fixed_effects)) {
    fixed_effects <- colnames(filtered_metadata)
  } else {
    # remove any fixed effects not found in metadata names
    to_remove <- setdiff(fixed_effects, colnames(filtered_metadata))
    if (length(to_remove) > 0)
      logging::logwarn(
        paste(
          "Feature name not found in metadata",
          "so not applied to formula as fixed effect: %s"
        ),
        paste(to_remove, collapse = " , ")
      )
    fixed_effects <- setdiff(fixed_effects, to_remove)
    if (length(fixed_effects) == 0) {
      logging::logerror("No fixed effects included in formula.")
      stop()
    }
  }
  
  if (!is.null(random_effects)) {
    # subtract random effects from fixed effects
    fixed_effects <- setdiff(fixed_effects, random_effects)
    # remove any random effects not found in metadata
    to_remove <-
      setdiff(random_effects, colnames(filtered_metadata))
    if (length(to_remove) > 0)
      logging::logwarn(
        paste(
          "Feature name not found in metadata",
          "so not applied to formula as random effect: %s"
        ),
        paste(to_remove, collapse = " , ")
      )
    random_effects <- setdiff(random_effects, to_remove)
    
    # create formula
    if (length(random_effects) > 0) {
      random_effects_formula_text <-
        paste("expr ~ (1 | ",
              paste(
                random_effects,
                ")",
                sep = '',
                collapse = " + (1 | "
              ),
              sep = '')
      logging::loginfo("Formula for random effects: %s",
                       random_effects_formula_text)
      random_effects_formula <-
        tryCatch(
          as.formula(random_effects_formula_text),
          error = function(e)
            stop(
              paste(
                "Invalid formula for random effects: ",
                random_effects_formula_text
              )
            )
        )
    }
  }
  
  # Reduce metadata to only include fixed/random effects in formula
  effects_names <- union(fixed_effects, random_effects)
  filtered_metadata <-
    filtered_metadata[, effects_names, drop = FALSE]
  
  # Create the fixed effects formula text
  formula_text <-
    paste("expr ~ ", paste(fixed_effects, collapse = " + "))
  logging::loginfo("Formula for fixed effects: %s", formula_text)
  formula <-
    tryCatch(
      as.formula(formula_text),
      error = function(e)
        stop(
          paste(
            "Invalid formula.",
            "Please provide a different formula: ",
            formula_text
          )
        )
    )

  #############################################################
  # Standardize metadata (excpet the offset variable), if set #
  #############################################################
  
  if (standardize) {
    logging::loginfo("Applying z-score to standardize continuous metadata")
    filtered_metadata <-
      filtered_metadata %>% dplyr::mutate_if(is.numeric, scale)
  } else {
    logging::loginfo("Bypass z-score application to metadata")
  }
  
  ##################################
  # Merge metadata and offset back #
  ##################################
  
  final_metadata <- as.data.frame(filtered_metadata)
  if (adjust_offset) {
    final_metadata$offset <- offset
  }

  analysis_features <- final_features
  transformed_features <- NULL
  if (!is.null(tweedie_p) && tweedie_p == 0) {
    analysis_features <- p0_transform_features(
      features = final_features,
      transform = p0_transform,
      pseudocount = p0_transform_pseudocount
    )
    if (p0_transform != "NONE") {
      transformed_features <- analysis_features
      logging::loginfo("Using p = 0 transformed feature copy with transform: %s", p0_transform)
      if (!no_output) {
        transformed_features_file <- file.path(output, "transformed_features.tsv")
        logging::loginfo("Writing p = 0 transformed feature copy to file: %s", transformed_features_file)
        write.table(
          transformed_features,
          file = transformed_features_file,
          sep = "\t",
          quote = FALSE,
          row.names = TRUE
        )
      }
    }
  } else if (p0_transform != "NONE") {
    message("p0_transform is only applied when the resolved tweedie_p is 0; ignoring p0_transform.")
  }
  
  ##############################################################
  # Apply the base model to the filtered data with user inputs #
  ##############################################################
  
  logging::loginfo("Running selected analysis method: %s", base_model)
  
  fit_data <- fit.Tweedieverse(
    features = analysis_features,
    metadata = final_metadata,
    base_model = base_model,
    link = link,
    tweedie_p = tweedie_p,
    Maaslin2_run = Maaslin2_run,
    method_args = method_args,
    formula = formula,
    random_effects_formula = random_effects_formula,
    adjust_offset = adjust_offset,
    correction = correction,
    cores = cores,
    optimizer = optimizer,
    na.action = na.action
  )
  
  
  ###################################################
  # Count the N and Zero-inflation for each feature #
  ###################################################
  
  logging::loginfo("Counting prevalence for each feature")
  try(fit_data$results$N <-
        apply(
          fit_data$results,
          1,
          FUN = function(x)
            length(final_features[, x[1]])
        ))
  try(fit_data$results$N.not.zero <-
        apply(
          fit_data$results,
          1,
          FUN = function(x)
            length(which(final_features[, x[1]] > 0))
        ))
  try(fit_data$results$percent.zero <-
        apply(
          fit_data$results,
          1,
          FUN = function(x)
            round(mean(final_features[, x[1]] == 0, na.rm = TRUE), 2) *
            100
        ))
  
  #########################
  # Write out the results #
  #########################
  
  ordered_results <-
    fit_data$results[order(fit_data$results$qval),]
  ordered_results <-
    ordered_results[!is.na(ordered_results$qval),] # Remove NA's
  
  if (median_comparison) {

    mc_input <- ordered_results
    names(mc_input)[names(mc_input) == "feature"] <- "taxon"
    names(mc_input)[names(mc_input) == "coef"] <- "effect_size"
    
    mc_out <- median_comparison_tweedie(mc_input,
                                        p_cutoff = 0.95,  # ignore p>=0.95
                                        subtract_median = median_subtraction,
                                        n_sims = 10000,
                                        median_threshold = 0)
    
    ## Replace the classical columns with median-based ones
    ordered_results$coef <- mc_out$coef_median
    ordered_results$pval <- mc_out$pval_median
    ordered_results$qval <- p.adjust(mc_out$pval_median,
                                     method = correction)
    
    ordered_results <- ordered_results[order(ordered_results$qval),]
    rownames(ordered_results) <- NULL 
  }

  if (run_presence_absence_model) {
    presence_results <- fit_presence_absence_model(
      features = final_features,
      metadata = final_metadata,
      formula = formula,
      random_effects_formula = random_effects_formula,
      correction = correction,
      cores = cores
    )
    ordered_results <- combine_abundance_presence_results(
      abundance_results = ordered_results,
      presence_results = presence_results,
      correction = correction
    )
    ordered_results <- ordered_results[order(ordered_results$qval),]
    rownames(ordered_results) <- NULL
  }
  
  ordered_results <-
    dplyr::select(
      ordered_results,
      c(
        'feature',
        'metadata',
        'value',
        "coef",
        "stderr",
        "pval",
        "qval"
      ),
      everything()
    )
  
  
  if (!no_output) {
    results_file <- file.path(output, "all_results.tsv")
    logging::loginfo("Writing all results to file (ordered by increasing q-values): %s",
                     results_file)
    write.table(
      ordered_results,
      file = results_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  # Write results passing threshold to file
  # (removing any that are NA for the q-value)
  significant_results <-
    ordered_results[ordered_results$qval <= max_significance,]
  
  if (!no_output) {
  significant_results_file <-
    file.path(output, "significant_results.tsv")
  logging::loginfo(
    paste(
      "Writing the significant results",
      "(those which are less than or equal to the threshold",
      "of %f ) to file (ordered by increasing q-values): %s"
    ),
    max_significance,
    significant_results_file
  )
  write.table(
    significant_results,
    file = significant_results_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  #######################################################
  # Create visualizations for results passing threshold #
  #######################################################
  
  logging::loginfo("Writing Tweedie inxed plot to file: %s",
                   output)
  tryCatch({
    tweedie_index_plot(ordered_results, figures_folder)
    
  }, error = function(err) {
    logging::logerror("Unable to do make a Tweedie inxed plot of results!!!")
    logging::logerror(err)
    # dev.off()
  })
  
  if (plot_heatmap &&
      nrow(significant_results) > 0 &&
      length(unique(significant_results$metadata)) > 1) {
    heatmap_file <- file.path(output, "Tweedieverse_Heatmap.pdf")
    logging::loginfo("Writing heatmap of significant results to file: %s",
                     heatmap_file)
    tryCatch({
      save_heatmap(significant_results_file,
                   heatmap_file,
                   figures_folder,
                   first_n = heatmap_first_n)
    }, error = function(err) {
      logging::logerror("Unable to do make a hetamp of results!!!")
      logging::logerror(err)
      # dev.off()
    })
  }
  
  if (plot_scatter) {
    logging::loginfo(
      paste(
        "Writing association plots",
        "(one for each significant association)",
        "to output folder: %s"
      ),
      output
    )
    association_plots(
      metadata,
      final_features/offset,
      significant_results_file,
      output,
      figures_folder
    )
  }
  }
  return(significant_results)
}


option_not_valid_error <- function(message, valid_options) {
  logging::logerror(paste(message, ": %s"), toString(valid_options))
  stop("Option not valid", call. = FALSE)
}

## Quiets concerns of R CMD check
utils::globalVariables(c(
  "base.model",
  "base.model_abundance",
  "base.model_presence",
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
