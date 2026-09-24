get_AICtab<-function(fit){
  
  ########################
  # Flag invalid options #
  ########################
  
  if (!inherits(fit, c("cpglm", "glmmTMB"))) {
    stop('Not supported. Valid options are cpglm and glmmTMB')
  }
  
  ######################
  #  Initialize AICtab #
  ######################
  
  AICtab<-rep(NA, 5)
  
  ###########################
  # Case-by-Case Extraction #
  ###########################
  
  if (inherits(fit, "cpglm")) {
    
    ##########################################
    # Back calculate logLik and BIC from AIC #
    ##########################################
    
    AIC<-fit$aic
    AIC_multiplier<-length(fit$y) - fit$df.residual
    logLik<-(AIC - 2*AIC_multiplier)/2
    BIC_multiplier<-AIC_multiplier*log(length(fit$y))
    BIC<-BIC_multiplier + 2*logLik
    deviance<-fit$deviance
    df.resid<-fit$df.residual
    
    # Coherent output
    AICtab<-c(AIC, BIC, logLik, deviance, df.resid)
  }
  
  if (inherits(fit, "glmmTMB")) {
    
    #######################################
    # Extract AICtab from glmmTMB objects #
    #######################################
    
    AICtab<-summary(fit)["AICtab"]$AICtab
    
  }
  
  ##########
  # Return #
  ##########
  
  names(AICtab)<-c('AIC', 'BIC', 'logLik', 'deviance', 'df.resid')
  return(AICtab)
}




# Adapted form: https://rstudio-pubs-static.s3.amazonaws.com/455435_30729e265f7a4d049400d03a18e218db.html

#' Entropy of a Vector
#'
#' Compute Shannon entropy for a vector.
#'
#' @param target A vector of values.
#' @return A numeric entropy value.
#' @examples
#' entropy(c("A", "A", "B", "B", "C"))
#' @export
entropy <- function(target) {
  #if(all(is.na(target)))  0 
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}

IG_numeric<-function(data, feature, target, bins=4) {
  #Strip out rows where feature is NA
  data<-data[!is.na(data[,feature]),]
  #compute entropy for the parent
  e0<-entropy(data[,target])
  
  data$cat<-cut(data[,feature], breaks = bins, labels = seq_len(bins))
  
  #use dplyr to compute e and p for each value of the feature
  dd_data <- data %>% dplyr::group_by(cat) %>% dplyr::summarise(
    e = entropy(get(target)),
    n = length(get(target)),
    min = min(get(feature)),
    max = max(get(feature))
  )
  
  #calculate p for each value of feature
  dd_data$p<-dd_data$n/nrow(data)
  #compute IG
  IG<-e0-sum(dd_data$p*dd_data$e)
  
  return(IG)
}



#returns IG for categorical variables.
IG_cat<-function(data,feature,target){
  #Strip out rows where feature is NA
  data<-data[!is.na(data[,feature]),] 
  #use dplyr to compute e and p for each value of the feature
  dd_data <- data %>% dplyr::group_by_at(feature) %>% dplyr::summarise(
    e = entropy(get(target)),
    n = length(get(target))
  )
  
  #compute entropy for the parent
  e0<-entropy(data[,target])
  #calculate p for each value of feature
  dd_data$p<-dd_data$n/nrow(data)
  #compute IG
  IG<-e0-sum(dd_data$p*dd_data$e)
  
  return(IG)
}

# entropy (c("A", "A", "A", "A", "A", "B", "B"))
# 0.8631206

#entropy (c("A", "A", "A", "A"))
# 0

#entropy (c("A", "A", "A", "A", "B", "B", "B", "B"))
#1

#entropy (c("C", "A", "A", "A", "B", "B", "B", "B"))
# 1.405639

#entropy (c("C", "A", "D", "A", "B", "B", "B", "B"))
# 1.75

# entropy (c(1, 1, 2, 1, 1, 1, 2, 1))
#0.8112781


# Written by Grace
read_input_table <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    comment.char = "",
    check.names = FALSE,
    row.names = 1
  )
}

merge_method_args <- function(defaults, overrides) {
  if (is.null(overrides)) {
    return(defaults)
  }
  if (!is.list(overrides)) {
    stop("Method-specific arguments must be provided as a list.")
  }
  utils::modifyList(defaults, overrides, keep.null = TRUE)
}

extract_method_args <- function(method_args, method) {
  if (is.null(method_args)) {
    return(list())
  }
  if (!is.list(method_args)) {
    stop("method_args must be NULL or a named list.")
  }

  candidates <- c(method, tolower(method))
  for (candidate in candidates) {
    if (!is.null(method_args[[candidate]])) {
      if (!is.list(method_args[[candidate]])) {
        stop(sprintf("method_args$%s must be a list.", candidate))
      }
      return(method_args[[candidate]])
    }
  }

  list()
}

p0_transform_features <- function(features,
                                  transform = "NONE",
                                  pseudocount = 1) {
  transform <- toupper(transform)
  if (transform == "NONE") {
    return(as.data.frame(features))
  }

  x <- as.matrix(features)
  storage.mode(x) <- "numeric"
  if (any(!is.finite(x), na.rm = TRUE)) {
    stop("p0_transform requires finite feature values.")
  }

  if (transform == "CLR") {
    if (any(x < 0, na.rm = TRUE)) {
      stop("CLR transformation requires non-negative feature values.")
    }
    if (pseudocount <= 0) {
      stop("CLR transformation requires p0_transform_pseudocount > 0.")
    }
    z <- log(x + pseudocount)
    z <- sweep(z, 1, rowMeans(z, na.rm = TRUE), "-")
  } else if (transform == "RCLR") {
    if (any(x < 0, na.rm = TRUE)) {
      stop("RCLR transformation requires non-negative feature values.")
    }
    z <- t(apply(x, 1, function(row) {
      positive <- row > 0
      if (!any(positive)) {
        stop("RCLR transformation requires at least one positive feature per sample.")
      }
      out <- rep(0, length(row))
      logged <- log(row[positive])
      out[positive] <- logged - mean(logged)
      out
    }))
  } else if (transform == "LOG") {
    if (any(x + pseudocount <= 0, na.rm = TRUE)) {
      stop("LOG transformation requires feature values + p0_transform_pseudocount to be strictly positive.")
    }
    z <- log(x + pseudocount)
  } else if (transform == "ARC_SIGNED_SQRT") {
    z <- sign(x) * sqrt(abs(x))
  } else {
    stop(paste("Unsupported p0_transform:", transform))
  }

  z <- as.data.frame(z)
  rownames(z) <- rownames(features)
  colnames(z) <- colnames(features)
  z
}

domain_bioc_container_map <- function() {
  list(
    microbiome = c("TreeSummarizedExperiment", "SummarizedExperiment"),
    single_cell = c("SingleCellExperiment"),
    bulk_rnaseq = c("SummarizedExperiment"),
    custom = c(
      "SummarizedExperiment",
      "SingleCellExperiment",
      "TreeSummarizedExperiment"
    )
  )
}

validate_domain_bioc_container <- function(input, domain) {
  allowed <- domain_bioc_container_map()[[domain]]
  if (is.null(allowed)) {
    stop(paste("Unknown domain:", domain))
  }

  if (!inherits(input, allowed)) {
    stop(
      paste(
        "For domain =", domain,
        "Bioconductor input must inherit from one of:",
        paste(allowed, collapse = ", "),
        ". Use a data.frame/file path for domain-agnostic input."
      )
    )
  }
  invisible(TRUE)
}

extractAssay <- function(input, assay_name = "counts") {
  
  # Extract assay name based on the user input
  if (assay_name %in% SummarizedExperiment::assayNames(input)) {
    counts_data <- SummarizedExperiment::assay(input, assay_name)
    cat("The specified assay has been extracted\n")
    return(as.data.frame(as.matrix(counts_data)))
  } else {
    cat("The specified assay was not found\n")
    return(NULL)
  }
}

is_multiassay_experiment <- function(input) {
  inherits(input, "MultiAssayExperiment")
}

extract_multiassay_experiments <- function(input) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop(
      paste(
        "MultiAssayExperiment input requires the MultiAssayExperiment package.",
        "Please install it or provide a single-omics input."
      )
    )
  }
  experiments <- MultiAssayExperiment::experiments(input)
  if (length(experiments) == 0L) {
    stop("MultiAssayExperiment input must contain at least one experiment.")
  }
  experiments
}

extract_multiassay_metadata <- function(input, input_metadata = NULL) {
  if (!is.null(input_metadata)) {
    if (is.character(input_metadata)) {
      return(read_input_table(input_metadata))
    }
    return(input_metadata)
  }
  metadata <- as.data.frame(SummarizedExperiment::colData(input))
  if (nrow(metadata) == 0L) {
    stop(
      paste(
        "MultiAssayExperiment input needs sample metadata.",
        "Provide input_metadata or populate colData(input_features)."
      )
    )
  }
  metadata
}

coerce_multiassay_experiment_input <- function(experiment) {
  valid_classes <- unique(unlist(domain_bioc_container_map()))
  if (inherits(experiment, valid_classes) || is.data.frame(experiment)) {
    return(experiment)
  }
  if (is.matrix(experiment)) {
    return(as.data.frame(experiment))
  }
  stop(
    sprintf(
      paste(
        "MultiAssayExperiment experiment of class <%s> is not supported.",
        "Use SummarizedExperiment-like experiments, matrices, or data.frames."
      ),
      class(experiment)[1]
    )
  )
}

resolve_multiomics_arg <- function(value,
                                   omics_name,
                                   omics_names,
                                   arg_name,
                                   allow_null_element = TRUE) {
  if (is.null(value)) {
    return(NULL)
  }

  if (is.list(value)) {
    value_names <- names(value)
    if (!is.null(value_names) && omics_name %in% value_names) {
      return(value[[omics_name]])
    }
    if (!is.null(value_names) && any(value_names %in% omics_names)) {
      if (allow_null_element) {
        return(NULL)
      }
      stop(sprintf("No value supplied for %s in omics layer %s.", arg_name, omics_name))
    }
    if (is.null(value_names) && length(value) == length(omics_names)) {
      return(value[[match(omics_name, omics_names)]])
    }
    return(value)
  }

  value_names <- names(value)
  if (!is.null(value_names) && omics_name %in% value_names) {
    return(unname(value[[omics_name]]))
  }
  if (!is.null(value_names) && any(value_names %in% omics_names)) {
    stop(sprintf("No value supplied for %s in omics layer %s.", arg_name, omics_name))
  }

  value
}

run_multiassay_tweedieverse <- function(input_features,
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
  if (!is.null(method_args) && !is.null(method.args)) {
    stop("Please provide only one of method_args or method.args.")
  }
  if (is.null(method_args)) {
    method_args <- method.args
  }

  experiments <- extract_multiassay_experiments(input_features)
  omics_names <- names(experiments)
  if (is.null(omics_names) || any(omics_names == "")) {
    omics_names <- paste0("omics", seq_along(experiments))
    names(experiments) <- omics_names
  }
  metadata <- extract_multiassay_metadata(input_features, input_metadata)

  results <- lapply(omics_names, function(omics_name) {
    experiment <- coerce_multiassay_experiment_input(experiments[[omics_name]])
    omics_output <- NULL
    if (!is.null(output)) {
      omics_output <- file.path(output, make.names(omics_name))
    }

    fit <- .Tweedieverse_single(
      input_features = experiment,
      input_metadata = metadata,
      output = omics_output,
      assay_name = resolve_multiomics_arg(assay_name, omics_name, omics_names, "assay_name"),
      abd_threshold = resolve_multiomics_arg(abd_threshold, omics_name, omics_names, "abd_threshold"),
      prev_threshold = resolve_multiomics_arg(prev_threshold, omics_name, omics_names, "prev_threshold"),
      var_threshold = resolve_multiomics_arg(var_threshold, omics_name, omics_names, "var_threshold"),
      entropy_threshold = resolve_multiomics_arg(entropy_threshold, omics_name, omics_names, "entropy_threshold"),
      base_model = resolve_multiomics_arg(base_model, omics_name, omics_names, "base_model"),
      link = resolve_multiomics_arg(link, omics_name, omics_names, "link"),
      tweedie_p = resolve_multiomics_arg(tweedie_p, omics_name, omics_names, "tweedie_p"),
      p0_transform = resolve_multiomics_arg(p0_transform, omics_name, omics_names, "p0_transform"),
      p0_transform_pseudocount = resolve_multiomics_arg(p0_transform_pseudocount, omics_name, omics_names, "p0_transform_pseudocount"),
      fixed_effects = resolve_multiomics_arg(fixed_effects, omics_name, omics_names, "fixed_effects"),
      random_effects = resolve_multiomics_arg(random_effects, omics_name, omics_names, "random_effects"),
      domain = resolve_multiomics_arg(domain, omics_name, omics_names, "domain"),
      normalization = resolve_multiomics_arg(normalization, omics_name, omics_names, "normalization"),
      adjust_offset = resolve_multiomics_arg(adjust_offset, omics_name, omics_names, "adjust_offset"),
      scale_factor = resolve_multiomics_arg(scale_factor, omics_name, omics_names, "scale_factor"),
      max_significance = resolve_multiomics_arg(max_significance, omics_name, omics_names, "max_significance"),
      correction = resolve_multiomics_arg(correction, omics_name, omics_names, "correction"),
      Maaslin2_run = resolve_multiomics_arg(Maaslin2_run, omics_name, omics_names, "Maaslin2_run"),
      median_comparison = resolve_multiomics_arg(median_comparison, omics_name, omics_names, "median_comparison"),
      median_subtraction = resolve_multiomics_arg(median_subtraction, omics_name, omics_names, "median_subtraction"),
      run_presence_absence_model = resolve_multiomics_arg(run_presence_absence_model, omics_name, omics_names, "run_presence_absence_model"),
      method_args = resolve_multiomics_arg(method_args, omics_name, omics_names, "method_args"),
      standardize = resolve_multiomics_arg(standardize, omics_name, omics_names, "standardize"),
      cores = resolve_multiomics_arg(cores, omics_name, omics_names, "cores"),
      optimizer = resolve_multiomics_arg(optimizer, omics_name, omics_names, "optimizer"),
      na.action = na.action,
      plot_heatmap = resolve_multiomics_arg(plot_heatmap, omics_name, omics_names, "plot_heatmap"),
      plot_scatter = resolve_multiomics_arg(plot_scatter, omics_name, omics_names, "plot_scatter"),
      heatmap_first_n = resolve_multiomics_arg(heatmap_first_n, omics_name, omics_names, "heatmap_first_n"),
      reference = resolve_multiomics_arg(reference, omics_name, omics_names, "reference")
    )
    fit$omics <- omics_name
    fit <- dplyr::select(fit, dplyr::all_of("omics"), dplyr::everything())
    fit
  })
  names(results) <- omics_names
  class(results) <- c("TweedieverseMultiAssayResult", class(results))
  results
}

resolve_tweedieverse_normalization <- function(domain,
                                               normalization = NULL,
                                               scale_factor = NULL) {
  if (!is.null(scale_factor)) {
    return("USER")
  }

  defaults <- c(
    microbiome = "TSS",
    single_cell = "SCRAN",
    bulk_rnaseq = "TMM",
    custom = "TSS"
  )

  if (is.null(normalization)) {
    normalization <- defaults[[domain]]
  }

  normalization <- toupper(normalization)
  if (normalization == "DESEQ2") {
    normalization <- "RLE"
  }

  all_choices <- c("TSS", "GMPR", "CSS", "SCRAN", "TMM", "RLE", "CPM", "MEDIAN", "NONE")
  if (!normalization %in% all_choices) {
    stop(
      paste(
        "normalization must be one of:",
        paste(all_choices, collapse = ", ")
      )
    )
  }

  domain_choices <- switch(
    domain,
    microbiome = c("TSS", "GMPR", "CSS", "MEDIAN", "NONE"),
    single_cell = c("SCRAN", "MEDIAN", "NONE"),
    bulk_rnaseq = c("TMM", "RLE", "CPM", "MEDIAN", "NONE"),
    custom = all_choices
  )
  if (!normalization %in% domain_choices) {
    stop(
      paste(
        "normalization =", normalization,
        "is not supported for domain =", domain,
        ". Supported choices are:",
        paste(domain_choices, collapse = ", ")
      )
    )
  }

  normalization
}

validate_tweedieverse_size_factor <- function(size_factor, label = "size factor") {
  size_factor <- as.numeric(size_factor)
  if (any(!is.finite(size_factor)) || any(is.na(size_factor))) {
    stop(paste(label, "must contain only finite numeric values."))
  }
  if (any(size_factor <= 0)) {
    stop(paste(label, "must be strictly positive because Tweedieverse uses log(size_factor) as an offset."))
  }
  size_factor
}

center_tweedieverse_size_factor <- function(size_factor) {
  size_factor <- validate_tweedieverse_size_factor(size_factor)
  gm <- exp(mean(log(size_factor)))
  size_factor / gm
}

compute_tweedieverse_size_factor <- function(features,
                                             metadata,
                                             normalization,
                                             scale_factor = NULL) {
  if (!is.null(scale_factor)) {
    if (!scale_factor %in% colnames(metadata)) {
      stop(
        paste(
          "The specified scale_factor variable is not present in the metadata table:\n",
          scale_factor
        )
      )
    }
    return(validate_tweedieverse_size_factor(metadata[, scale_factor], "scale_factor"))
  }

  normalization <- toupper(normalization)
  if (normalization == "NONE") {
    return(rep(1, nrow(features)))
  }

  x <- as.matrix(features)
  storage.mode(x) <- "numeric"

  if (any(!is.finite(x), na.rm = TRUE)) {
    stop("Features must contain only finite values for size-factor normalization.")
  }

  if (normalization %in% c("TSS", "CPM")) {
    if (any(x < 0, na.rm = TRUE)) {
      stop(paste(normalization, "size factors require non-negative features."))
    }
    return(center_tweedieverse_size_factor(rowSums(x, na.rm = TRUE)))
  }

  if (normalization == "MEDIAN") {
    size_factor <- apply(x, 1, function(row) {
      vals <- row[is.finite(row) & !is.na(row) & row > 0]
      if (length(vals) == 0L) {
        return(NA_real_)
      }
      stats::median(vals)
    })
    return(center_tweedieverse_size_factor(size_factor))
  }

  if (any(x < 0, na.rm = TRUE)) {
    stop(paste(normalization, "size factors require non-negative features."))
  }

  if (normalization == "RLE") {
    return(compute_rle_size_factor(x))
  }

  if (normalization == "GMPR") {
    return(compute_gmpr_size_factor(x))
  }

  if (normalization == "CSS") {
    return(compute_css_size_factor(x))
  }

  if (normalization == "TMM") {
    return(compute_tmm_size_factor(x))
  }

  if (normalization == "SCRAN") {
    return(compute_scran_size_factor(x))
  }

  stop(paste("Unsupported normalization:", normalization))
}

compute_rle_size_factor <- function(x) {
  geom_means <- apply(x, 2, function(col) {
    vals <- col[col > 0]
    if (length(vals) == 0L) {
      return(NA_real_)
    }
    exp(mean(log(vals)))
  })
  keep <- is.finite(geom_means) & geom_means > 0
  if (!any(keep)) {
    stop("RLE size factors could not be computed because no feature has positive values.")
  }
  ratios <- sweep(x[, keep, drop = FALSE], 2, geom_means[keep], "/")
  ratios[ratios <= 0] <- NA_real_
  size_factor <- apply(ratios, 1, stats::median, na.rm = TRUE)
  center_tweedieverse_size_factor(size_factor)
}

compute_gmpr_size_factor <- function(x) {
  n_samples <- nrow(x)
  ratios <- matrix(NA_real_, n_samples, n_samples)
  for (i in seq_len(n_samples)) {
    for (j in seq_len(n_samples)) {
      if (i == j) {
        next
      }
      shared <- x[i, ] > 0 & x[j, ] > 0
      if (any(shared)) {
        ratios[i, j] <- stats::median(x[i, shared] / x[j, shared])
      }
    }
  }
  size_factor <- apply(ratios, 1, function(row) {
    vals <- row[is.finite(row) & row > 0]
    if (length(vals) == 0L) {
      return(NA_real_)
    }
    exp(mean(log(vals)))
  })
  center_tweedieverse_size_factor(size_factor)
}

compute_css_size_factor <- function(x, percentile = 0.75) {
  positive_values <- x[x > 0]
  if (length(positive_values) == 0L) {
    stop("CSS size factors require at least one positive feature value.")
  }
  threshold <- as.numeric(stats::quantile(positive_values, probs = percentile, names = FALSE))
  size_factor <- rowSums(ifelse(x <= threshold, x, 0), na.rm = TRUE)
  center_tweedieverse_size_factor(size_factor)
}

compute_tmm_size_factor <- function(x) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR is required for normalization = 'TMM'. Install edgeR or choose another normalization.")
  }
  counts <- t(x)
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  size_factor <- dge$samples$lib.size * dge$samples$norm.factors
  center_tweedieverse_size_factor(size_factor)
}

compute_scran_size_factor <- function(x) {
  if (!requireNamespace("scran", quietly = TRUE) ||
      !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("scran and SingleCellExperiment are required for normalization = 'SCRAN'. Install them or choose another normalization.")
  }
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = t(x)))
  sce <- scran::computeSumFactors(sce)
  size_factor <- SingleCellExperiment::sizeFactors(sce)
  center_tweedieverse_size_factor(size_factor)
}


#' Median Comparison for Compositionality Adjustment
#'
#' Adjust Tweedieverse(or any other differential analysis methods) coefficient estimates and p-values by testing each taxon
#' against the *median* effect for the same metadata variable - a simple
#' post-hoc strategy to curb false discoveries driven by the compositional
#' nature of microbiome count data (after the approach adopted in **MaAsLin 3**).
#'
#' @param df A `data.frame` returned by **Tweedieverse** or any other differential analysis methods containing (at
#'   minimum) the columns  
#'   `taxon`, `metadata`, `effect_size`, `pval`, `stderr`, and `qval`.
#' @param p_cutoff Numeric.  Upper bound on the original p-value to include an
#'   effect in the median calculation (default `0.95` = all non-missing).
#' @param subtract_median Logical.  If `TRUE`, subtracts the group
#'   median from every coefficient before returning it in `coef_median`;
#'   otherwise the original `effect_size` is copied unchanged
#'   (default `FALSE`).
#' @param n_sims Integer.  Number of Monte-Carlo simulations used to estimate
#'   the covariance between each coefficient and the group median
#'   (default `10 000`).
#' @param median_threshold Numeric.  Absolute difference below which a
#'   coefficient is considered effectively equal to the median and assigned
#'   `pval_median = 1` (default `0`).
#'
#' @details
#' For every distinct value in `metadata` the algorithm
#' \enumerate{
#'   \item keeps coefficients with `pval < p_cutoff` and computes their median;
#'   \item simulates `n_sims` draws of coefficients using a normal
#'         approximation (`N(effect_size, stderr^2)`) and records the empirical
#'         distribution of the simulated medians;
#'   \item derives a variance-inflated *offset* that accounts for the
#'         covariance between each coefficient and the group median;
#'   \item performs a two-sided Z-test of `H0 : beta_i = offset_i`, returning the
#'         resulting p-value in `pval_median`.
#' }
#'
#' @return
#' The input `df` with two new columns:
#' \describe{
#'   \item{`coef_median`}{Median-centred coefficient (or the original
#'   `effect_size` if `subtract_median = FALSE`).}
#'   \item{`pval_median`}{Two-sided p-value from the median-comparison test.}
#' }
#'
#' @section Warning:
#' The procedure is heuristic and relies on normal approximations.  Results
#' may be unstable for very small sample sizes or when `stderr` values are
#' zero or missing.
#'
#' @seealso [MaAsLin 3 GitHub](https://github.com/biobakery/maaslin3)
#'
#' @examples
#' toy <- data.frame(
#'   taxon = c("tax1", "tax2"),
#'   metadata = c("grp", "grp"),
#'   effect_size = c(0.4, 0.2),
#'   pval = c(0.01, 0.2),
#'   stderr = c(0.1, 0.1),
#'   qval = c(0.02, 0.25)
#' )
#' median_comparison_tweedie(toy, n_sims = 100)
#'
#' \dontrun{
#' 
#' ######################
#' # HMP2 input_features Analysis #
#' ######################
#'
#' #############
#' # Load input_features #
#' #############
#' 
#' library(data.table)
#' input_features <- fread("https://raw.githubusercontent.com/biobakery/Maaslin2/master/inst/extdata/HMP2_taxonomy.tsv", sep ="\t")
#' input_metadata <-fread("https://raw.githubusercontent.com/biobakery/Maaslin2/master/inst/extdata/HMP2_metadata.tsv", sep ="\t")
#'
#' ###############
#' # Format data #
#' ###############
#'
#' library(tibble)
#' features<- column_to_rownames(input_features, 'ID')
#' metadata<- column_to_rownames(input_metadata, 'ID')
#'
#' #############
#' # Fit Model #
#' #############
#'
#' library(Tweedieverse)
#' HMP2 <- Tweedieverse(
#' features,
#' metadata,
#' output = './demo_output/HMP2', # Assuming demo_output exists
#' fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'),
#' random_effects = c('site', 'subject'),
#' base_model = 'CPLM',
#' adjust_offset = FALSE, # No offset as the values are relative abundances
#' cores = 8, # Make sure your computer has the capability
#' median_comparison = TRUE,
#' median_subtraction = TRUE,
#' standardize = FALSE,
#' reference = c('diagnosis,nonIBD'))
#' 
#' HMP2_adj <- median_comparison_tweedie(HMP2,
#'                                         p_cutoff = 0.95,
#'                                         subtract_median = TRUE,
#'                                         n_sims = 10000,
#'                                         median_threshold = 0)
#'
#' head(HMP2_adj[, c("taxon", "metadata", "coef_median", "pval_median")])
#' 
#' }
#'
#' @export
median_comparison_tweedie <- function(df,
                                      p_cutoff = 0.95,
                                      subtract_median = FALSE,
                                      n_sims = 10000,
                                      median_threshold = 0) {
  # df is your Tweedieverse output data.frame with columns:
  #   taxon, metadata, effect_size, pval, stderr, qval
  #
  # We'll store results here:
  df$pval_median <- NA_real_
  df$coef_median <- df$effect_size  # By default, same as effect_size
  
  # Process each metadata variable separately
  for (md in unique(df$metadata)) {
    # 1) Subset to just this metadata predictor
    sub_idx <- which(df$metadata == md)
    sub_df  <- df[sub_idx, ]
    
    # 2) Filter out obviously "bad" or huge p-values before computing the median
    use_idx <- which(!is.na(sub_df$pval) & sub_df$pval < p_cutoff)
    if (length(use_idx) == 0) {
      # If none are usable, move on
      next
    }
    
    # 3) Compute the "group-wide" median of the usable coefficients
    cur_median <- stats::median(sub_df$effect_size[use_idx], na.rm = TRUE)
    if (is.na(cur_median)) {
      # If no valid median, skip
      next
    }
    
    # 4) Optionally shift each coefficient by the median
    if (subtract_median) {
      sub_df$coef_median <- sub_df$effect_size - cur_median
    } else {
      sub_df$coef_median <- sub_df$effect_size
    }
    
    coefs    <- sub_df$effect_size
    ses      <- sub_df$stderr
    n_coefs  <- length(coefs)
    
    # Identify which coefficients were used for the median
    use_bool <- rep(FALSE, n_coefs)
    use_bool[use_idx] <- TRUE
    
    # 5) Simulate draws to approximate correlation of each coefficient w/ median
    #    sim_results has columns = draws, row 1 = simulated median, next rows = coefs
    sim_results <- replicate(n_sims, {
      sim_coefs   <- stats::rnorm(n_coefs, mean = coefs, sd = ses)
      sim_median  <- stats::median(sim_coefs[use_bool])
      c(sim_median, sim_coefs)
    })
    
    sim_medians <- sim_results[1, ]
    all_sims    <- sim_results[-1, , drop = FALSE]  # row i => draws for coef i
    
    # Covariance of each coefficient with the median, across draws
    cov_adjust <- apply(all_sims, 1, function(x) stats::cov(x, sim_medians))
    
    # 6) "offset to test" for each coefficient, per the MaAsLin 3 approach:
    #    offset_i = coefs[i] +/- ... depends on difference from median & correlation
    median_sd <- sd(sim_medians)  # the empirical SD of the simulated median
    offsets_to_test <- abs(cur_median - coefs) *
      sqrt( (ses^2) / ( ses^2 + median_sd^2 - 2*cov_adjust ) ) + coefs
    
    # 7) For each coefficient, finalize the p-value vs. the median
    #    - If difference from median < threshold => p=1
    #    - Else do test:  H0: coef[i] = offsets_to_test[i]
    #        =>  z = (coef - offset) / SE
    pvals_median <- numeric(n_coefs)
    
    for (i in seq_len(n_coefs)) {
      if (abs(coefs[i] - cur_median) < median_threshold) {
        # If difference is trivially small => p=1
        pvals_median[i] <- 1
      } else if (is.na(coefs[i]) || is.na(ses[i]) || ses[i] == 0) {
        pvals_median[i] <- NA
      } else {
        # Normal approx. test for H0: coefs[i] == offsets_to_test[i]
        z_stat <- (coefs[i] - offsets_to_test[i]) / ses[i]
        pvals_median[i] <- 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
      }
    }
    
    # Save the results back into the subset
    sub_df$pval_median <- pvals_median
    df[sub_idx, ]      <- sub_df
  }
  
  # Return the augmented data
  return(df)
}

tweedieverse_cct <- function(pvals, weights = NULL) {
  pvals <- ifelse(is.na(pvals), 1, pvals)

  if (any(pvals < 0 | pvals > 1)) {
    stop("All p-values must be between 0 and 1.")
  }

  if (any(pvals == 0) && any(pvals == 1)) {
    stop("Cannot combine exact 0 and exact 1 p-values.")
  }
  if (any(pvals == 0)) {
    return(0)
  }
  if (any(pvals == 1)) {
    return(1)
  }

  if (is.null(weights)) {
    weights <- rep(1 / length(pvals), length(pvals))
  } else if (length(weights) != length(pvals)) {
    stop("weights must have the same length as pvals.")
  } else if (any(weights < 0)) {
    stop("weights must be non-negative.")
  } else {
    weights <- weights / sum(weights)
  }

  is_small <- pvals < 1e-16
  if (!any(is_small)) {
    cct_stat <- sum(weights * tan((0.5 - pvals) * pi))
  } else {
    cct_stat <- sum((weights[is_small] / pvals[is_small]) / pi)
    cct_stat <- cct_stat +
      sum(weights[!is_small] * tan((0.5 - pvals[!is_small]) * pi))
  }

  if (cct_stat > 1e15) {
    return((1 / cct_stat) / pi)
  }
  1 - stats::pcauchy(cct_stat)
}

tweedieverse_cct_rows <- function(mat) {
  if (ncol(mat) == 1L) {
    return(as.numeric(mat[, 1L]))
  }
  apply(mat, 1L, function(pv) tweedieverse_cct(as.numeric(pv)))
}

presence_augmentation_weight <- function(formula, data) {
  rhs_formula <- stats::delete.response(stats::terms(formula))
  n_predictors <- ncol(stats::model.matrix(rhs_formula, data = data))
  n_predictors / (2 * nrow(data))
}

augment_presence_data <- function(formula,
                                  data,
                                  response = "expr",
                                  weights = NULL) {
  n_samples <- nrow(data)
  if (is.null(weights)) {
    weights <- rep(1, n_samples)
  }
  if (length(weights) != n_samples) {
    stop("weights must have length equal to the number of samples.")
  }

  augmentation_weight <- presence_augmentation_weight(formula, data)
  data_zero <- data
  data_one <- data
  data_zero[[response]] <- 0L
  data_one[[response]] <- 1L

  list(
    data = rbind(data, data_zero, data_one),
    weights = c(
      weights,
      rep(augmentation_weight, n_samples),
      rep(augmentation_weight, n_samples)
    )
  )
}

fit_augmented_presence_model <- function(formula,
                                         data,
                                         has_random_effects = FALSE,
                                         offset = NULL,
                                         na.action = na.exclude) {
  augmented <- augment_presence_data(formula = formula, data = data)
  model_offset <- NULL
  if (!is.null(offset)) {
    if (length(offset) != nrow(data)) {
      stop("offset must have length equal to the number of samples.")
    }
    model_offset <- rep(offset, 3L)
  }

  fit_fun <- if (has_random_effects) glmmTMB::glmmTMB else stats::glm
  args <- list(
    formula = formula,
    family = stats::binomial(),
    data = augmented$data,
    weights = augmented$weights,
    na.action = na.action
  )
  if (!is.null(model_offset)) {
    args$offset <- model_offset
  }

  withCallingHandlers(
    do.call(fit_fun, args),
    warning = function(w) {
      if (grepl("non-integer #successes in a binomial glm!",
                conditionMessage(w),
                fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

fit_presence_absence_model <- function(features,
                                       metadata,
                                       formula,
                                       random_effects_formula = NULL,
                                       correction = "BH",
                                       cores = 1,
                                       na.action = na.exclude) {
  if (!is.null(random_effects_formula) &&
      !requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("glmmTMB is required for presence-absence models with random_effects.")
  }

  presence_formula <- formula
  has_random_effects <- !is.null(random_effects_formula)
  if (has_random_effects) {
    fixed_terms <- setdiff(all.vars(formula)[-1], "offset")
    formula_text <-
      paste(". ~", paste(fixed_terms, collapse = " + "), ".", sep = " + ")
    presence_formula <- stats::update(random_effects_formula, formula_text)
  }

  log_offset <- NULL
  if ("offset" %in% colnames(metadata)) {
    log_offset <- log(metadata$offset)
  }
  metadata_names <- setdiff(colnames(metadata), "offset")

  cluster <- NULL
  if (cores > 1) {
    logging::loginfo("Creating cluster of %s R processes for presence-absence models", cores)
    cluster <- parallel::makeCluster(cores)
    parallel::clusterExport(
      cluster,
      c(
        "features",
        "metadata",
        "presence_formula",
        "has_random_effects",
        "log_offset",
        "na.action",
        "fit_augmented_presence_model",
        "augment_presence_data",
        "presence_augmentation_weight"
      ),
      envir = environment()
    )
  }

  outputs <- pbapply::pblapply(seq_len(ncol(features)), cl = cluster, function(x) {
    expr <- as.integer(features[, x] > 0)
    data_sub <- data.frame(metadata, expr = expr)

    if (length(unique(expr)) < 2L) {
      para <- as.data.frame(matrix(NA_real_, nrow = length(metadata_names), ncol = 5))
      para$name <- metadata_names
    } else {
      fit <- try(
        fit_augmented_presence_model(
          formula = presence_formula,
          data = data_sub,
          has_random_effects = has_random_effects,
          offset = log_offset,
          na.action = na.action
        ),
        silent = TRUE
      )

      if (!inherits(fit, "try-error")) {
        summary_matrix <- if (has_random_effects) {
          summary(fit)$coefficients$cond
        } else {
          stats::coef(summary(fit))
        }
        p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), colnames(summary_matrix))[1]
        if (!is.na(p_col)) {
          para <- as.data.frame(summary_matrix)[-1, c("Estimate", "Std. Error", p_col), drop = FALSE]
          para$base.model <- "Presence-absence LR"
          para$tweedie.index <- NA_real_
          para$name <- rownames(summary_matrix)[-1]
        } else {
          para <- as.data.frame(matrix(NA_real_, nrow = length(metadata_names), ncol = 5))
          para$name <- metadata_names
        }
      } else {
        para <- as.data.frame(matrix(NA_real_, nrow = length(metadata_names), ncol = 5))
        para$name <- metadata_names
      }
    }

    colnames(para) <- c("coef", "stderr", "pval", "base.model", "tweedie.index", "name")
    para$feature <- colnames(features)[x]
    para
  })

  if (!is.null(cluster)) {
    parallel::stopCluster(cluster)
  }

  paras <- do.call(rbind, outputs)
  paras$qval <- as.numeric(stats::p.adjust(paras$pval, method = correction))

  metadata_names_ordered <- metadata_names[order(nchar(metadata_names), decreasing = TRUE)]
  extract_metadata_name <- function(name) {
    hit <- metadata_names_ordered[mapply(startsWith, name, metadata_names_ordered)][1]
    if (is.na(hit)) {
      return(name)
    }
    hit
  }
  paras$metadata <- unlist(lapply(paras$name, extract_metadata_name))
  paras$value <- mapply(function(x, y) {
    if (is.na(x) || is.na(y) || x == y) {
      x
    } else {
      gsub(x, "", y)
    }
  }, paras$metadata, paras$name)
  paras <- paras[order(paras$qval, decreasing = FALSE),]
  paras <- dplyr::select(paras, c("feature", "metadata", "value"), dplyr::everything())
  paras <- dplyr::select(paras, -name)
  rownames(paras) <- NULL
  paras
}

combine_abundance_presence_results <- function(abundance_results,
                                               presence_results,
                                               correction = "BH") {
  abundance_keep <- dplyr::select(
    abundance_results,
    feature,
    metadata,
    value,
    coef_abundance = coef,
    stderr_abundance = stderr,
    pval_abundance = pval,
    qval_abundance = qval,
    base.model_abundance = base.model,
    tweedie.index
  )
  presence_keep <- dplyr::select(
    presence_results,
    feature,
    metadata,
    value,
    coef_presence = coef,
    stderr_presence = stderr,
    pval_presence = pval,
    qval_presence = qval,
    base.model_presence = base.model
  )

  combined <- dplyr::left_join(
    abundance_keep,
    presence_keep,
    by = c("feature", "metadata", "value")
  )
  pmat <- as.matrix(combined[, c("pval_abundance", "pval_presence"), drop = FALSE])
  combined$pval <- tweedieverse_cct_rows(pmat)
  combined$qval <- as.numeric(stats::p.adjust(combined$pval, method = correction))
  combined$coef <- combined$coef_abundance
  combined$stderr <- combined$stderr_abundance
  combined$base.model <- "CCT"
  combined <- dplyr::select(
    combined,
    feature,
    metadata,
    value,
    coef,
    stderr,
    pval,
    qval,
    coef_abundance,
    stderr_abundance,
    pval_abundance,
    qval_abundance,
    coef_presence,
    stderr_presence,
    pval_presence,
    qval_presence,
    base.model,
    base.model_abundance,
    base.model_presence,
    tweedie.index
  )
  rownames(combined) <- NULL
  combined
}
