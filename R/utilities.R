get_AICtab<-function(fit){
  
  ########################
  # Flag invalid options #
  ########################
  
  if (!class(fit) %in% c('cpglm', 'zcpglm', 'glmmTMB')){
    stop('Not supported. Valid options are cplm , zcpglm, and glmmTMB')
  }
  
  ######################
  #  Initialize AICtab #
  ######################
  
  AICtab<-rep(NA, 5)
  
  ###########################
  # Case-by-Case Extraction #
  ###########################
  
  if (class(fit)=='cpglm'){
    
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
  
  if (class(fit)=='zcpglm'){
    
    ##########################################
    # Back calculate AIC and BIC from logLik #
    ##########################################
    
    logLik<--fit$llik
    AIC_multiplier<-length(fit$y) - fit$df.residual
    BIC_multiplier<-AIC_multiplier*log(length(fit$y))
    AIC<-2*AIC_multiplier + 2*logLik
    BIC<-BIC_multiplier + 2*logLik
    deviance<-NA
    df.resid<-fit$df.residual
    
    # Coherent output
    AICtab<-c(AIC, BIC, logLik, deviance, df.resid)
    
  }
  
  if (class(fit)=='glmmTMB'){
    
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
  
  data$cat<-cut(data[,feature], breaks=bins, labels=c(1:bins))
  
  #use dplyr to compute e and p for each value of the feature
  dd_data <- data %>% group_by(cat) %>% summarise(e=entropy(get(target)), 
                                                  n=length(get(target)),
                                                  min=min(get(feature)),
                                                  max=max(get(feature))
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
  dd_data <- data %>% group_by_at(feature) %>% summarise(e=entropy(get(target)), 
                                                         n=length(get(target))
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
extractAssay <- function(input, assay_name = "counts") {
  
  # Extract assay name based on the user input
  if ("counts" %in% assayNames(input)) {
    counts_data <- assay(input, assay_name)
    cat("The specified assay has been extracted\n")
    return(as.data.frame(as.matrix(counts_data)))
  } else {
    cat("The specified assay was not found\n")
    return(NULL)
  }
}


#' Median Comparison for Compositionality Adjustment
#'
#' Adjust Tweedieverse(or any other differential analysis methods) coefficient estimates and p-values by testing each taxon
#' against the *median* effect for the same metadata variable — a simple
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
#'         approximation (`N(effect_size, stderr²)`) and records the empirical
#'         distribution of the simulated medians;
#'   \item derives a variance-inflated *offset* that accounts for the
#'         covariance between each coefficient and the group median;
#'   \item performs a two-sided Z-test of `H₀ : βᵢ = offsetᵢ`, returning the
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
    cur_median <- median(sub_df$effect_size[use_idx], na.rm = TRUE)
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
      sim_coefs   <- rnorm(n_coefs, mean = coefs, sd = ses)
      sim_median  <- median(sim_coefs[use_bool])
      c(sim_median, sim_coefs)
    })
    
    sim_medians <- sim_results[1, ]
    all_sims    <- sim_results[-1, , drop = FALSE]  # row i => draws for coef i
    
    # Covariance of each coefficient with the median, across draws
    cov_adjust <- apply(all_sims, 1, function(x) cov(x, sim_medians))
    
    # 6) "offset to test" for each coefficient, per the MaAsLin 3 approach:
    #    offset_i = coefs[i] ± ... depends on difference from median & correlation
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
        pvals_median[i] <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
      }
    }
    
    # Save the results back into the subset
    sub_df$pval_median <- pvals_median
    df[sub_idx, ]      <- sub_df
  }
  
  # Return the augmented data
  return(df)
}
