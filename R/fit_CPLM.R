fit.CPLM <- function(features,
                     metadata,
                     link = "log",
                     tweedie_p = NULL,
                     Maaslin2_run = TRUE,
                     method_args = NULL,
                     formula = NULL,
                     random_effects_formula = NULL,
                     correction = 'BH',
                     cores = 4,
                     optimizer = 'nlminb',
                     na.action = na.exclude) {

  tweedie_link_power <- function(link) {
    switch(link,
           log = 0,
           identity = 1,
           sqrt = 0.5,
           inverse = -1,
           stop("Unsupported Tweedie link."))
  }

  fit_maaslin2_lm <- function(features,
                              metadata,
                              formula,
                              random_effects_formula,
                              correction,
                              cores,
                              method_args = NULL) {
    if (!requireNamespace("Maaslin2", quietly = TRUE)) {
      stop("Maaslin2 is required when tweedie_p = 0.")
    }

    maaslin_metadata <- metadata
    if ("offset" %in% colnames(maaslin_metadata)) {
      maaslin_metadata <- dplyr::select(maaslin_metadata, -offset)
    }

    fixed_effects <- setdiff(all.vars(formula)[-1], "offset")
    random_effects <- NULL
    if (!is.null(random_effects_formula)) {
      random_effects <- setdiff(all.vars(random_effects_formula)[-1], "offset")
    }

    tmp <- file.path(tempdir(), paste0("maaslin2_", sample(1e8, 1)))
    dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
    on.exit({
      if (requireNamespace("logging", quietly = TRUE)) {
        try(logging::removeHandler("logging::writeToFile"), silent = TRUE)
      }
      unlink(tmp, recursive = TRUE)
    }, add = TRUE)

    maaslin2_args <- extract_method_args(method_args, "Maaslin2")
    fit <- do.call(
      Maaslin2::Maaslin2,
      merge_method_args(
        list(
          input_data = features,
          input_metadata = maaslin_metadata,
          output = tmp,
          fixed_effects = fixed_effects,
          random_effects = random_effects,
          min_abundance = -Inf,
          normalization = "NONE",
          transform = "NONE",
          save_scatter = FALSE,
          save_models = FALSE,
          plot_heatmap = FALSE,
          plot_scatter = FALSE,
          max_significance = 1,
          correction = correction,
          standardize = FALSE,
          cores = cores
        ),
        maaslin2_args
      )
    )

    paras <- fit$results
    paras$base.model <- "Maaslin2"
    paras$tweedie.index <- 0
    paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
    if (!"value" %in% colnames(paras)) {
      paras$value <- paras$metadata
    }
    paras <- paras[order(paras$qval, decreasing = FALSE),]
    paras <- dplyr::select(
      paras,
      c('feature', 'metadata', 'value'),
      dplyr::everything()
    )
    rownames(paras) <- NULL
    return(list("results" = paras))
  }

  if (!is.null(tweedie_p) && tweedie_p == 0 && isTRUE(Maaslin2_run)) {
    return(fit_maaslin2_lm(
      features = features,
      metadata = metadata,
      formula = formula,
      random_effects_formula = random_effects_formula,
      correction = correction,
      cores = cores,
      method_args = method_args
    ))
  }

  if (!is.null(tweedie_p) && !is.null(random_effects_formula) &&
      !(tweedie_p %in% c(0, 1, 2, 3) || (tweedie_p > 1 && tweedie_p < 2))) {
    stop("With random_effects, fixed tweedie_p values are supported only for p = 0, p = 1, 1 < p < 2, p = 2, or p = 3.")
  }
  
  ######################################
  # Fit and summary functions for CPLM #
  ######################################
  
  if (is.null(random_effects_formula)) {
    
    ##########################
    # Fixed effects modeling #
    ##########################
    
    model_function <- function(formula,
                               data,
                               link,
                               tweedie_p,
                               optimizer,
                               na.action) {
      if (is.null(tweedie_p)) {
        return(
          cplm::cpglm(
            formula = formula,
            data = data,
            link = link,
            optimizer = optimizer,
            na.action = na.action
          )
        )
      }

      return(
        stats::glm(
          formula = formula,
          data = data,
          family = statmod::tweedie(
            var.power = tweedie_p,
            link.power = tweedie_link_power(link)
          ),
          na.action = na.action
        )
      )
    }
    
    summary_function <- function(fit) {
      if (inherits(fit, "cpglm")) {
        cplm_out <-
          capture.output(cplm_summary <- cplm::summary(fit)$coefficients)
        para <- as.data.frame(cplm_summary)[-1,-3]
        para$base.model <- 'CPLM'
        para$tweedie.index <- round(fit$p, 3)
        para$name <- rownames(cplm_summary)[-1]
      } else {
        glm_summary <- stats::coef(summary(fit))
        para <- as.data.frame(glm_summary)[-1, c(1, 2, 4), drop = FALSE]
        para$base.model <- 'Tweedie GLM'
        para$tweedie.index <- tweedie_p
        para$name <- rownames(glm_summary)[-1]
      }
      return(para)
    }
    
  } else{
    
    ###########################
    # Random effects modeling #
    ###########################
    
    fixed_terms <- setdiff(all.vars(formula)[-1], "offset")
    formula <-
      paste('. ~', paste(fixed_terms, collapse = ' + '), '.', sep = ' + ')
    formula <- update(random_effects_formula, formula)
    
    model_function <- function(formula,
                               data,
                               link,
                               tweedie_p,
                               optimizer,
                               na.action) {
      family <- glmmTMB::tweedie(link = link)
      start <- NULL
      map <- NULL
      if (!is.null(tweedie_p)) {
        if (tweedie_p > 1 && tweedie_p < 2) {
          start <- list(thetaf = stats::qlogis(tweedie_p - 1))
          map <- list(thetaf = factor(NA))
        } else {
          family <- switch(as.character(tweedie_p),
                           "0" = stats::gaussian(link = link),
                           "1" = stats::poisson(link = link),
                           "2" = stats::Gamma(link = link),
                           "3" = stats::inverse.gaussian(link = link),
                           stop("Unsupported random-effects Tweedie variance power."))
        }
      }
      args <- list(
          formula = formula,
          data = data,
          family = family,
          ziformula = ~ 0,
          na.action = na.action
      )
      if (!is.null(start)) {
        args$start <- start
      }
      if (!is.null(map)) {
        args$map <- map
      }
      return(do.call(glmmTMB::glmmTMB, args))
      
    }
    
    summary_function <- function(fit) {
      glmmTMB_summary <- coef(summary(fit))
      para <- as.data.frame(glmmTMB_summary$cond)[-1,-3]
      para$base.model <- ifelse(is.null(tweedie_p), 'CPLM', 'Tweedie GLMM')
      para$tweedie.index <- if (is.null(tweedie_p)) {
        round(unname(plogis(fit$fit$par["thetaf"]) + 1), 3)
      } else {
        tweedie_p
      }
      para$name <- rownames(glmmTMB_summary$cond)[-1]
      return(para)
    }
  }
  
  
  #######################################
  # Init cluster for parallel computing #
  #######################################
  
  cluster <- NULL
  if (cores > 1)
  {
    logging::loginfo("Creating cluster of %s R processes", cores)
    cluster <- parallel::makeCluster(cores)
    clusterExport(
      cluster,
      c(
        "features",
        "metadata",
        "formula",
        "link",
        "tweedie_p",
        "tweedie_link_power",
        "optimizer",
        "na.action",
        "model_function",
        "summary_function"
      ),
      envir = environment()
    )
  }
  
  ##############################
  # Apply per-feature modeling #
  ##############################
  
  outputs <-
    pbapply::pblapply(seq_len(ncol(features)), cl = cluster, function(x) {
      metadata_names <- setdiff(colnames(metadata), "offset")
      
      #################################
      # Create per-feature data frame #
      #################################
      
      featuresVector <- features[, x]
      logging::loginfo("Fitting model to feature number %d, %s",
                       x,
                       colnames(features)[x])
      dat_sub <-
        data.frame(expr = as.numeric(featuresVector), metadata)
      
      #############
      # Fit model #
      #############
      
      fit <- tryCatch({
        fit1 <-
          model_function(
            formula = formula,
            data = dat_sub,
            link = link,
            tweedie_p = tweedie_p,
            optimizer = optimizer,
            na.action = na.action
          )
      }, error = function(err) {
        fit1 <-
          try({
            model_function(
              formula = formula,
              data = dat_sub,
              link = link,
              tweedie_p = tweedie_p,
              optimizer = optimizer,
              na.action = na.action
            )
          })
        return(fit1)
      })
      
      #################
      # Gather Output #
      #################
      
      output <- list()
      if (all(!inherits(fit, "try-error"))) {
        output$para <- summary_function(fit)
      }
      else{
        logging::logwarn("Fitting problem for feature %s returning NA", x)
        output$para <-
          as.data.frame(matrix(NA,  nrow = length(metadata_names), ncol = 5))
        output$para$name <-
          metadata_names
      }
      colnames(output$para) <-
        c('coef',
          'stderr' ,
          'pval',
          'base.model',
          'tweedie.index',
          'name')
      output$para$feature <- colnames(features)[x]
      return(output)
    })
  
  ####################
  # Stop the cluster #
  ####################
  
  if (!is.null(cluster))
    parallel::stopCluster(cluster)
  
  #####################################
  # Bind the results for each feature #
  #####################################
  
  paras <-
    do.call(rbind, lapply(outputs, function(x) {
      return(x$para)
    }))
  
  ################################
  # Apply correction to p-values #
  ################################
  
  paras$qval <-
    as.numeric(p.adjust(paras$pval, method = correction))
  
  #####################################################
  # Determine the metadata names from the model names #
  #####################################################
  
  metadata_names <- setdiff(colnames(metadata), "offset")
  # order the metadata names by decreasing length
  metadata_names_ordered <-
    metadata_names[order(nchar(metadata_names), decreasing = TRUE)]
  # find the metadata name based on the match
  # to the beginning of the string
  extract_metadata_name <- function(name) {
    return(metadata_names_ordered[mapply(startsWith,
                                         name,
                                         metadata_names_ordered)][1])
  }
  paras$metadata <-
    unlist(lapply(paras$name, extract_metadata_name))
  # compute the value as the model contrast minus metadata
  paras$value <-
    mapply(function(x, y) {
      if (x == y)
        x
      else
        gsub(x, "", y)
    }, paras$metadata, paras$name)
  
  ##############################
  # Sort by decreasing q-value #
  ##############################
  
  paras <- paras[order(paras$qval, decreasing = FALSE),]
  paras <-
    dplyr::select(paras,
                  c('feature', 'metadata', 'value'),
                  dplyr::everything())
  paras <- dplyr::select(paras,-name)
  rownames(paras) <- NULL
  return(list("results" = paras))
}
