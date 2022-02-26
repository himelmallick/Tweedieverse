fit.ZACP <- function(features,
                     metadata,
                     link = "log",
                     formula = NULL,
                     random_effects_formula = NULL,
                     criteria_ZACP = 'BIC',
                     correction = 'BH',
                     cores = 4,
                     optimizer = 'nlminb',
                     na.action = na.exclude) {
  
  ######################################
  # Fit and summary functions for ZACP #
  ######################################
  
  if (is.null(random_effects_formula)) {
    
    ##########################
    # Fixed effects modeling #
    ##########################
    
    model_function <- function(formula,
                               data,
                               link,
                               optimizer,
                               na.action,
                               criteria_ZACP) {
      
      
      ###################################
      # Determing best fitting ZI model #
      ###################################
      
      fitA <- tryCatch({
        fit1 <- cplm::cpglm(
          formula = formula,
          data = data,
          link = link,
          optimizer = optimizer,
          na.action = na.action
        )
      }, error = function(err) {
        fit1 <- try({
          cplm::cpglm(
            formula = formula,
            data = data,
            link = link,
            optimizer = optimizer,
            na.action = na.action
          )
        })
        return(fit1)
      })
      
      fitB <- tryCatch({
        fit1 <- cplm::zcpglm(
          formula = formula,
          data = data,
          link = link,
          optimizer = optimizer,
          na.action = na.action
        )
      }, error = function(err) {
        fit1 <- try({
          cplm::zcpglm(
            formula = formula,
            data = data,
            link = link,
            optimizer = optimizer,
            na.action = na.action
          )
        })
        return(fit1)
      })
      
      ####################
      # Winner Selection #
      ####################
      
      # If only one model has failed the other becomes winner by default
      # If both have failed
      if (class(fitA) == "try-error" && class(fitB) != "try-error") {
        fit_indicator = "ZICP"
      } else if (class(fitA) != "try-error"  && class(fitB) == "try-error") {
        fit_indicator = "CPLM"
      } else if (class(fitA) == "try-error"  && class(fitB) == "try-error") {
        fit_indicator = "CPLM"
      } else {
        # If both converged, pick the best based on ZACP criteria
        if (!criteria_ZACP %in% c('AIC', 'BIC')) {
          stop('The selected ZACP criteria is not supported. Use either AIC or BIC.')
        } 
        else{
          AICtab<-c(get_AICtab(fitA)[criteria_ZACP], get_AICtab(fitB)[criteria_ZACP])
          best.fit <- which.min(AICtab)
          fit_indicator <- ifelse(best.fit == 1, "CPLM", "ZICP")
        }
      }
      
      
      ###############################
      # Optimized per-feature model #
      ###############################
      
      ifelse(fit_indicator == 'CPLM', fitA, fitB)
    }
    
    
    summary_function <- function(fit) {
      if (class(fit) == "cpglm") {
        cplm_out <-
          capture.output(cplm_summary <- cplm::summary(fit)$coefficients)
        para <- as.data.frame(cplm_summary)[-1,-3]
        para$base.model <- 'CPLM'
        para$tweedie.index <- round(fit$p, 3)
        para$name <- rownames(cplm_summary)[-1]
        return(para)
      }
      if (class(fit) == "zcpglm") {
        zicp_out <-
          capture.output(zicp_summary <-
                           cplm::summary(fit)$coefficients$tweedie)
        para <- as.data.frame(zicp_summary)[-1,-3]
        para$base.model <- 'ZICP'
        para$tweedie.index <- round(fit$p, 3)
        para$name <- rownames(zicp_summary)[-1]
        return(para)
      }
    }
    
  } else{
    
    ###########################
    # Random effects modeling #
    ###########################
    
    formula <-paste('. ~', paste(all.vars(formula)[-1], collapse = ' + '), '.', sep = ' + ')
    formula <- update(random_effects_formula, formula)
    
    model_function <- function(formula,
                               data,
                               link,
                               optimizer,
                               na.action,
                               criteria_ZACP) {
      
      ###################################
      # Determing best fitting ZI model #
      ###################################
      
      fitA <- tryCatch({
        fit1 <- glmmTMB::glmmTMB(
          formula = formula,
          data = data,
          family = glmmTMB::tweedie(link = link),
          ziformula = ~ 0,
          na.action = na.action
        )
      }, error = function(err) {
        fit1 <- try({
          glmmTMB::glmmTMB(
            formula = formula,
            data = data,
            family = glmmTMB::tweedie(link = link),
            ziformula = ~ 0,
            na.action = na.action
          )
        })
        return(fit1)
      })
      
      fitB <- tryCatch({
        fit1 <- glmmTMB::glmmTMB(
          formula = formula,
          data = data,
          family = glmmTMB::tweedie(link = link),
          ziformula = ~ 1,
          na.action = na.action
        )
      }, error = function(err) {
        fit1 <- try({
          glmmTMB::glmmTMB(
            formula = formula,
            data = data,
            family = glmmTMB::tweedie(link = link),
            ziformula = ~ 1,
            na.action = na.action
          )
        })
        return(fit1)
      })
      
      ####################
      # Winner Selection #
      ####################
      
      # If only one model has failed the other becomes winner by default
      # If both have failed
      if (class(fitA) == "try-error" && class(fitB) != "try-error") {
        fit_indicator = "ZICP"
      } else if (class(fitA) != "try-error"  && class(fitB) == "try-error") {
        fit_indicator = "CPLM"
      } else if (class(fitA) == "try-error"  && class(fitB) == "try-error") {
        fit_indicator = "CPLM"
      } else {
        # If both converged, pick the best based on ZACP criteria
        if (!criteria_ZACP %in% c('AIC', 'BIC')) {
          stop('The selected ZACP criteria is not supported. Use either AIC or BIC.')
        } 
        else{
          AICtab<-c(get_AICtab(fitA)[criteria_ZACP], get_AICtab(fitB)[criteria_ZACP])
          best.fit <- which.min(AICtab)
          fit_indicator <- ifelse(best.fit == 1, "CPLM", "ZICP")
        }
      }
      
      
      ###############################
      # Optimized per-feature model #
      ###############################
      
      ifelse(fit_indicator == 'CPLM', fitA, fitB)
    }
    
    
    summary_function <- function(fit) {
      glmmTMB_summary <- coef(summary(fit))
      para <- as.data.frame(glmmTMB_summary$cond)[-1,-3]
      para$base.model <-ifelse(is.null(glmmTMB_summary$zi), 'CPLM', 'ZICP')
      para$tweedie.index <-round(unname(plogis(fit$fit$par["thetaf"]) + 1), 3)
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
        "optimizer",
        "na.action",
        "model_function",
        "summary_function",
        "get_AICtab",
        "criteria_ZACP"
      ),
      envir = environment()
    )
  }
  
  ##############################
  # Apply per-feature modeling #
  ##############################
  
  outputs <-
    pbapply::pblapply(seq_len(ncol(features)), cl = cluster, function(x) {
      
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
            optimizer = optimizer,
            na.action = na.action,
            criteria_ZACP = criteria_ZACP
          )
      }, error = function(err) {
        fit1 <-
          try({
            model_function(
              formula = formula,
              data = dat_sub,
              link = link,
              optimizer = optimizer,
              na.action = na.action,
              criteria_ZACP = criteria_ZACP
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
        logging::logwarn(paste("Fitting problem for feature", x, "returning NA"))
        output$para <-
          as.data.frame(matrix(NA,  nrow = ncol(metadata) - 1, ncol = 5)) # Everything except offset
        output$para$name <-
          colnames(metadata)[-ncol(metadata)] # Everything except offset
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
