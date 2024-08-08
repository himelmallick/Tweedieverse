fit.Tweedieverse <- function(features,
                             metadata,
                             base_model = 'CPLM',
                             link = "log",
                             formula = NULL,
                             random_effects_formula = NULL,
                             cutoff_ZSCP = 0.3,
                             criteria_ZACP = 'BIC',
                             adjust_offset = TRUE,
                             correction = 'BH',
                             cores = 4,
                             optimizer = 'nlminb',
                             na.action = na.exclude) {
  
  ################################################################
  # Set the formula default to all fixed effects if not provided #
  ################################################################
  
  if ("offset" %in% colnames(metadata)) {
    all_available_metadata <- setdiff(colnames(metadata), "offset")
    if (is.null(formula))
      formula <-
        as.formula(paste("expr ~ ", paste(all_available_metadata, collapse = "+")))
    if (adjust_offset)
      formula <- update(formula, . ~ . - offset(log(offset)))
  } else{
    if (is.null(formula))
      formula <-
        as.formula(paste("expr ~ ", paste(colnames(metadata), collapse = "+")))
  }
  
  
  ##############################################################
  # Call per-feature base models and return results for output #
  ##############################################################
  
  if (base_model == 'CPLM') {
    fit <- fit.CPLM(
      features = features,
      metadata = metadata,
      link = link,
      formula = formula,
      random_effects_formula = random_effects_formula,
      correction = correction,
      cores = cores,
      optimizer = optimizer,
      na.action = na.action
    )
    
  }
  
  if (base_model == 'ZICP') {
    fit <- fit.ZICP(
      features = features,
      metadata = metadata,
      link = link,
      formula = formula,
      random_effects_formula = random_effects_formula,
      correction = correction,
      cores = cores,
      optimizer = optimizer,
      na.action = na.action
    )
  }
  
  
  if (base_model == 'ZSCP') {
    fit <- fit.ZSCP(
      features = features,
      metadata = metadata,
      link = link,
      formula = formula,
      random_effects_formula = random_effects_formula,
      cutoff_ZSCP = cutoff_ZSCP,
      correction = correction,
      cores = cores,
      optimizer = optimizer,
      na.action = na.action
    )
  }
  
  if (base_model == 'ZACP') {
    fit <- fit.ZACP(
      features = features,
      metadata = metadata,
      link = link,
      formula = formula,
      random_effects_formula = random_effects_formula,
      criteria_ZACP = criteria_ZACP,
      correction = correction,
      cores = cores,
      optimizer = optimizer,
      na.action = na.action
    )
  }
  return(fit)
}
