fit.Tweedieverse <- function(features,
                             metadata,
                             base_model = 'CPLM',
                             link = "log",
                             tweedie_p = NULL,
                             Maaslin2_run = TRUE,
                             method_args = NULL,
                             formula = NULL,
                             random_effects_formula = NULL,
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
  
  if (base_model != "CPLM") {
    stop("Only base_model = 'CPLM' is supported.")
  }

  fit <- fit.CPLM(
    features = features,
    metadata = metadata,
    link = link,
    tweedie_p = tweedie_p,
    Maaslin2_run = Maaslin2_run,
    method_args = method_args,
    formula = formula,
    random_effects_formula = random_effects_formula,
    correction = correction,
    cores = cores,
    optimizer = optimizer,
    na.action = na.action
  )
  return(fit)
}
