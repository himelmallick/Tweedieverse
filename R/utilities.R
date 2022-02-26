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
