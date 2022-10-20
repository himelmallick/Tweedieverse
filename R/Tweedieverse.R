#' Differential analysis of multi-omics data using Tweedie GLMs
#'
#' Fit a per-feature Tweedie generalized linear model to omics features.

#' @param input_features A tab-delimited input file or an R data frame of features (can be in rows/columns)
#' and samples (or cells). Samples are expected to have matching names with \code{input_metadata}. 
#' \code{input_features} can also be an object of class \code{SummarizedExperiment} or \code{SingleCellExperiment} 
#' that contains the expression or abundance matrix and other metadata; the \code{assays} 
#' slot contains the expression or abundance matrix and is named \code{"counts"}.  
#' This matrix should have one row for each feature and one sample for each column.  
#' The \code{colData} slot should contain a data frame with one row per 
#' sample and columns that contain metadata for each sample. 
#' Additional information about the experiment can be contained in the
#' \code{metadata} slot as a list.
#' @param input_metadata A tab-delimited input file or an R data frame of metadata (rows/columns).
#' Samples are expected to have matching sample names with \code{input_features}. 
#' This file is ignored when \code{input_features} is a \code{SummarizedExperiment} 
#' or a \code{SingleCellExperiment} object with \code{colData} containing the same information.
#' @param output The output folder to write results.
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
#' @param base_model The per-feature base model. Default is "CPLM". Must be one of "CPLM", "ZICP", "ZSCP", or "ZACP".
#' @param link A specification of the GLM link function. Default is "log". Must be one of "log", "identity", "sqrt", or "inverse".
#' @param fixed_effects Metadata variable(s) describing the fixed effects coefficients.
#' @param random_effects Metadata variable(s) describing the random effects part of the model.
#' @param cutoff_ZSCP For \code{base_model = "ZSCP"}, the cutoff to stratify features for
#' adaptive ZI modeling based on sparsity (zero-inflation proportion). Default is 0.3. Must be between 0 and 1.
#' @param criteria_ZACP For \code{base_model = "ZACP"}, the criteria to select the
#' best fitting model per feature.  The possible options are 'AIC' and BIC' (default).
#' More criteria will be supported in a future release.
#' @param adjust_offset If TRUE (default), an offset term will be included as the logarithm of \code{scale_factor}.
#' @param scale_factor Name of the numerical variable containing library size (for non-normalized data) or scale factor
#' (for normalized data) across samples to be included as an offset in the base model (when \code{adjust_offset = TRUE}).
#' If not found in metadata, defaults to the sample-wise total sums, unless \code{adjust_offset = FALSE}.
#' @param max_significance The q-value threshold for significance. Default is 0.05.
#' @param correction The correction method for computing the q-value (see \code{\link[stats]{p.adjust}} for options, default is 'BH').
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
#' @importFrom stats coef fitted as.formula na.exclude p.adjust plogis sd update
#' @importFrom utils capture.output read.table write.table
#' @importFrom dplyr %>% everything
#' @importFrom parallel clusterExport
#' @return A data frame containing coefficient estimates, p-values, and q-values (multiplicity-adjusted p-values) are returned.
#'
#' @author Himel Mallick, \email{himel.mallick@@merck.com}
#'
#' @examples
#'
#' \dontrun{
#'
#' ##############################################################################
#' # Example 1 - Differential Abundance Analysis of Synthetic Microbiome Counts #
#' ##############################################################################
#'
#' #######################################
#' # Install and Load Required Libraries #
#' #######################################
#'
#' library(devtools)
#' devtools::install_github('biobakery/sparseDOSSA@@varyLibSize')
#' library(sparseDOSSA)
#' library(stringi)
#'
#' ######################
#' # Specify Parameters #
#' ######################
#'
#' n.microbes <- 200 # Number of Features
#' n.samples <- 100 # Number of Samples
#' spike.perc <- 0.02 # Percentage of Spiked-in Bugs
#' spikeStrength<-"20" # Effect Size
#'
#' ###########################
#' # Specify Binary Metadata #
#' ###########################
#'
#' n.metadata <- 1
#' UserMetadata<-as.matrix(rep(c(0,1), each=n.samples/2))
#' UserMetadata<-t(UserMetadata) # Transpose
#'
#' ###################################################
#' # Spiked-in Metadata (Which Metadata to Spike-in) #
#' ###################################################
#'
#' Metadatafrozenidx<-1
#' spikeCount<-as.character(length(Metadatafrozenidx))
#' significant_metadata<-paste('Metadata', Metadatafrozenidx, sep='')
#'
#' #############################################
#' # Generate SparseDOSSA Synthetic Abundances #
#' #############################################
#'
#' DD<-sparseDOSSA::sparseDOSSA(number_features = n.microbes,
#' number_samples = n.samples,
#' UserMetadata=UserMetadata,
#' Metadatafrozenidx=Metadatafrozenidx,
#' datasetCount = 1,
#' spikeCount = spikeCount,
#' spikeStrength = spikeStrength,
#' noZeroInflate=TRUE,
#' percent_spiked=spike.perc,
#' seed = 1234)
#'
#' ##############################
#' # Gather SparseDOSSA Outputs #
#' ##############################
#'
#' sparsedossa_results <- as.data.frame(DD$OTU_count)
#' rownames(sparsedossa_results)<-sparsedossa_results$X1
#' sparsedossa_results<-sparsedossa_results[-1,-1]
#' colnames(sparsedossa_results)<-paste('Sample', 1:ncol(sparsedossa_results), sep='')
#' data<-as.matrix(sparsedossa_results[-c((n.metadata+1):(2*n.microbes+n.metadata)),])
#' data<-data.matrix(data)
#' class(data) <- "numeric"
#' truth<-c(unlist(DD$truth))
#' truth<-truth[!stri_detect_fixed(truth,":")]
#' truth<-truth[(5+n.metadata):length(truth)]
#' truth<-as.data.frame(truth)
#' significant_features<-truth[seq(1,
#' (as.numeric(spikeCount)+1)*(n.microbes*spike.perc), (as.numeric(spikeCount)+1)),]
#' significant_features<-as.vector(significant_features)
#'
#' ####################
#' # Extract Features #
#' ####################
#'
#' features<-as.data.frame(t(data[-c(1:n.metadata),]))
#'
#' ####################
#' # Extract Metadata #
#' ####################
#'
#' metadata<-as.data.frame(data[1,])
#' colnames(metadata)<-rownames(data)[1]
#'
#' ###############################
#' # Mark True Positive Features #
#' ###############################
#'
#' wh.TP = colnames(features) %in% significant_features
#' colnames(features)<-paste("Feature", 1:n.microbes, sep = "")
#' newname = paste0(colnames(features)[wh.TP], "_TP")
#' colnames(features)[wh.TP] <- newname;
#' colnames(features)[grep('TP', colnames(features))]
#'
#' ####################
#' # Run Tweedieverse #
#' ###################
#'
#' ###################
#' # Default options #
#' ###################
#'
#' CPLM <-Tweedieverse(
#' features,
#' metadata,
#' output = './demo_output/CPLM') # Assuming demo_output exists
#'
#' ###############################################
#' # User-defined prevalence-abundance filtering #
#' ###############################################
#'
#' ZICP<-Tweedieverse(
#' features,
#' metadata,
#' output = './demo_output/ZICP', # Assuming demo_output exists
#' base_model = 'ZICP',
#' abd_threshold = 0.0,
#' prev_threshold = 0.2)
#'
#' ####################################
#' # User-defined variance filtering  #
#' ####################################
#'
#' sds<-apply(features, 2, sd)
#' var_threshold = median(sds)/2
#' ZSCP<-Tweedieverse(
#' features,
#' metadata,
#' output = './demo_output/ZSCP', # Assuming demo_output exists
#' base_model = 'ZSCP',
#' var_threshold = var_threshold)
#'
#' ##################
#' # Multiple cores #
#' ##################
#'
#' ZACP<-Tweedieverse(
#' features,
#' metadata,
#' output = './demo_output/ZACP', # Assuming demo_output exists
#' base_model = 'ZACP',
#' cores = 4)
#'
#' ##########################################################################
#' # Example 2 - Multivariable Association on HMP2 Longitudinal Microbiomes #
#' ##########################################################################
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
#' standardize = FALSE,
#' reference = c('diagnosis,nonIBD'))
#'
#' }
#' @keywords microbiome, metagenomics, multiomics, scRNASeq, tweedie, singlecell
#' @export
Tweedieverse <- function(input_features,
                         input_metadata,
                         output,
                         abd_threshold = 0.0,
                         prev_threshold = 0.1,
                         var_threshold = 0.0,
                         entropy_threshold = 0.0,
                         base_model = "CPLM",
                         link = "log",
                         fixed_effects = NULL,
                         random_effects = NULL,
                         cutoff_ZSCP = 0.3,
                         criteria_ZACP = "BIC",
                         adjust_offset = TRUE,
                         scale_factor = NULL,
                         max_significance = 0.05,
                         correction = "BH",
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
  
  model_choices <- c("CPLM", "ZICP", "ZACP", "ZSCP")
  link_choices <- c("log", "identity", "sqrt", "inverse")
  criteria_ZACP_choices <- c("AIC", "BIC")
  correction_choices <-
    c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY")
  optimizer_choices <- c("nlminb", "bobyqa", "L-BFGS-B")
  
  #######################################################################
  #=====================================================================#
  # Read in the data and metadata, create output folder, initialize log #
  #=====================================================================#
  #######################################################################
  
  ##############################################################
  # Extract features and metadata based on user-provided input #
  ##############################################################
  
  if (!(methods::is(input_features, "SummarizedExperiment")) & !(methods::is(input_features, "SingleCellExperiment")) & !(is.character(input_features)) & !(is.data.frame(input_features))) {
    stop(cat(paste('Input data of class <',class(input_features), '> not supported. Please use either SummarizedExperiment or SingleCellExperiment or data.frame')))
  } else if (methods::is(input_features, "SummarizedExperiment") | methods::is(input_features, "SingleCellExperiment")) {
    SumExp <- methods::as(input_features, "SummarizedExperiment")
    data <- assay(SumExp); metadata <- as(colData(SumExp),"data.frame")
      if (is.null(SummarizedExperiment::assayNames(SumExp)) || SummarizedExperiment::assayNames(SumExp)[1] != "counts") {
        message("Renaming the first element in assays(input_features) to 'counts'")
        SummarizedExperiment::assayNames(SumExp)[1] <- "counts"
        if (is.null(colnames(counts(SumExp)))) {stop("Must supply sample/cell names!")}
        }
  } else{
    
    # if a character string then this is a file name, else it
    # is a data frame
    if (is.character(input_features)) {
      data <-
        data.frame(
          read.delim::fread(input_features, header = TRUE, sep = '\t'),
          header = TRUE,
          fill = T,
          comment.char = "" ,
          check.names = F,
          row.names = 1
        )
      if (nrow(input_features) == 1) {
        # read again to get row name
        data <- read.delim(
          input_features,
          header = TRUE,
          fill = T,
          comment.char = "" ,
          check.names = F,
          row.names = 1
        )
      }
    } else {
      data <- input_features
    }
    if (is.character(input_metadata)) {
      input_metadata <-
        data.frame(
          read.delim::fread(input_metadata, header = TRUE, sep = '\t'),
          header = TRUE,
          fill = T,
          comment.char = "" ,
          check.names = F,
          row.names = 1
        )
      if (nrow(input_metadata) == 1) {
        input_metadata <- read.delim(
          input_metadata,
          header = TRUE,
          fill = T,
          comment.char = "" ,
          check.names = F,
          row.names = 1
        )
      }
    } else {
      metadata <- input_metadata
    }
  }
    
  # create an output folder and figures folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  if (plot_heatmap || plot_scatter) {
    figures_folder <- file.path(output, "figures")
    if (!file.exists(figures_folder)) {
      print("Creating output figures folder")
      dir.create(figures_folder)
    }
  }
 
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
  logging::logdebug("Fixed effects: %s", fixed_effects)
  logging::logdebug("Random effects: %s", random_effects)
  logging::logdebug("ZSCP cutoff: %f", cutoff_ZSCP)
  logging::logdebug("ZACP criteria: %s", criteria_ZACP)
  logging::logdebug("Offset adjustment: %s", adjust_offset)
  logging::logdebug("Scale factor: %s", scale_factor)
  logging::logdebug("Max significance: %f", max_significance)
  logging::logdebug("Correction method: %s", correction)
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
  
  # Check if the selected criteria_ZACP is valid
    if (!criteria_ZACP %in% criteria_ZACP_choices) {
      option_not_valid_error(
        paste(
          "Please select a criteria",
          "from the list of available options"
        ),
        toString(criteria_ZACP_choices)
      )
    }
  
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
  
  ############################################################
  # Check if the selected numerical options are within range #
  ############################################################
  
  prop_options <- c(prev_threshold, cutoff_ZSCP, max_significance)
  if (any(prop_options < 0) || any(prop_options > 1)) {
    stop(
      paste(
        "One of the following is outside [0, 1]:",
        "prev_threshold, cutoff_ZSCP, max_significance"
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
      data <- type.convert(as.data.frame(t(data)))
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
        data <- type.convert(as.data.frame(t(data)))
        metadata <- type.convert(as.data.frame(t(metadata)))
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
          metadata <- type.convert(as.data.frame(t(metadata)))
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
        metadata[,i] = relevel(metadata[,i], ref = ref)
      }
    } else if (length(mlevels) > 2) {
      if (!is.na(ref)) {
        metadata[,i] = relevel(metadata[,i], ref = ref)
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
  
  sds <- apply(filtered_data, 2, na.rm = T, sd)
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
  
  
  ########################################################################
  # Set the scale factor to rowsum if not provided and create a modified #
  # metadata table without the scale factor variable (if present) ########
  ########################################################################
  
  if (is.null(scale_factor)) {
    offset <- rowSums(final_features)
  } else {
    if (!scale_factor %in% colnames(unfiltered_metadata)) {
      stop(
        paste(
          "The specified scale_factor variable is not present in the metadata table:\n",
          scale_factor
        )
      )
    } else {
      unfiltered_metadata <-
        dplyr::select(unfiltered_metadata,-scale_factor)
      offset <- unfiltered_metadata[, scale_factor]
    }
  }
  
  
  ####################################
  # Filter metadata based on entropy #
  ####################################
  
  # Reduce metadata to only include those pass entropy threshold
  temp_filtered_metadata <- unfiltered_metadata[, apply(unfiltered_metadata, 2, entropy) > entropy_threshold, drop=F]
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
    fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
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
    random_effects <-
      unlist(strsplit(random_effects, ",", fixed = TRUE))
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
  
  final_metadata <- as.data.frame(cbind.data.frame(filtered_metadata, offset))
  
  ##############################################################
  # Apply the base model to the filtered data with user inputs #
  ##############################################################
  
  logging::loginfo("Running selected analysis method: %s", base_model)
  
  fit_data <- fit.Tweedieverse(
    features = final_features,
    metadata = final_metadata,
    base_model = base_model,
    link = link,
    formula = formula,
    random_effects_formula = random_effects_formula,
    cutoff_ZSCP = cutoff_ZSCP,
    criteria_ZACP = criteria_ZACP,
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
  
  results_file <- file.path(output, "all_results.tsv")
  logging::loginfo("Writing all results to file (ordered by increasing q-values): %s",
                   results_file)
  ordered_results <-
    fit_data$results[order(fit_data$results$qval),]
  ordered_results <-
    ordered_results[!is.na(ordered_results$qval),] # Remove NA's
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
  write.table(
    ordered_results,
    file = results_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # Write results passing threshold to file
  # (removing any that are NA for the q-value)
  significant_results <-
    ordered_results[ordered_results$qval <= max_significance,]
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
  
  if (plot_heatmap & length(unique(significant_results_file["metdata"])) > 1) {
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
  logging::loginfo("Writing Tweedie inxed plot to file: %s",
                   output)
  tryCatch({
    tweedie_index_plot(ordered_results, figures_folder)
    
  }, error = function(err) {
    logging::logerror("Unable to do make a Tweedie inxed plot of results!!!")
    logging::logerror(err)
    # dev.off()
  })
  
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
  
  return(significant_results)
}


option_not_valid_error <- function(message, valid_options) {
  logging::logerror(paste(message, ": %s"), toString(valid_options))
  stop("Option not valid", call. = FALSE)
}

## Quiets concerns of R CMD check
utils::globalVariables(c("name"))
