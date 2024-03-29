###############################################################################
# Tweedieverse visualizations

# Copyright (c) 2020 Himel Mallick and Ali Rahnavard

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################

# Author: Ali Rahnavard
# Email: gholamali.rahnavard@gmail.com
# This script includes functions for visualizing overall output of Tweedieverse and
# individual associations as scatterplot and boxplot
#  Visualization results are provided as pdf and RDS format to be used with manuscript quality.

# Load libraries
for (lib in c('ggplot2', "grid", 'pheatmap', 'cowplot')) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}


#' @export
theme_omicsEye <- function() list(
  cowplot::theme_cowplot(),
  ggplot2::theme(
    text               = ggplot2::element_text(size = 6),
    axis.text          = ggplot2::element_text(size = 5),
    axis.title.x       = ggplot2::element_text(margin=ggplot2::margin(1, 0, 0.5, 0)),
    axis.title.x.top   = ggplot2::element_text(margin=ggplot2::margin(0, 0, 2, 0)),
    axis.title.y       = ggplot2::element_text(margin=ggplot2::margin(0, 1, 0, 0.5)),
    axis.title.y.right = ggplot2::element_text(margin=ggplot2::margin(0, 0, 0, 2)),
    axis.text.x        = ggplot2::element_text(margin=ggplot2::margin(1, 0, 0, 0)),
    axis.text.x.top    = ggplot2::element_text(margin=ggplot2::margin(0, 0, 1, 0)),
    axis.text.y        = ggplot2::element_text(margin=ggplot2::margin(0, 1, 0, 0)),
    axis.text.y.right  = ggplot2::element_text(margin=ggplot2::margin(0, 0, 0, 1)),
    axis.ticks         = ggplot2::element_line(size=0.3),
    axis.ticks.length  = ggplot2::unit(2, "pt"),
    axis.line          = ggplot2::element_line(size=0.3),
    axis.line.x        = ggplot2::element_line(size=0.3),
    axis.line.y        = ggplot2::element_line(size=0.3),
    line               = ggplot2::element_line(size=0.3),
    legend.margin      = ggplot2::margin(4, 4, 4, 4),
    legend.key.size    = ggplot2::unit(8, "pt"),
    legend.box.spacing = ggplot2::unit(4, "pt"),
    panel.spacing      = ggplot2::unit(1.5, "pt"),
    plot.title         = ggplot2::element_text(size=8),
    plot.margin        = ggplot2::margin(1, 1, 1, 1),
    strip.background   = ggplot2::element_blank(),
    strip.text         = ggplot2::element_text(size=6),
    strip.text.x       = ggplot2::element_text(margin= ggplot2::margin(3, 0, 3, 0)),
    strip.text.y       = ggplot2::element_text(margin= ggplot2::margin(0, 3, 0, 3))
  )
)



# Tweedieverse heatmap function for overall view of associations
#' @export
Tweedieverse_heatmap <-
  function(output_results,
           title = NA,
           cell_value = 'qval',
           data_label = 'data',
           metadata_label = 'metadata',
           border_color = 'grey93',
           color = colorRampPalette(c("darkblue", "grey90", "darkred")),
           col_rotate = 90,
           first_n = 50,
           write_to = NA) {
    # read Tweedieverse output
    if (is.character(output_results)) {
      df <- read.table(
        output_results,
        header = TRUE,
        sep = "\t",
        fill = TRUE,
        comment.char = "" ,
        check.names = FALSE
      )
    } else {
      data <- output_results
    }

    title_additional <- ""

    title_additional <- ""
    if (!is.na(first_n) & first_n > 0 & first_n < dim(df)[1]) {
      if (cell_value == 'coef') {
        df <- df[order(-abs(df[cell_value])) ,]
      } else{
        df <- df[order(df[cell_value]),]
      }
      # get the top n features with significant associations
      df_sub <- df[1:first_n, ]
      for (first_n_index in seq(first_n, dim(df)[1]))
      {
        if (length(unique(df_sub$feature)) == first_n)
        {
          break
        }
        df_sub <- df[1:first_n_index, ]
      }
      # get all rows that have the top N features
      df <- df[which(df$feature %in% df_sub$feature), ]
      title_additional <- paste("Top", first_n, sep = " ")
    }

    if (dim(df)[1] < 2) {
      print('There are no associations to plot!')
      return(NULL)
    }

    metadata <- df$metadata
    data <- df$feature
    dfvalue <- df$value
    value <- NA

    # values to use for coloring the heatmap
    # and set the colorbar boundaries
    if (cell_value == "pval") {
      value <- -log(df$pval) * sign(df$coef)
      value <- pmax(-20, pmin(20, value))
      if (is.null(title))
        title <- "(-log(pval)*sign(coeff))"
    } else if (cell_value == "qval") {
      value <- -log(df$qval) * sign(df$coef)
      value <- pmax(-20, pmin(20, value))
      if (is.null(title))
        title <- "(-log(qval)*sign(coeff))"
    } else if (cell_value == "coef") {
      value <- df$coef
      if (is.null(title))
        title <- "(coeff)"
    }

    if (title_additional != "") {
      title <-
        paste(title_additional,
              "features with significant associations",
              title,
              sep = " ")
    } else {
      title <- paste("Significant associations", title, sep = " ")
    }

    # identify variables with more than one level present
    verbose_metadata <- c()
    metadata_multi_level <- c()
    for (i in unique(metadata)) {
      levels <- unique(df$value[df$metadata == i])
      if (length(levels) > 1) {
        metadata_multi_level <- c(metadata_multi_level, i)
        for (j in levels) {
          verbose_metadata <- c(verbose_metadata, paste(i, j))
        }
      } else {
        verbose_metadata <- c(verbose_metadata, i)
      }
    }

    n <- length(unique(data))
    m <- length(unique(verbose_metadata))

    if (n < 2) {
      print(
        paste(
          "There is not enough features in the associations",
          "to create a heatmap plot.",
          "Please review the associations in text output file."
        )
      )
      return(NULL)
    }

    if (m < 2) {
      print(
        paste(
          "There is not enough metadata in the associations",
          "to create a heatmap plot.",
          "Please review the associations in text output file."
        )
      )
      return(NULL)
    }

    a = matrix(0, nrow = n, ncol = m)
    a <- as.data.frame(a)

    rownames(a) <- unique(data)
    colnames(a) <- unique(verbose_metadata)

    for (i in seq_len(dim(df)[1])) {
      current_metadata <- metadata[i]
      if (current_metadata %in% metadata_multi_level) {
        current_metadata <- paste(metadata[i], dfvalue[i])
      }
      if (abs(a[as.character(data[i]),
                as.character(current_metadata)]) > abs(value[i]))
        next
      a[as.character(data[i]), as.character(current_metadata)] <-
        value[i]
    }

    # get the range for the colorbar
    max_value <- ceiling(max(a))
    min_value <- ceiling(min(a))
    range_value <- max(c(abs(max_value), abs(min_value)))
    breaks <- seq(-1 * range_value, range_value, by = 1)

    p <- NULL
    tryCatch({
      p <-
        pheatmap::pheatmap(
          a,
          cellwidth = 5,
          cellheight = 5,
          # changed to 3
          main = title,
          fontsize = 6,
          kmeans_k = NA,
          border = TRUE,
          show_rownames = TRUE,
          show_colnames = TRUE,
          scale = "none",
          cluster_rows = FALSE,
          cluster_cols = TRUE,
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          legend = TRUE,
          border_color = border_color,
          color = color(range_value * 2),
          breaks = breaks,
          treeheight_row = 0,
          treeheight_col = 0,
          display_numbers = matrix(ifelse(a > 0.0, "+", ifelse(a < 0.0, "-", "")), nrow(a))
        )
    }, error = function(err) {
      logging::logerror("Unable to plot heatmap")
      logging::logerror(err)
    })
    if (!is.na(write_to))
      saveRDS(p, file = paste(write_to, "/gg_heatmap.RDS", sep = ""))
    return(p)
  }

save_heatmap <-
  function(results_file,
           heatmap_file,
           figures_folder,
           title = NULL,
           cell_value = "qval",
           data_label = 'data',
           metadata_label = 'metadata',
           border_color = "grey93",
           color = colorRampPalette(c("blue", "grey90", "red")),
           first_n = 50) {
    # generate a heatmap and save it to a pdf and as a jpg
    heatmap <-
      Tweedieverse_heatmap(
        results_file,
        title,
        cell_value,
        data_label,
        metadata_label,
        border_color,
        color,
        first_n,
        figures_folder
      )

    if (!is.null(heatmap)) {
      pdf(heatmap_file)
      print(heatmap)
      dev.off()

      jpg_file <- file.path(figures_folder, "heatmap.jpg")
      jpeg(jpg_file,
           res = 150,
           height = 800,
           width = 1100)
      print(heatmap)
      dev.off()
    }

  }

association_plots <-
  function(metadata,
           features,
           output_results,
           write_to = './',
           figures_folder = './figures/',
           max_jpgs = 3)
  {
    #Tweedieverse scatter plot function and theme

    # combine the data and metadata to one datframe using common rows
    # read Tweedieverse output
    if (is.character(features)) {
      features <-
        data.frame(
          read.delim::fread(features, header = TRUE, sep = '\t'),
          header = TRUE,
          fill = T,
          comment.char = "" ,
          check.names = F,
          row.names = 1
        )
      if (nrow(data) == 1) {
        # read again to get row name
        features <- read.delim(
          features,
          header = TRUE,
          fill = T,
          comment.char = "" ,
          check.names = F,
          row.names = 1
        )
      }
    } else {
      features <- features
    }
    if (is.character(metadata)) {
      metadata <-
        data.frame(
          read.delim::fread(input_metadata, header = TRUE, sep = '\t'),
          header = TRUE,
          fill = T,
          comment.char = "" ,
          check.names = F,
          row.names = 1
        )
      if (nrow(metadata) == 1) {
        metadata <- read.delim(
          input_metadata,
          header = TRUE,
          fill = T,
          comment.char = "" ,
          check.names = F,
          row.names = 1
        )
      }
    } else {
      metadata <- metadata
    }
    common_rows <-
      intersect(rownames(features), rownames(metadata))
    input_df_all <-
      cbind(features[common_rows, , drop = FALSE],
            metadata[common_rows, , drop = FALSE])

    # read Tweedieverse output
    if (is.character(output_results)) {
      output_df_all <- read.table(
        output_results,
        header = TRUE,
        row.names = NULL,
        sep = "\t",
        fill = FALSE,
        comment.char = "" ,
        check.names = FALSE
      )
    } else {
      output_df_all <- output_results
    }

    if (dim(output_df_all)[1] < 1) {
      print('There are no associations to plot!')
      return(NULL)
    }

    logging::loginfo(
      paste(
        "Plotting associations from most",
        "to least significant,",
        "grouped by metadata"
      )
    )
    metadata_types <- unlist(output_df_all[, 'metadata'])
    metadata_labels <-
      unlist(metadata_types[!duplicated(metadata_types)])
    metadata_number <- 1
    if (is.na(max_jpgs))
      max_jpgs <- dim(output_df_all)[1]
    saved_plots <- vector('list', max_jpgs)
    for (label in metadata_labels) {
      # for file name replace any non alphanumeric with underscore
      plot_file <-
        paste(write_to,
              "/",
              gsub("[^[:alnum:]_]", "_", label),
              ".pdf",
              sep = "")
      data_index <- which(label == metadata_types)
      saved_ggs <- vector('list', length(data_index))
      logging::loginfo("Plotting data for metadata number %s, %s",
                       metadata_number,
                       label)
      pdf(plot_file,
          width = 2.65,
          height = 2.5,
          onefile = TRUE)

      x <- NULL
      y <- NULL
      count <- 1
      for (i in data_index) {
        x_label <- as.character(output_df_all[i, 'metadata'])
        y_label <- as.character(output_df_all[i, 'feature'])
        results_value <-
          as.character(output_df_all[i, 'value'])
        qval <- as.numeric(output_df_all[i, 'qval'])
        coef_val <- as.numeric(output_df_all[i, 'coef'])
        input_df <- input_df_all[c(x_label, y_label)]
        colnames(input_df) <- c("x", "y")

        # if Metadata is continuous generate a scatter plot
        # Continuous is defined as numerical with more than
        # 2 values (to exclude binary data)
        temp_plot <- NULL
        if (is.numeric(input_df[1, 'x']) &
            length(unique(input_df[['x']])) > 2) {
          logging::loginfo("Creating scatter plot for continuous data, %s vs %s",
                           x_label,
                           y_label)
          temp_plot <- ggplot2::ggplot(data = input_df,
                                       ggplot2::aes(
                                         as.numeric(as.character(x)),
                                         as.numeric(as.character(y))
                                       )) +
            ggplot2::geom_point(
              fill = 'darkolivegreen4',
              color = 'black',
              alpha = .5,
              shape = 21,
              size = 1,
              stroke = 0.15

            ) +
            ggplot2::scale_x_continuous(limits = c(min(input_df['x']), max(input_df['x']))) +
            ggplot2::scale_y_continuous(limits = c(min(input_df['y']), max(input_df['y']))) +
            ggplot2::stat_smooth(
              method = "glm",
              size = 0.25,
              color = 'blue',
              na.rm = TRUE
            ) +
            ggplot2::guides(alpha = 'none') +
            ggplot2::labs("") +
            ggplot2::xlab(x_label) +
            ggplot2::ylab(y_label) +
            theme_omicsEye() +
            ggplot2::annotate(
              geom = "text",
              x = Inf,
              y = Inf,
              hjust = 1,
              vjust = 1,
              label = sprintf(
                "FDR: %s\nCoefficient: %s\nN: %s",
                formatC(qval, format = "e", digits = 3),
                formatC(coef_val,
                        format = "e",
                        digits = 2),
                formatC(length(input_df[, 'x']))
              ) ,
              color = "black",
              size = 2,
              fontface = "italic"
            ) + ggplot2::scale_y_log10()
        } else{
          # if Metadata is categorical generate a boxplot
          ### check if the variable is categorical

          logging::loginfo("Creating boxplot for categorical data, %s vs %s",
                           x_label,
                           y_label)
          input_df['x'] <-
            lapply(input_df['x'], as.character)

          # count the Ns for each group
          x_axis_label_names <- unique(input_df[['x']])
          renamed_levels <-
            as.character(levels(metadata[, x_label]))
          if (length(renamed_levels) == 0) {
            renamed_levels <- x_axis_label_names
          }
          for (name in x_axis_label_names) {
            total <- length(which(input_df[['x']] == name))
            new_n <-
              paste(name, " (n=", total, ")", sep = "")
            input_df[which(input_df[['x']] == name), 'x'] <-
              new_n
            renamed_levels <-
              replace(renamed_levels, renamed_levels == name, new_n)
          }
          input_df$xnames <-
            factor(input_df[['x']], levels = renamed_levels)
          temp_plot <-
            ggplot2::ggplot(data = input_df, ggplot2::aes(xnames, y)) +
            ggplot2::geom_boxplot(
              ggplot2::aes(fill = xnames),
              outlier.alpha = 0.0,
              na.rm = TRUE,
              alpha = .5,
              show.legend = FALSE
            ) +
            ggplot2::geom_point(
              ggplot2::aes(fill = xnames),
              alpha = 0.75 ,
              size = 1,
              shape = 21,
              stroke = 0.15,
              color = 'black',
              show.legend = FALSE,
              position = ggplot2::position_jitterdodge()
            ) +
            ggplot2::scale_fill_brewer(palette = "Spectral", direction = -1)

          # format the figure to default nature format
          # remove legend, add x/y labels
          temp_plot <- temp_plot +
            theme_omicsEye() +
            ggplot2::theme(
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_blank(),
              axis.line = ggplot2::element_line(colour = "black")
            ) +
            ggplot2::xlab("") + #x_label
            ggplot2::ylab(y_label) +
            ggplot2::annotate(
              geom = "text",
              x = Inf,
              y = Inf,
              hjust = 1,
              vjust = 1,
              label = sprintf(
                "FDR: %s\nCoefficient: %s\nValue: %s",
                formatC(qval, format = "e", digits = 3),
                formatC(coef_val,
                        format = "e",
                        digits = 2),
                results_value
              ) ,
              color = "black",
              size = 2,
              fontface = "italic"
            ) + ggplot2::scale_y_log10()
        }
        stdout <-
          capture.output(print(temp_plot), type = "message")
        if (length(stdout) > 0)
          logging::logdebug(stdout)
        if (count < max_jpgs + 1)
          saved_plots[[count]] <- temp_plot
        saved_ggs[[count]] <- temp_plot
        count <- count + 1
      }

      dev.off()
      # print the saved figures
      for (plot_number in seq(1, max_jpgs)) {
        jpg_file <- file.path(figures_folder,
                              paste0(substr(
                                basename(plot_file), 1, nchar(basename(plot_file)) - 4
                              ),
                              "_",
                              plot_number,
                              ".jpg"))
        jpeg(jpg_file,
             res = 300,
             width = 960,
             height = 960)
        stdout <-
          capture.output(print(saved_plots[[plot_number]]))
        dev.off()
      }
      saveRDS(saved_ggs,
              file = paste(figures_folder,
                           "/" ,
                           label,
                           "_gg_associations.RDS",
                           sep = ""))
      metadata_number <- metadata_number + 1
    }

  }

tweedie_index_plot <-
  function(output_results,
           figures_folder = './figures/')
  {
    # read Tweedieverse output
    if (is.character(output_results)) {
      output_df_all <- read.table(
        output_results,
        header = TRUE,
        row.names = NULL,
        sep = "\t",
        fill = FALSE,
        comment.char = "" ,
        check.names = FALSE
      )
    } else {
      output_df_all <- output_results
    }

    if (dim(output_df_all)[1] < 1) {
      print('There are no associations to plot!')
      return(NULL)
    }

    logging::loginfo(paste("Plotting tweedie.index ",
                           "colored by metadata"))
    plot_file <-
      paste(figures_folder,
            "/tweedie_index_plot.pdf",
            sep = "")

    logging::loginfo("Plotting tweedie.index")
    pdf(plot_file,
        width = 2.65,
        height = 1.5,
        onefile = TRUE)

    temp_plot <- NULL

    temp_plot <- ggplot2::ggplot(output_df_all,
      ggplot2::aes(x=tweedie.index))+
      ggplot2::geom_histogram(position = "identity", alpha = 0.8)  +
      ggplot2::xlab("Tweedie index") + ggplot2::ylab("Number of omics features") +
      theme_omicsEye() +
      ggplot2::theme(legend.justification = c(0, 0),
                     legend.position = c(.3, .5))

    # print the saved figures
    tdout <-
      capture.output(print(temp_plot), type = "message")
    dev.off()
    jpg_file <- paste(figures_folder,
                      "/tweedie_index_plot.jpg",
                      sep = "")
    jpeg(jpg_file,
         res = 300,
         width = 960,
         height = 600)
    stdout <-
      capture.output(print(temp_plot))
    dev.off()
    saveRDS(temp_plot,
            file = paste(figures_folder,
                         "/gg_tweedie_index_plot.RDS",
                         sep = ""))

  }
