utils::globalVariables(c("gene_id", "query", "count", "score"))

#' @title Base class for TranscriptionRates objects
#'
#' @description
#' Virtual base class that defines the common interface for transcription rate
#' objects. Both \code{ExperimentTranscriptionRates} and
#' \code{SimulationTranscriptionRates} inherit from this class.
#'
#' @slot rates a \code{\link[tibble]{tbl_df}} containing the estimated rates
#' @slot name a character value for the name of the experiment
#' @slot stericHindrance a logical value indicating whether steric hindrance
#' was modeled
#'
#' @name TranscriptionRates-class
#' @rdname TranscriptionRates-class
#' @importClassesFrom tibble tbl_df
#' @importFrom methods slot
#' @importFrom ggplot2 ggsave ggplot geom_histogram labs theme_minimal expansion
#' @importFrom ggplot2 geom_point geom_abline annotate theme_bw theme
#' @importFrom ggplot2 geom_density geom_smooth element_text geom_density_2d aes
#' @importFrom ggplot2 geom_vline geom_tile scale_fill_manual element_blank
#' @importFrom ggplot2 scale_color_gradient geom_segment scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous geom_line xlim ylim geom_violin
#' @importFrom ggplot2 geom_boxplot scale_color_manual geom_hline element_line
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom grid unit
#' @importFrom rlang .data sym
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame coverage strand<-
#' @importFrom GenomicRanges promoters findOverlaps start end width strand
#' @importFrom grDevices nclass.Sturges
#' @importFrom stats cor
#' @importFrom MASS kde2d
#' @exportClass TranscriptionRates
methods::setClass("TranscriptionRates",
    slots = c(
        rates = "tbl_df",
        name = "character",
        stericHindrance = "logical"
    ),
    contains = "VIRTUAL"
)

#' @title Generic function for estimating transcription rates
#'
#' @description
#' Generic function that estimates transcription rates from either simulation
#' data (SimulatePolymerase object) or experimental data (bigwig files and
#' genomic regions).
#'
#' @param x The input data (either a SimulatePolymerase object or bigwig files)
#' @param ... Additional arguments passed to the specific methods
#' @return An object containing estimated transcription rates
#' 
#' @examples
#' load(system.file("extdata", "granges_for_read_counting_DLD1_chr21.RData",
#' package = "STADyUM"))
#' expRates <- estimateTranscriptionRates(system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
#' bigwigMinus = system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw", package = "STADyUM"),
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control",
#'     stericHindrance = TRUE,
#'     omegaScale = 1000
#' )
#' @export
setGeneric("estimateTranscriptionRates", function(x, ...) {
    standardGeneric("estimateTranscriptionRates")
})

#' @title Accessor for estimated rates
#'
#' @description
#' Generic accessor for the estimated rates from any TranscriptionRates object
#'
#' @param object A TranscriptionRates object
#' @return A tibble containing the estimated rates
#' 
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load(system.file("extdata", 
#' "granges_for_read_counting_DLD1_chr21.RData", package = "STADyUM"))
#' expRates <- estimateTranscriptionRates(system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
#' bigwigMinus = system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw", package = "STADyUM"),
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control",
#'     stericHindrance = TRUE,
#'     omegaScale = 1000
#' )
#' rates(expRates)
#' 
#' @export
setGeneric("rates", function(object) standardGeneric("rates"))

#' @title Accessor for steric hindrance flag
#'
#' @description
#' Generic accessor for the steric hindrance flag from any TranscriptionRates
#' object
#'
#' @param object A TranscriptionRates object
#' @return A logical value indicating whether steric hindrance was modeled
#' 
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load(system.file("extdata", "granges_for_read_counting_DLD1_chr21.RData",
#' package = "STADyUM"))
#' expRates <- estimateTranscriptionRates(system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
#' bigwigMinus = system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw", package = "STADyUM"),
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control",
#'     stericHindrance = TRUE,
#'     omegaScale = 1000
#' )
#' stericHindrance(expRates)
#' 
#' @export
setGeneric("stericHindrance",
    function(object) standardGeneric("stericHindrance"))

## Plotting utilities

#' @title Plot Mean Pause Site Distribution
#'
#' @description
#' Creates a histogram plot showing the distribution of observed mean pause site
#' positions across all genes. This visualization helps identify the range and
#' shape of pause site positions
#'
#'
#' @param object an \code{\linkS4class{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotMeanPauseDistrib", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotMeanPauseDistrib")
})
#' @rdname TranscriptionRates-class
setMethod(
    "plotMeanPauseDistrib", "TranscriptionRates",
    function(object, file = NULL, width = 8, height = 6, dpi = 300) {
        cr <- rates(object)
        p <- ggplot(cr, aes(x = .data$fkMean)) +
            geom_histogram(
                bins = nclass.Sturges(cr$fkMean),
                fill = "#56B4E9", alpha = 0.8,
                color = "white", linewidth = 0.1
            ) +
            labs(
                x = "Mean Pause Site Position (bp)",
                y = "Count",
                title = "Distribution of Mean Pause Site Positions",
                subtitle = paste("n =", nrow(cr), "genes")
            ) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(size = 10, color = "gray50", hjust = 0.5),
                axis.title = element_text(size = 11, face = "bold"),
                axis.text = element_text(size = 10),
                axis.line = element_line(color = "black", linewidth = 0.5)
            )

        if (!is.null(file)) {
            ggsave(file, p, width = width, height = height, dpi = dpi)
        }

        return(p)
    }
)

#' @title Plot Expected vs Actual Pause Site Counts
#'
#' @description
#' Creates a scatter plot comparing actual pause site counts against
#' expected pause site counts from the EM algorithm. It is comparing the number
#' of polymerase at the pause site in the simulated/experimental data against
#' the number of polymerase at the pause site expected by the model. This
#' visualization assesses the goodness-of-fit of the pause site model by
#' showing how well the model predictions align with the actual data. A perfect
#' fit would show all points on the diagonal line (y=x). The R^2 value is
#' calculated and displayed on the plot to quantify the model fit quality. This
#' plot is useful for validating the accuracy of the pause site estimation and
#' identifying any systematic biases in the model predictions.
#'
#'
#' @param object an \code{\linkS4class{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotExpectedVsActualPauseSiteCounts", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotExpectedVsActualPauseSiteCounts")
})
#' @rdname TranscriptionRates-class
setMethod(
    "plotExpectedVsActualPauseSiteCounts", "TranscriptionRates",
    function(object, file = NULL, width = 8, height = 6, dpi = 300) {
        cr <- rates(object)

        allData <- data.frame(
            actual = unlist(cr$actualPauseSiteCounts),
            expected = unlist(cr$expectedPauseSiteCounts)
        )

        rSquared <- cor(allData$actual, allData$expected)^2
        r2Text <- paste("R^2 =", round(rSquared, 3))

        p <- ggplot(allData, aes(x = .data$actual, y = .data$expected)) +
            geom_point(alpha = 0.6, size = 0.8) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
            color = "red") +
            annotate("text",
                x = max(allData$actual) * 0.05,
                y = max(allData$expected) * 0.95,
                label = r2Text,
                size = 4, fontface = "bold", hjust = 0
            ) +
            labs(
                x = "Actual Pause Site Counts",
                y = "Expected Pause Site Counts",
                title = "Model Fit: Actual vs Expected Pause Site Counts"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
            )

        if (!is.null(file)) {
            ggsave(file, p,
                width = width, height = height,
                dpi = dpi
            )
        }
        return(p)
    }
)

#' @title Plot Chi Distribution
#'
#' @description
#' Creates a density plot showing the distribution of gene body RNAP density
#' (chi) across all genes. This visualization helps identify the range and
#' shape of RNA polymerase density in gene bodies, which can reveal patterns in
#' transcriptional activity.
#'
#'
#' @param object an \code{\linkS4class{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotChiDistrib", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotChiDistrib")
})
#' @rdname TranscriptionRates-class
setMethod(
    "plotChiDistrib", "TranscriptionRates",
    function(object, file = NULL, width = 8, height = 6, dpi = 300) {
        cr <- rates(object)

        p <- ggplot(cr, aes(x = .data$chi)) +
            geom_density(fill = "#56B4E9", alpha = 0.7) +
            labs(
                x = "RNAP Density (chi)",
                y = "Density",
                title = "Distribution of Gene Body RNAP Density"
            ) +
            theme_minimal() +
            theme(
                plot.title = element_text(hjust = 0.5)
            )

        if (!is.null(file)) {
            ggsave(file, p,
                width = width, height = height,
                dpi = dpi
            )
        }
        return(p)
    }
)

#' @title Plot Beta vs Chi
#'
#' @description
#' Plot a scatter plot with gene body RNAP density on the x-axis and beta (ratio
#' of gene body RNAP density to pause region RNAP density) on the y-axis. Fits a
#' linear model to the data and plots the line. Can plot beta for either the
#' adapted model or the single pause site model.
#'
#' @param object an \code{\linkS4class{TranscriptionRates}} object
#' @param label_x the x-axis label
#' @param label_y the y-axis label
#' @param xlim the x-axis limits
#' @param ylim the y-axis limits
#' @param file the path to a file to save the plot to
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotBetaVsChi", function(
    object, label_x, label_y, xlim = NULL, ylim = NULL,
    file = NULL) {
    standardGeneric("plotBetaVsChi")
})
#' @rdname TranscriptionRates-class
setMethod("plotBetaVsChi", "TranscriptionRates",
    function(object, label_x, label_y, xlim = NULL, ylim = NULL, file = NULL) {
        cr <- rates(object)

       p <- ggplot(cr,
              aes(x = log10(.data$betaAdp), y = log10(.data$chi))) +
              geom_pointdensity(adjust = 0.3, size = 0.6) +
              labs(
                x     = expression(log[10] ~ beta),
                y     = expression(log[10] ~ chi),
                color = "Density"
              ) +
            theme(
            axis.title = element_text(size = 14),
            axis.text  = element_text(size = 14),
            strip.text = element_text(size = 14),
            legend.title =  element_text(size = 10),
            legend.text  = element_text(size = 8),
            legend.key.size = unit(0.4, "cm"),
            legend.spacing = unit(0.2, "cm")

            )

            if (!is.null(xlim) || !is.null(ylim)) {
                p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
            }

        if (!is.null(file)) {
            ggsave(file, p, width = 3.5, height = 3)
        }
        return(p)
    }
)


#' @title Plot Scatter with Point Density
#'
#' @description
#' Creates a scatter plot colored by local point density for any two columns in
#' the rates tibble, with optional log10 transformation on either axis. When a
#' discrete \code{color_var} is supplied the points are colored by that variable
#' instead of by density.
#'
#' @param object an \code{\linkS4class{TranscriptionRates}} object
#' @param xvar character; column name in \code{rates(object)} for the x-axis
#' @param yvar character; column name in \code{rates(object)} for the y-axis
#' @param xlab x-axis label (string or \code{expression}); defaults to
#'   \code{xvar} with a \code{log10()} prefix when \code{log_x = TRUE}
#' @param ylab y-axis label (string or \code{expression}); same rules as
#'   \code{xlab}
#' @param log_x logical; apply \code{log10} to x before plotting
#' @param log_y logical; apply \code{log10} to y before plotting
#' @param color_var optional character; column name for discrete color grouping;
#'   when supplied \code{geom_point} is used instead of
#'   \code{geom_pointdensity}
#' @param color_values optional named character vector passed to
#'   \code{scale_color_manual} when \code{color_var} is set
#' @param color_lab legend title for the color scale; defaults to
#'   \code{"Density"} or \code{color_var}
#' @param xlim optional numeric(2); forwarded to \code{coord_cartesian}
#' @param ylim optional numeric(2); forwarded to \code{coord_cartesian}
#' @param file optional file path to save the plot
#' @param width plot width in inches
#' @param height plot height in inches
#' @param dpi plot resolution
#'
#' @return a \code{\link[ggplot2]{ggplot}} object
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotScatterDensity", function(
    object, xvar, yvar,
    xlab = NULL, ylab = NULL,
    log_x = FALSE, log_y = FALSE,
    color_var = NULL, color_values = NULL, color_lab = NULL,
    xlim = NULL, ylim = NULL,
    file = NULL, width = 3.5, height = 3, dpi = 300) {
    standardGeneric("plotScatterDensity")
})
#' @rdname TranscriptionRates-class
setMethod("plotScatterDensity", "TranscriptionRates",
    function(object, xvar, yvar,
             xlab = NULL, ylab = NULL,
             log_x = FALSE, log_y = FALSE,
             color_var = NULL, color_values = NULL, color_lab = NULL,
             xlim = NULL, ylim = NULL,
             file = NULL, width = 3.5, height = 3, dpi = 300) {
        cr <- rates(object)

        x_vals <- if (log_x) log10(cr[[xvar]]) else cr[[xvar]]
        y_vals <- if (log_y) log10(cr[[yvar]]) else cr[[yvar]]
        plot_df <- data.frame(x = x_vals, y = y_vals)

        if (is.null(xlab)) xlab <- if (log_x) paste0("log10(", xvar, ")") else xvar
        if (is.null(ylab)) ylab <- if (log_y) paste0("log10(", yvar, ")") else yvar

        if (!is.null(color_var)) {
            plot_df$color_group <- cr[[color_var]]
            if (is.null(color_lab)) color_lab <- color_var
            p <- ggplot(plot_df, aes(x = .data$x, y = .data$y,
                                     color = .data$color_group)) +
                geom_point(size = 0.8, alpha = 0.7)
            if (!is.null(color_values)) {
                p <- p + scale_color_manual(values = color_values)
            }
        } else {
            if (is.null(color_lab)) color_lab <- "Density"
            p <- ggplot(plot_df, aes(x = .data$x, y = .data$y)) +
                geom_pointdensity(adjust = 0.3, size = 0.6)
        }

        p <- p +
            labs(x = xlab, y = ylab, color = color_lab) +
            theme(
                axis.title      = element_text(size = 14),
                axis.text       = element_text(size = 14),
                strip.text      = element_text(size = 14),
                legend.title    = element_text(size = 10),
                legend.text     = element_text(size = 8),
                legend.key.size = unit(0.4, "cm"),
                legend.spacing  = unit(0.2, "cm")
            )

        if (!is.null(xlim) || !is.null(ylim)) {
            p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
        }

        if (!is.null(file)) {
            ggsave(file, p, width = width, height = height, dpi = dpi)
        }
        return(p)
    }
)

#' @title Plot pause site contour map
#'
#' @description
#' Plot a contour map with mean pause site position on the x-axis and pause site
#' variance on the y-axis.
#'
#' @param object an \code{\linkS4class{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load(system.file("extdata", "granges_for_read_counting_DLD1_chr21.RData",
#' package = "STADyUM"))
#' expRates <- estimateTranscriptionRates(system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
#' bigwigMinus = system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw", package = "STADyUM"),
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control"
#' )
#' plotPauseSiteContourMap(expRates, file="pause_sites_contour_map.png")
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotPauseSiteContourMap", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotPauseSiteContourMap")
})
#' @rdname TranscriptionRates-class
setMethod(
    "plotPauseSiteContourMap", "TranscriptionRates",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {
        cr <- rates(object)

        p <- ggplot(cr, aes(x = .data$fkMean, y = .data$fkVar)) +
            geom_density_2d(color = "blue", linewidth = 0.8) +
            geom_point(alpha = 0.6, size = 1.5, color = "#E69F00") +
            labs(
                x = "Mean Pause Site Position (bp)",
                y = "Pause Site Variance (bp^2)",
                title = "Pause Site Mean vs Variance Distribution"
            ) +
            theme_bw()

        if (!is.null(file)) {
            ggsave(file, p,
                width = width, height = height,
                dpi = dpi
            )
        }
        return(p)
    }
)