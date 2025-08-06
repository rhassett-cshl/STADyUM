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
#' @import ggplot2
#' @importFrom grDevices nclass.Sturges
#' @import ggpubr
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
#' load("inst/extdata/granges_for_read_counting_DLD1_chr21.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control",
#'     stericHindrance = TRUE,
#'     omegaScale = 1000,
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
#' load("inst/extdata/granges_for_read_counting_DLD1_chr21.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control",
#'     stericHindrance = TRUE,
#'     omegaScale = 1000,
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
#' load("inst/extdata/granges_for_read_counting_DLD1_chr21.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control",
#'     stericHindrance = TRUE,
#'     omegaScale = 1000,
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
#' @param object an \code{\link{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_DLD1_chr21.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control"
#' )
#' plotMeanPauseDistrib(expRates, file="mean_pause_distrib.png")
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotMeanPauseDistrib", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotMeanPauseDistrib")
})

setMethod(
    "plotMeanPauseDistrib", "TranscriptionRates",
    function(object, file = NULL, width = 8, height = 6, dpi = 300) {
        cr <- rates(object)
        p <- gghistogram(
            cr, x = "fkMean",
            add = "mean",
            bins = nclass.Sturges(cr$fkMean),
            fill = "#56B4E9", alpha = 0.8,
            color = "white", size = 0.1
        ) +
            labs(
                x = "Mean Pause Site Position (bp)",
                y = "Count",
                title = "Distribution of Mean Pause Site Positions",
                subtitle = paste("n =", nrow(cr), "genes")
            ) +
            theme_pubr() +
            theme(
                plot.title = element_text(size = 14, face = "bold", 
                hjust = 0.5), plot.subtitle = element_text(size = 10, 
                color = "gray50", hjust = 0.5),
                axis.title = element_text(size = 11, face = "bold"),
                axis.text = element_text(size = 10),
                axis.line = element_line(color = "black", size = 0.5)
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
#' fit would show all points on the diagonal line (y=x). The R² value is
#' calculated and displayed on the plot to quantify the model fit quality. This
#' plot is useful for validating the accuracy of the pause site estimation and
#' identifying any systematic biases in the model predictions.
#'
#'
#' @param object an \code{\link{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_DLD1_chr21.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control"
#' )
#' plotExpectedVsActualPauseSiteCounts(expRates,
#' file="expected_vs_actual_pause_site_counts.png")
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotExpectedVsActualPauseSiteCounts", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotExpectedVsActualPauseSiteCounts")
})

setMethod(
    "plotExpectedVsActualPauseSiteCounts", "TranscriptionRates",
    function(object, file = NULL, width = 8, height = 6, dpi = 300) {
        cr <- rates(object)

        allData <- data.frame(
            actual = unlist(cr$actualPauseSiteCounts),
            expected = unlist(cr$expectedPauseSiteCounts)
        )

        rSquared <- cor(allData$actual, allData$expected)^2
        r2Text <- paste("R² =", round(rSquared, 3))

        p <- ggplot(allData, aes(x = actual, y = expected)) +
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
#' @param object an \code{\link{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_DLD1_chr21.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control"
#' )
#' plotChiDistrib(expRates, file="chi_distrib.png")
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotChiDistrib", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotChiDistrib")
})

setMethod(
    "plotChiDistrib", "TranscriptionRates",
    function(object, file = NULL, width = 8, height = 6, dpi = 300) {
        cr <- rates(object)

        p <- ggdensity(
            cr, x = "chi",
            fill = "#56B4E9", alpha = 0.7
        ) +
            labs(
                x = "RNAP Density (chi)",
                y = "Density",
                title = "Distribution of Gene Body RNAP Density"
            ) +
            theme_pubr() +
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
#' @param object an \code{\link{TranscriptionRates}} object
#' @param betaType the type of beta to plot. Can be "betaAdp" for the adapted
#' model or "betaOrg" for the single pause site model. Defaults to "betaAdp".
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_DLD1_chr21.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control"
#' )
#' plotBetaVsChi(expRates, betaType = "betaAdp", file="beta_vs_chi.png")
#'
#' @rdname TranscriptionRates-class
#' @export
setGeneric("plotBetaVsChi", function(
    object, betaType = "betaAdp",
    file = NULL, width = 8, height = 6, dpi = 300, ...) {
    standardGeneric("plotBetaVsChi")
})

setMethod("plotBetaVsChi", "TranscriptionRates",
    function(object, betaType = "betaAdp", file = NULL,
        width = 8, height = 6, dpi = 300) {
        cr <- rates(object)

        if (!betaType %in% c("betaAdp", "betaOrg")) {
            stop("betaType must be either 'betaAdp' or 'betaOrg'")
        }

        # Set y-axis label based on beta type
        yLabel <- if (betaType == "betaAdp") {
            "Pause Escape Rate (betaAdp)"
        } else {
            "Pause Escape Rate (betaOrg)"
        }

        titleText <- if (betaType == "betaAdp") {
            "Gene Activity vs Pause Escape Rate (Adapted Model)"
        } else {
            "Gene Activity vs Pause Escape Rate (Single Pause Site)"
        }

        p <- ggplot(cr, aes(x = chi, y = !!sym(betaType))) +
            geom_point(alpha = 0.7, color = "#CC79A7") +
            geom_smooth(method = "loess", se = TRUE, color = "red") +
            labs(
                x = "Gene Body RNAP Density (chi)",
                y = yLabel,
                title = titleText
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
            )

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
#' @param object an \code{\link{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_DLD1_chr21.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
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

setMethod(
    "plotPauseSiteContourMap", "TranscriptionRates",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {
        cr <- rates(object)

        p <- ggplot(cr, aes(x = fkMean, y = fkVar)) +
            geom_density_2d(color = "blue", size = 0.8) +
            geom_point(alpha = 0.6, size = 1.5, color = "#E69F00") +
            labs(
                x = "Mean Pause Site Position (bp)",
                y = "Pause Site Variance (bp²)",
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