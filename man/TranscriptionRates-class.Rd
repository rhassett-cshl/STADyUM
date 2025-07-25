% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base-transcription-rates-class.R
\docType{class}
\name{TranscriptionRates-class}
\alias{TranscriptionRates-class}
\alias{plotMeanPauseDistrib}
\alias{plotExpectedVsActualPauseSiteCounts}
\alias{plotChiDistrib}
\alias{plotBetaVsChi}
\alias{plotPauseSiteContourMap}
\title{Base class for TranscriptionRates objects}
\usage{
plotMeanPauseDistrib(object, file = NULL, width = 8, height = 6, dpi = 300)

plotExpectedVsActualPauseSiteCounts(
  object,
  file = NULL,
  width = 8,
  height = 6,
  dpi = 300
)

plotChiDistrib(object, file = NULL, width = 8, height = 6, dpi = 300)

plotBetaVsChi(
  object,
  betaType = "betaAdp",
  file = NULL,
  width = 8,
  height = 6,
  dpi = 300,
  ...
)

plotPauseSiteContourMap(object, file = NULL, width = 8, height = 6, dpi = 300)
}
\arguments{
\item{object}{an \code{\link{TranscriptionRates}} object}

\item{file}{the path to a file to save the plot to}

\item{width}{the width of the plot in inches}

\item{height}{the height of the plot in inches}

\item{dpi}{the resolution of the plot in dpi}

\item{betaType}{the type of beta to plot. Can be "betaAdp" for the adapted
model or "betaOrg" for the single pause site model. Defaults to "betaAdp".}
}
\value{
an \code{\link{ggplot2}} object

an \code{\link{ggplot2}} object

an \code{\link{ggplot2}} object

an \code{\link{ggplot2}} object

an \code{\link{ggplot2}} object
}
\description{
Virtual base class that defines the common interface for transcription rate
objects. Both \code{ExperimentTranscriptionRates} and
\code{SimulationTranscriptionRates} inherit from this class.

Creates a histogram plot showing the distribution of observed mean pause site
positions across all genes. This visualization helps identify the range and
shape of pause site positions

Creates a scatter plot comparing actual pause site counts against
expected pause site counts from the EM algorithm. It is comparing the number
of polymerase at the pause site in the simulated/experimental data against
the number of polymerase at the pause site expected by the model. This
visualization assesses the goodness-of-fit of the pause site model by
showing how well the model predictions align with the actual data. A perfect
fit would show all points on the diagonal line (y=x). The R² value is
calculated and displayed on the plot to quantify the model fit quality. This
plot is useful for validating the accuracy of the pause site estimation and
identifying any systematic biases in the model predictions.

Creates a density plot showing the distribution of gene body RNAP density
(chi) across all genes. This visualization helps identify the range and
shape of RNA polymerase density in gene bodies, which can reveal patterns in
transcriptional activity.

Plot a scatter plot with gene body RNAP density on the x-axis and beta (ratio
of gene body RNAP density to pause region RNAP density) on the y-axis. Fits a
linear model to the data and plots the line. Can plot beta for either the
adapted model or the single pause site model.

Plot a contour map with mean pause site position on the x-axis and pause site
variance on the y-axis.
}
\section{Slots}{

\describe{
\item{\code{rates}}{a \code{\link[tibble]{tbl_df}} containing the estimated rates}

\item{\code{stericHindrance}}{a logical value indicating whether steric hindrance
was modeled}
}}

\examples{
# Create an ExperimentTranscriptionRates object
load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
expRates <- estimateTranscriptionRates(
    "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
    bigwigMinus = 
     "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
    pauseRegions = bw_pause_21_subset,
    geneBodyRegions = bw_gene_body_21_subset,
    stericHindrance = TRUE,
    omegaScale = 1000,
)
plotMeanPauseDistrib(expRates, file="mean_pause_distrib.png")

# Create an ExperimentTranscriptionRates object
load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
expRates <- estimateTranscriptionRates(
    "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
    bigwigMinus = 
     "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
    pauseRegions = bw_pause_21_subset,
    geneBodyRegions = bw_gene_body_21_subset,
    stericHindrance = TRUE,
    omegaScale = 1000,
)
plotExpectedVsActualPauseSiteCounts(expRates,
file="expected_vs_actual_pause_site_counts.png")

# Create an ExperimentTranscriptionRates object
load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
expRates <- estimateTranscriptionRates(
    "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
    bigwigMinus = 
     "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
    pauseRegions = bw_pause_21_subset,
    geneBodyRegions = bw_gene_body_21_subset,
    stericHindrance = TRUE,
    omegaScale = 1000,
)
plotChiDistrib(expRates, file="chi_distrib.png")

# Create an ExperimentTranscriptionRates object
load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
expRates <- estimateTranscriptionRates(
    "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
    bigwigMinus = 
     "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
    pauseRegions = bw_pause_21_subset,
    geneBodyRegions = bw_gene_body_21_subset,
    stericHindrance = TRUE,
    omegaScale = 1000,
)
plotBetaVsChi(expRates, betaType = "betaAdp", file="beta_vs_chi.png")

# Create an ExperimentTranscriptionRates object
load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
expRates <- estimateTranscriptionRates(
    "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
    bigwigMinus = 
     "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
    pauseRegions = bw_pause_21_subset,
    geneBodyRegions = bw_gene_body_21_subset,
    stericHindrance = TRUE,
    omegaScale = 1000,
)
plotPauseSiteContourMap(expRates, file="pause_sites_contour_map.png")

}
