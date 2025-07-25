% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base-transcription-rates-class.R
\name{rates}
\alias{rates}
\title{Accessor for estimated rates}
\usage{
rates(object)
}
\arguments{
\item{object}{A TranscriptionRates object}
}
\value{
A tibble containing the estimated rates
}
\description{
Generic accessor for the estimated rates from any TranscriptionRates object
}
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
rates(expRates)

}
