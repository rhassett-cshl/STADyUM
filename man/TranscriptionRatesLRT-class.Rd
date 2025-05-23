% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transcription_rates_lrt-class.R
\docType{class}
\name{TranscriptionRatesLRT-class}
\alias{TranscriptionRatesLRT-class}
\alias{expData1}
\alias{expData2}
\alias{spikeInScalingFactor}
\title{Class TranscriptionRatesLRT for comparing aspects of transcriptional
dynamics of the same TU under different conditions using Likelihood Ratio
Test Statistics}
\usage{
expData1(object)

expData2(object)

spikeInScalingFactor(object)
}
\description{
Class \code{TranscriptionRatesLRT}
}
\examples{
expData1 <- estimateExperimentTranscriptionRates(
    bigwigPlus = "path/to/plus.bw",
    bigwigMinus = "path/to/minus.bw",
    pauseRegions = GRanges("chr1:1-1000"),
    geneBodyRegions = GRanges("chr1:1-2000"),
    geneNameColumn = "gene_id"
)
expData2 <- estimateExperimentTranscriptionRates(
    bigwigPlus = "path/to/plus.bw",
    bigwigMinus = "path/to/minus.bw",
    pauseRegions = GRanges("chr1:1-1000"),
    geneBodyRegions = GRanges("chr1:1-2000"),
    geneNameColumn = "gene_id"
)
lrts <- likelihoodRatioTest(expData1, expData2, "spikeInScalingFactor.csv")
expData1(lrts)
expData1 <- estimateExperimentTranscriptionRates(
    bigwigPlus = "path/to/plus.bw",
    bigwigMinus = "path/to/minus.bw",
    pauseRegions = GRanges("chr1:1-1000"),
    geneBodyRegions = GRanges("chr1:1-2000"),
    geneNameColumn = "gene_id"
)
expData2 <- estimateExperimentTranscriptionRates(
    bigwigPlus = "path/to/plus.bw",
    bigwigMinus = "path/to/minus.bw",
    pauseRegions = GRanges("chr1:1-1000"),
    geneBodyRegions = GRanges("chr1:1-2000"),
    geneNameColumn = "gene_id"
)
lrts <- likelihoodRatioTest(expData1, expData2, "spikeInScalingFactor.csv")
expData2(lrts)
expData1 <- estimateExperimentTranscriptionRates(
    bigwigPlus = "path/to/plus.bw",
    bigwigMinus = "path/to/minus.bw",
    pauseRegions = GRanges("chr1:1-1000"),
    geneBodyRegions = GRanges("chr1:1-2000"),
    geneNameColumn = "gene_id"
)
expData2 <- estimateExperimentTranscriptionRates(
    bigwigPlus = "path/to/plus.bw",
    bigwigMinus = "path/to/minus.bw",
    pauseRegions = GRanges("chr1:1-1000"),
    geneBodyRegions = GRanges("chr1:1-2000"),
    geneNameColumn = "gene_id"
)
lrts <- likelihoodRatioTest(expData1, expData2, "spikeInScalingFactor.csv")
spikeInScalingFactor(lrts)
}
