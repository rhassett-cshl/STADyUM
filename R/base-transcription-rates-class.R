#' @title Base class for TranscriptionRates objects
#'
#' @description
#' Virtual base class that defines the common interface for transcription rate
#' objects. Both \code{ExperimentTranscriptionRates} and
#' \code{SimulationTranscriptionRates} inherit from this class.
#'
#' @slot rates a \code{\link[tibble]{tbl_df}} containing the estimated rates
#' @slot stericHindrance a logical value indicating whether steric hindrance
#' was modeled
#'
#' @name TranscriptionRates-class
#' @rdname TranscriptionRates-class
#' @importClassesFrom tibble tbl_df
#' @importFrom methods slot
#' @exportClass TranscriptionRates
methods::setClass("TranscriptionRates",
    slots = c(
        rates = "tbl_df",
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
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
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
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
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
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     stericHindrance = TRUE,
#'     omegaScale = 1000,
#' )
#' stericHindrance(expRates)
#' 
#' @export
setGeneric("stericHindrance",
    function(object) standardGeneric("stericHindrance"))
