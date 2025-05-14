#' Generic function for stericHindrance
#' @param object An object
#' @return The stericHindrance value
#' @examples
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
#' )
#' stericHindrance(expRates)
#' @export
setGeneric("stericHindrance", function(object)
    standardGeneric("stericHindrance"))