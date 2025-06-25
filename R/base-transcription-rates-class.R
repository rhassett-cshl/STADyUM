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

#' @title Accessor for estimated rates
#'
#' @description
#' Generic accessor for the estimated rates from any TranscriptionRates object
#'
#' @param object A TranscriptionRates object
#' @return A tibble containing the estimated rates
#' @export
#' @examples
#' # Works with both ExperimentTranscriptionRates and SimulationTranscriptionRates
#' rates(expRates)
#' rates(simRates)
setGeneric("rates", function(object) standardGeneric("rates"))

#' @title Accessor for steric hindrance flag
#'
#' @description
#' Generic accessor for the steric hindrance flag from any TranscriptionRates object
#'
#' @param object A TranscriptionRates object
#' @return A logical value indicating whether steric hindrance was modeled
#' @export
#' @examples
#' # Works with both ExperimentTranscriptionRates and SimulationTranscriptionRates
#' stericHindrance(expRates)
#' stericHindrance(simRates)
setGeneric("stericHindrance", function(object) standardGeneric("stericHindrance"))