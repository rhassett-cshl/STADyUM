

#' Class simulate_polymerase
#'
#' Class \code{simulate_polymerase} has
#'
#' @slot counts a \code{data.frame} with five columns gene_id, summarized_pause_counts, pause_length
#' summarized_gb_counts, gb_length

#' @name simulate_polymerase-class
#' @rdname simulate_polymerase-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @exportClass simulate_polymerase
methods::setClass("simulate_polymerase",
                  slots = c(counts = "data.frame"),
                  validity = experiment_transcription_rates_valid
)