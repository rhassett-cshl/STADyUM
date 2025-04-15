

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


#' simulate_polymerase
#'
#' Estimates 
#' @param k a numeric (int) value for the mean of pause sites across cells. Defaults to 50bp.
#' @param kMin a numeric (int) for the lower bound of allowed pause site. Defaults to 17bp. 
#' @param kMax a numeric (int) for the lower bound of allowed pause site. Defaults to 200bp. 
#' @param kSd a numeric (double) value for std dev of pause sites across cells. Defaults to 0bp.
#' @param gene_len a numeric (int) for the length of the gene. Defaults to 1000bp.
#' @param alpha a numeric (double) for the rate of transcription initiation. Defaults to 1 event per min.
#' @param beta a numeric (double) for the pause release rate. Defaults to 1 event per min.
#' @param zeta a numeric (double) for the mean of elongation rates across sites. Defaults to 2000bp per min.
#' @param zeta_sd a numeric (double) value for std dev of elongation rates across sites. Defaults to 1000.
#' @param zeta_max a numeric (double) for the max elongation rate allowed. Defaults to 2500bp per min.
#' @param zeta_min a numeric (double) for the minimum of elongation rate allowed. Defaults to 1500bp per min.
#' @param zeta_vec TODO
#' @param total_cells a numeric (int) for the total number of cells to simulate. Defaults to 10.
#' @param s a numeric (int) for the polymerase II size. Defaults to 33bp.
#' @param addSpace a numeric (int) for additional space in addition to RNAP size. Defaults to 17bp.
#' @param time a numeric (double) for the total time to simulate. Defaults to 0.1 min.
#' @param csv a numeric (int) for the number of steps to record. Defaults to 1 step.
#' @param output_dir a character (string) for the output directory. Defaults to TODO.
#' 
#' @return an \code{\link{simulate_polymerase-class}} object
#'
#' @export
simulate_polymerase <- function(k, k_min, k_max, ksd, gene_len, alpha, beta, zeta, zeta_sd, zeta_max, zeta_min, zeta_vec, total_cells, s, h, time, delta_t, csv_steps_to_record, output_dir) {
  # Input validation
  if (!is.numeric(k) || k <= 0) stop("k must be a positive integer")
  if (!is.numeric(k_min) || k_min <= 0) stop("k_min must be a positive integer")
  if (!is.numeric(k_max) || k_max <= 0) stop("k_max must be a positive integer")
  if (!is.numeric(ksd) || ksd < 0) stop("ksd must be non-negative")
  if (!is.numeric(gene_len) || gene_len <= 0) stop("gene_len must be positive")
  if (!is.numeric(alpha) || alpha <= 0) stop("alpha must be positive")
  if (!is.numeric(beta) || beta <= 0) stop("beta must be positive")
  if (!is.numeric(zeta) || zeta <= 0) stop("zeta must be positive")
  if (!is.numeric(zeta_sd) || zeta_sd < 0) stop("zeta_sd must be non-negative")
  if (!is.numeric(zeta_max) || zeta_max <= 0) stop("zeta_max must be positive")
  if (!is.numeric(zeta_min) || zeta_min <= 0) stop("zeta_min must be positive")
  if (!is.numeric(total_cells) || total_cells <= 0) stop("total_cells must be positive")
  if (!is.numeric(s) || s < 0) stop("s must be non-negative")
  if (!is.numeric(h) || h < 0) stop("h must be non-negative")
  if (!is.numeric(time) || time <= 0) stop("time must be positive")
  if (!is.numeric(delta_t) || delta_t <= 0) stop("delta_t must be positive")
  if (!is.numeric(csv_steps_to_record) || csv_steps_to_record < 0) stop("csv_steps_to_record must be non-negative")
# output dir? zeta_vec?

  # Call C++ function
  result <- .Call('_STADyUM_simulate_polymerase_cpp',
                 as.integer(k),
                 as.integer(k_min),
                 as.integer(k_max),
                 as.double(ksd),
                 as.integer(gene_len),
                 as.double(alpha),
                 as.double(beta),
                 as.integer(zeta),
                 as.double(zeta_sd),
                 as.double(zeta_max),
                 as.double(zeta_min),
                 as.integer(total_cells),
                 as.integer(s),
                 as.integer(h),
                 as.double(time),
                 as.double(delta_t),
                 as.integer(csv_steps_to_record),
                 PACKAGE = 'STADyUM')

  return(result)
}