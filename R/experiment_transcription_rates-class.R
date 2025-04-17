#' @importFrom rtracklayer import.bw
#' @import GenomicRanges
#' @import tibble
#' @importFrom methods new
#' @importFrom S4Vectors elementMetadata DataFrame splitAsList
#' @importFrom dplyr mutate select %>% left_join
#' @importFrom purrr map map_dbl
#' @importFrom stats dnorm
#' @importFrom utils write.csv
#' @import ggplot2
#' @import progress
#' 
#' @title Check validity of experiment_transcription_rates object
#' @description Validates an experiment_transcription_rates object
#' @param object A \code{experiment_transcription_rates} object
#' @return TRUE if valid, else errors
#' @keywords internal
experiment_transcription_rates_valid <- function(object) {
  errors <- c()

  # Check counts data.frame structure
  if (!is.data.frame(object@counts)) {
    errors <- c(errors, "counts must be a data.frame")
  } else {
    required_cols <- c("gene_id", "summarized_pause_counts", "pause_length",
                      "summarized_gb_counts", "gb_length")
    missing_cols <- setdiff(required_cols, colnames(object@counts))
    if (length(missing_cols) > 0) {
      errors <- c(errors, paste("counts missing required columns:", 
                              paste(missing_cols, collapse = ", ")))
    }
  }

  # Check bigwig file paths
  if (!file.exists(object@bigwig_plus)) {
    errors <- c(errors, "bigwig_plus file does not exist")
  }
  if (!file.exists(object@bigwig_minus)) {
    errors <- c(errors, "bigwig_minus file does not exist")
  }

  # Check GRanges objects
  if (!inherits(object@pause_regions, "GRanges")) {
    errors <- c(errors, "pause_regions must be a GRanges object")
  }
  if (!inherits(object@gene_body_regions, "GRanges")) {
    errors <- c(errors, "gene_body_regions must be a GRanges object")
  }

  # Check gene_name_column
  if (!is.character(object@gene_name_column) || length(object@gene_name_column) != 1) {
    errors <- c(errors, "gene_name_column must be a single character string")
  }

  # Check steric_hindrance
  if (!is.logical(object@steric_hindrance) || length(object@steric_hindrance) != 1) {
    errors <- c(errors, "steric_hindrance must be a single logical value")
  }

  # Check omega_scale
  if (object@steric_hindrance) {
    if (!is.numeric(object@omega_scale) || length(object@omega_scale) != 1 || object@omega_scale <= 0) {
      errors <- c(errors, "omega_scale must be a positive numeric value when steric_hindrance is TRUE")
    }
  }

  # Check rates tbl_df
  if (!inherits(object@rates, "tbl_df")) {
    errors <- c(errors, "rates must be a tbl_df")
  } else {
    required_cols <- c("gene_id", "chi", "beta_org", "beta_adp", "fk_mean", "fk_var")
    if (object@steric_hindrance) {
      required_cols <- c(required_cols, "phi", "omega_zeta", "beta_zeta", "alpha_zeta")
    }
    missing_cols <- setdiff(required_cols, colnames(object@rates))
    if (length(missing_cols) > 0) {
      errors <- c(errors, paste("rates missing required columns:", 
                              paste(missing_cols, collapse = ", ")))
    }
  }

  if (length(errors) == 0) TRUE else errors
}

#' Class experiment_transcription_rates
#'
#' Class \code{experiment_transcription_rates} has read counts, pause and gene body genomic region
#' coordinates, steric hindrance and omega scale factor parameters used to estimate the transcription 
#' rates
#'
#' @slot counts a \code{data.frame} with five columns gene_id, summarized_pause_counts, pause_length
#' summarized_gb_counts, gb_length
#' @slot bigwig_plus a path to bigwig for plus strand
#' @slot bigwig_minus a path to bigwig for minus strand
#' @slot pause_regions a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the pause region coordinates
#' @slot gene_body_regions a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the gene body region coordinates
#' @slot gene_name_column a string for the gene name column in the GRanges objects
#' @slot steric_hindrance a logical value representing whether landing-pad occupancy was inferred when
#' estimating the rates
#' @slot omega_scale a numeric value for the scale factor used to calculate omega
#' @slot rates a \code{tbl_df} containing estimated transcription rates such as chi estimates, 
#' beta_org estimates from the initial model, beta_adp estimates from the model with varying 
#' pause sites, fk_mean giving the mean position of pause sites, fk_var for variance of pause sites,
#' phi estimates for landing-pad occupancy, omega_zeta for the effective initiation rate, beta_zeta
#' for the pause-escape rate, alpha_zeta for the potential initiation rate, and likelihoods
#'
#' @name experiment_transcription_rates-class
#' @rdname experiment_transcription_rates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @exportClass experiment_transcription_rates
methods::setClass("experiment_transcription_rates",
                  slots = c(counts = "data.frame",
                            bigwig_plus = "character",
                            bigwig_minus = "character",
                            pause_regions = "GRanges",
                            gene_body_regions = "GRanges",
                            gene_name_column = "character",
                            steric_hindrance = "logical",
                            omega_scale = "numeric",
                            rates = "tbl_df"),
                  validity = experiment_transcription_rates_valid
)

#' estimate_experiment_transcription_rates
#'
#' Estimates the transcription rates from experimental data and contructs an object that holds
#' these rates
#' @param bigwig_plus the path to a bigwig file from the plus strand recording PRO-seq read counts
#' @param bigwig_minus the path to a bigwig file from the minus strand recording PRO-seq read counts
#' @param pause_regions a \link[GenomicRanges]{GRanges-class} object that must contain a gene_id
#' metadata column
#' @param gene_body_regions a \link[GenomicRanges]{GRanges-class} object that must contain a gene_id
#' metadata column
#' @param gene_name_column a string that indicates which column in the GRanges object contains the
#' gene names
#' @param steric_hindrance a logical value to determine whether to infer landing-pad occupancy or not.
#' Defaults to FALSE.
#' @param omega_scale a numeric value for scaling omega. Defaults to NULL.
#'
#' @return an \code{\link{experiment_transcription_rates-class}} object
#'
#' @export
estimate_experiment_transcription_rates <- function(bigwig_plus, bigwig_minus,
                                                pause_regions, gene_body_regions,
                                                gene_name_column, steric_hindrance = FALSE, 
                                                omega_scale = NULL) {
  # 1. Input validation for bigwig files
  if (!file.exists(bigwig_plus)) {
    stop("bigwig_plus file does not exist")
  }
  if (!file.exists(bigwig_minus)) {
    stop("bigwig_minus file does not exist")
  }
  if (!file.access(bigwig_plus, 4) == 0) {
    stop("bigwig_plus file is not readable")
  }
  if (!file.access(bigwig_minus, 4) == 0) {
    stop("bigwig_minus file is not readable")
  }

  # 2. Check for empty or invalid regions
  if (length(pause_regions) == 0) {
    stop("pause_regions is empty")
  }
  if (length(gene_body_regions) == 0) {
    stop("gene_body_regions is empty")
  }

  # Check regions have valid ranges
  if (any(width(pause_regions) <= 0)) {
    stop("pause_regions contains invalid (non-positive) widths")
  }
  if (any(width(gene_body_regions) <= 0)) {
    stop("gene_body_regions contains invalid (non-positive) widths")
  }

  # **Some checks prior to beginning construction**
  # Check for correct GRanges object metadata
  if (!gene_name_column %in%
      colnames(S4Vectors::elementMetadata(pause_regions))) {
    stop(paste("pause_regions does not have a column matching",
               gene_name_column))
  }
  if (!gene_name_column %in%
      colnames(S4Vectors::elementMetadata(gene_body_regions))) {
    stop(paste("gene_body_regions does not have a column matching",
               gene_name_column))
  }
  if (!is.character(
    S4Vectors::elementMetadata(pause_regions)[, gene_name_column])) {
    stop(paste("gene name column", gene_name_column,
               "must be of class 'character' in pause_regions object"))
  }
  if (!is.character(
    S4Vectors::elementMetadata(gene_body_regions)[, gene_name_column])) {
    stop(paste("gene name column", gene_name_column,
               "must be of class 'character' in gene_body_regions object"))
  }
  duplicated_pause_region_gene_names <-
    any(duplicated(S4Vectors::elementMetadata(pause_regions)[, gene_name_column]))
  if (duplicated_pause_region_gene_names) {
    stop("One or more gene names are duplicated in pause region, gene names must be unique")
  }
  duplicated_gene_body_region_gene_names <-
    any(duplicated(S4Vectors::elementMetadata(gene_body_regions)[, gene_name_column]))
  if (duplicated_gene_body_region_gene_names) {
    stop("One or more gene names are duplicated in gene body region, gene names must be unique")
  }
  if (steric_hindrance && (is.null(omega_scale) && !is.numeric(omega_scale) && omega_scale <= 0) ) {
    stop("For steric hindrance case, omega_scale parameter must be set to numeric greater than 0")
  }
  
  # **End checks**

  # Force copy of object underlying GRanges to prevent any weird side effects if
  # GRanges is using a data.table or something else that can modify in place
  pause_regions <- GenomicRanges::makeGRangesFromDataFrame(
    data.table::copy(data.table::as.data.table(pause_regions)),
    keep.extra.columns = T)

  gene_body_regions <- GenomicRanges::makeGRangesFromDataFrame(
    data.table::copy(data.table::as.data.table(gene_body_regions)),
    keep.extra.columns = T)


  # set up parameters
  k <- 50
  kmin <- 1
  kmax <- 200 # also used as k on the poisson case

  rnap_size <- 50
  zeta <- 2000

  message("processing bigwigs...")

  rc_cutoff <- 20 # read count cut-off for both gene body and pause peak

  # Add progress bar for bigwig import
  pb <- progress::progress_bar$new(
    format = "Importing bigwigs [:bar] :percent eta: :eta",
    total = 3,
    clear = FALSE
  )
  pb$tick(0)
  
  pb$tick()
  bwp1_p3 <- import.bw(bigwig_plus)
  pb$tick()
  bwm1_p3 <- import.bw(bigwig_minus)
  pb$tick()

  # 5. Check for reasonable read counts
  if (sum(bwp1_p3$score) == 0) {
    stop("No reads found in plus strand bigwig file")
  }
  if (sum(bwm1_p3$score) == 0) {
    stop("No reads found in minus strand bigwig file")
  }

  # Add progress bar for processing
  pb <- progress::progress_bar$new(
    format = "Processing bigwigs [:bar] :percent eta: :eta",
    total = 4,
    clear = FALSE
  )
  pb$tick(0)
  
  bwp1_p3 <- process_bw(bw = bwp1_p3, strand = "+")
  pb$tick()
  bwm1_p3 <- process_bw(bw = bwm1_p3, strand = "-")
  pb$tick()
  bw1_p3 <- c(bwp1_p3, bwm1_p3)
  rm(bwp1_p3, bwm1_p3)

  # make sure pause region is the same as kmax used in EM
  pause_regions <- promoters(pause_regions, upstream = 0, downstream = kmax)

  # summarize read counts
  rc1_pause <- summarise_bw(bw1_p3, pause_regions, "summarized_pause_counts")
  rc1_pause$pause_length <- kmax
  pb$tick()

  rc1_gb <- summarise_bw(bw1_p3, gene_body_regions, "summarized_gb_counts")
  rc1_gb$gb_length <-
  width(gene_body_regions)[match(rc1_gb$gene_id, gene_body_regions$gene_id)]
  pb$tick()

  # prepare read count table
  rc1 <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
              list(rc1_pause, rc1_gb))

  # clean up some genes with missing values in tss length or gene body length
  rc1 <- rc1[!(is.na(rc1$pause_length) | is.na(rc1$gb_length)), ]
  rc1 <- rc1[(rc1$summarized_pause_counts > rc_cutoff) & (rc1$summarized_gb_counts > rc_cutoff), ]

  message("estimating rates...")

  #### Initial model: Poisson-based Maximum Likelihood Estimation ####
  analytical_rate_tbl <-
    tibble::tibble(gene_id = rc1$gene_id,
          beta_org =  (rc1$summarized_gb_counts / rc1$gb_length) / (rc1$summarized_pause_counts / rc1$pause_length))

  #### Adapted model: allow uncertainty in the pause site and steric hindrance ####
  # prepare data for running EM
  em_rate <-
    DataFrame(gene_id = rc1$gene_id,
              s = rc1$summarized_gb_counts,
              N = rc1$gb_length)


  # use read count within gene body to pre-estimate chi hat
  em_rate$chi = em_rate$s / em_rate$N

  # get read counts on each position within pause peak (from kmin to kmax)
  Xk <-
    BRGenomics::getCountsByPositions(bw1_p3, pause_regions, melt = TRUE)
  Xk <- splitAsList(Xk$signal, Xk$region)
  names(Xk) <- pause_regions$gene_id

  em_rate$Xk <- Xk[em_rate$gene_id]

  # initialize beta using sum of read counts within pause peak
  em_rate$Xk_sum <- sapply(em_rate$Xk, sum)
  em_rate$beta_int <- em_rate$chi / em_rate$Xk_sum

  # initialize fk with some reasonable values based on heuristic
  fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
  fk_int <- fk_int / sum(fk_int)

  # estimate rates using EM
  em_ls <- list()

  if(steric_hindrance) {
    em_rate$omega_zeta <- em_rate$chi * omega_scale
    em_rate$omega <- em_rate$omega_zeta / zeta

    # compute a scaling factor lambda for the purpose of using the same parameters as simulations
    lambda <- zeta ^ 2 / omega_scale
  }

  for (i in 1:NROW(em_rate)) {
    rc <- em_rate[i, ]

    if(!steric_hindrance) {
      em_ls[[i]] <- pause_escape_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                            fk_int = fk_int, beta_int = rc$beta_int[[1]],
                            chi_hat = rc$chi, max_itr = 500, tor = 1e-4)
    } else {
      em_ls[[i]] <- steric_hindrance_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax, f1 = 0.517, f2 = 0.024,
                            fk_int = fk_int, beta_int = rc$beta_int[[1]], phi_int = 0.5,
                            chi_hat = rc$chi, lambda = lambda, zeta = zeta,
                            max_itr = 500, tor = 1e-4)
    }
  }

  names(em_ls) <- em_rate$gene_id

  # get rate estimates and posterior distribution of pause sites
  em_rate$beta_adp <- map_dbl(em_ls, "beta", .default = NA)
  em_rate$Yk <- map(em_ls, "Yk", .default = NA)
  em_rate$fk <- map(em_ls, "fk", .default = NA)
  em_rate$fk_mean <- map_dbl(em_ls, "fk_mean", .default = NA)
  em_rate$fk_var <- map_dbl(em_ls, "fk_var", .default = NA)
  # calculate Yk / Xk
  em_rate$t <- sapply(em_rate$Yk, sum)
  em_rate$proportion_Yk <- em_rate$t / sapply(em_rate$Xk, sum)
  em_rate$likelihood <- map_dbl(em_ls, ~ .x$likelihoods[[length(.x$likelihoods)]])

  em_rate <- em_rate %>% as_tibble()

  if (steric_hindrance) {
    em_rate$phi <- map_dbl(em_ls, "phi", .default = NA)
    em_rate <- em_rate %>%
      mutate(alpha_zeta = omega_zeta / (1 - phi))
  }

  em_rate <- em_rate %>% left_join(analytical_rate_tbl, by = "gene_id")

  if (!steric_hindrance) {
      em_rate <- em_rate %>%
        select(gene_id, chi, beta_org, beta_adp, fk_mean, fk_var)
  } else {
    em_rate <- em_rate %>%
      select(gene_id, chi, beta_org, beta_adp, fk_mean, fk_var, phi, omega_zeta) %>%
      mutate(beta_zeta = beta_adp * zeta,
            alpha_zeta = omega_zeta / (1 - phi))
  }

  # Return experiment transcription rates object
  return(methods::new(Class = "experiment_transcription_rates",
                      counts = as.data.frame(rc1),
                      bigwig_plus = bigwig_plus,
                      bigwig_minus = bigwig_minus,
                      pause_regions = pause_regions,
                      gene_body_regions = gene_body_regions,
                      gene_name_column = gene_name_column,
                      steric_hindrance = steric_hindrance,
                      omega_scale = omega_scale,
                      rates = em_rate
  ))
}

#' @inherit methods::show
methods::setMethod("show", signature = "experiment_transcription_rates", function(object) {
  num_genes <- length(unique(object@counts$gene_id))
  gene_string <- paste("Number of Genes:", num_genes)

  write("A experiment_transcription_rates object with:", file = stdout())
  write(gene_string, file = stdout())
  write(paste("Estimated Transcription Rates:", nrow(object@rates)), file = stdout())
  write(paste("Steric Hindrance:", object@steric_hindrance), file = stdout())
  write(paste("Omega Scale:", object@omega_scale), file = stdout())
  write(paste("Summarized Read Counts:", nrow(object@counts)), file = stdout())
})

#' Get transcription rates
#'
#' @param object An experiment_transcription_rates object
#' @return A tibble containing the transcription rates
#' @export
setGeneric("getRates", function(object) standardGeneric("getRates"))
setMethod("getRates", "experiment_transcription_rates", function(object) {
  object@rates
})

#' Get read counts
#'
#' @param object An experiment_transcription_rates object
#' @return A data.frame containing the read counts
#' @export
setGeneric("getCounts", function(object) standardGeneric("getCounts"))
setMethod("getCounts", "experiment_transcription_rates", function(object) {
  object@counts
})

#' Get genomic regions
#'
#' @param object An experiment_transcription_rates object
#' @param type Either "pause" or "gene_body" to specify which regions to return
#' @return A GRanges object containing the specified regions
#' @export
setGeneric("getRegions", function(object, type) standardGeneric("getRegions"))
setMethod("getRegions", "experiment_transcription_rates", function(object, type) {
  if (!type %in% c("pause", "gene_body")) {
    stop("type must be either 'pause' or 'gene_body'")
  }
  if (type == "pause") {
    object@pause_regions
  } else {
    object@gene_body_regions
  }
})

#' Export rates to CSV
#'
#' @param object An experiment_transcription_rates object
#' @param file Path to output CSV file
#' @export
setGeneric("exportRatesToCSV", function(object, file) standardGeneric("exportRatesToCSV"))
setMethod("exportRatesToCSV", "experiment_transcription_rates", function(object, file) {
  write.csv(object@rates, file = file, row.names = FALSE)
})

#' Plot transcription rates
#'
#' @param object An experiment_transcription_rates object
#' @param type Type of plot to create ("scatter", "histogram", or "density")
#' @param rate_type Which rate to plot ("beta_org", "beta_adp", "chi", etc.)
#' @param file Optional path to save the plot. If provided, the plot will be saved to this location.
#' @param width Width of the saved plot in inches. Default is 8.
#' @param height Height of the saved plot in inches. Default is 6.
#' @param dpi Resolution of the saved plot. Default is 300.
#' @param ... Additional arguments passed to the plotting function
#' @return A ggplot object
#' @export
setGeneric("plotRates", function(object, type = "scatter", rate_type = "beta_adp", file = NULL, width = 8, height = 6, dpi = 300, ...) standardGeneric("plotRates"))
setMethod("plotRates", "experiment_transcription_rates", function(object, type = "scatter", rate_type = "beta_adp", file = NULL, width = 8, height = 6, dpi = 300, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for plotting")
  }
  
  if (!rate_type %in% colnames(object@rates)) {
    stop(paste("rate_type", rate_type, "not found in rates data"))
  }
  
  data <- object@rates
  
  p <- switch(type,
         scatter = {
           ggplot2::ggplot(data, ggplot2::aes(x = .data$beta_org, y = .data[[rate_type]])) +
             ggplot2::geom_point(color = "#1E88E5", alpha = 0.7, size = 2) +  # Blue points
             ggplot2::labs(x = "Original Beta", y = rate_type) +
             ggplot2::theme_bw() +  # White background with grid
             ggplot2::theme(
               panel.grid.major = ggplot2::element_line(color = "gray90"),
               panel.grid.minor = ggplot2::element_line(color = "gray95"),
               axis.text = ggplot2::element_text(color = "black", size = 12),
               axis.title = ggplot2::element_text(color = "black", size = 14)
             )
         },
         histogram = {
           ggplot2::ggplot(data, ggplot2::aes(x = .data[[rate_type]])) +
             ggplot2::geom_histogram(bins = 30, fill = "#1E88E5", color = "white", alpha = 0.7) +
             ggplot2::labs(x = rate_type, y = "Count") +
             ggplot2::theme_bw() +
             ggplot2::theme(
               panel.grid.major = ggplot2::element_line(color = "gray90"),
               panel.grid.minor = ggplot2::element_line(color = "gray95"),
               axis.text = ggplot2::element_text(color = "black", size = 12),
               axis.title = ggplot2::element_text(color = "black", size = 14)
             )
         },
         density = {
           ggplot2::ggplot(data, ggplot2::aes(x = .data[[rate_type]])) +
             ggplot2::geom_density(fill = "#1E88E5", color = "#0D47A1", alpha = 0.7) +
             ggplot2::labs(x = rate_type, y = "Density") +
             ggplot2::theme_bw() +
             ggplot2::theme(
               panel.grid.major = ggplot2::element_line(color = "gray90"),
               panel.grid.minor = ggplot2::element_line(color = "gray95"),
               axis.text = ggplot2::element_text(color = "black", size = 12),
               axis.title = ggplot2::element_text(color = "black", size = 14)
             )
         },
         stop("Invalid plot type. Choose from 'scatter', 'histogram', or 'density'")
  )
  
  # Save the plot if file is provided
  if (!is.null(file)) {
    ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
})
