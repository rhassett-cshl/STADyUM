#' experiment_transcription_rates_valid
#'
#' Checks that a \code{experiment_transcription_rates} object is valid
#' @param object a \code{experiment_transcription_rates} object
#'
#' @return TRUE if valid, else errors
experiment_transcription_rates_valid <- function(object) {
  errors <- c()

  if (length(errors) == 0) TRUE else errors
}

#' Class experiment_transcription_rates
#'
#' Class \code{experiment_transcription_rates} has
#'
#' @slot transcripts a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the transcript coordinates
#'
#' @name experiment_transcription_rates-class
#' @rdname experiment_transcription_rates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @exportClass experiment_transcription_rates
methods::setClass("experiment_transcription_rates",
                  slots = c(counts = "list",
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
#' @param scale a csv file providing scaling factors for omega with columns sample_id, omega_scale_h,
#' omega_scale_l. Defaults to NULL.
#'
#' @return an \code{\link{experiment_transcription_rates-class}} object
#'
#' @export
estimate_experiment_transcription_rates <- function(bigwig_plus, bigwig_minus,
                                                pause_regions, gene_body_regions,
                                                gene_name_column, steric_hindrance = FALSE, 
                                                omega_scale = NULL) {
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
  if (duplicated_gene_body_regions_gene_names) {
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

  bwp1_p3 <- import.bw(bigwig_plus)
  bwm1_p3 <- import.bw(bigwig_minus)

  # TODO: bigwig validity checks

  bwp1_p3 <- process_bw(bw = bwp1_p3, strand = "+")
  bwm1_p3 <- process_bw(bw = bwm1_p3, strand = "-")
  bw1_p3 <- c(bwp1_p3, bwm1_p3)
  rm(bwp1_p3, bwm1_p3)

  # make sure pause region is the same as kmax used in EM
  pause_regions <- promoters(pause_regions, upstream = 0, downstream = kmax)

  # summarize read counts
  rc1_pause <- summarise_bw(bw1_p3, pause_regions, "sp1")
  rc1_pause$pause_length <- kmax

  rc1_gb <- summarise_bw(bw1_p3, gene_body_regions, "sb1")
  rc1_gb$gb_length <-
  width(gene_body_regions)[match(rc1_gb$gene_id, gene_body_regions$gene_id)]

  # prepare read count table
  rc1 <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
              list(rc1_pause, rc1_gb))

  # clean up some genes with missing values in tss length or gene body length
  rc1 <- rc1[!(is.na(rc1$pause_length) | is.na(rc1$gb_length)), ]
  rc1 <- rc1[(rc1$sp1 > rc_cutoff) & (rc1$sb1 > rc_cutoff), ]

  message("estimating rates...")

  #### Initial model: Poisson-based Maximum Likelihood Estimation ####
  analytical_rate_tbl <-
    tibble(gene_id = rc1$gene_id,
          beta_org =  (rc1$sb1 / rc1$gb_length) / (rc1$sp1 / rc1$pause_length))

  #### Adapted model: allow uncertainty in the pause site and steric hindrance ####
  # prepare data for running EM
  em_rate <-
    DataFrame(gene_id = rc1$gene_id,
              s = rc1$sb1,
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

  em_rate <- em_rate %>% as_tibble()

  if (steric_hindrance) {
    em_rate$phi <- map_dbl(em_ls, "phi", .default = NA)
    em_rate <- em_rate %>%
      mutate(alpha_zeta = omega_zeta / (1 - phi))
  }

  em_rate <- em_rate %>% left_join(analytical_rate_tbl, by = "gene_id")

  # Return experiment transcription rates object
  return(methods::new(Class = "experiment_transcription_rates",
                      counts = rc1,
                      pause_regions = pause_regions,
                      gene_body_regions = gene_body_regions,
                      gene_name_column = gene_name_column,
                      steric_hindrance = steric_hindrance,
                      omega_scale = omega_scale,
                      rates = em_rate
  ))
}

#' @inherit methods::show
methods::setMethod("show", signature = "transcript_quantifier", function(object) {
  num_transcripts <- length(object@transcripts)
  num_models <- sum(unlist(lapply(object@models, ncol)))
  num_loci <- length(object@models)
  bin_size <- object@bin_size
  bwp <- object@count_metadata$bigwig_plus
  bwm <- object@count_metadata$bigwig_minus

  if (!is.na(object@column_identifiers["gene_id"])) {
    num_genes <- length(unique(
      S4Vectors::mcols(object@transcripts)[[object@column_identifiers["gene_id"]]]))
    gene_string <- paste("Number of Genes:", num_genes)
  } else {
    gene_string <- "No gene id present"
  }

  write("A transcript_quantifier object with:", file = stdout())
  write(paste(num_transcripts, "transcripts converted to", num_models, "models",
              "grouped into", num_loci, "loci"), file = stdout())
  write(gene_string, file = stdout())
  write(paste("bin size:", bin_size), file = stdout())
  write(paste("Bigwig data (plus):", bwp), file = stdout())
  write(paste("Bigwig data (minus):", bwm), file = stdout())
})
