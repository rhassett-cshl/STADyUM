#' simulation_transcription_rates_valid
#'
#' Checks that a \code{simulation_transcription_rates} object is valid
#' @param object a \code{simulation_transcription_rates} object
#'
#' @return TRUE if valid, else errors
simulation_transcription_rates_valid <- function(object) {
  errors <- c()

  if (length(errors) == 0) TRUE else errors
}

#' Class simulation_transcription_rates
#'
#' Class \code{simulation_transcription_rates} 
#'
#' @name simulation_transcription_rates-class
#' @rdname simulation_transcription_rates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @importClassesFrom data.table data.table
#' @exportClass simulation_transcription_rates
methods::setClass("simulation_transcription_rates",
                  slots = c(simpol="simulate_polymerase",
                            steric_hindrance="logical"),
                  validity = simulation_transcription_rates_valid
)

#' estimate_simulation_transcription_rates
#'
#' Estimates the transcription rates from simulated data and contructs an object that holds
#' these rates
#' @param simpol a \code{\link{simulate_polymerase-class}} object
#' @param steric_hindrance a logical value to determine whether to infer landing-pad occupancy or not.
#' Defaults to FALSE.
#'
#' @return an \code{\link{_transcription_rates-class}} object
#'
#' @export
estimate_simulation_transcription_rates <- function() {

  sample_cell <- 5000 # number of cells sampled each time
  sample_n <- 50 # number of times to sample
  matched_len <- 2e4 # gene length to be matched

  kmin <- 1
  kmax <- 200
  matched_gb_len <- matched_len - kmax
  count_rnap <- FALSE

  total_spacing <- simpol@s + simpol@h

  k <- simpol@k

  start_point <- 0.99 * 1e6 # set a start coordinate for the simulated gene
  # set lambda according to the read coverage of PRO-seq in the control samples
  # from Dukler et al. 2017
  lambda <- 102.1 

  # last step of position matrix
  rnap_pos <- simpol@position_matrix
  total_cell <- NCOL(rnap_pos)
  gene_len <- NROW(rnap_pos) - 1

  # count number of RNAP before the pause site
  #if (count_rnap) {
    # get probability vector
  #  prob <- readRDS(str_replace(rds_in, ".RDS", "_prob_init.RDS"))
  #  alpha <- as.double(gsub(".*a([0-9].*)b.*", "\\1",sel_sample))
  #  beta <- as.double(gsub(".*b([0-9].*)g.*", "\\1",sel_sample))
  #  zeta <- as.double(gsub(".*z([0-9].*)zsd.*", "\\1",sel_sample))

    # calculate time slice first then get the corresponding probability for beta
  #  beta_prob <- prob[1, 1] / alpha * beta
    # get pause position for every cell
  #  idx <- which(prob == beta_prob, arr.ind = TRUE)
  #  idx <- idx[idx[, 1] != 1, ]
  #}
  
  #### initiation and pause release rate estimates ####
  # generate regions for read counting
  gn_rng <-
  GRanges(seqnames = rep("chr1", 3),
          IRanges(start = c(1, kmax + 1, 1),
                  end = c(kmax, gene_len, spacing)))
  
  gn_rng <- shift(gn_rng, shift = start_point)

  region_names <- c("tss", "gb", "landing")
  names(gn_rng) <- region_names

  len <- as.list(width(gn_rng))
  names(len) <- region_names

  # set seeds for random sampling
  seeds <- seq(from = 2013, by = 1, length.out = sample_n)
  # a list to Granges for rnap positions
  rnap_grng <- list()
  # a list recording number of RNAPs at or before the pause site
  if (count_rnap) rnap_n_ls <- list()

  for (i in 1:sample_n) {
    set.seed(seeds[[i]])
    sel_cells <- sample(1:total_cell, size = sample_cell, replace = TRUE)
    res_pos <- rnap_pos[, sel_cells]
    # get rid of position 1, which is always 1
    res_pos <- res_pos[-1, ]
    if (count_rnap) {
      # get pause sites
      pause_site <- idx[sel_cells, 1] - 1
      # generate data mask
      # inspired by https://stackoverflow.com/questions/47732085/sum-of-some-positions-in-a-row-r
      res_shape <- dim(res_pos)
      after_pause_len <- res_shape[1] - pause_site
      mask_mx <- map2(pause_site, after_pause_len,
          function(x, y) c(rep(TRUE, x), rep(FALSE, y)))
      mask_mx <- Matrix::Matrix(unlist(mask_mx), nrow = res_shape[1], ncol = res_shape[2])
      # calculate rnap number before pause site for every cell
      rnap_n_ls[[i]] <- colSums(res_pos * mask_mx)
    }
    # combine rnap positions across all cells
    res_all <- rowSums(res_pos)
    # generate bigwig for positive strand
    rnap_grng[[i]] <-
      GRanges(seqnames = "chr1",
              IRanges(start = (1 + start_point) : (gene_len + start_point),
                      width = 1),
              score = res_all,
              strand = "+",
              seqlengths = c("chr1" = gene_len * 10)  + start_point)

    rm(res_pos, res_all)
  }

  summarise_bw <-
  function(bw, grng) {
    rc <- grng %>%
      plyranges::group_by_overlaps(bw) %>%
      plyranges::group_by(query) %>%
      plyranges::summarise(score = sum(score))
    if (!1 %in% rc$query) {
      rc <- rbind(DataFrame(list(query = 1, score = 0)), rc)
    }
    rc <- as.list(rc$score)
    names(rc) <- region_names
    return(rc)
  }

  bw_dfs <- tibble(trial = 1:sample_n)
  bw_dfs$rc_region <- map(rnap_grng, ~ summarise_bw(.x, gn_rng))

  bw_dfs$rc_tss <- map_dbl(bw_dfs$rc_region, "tss")
  bw_dfs$rc_gb <-map_dbl(bw_dfs$rc_region, "gb")
  bw_dfs$rc_landing <-map_dbl(bw_dfs$rc_region, "landing")

  #### empirical way to calculate steric hindrance at pause site ####
  bw_dfs <- bw_dfs %>%
    mutate(# number of RNAPs per cell per gene
          R = (rc_tss + rc_gb) / sample_cell,
          # number of RNAPs in the pause peak per cell per gene
          R_pause = rc_tss / sample_cell,
          # proportion of landing pad being occupied by RNAP, i.e., empirical phi
          rnap_prop = rc_landing / sample_cell
          )
  
  # whether to match the simulated number of RNAPs to read coverage in experimental data or not
  # here match RNAP number within kmin to kmax, RNAP in gene body will be taken care
  # of afterwards
  if (!is.null(lambda)) {
    rnap_grng <- map(rnap_grng, function(grng) {
      grng$score[kmin:kmax] <-
        rpois(length(kmin:kmax), grng$score[kmin:kmax] / sample_cell * lambda)
      # first 20bp get removed because they are usually not seen in sequencing
      grng$score[1:20] <- 0
      return(grng)
    })
    bw_dfs$rc_region <- map(rnap_grng, ~ summarise_bw(.x, gn_rng))
    bw_dfs$rc_tss <- map_dbl(bw_dfs$rc_region, "tss")
  }

  # match RNAP number within gene bodies to desired read coverage
  if (!is.null(lambda)) {
    pois_mean <- (lambda * bw_dfs$rc_gb / sample_cell) * (matched_gb_len / len$gb)
    bw_dfs$rc_gb <- rpois(length(pois_mean), pois_mean)
    # assign matched gene body length as gene body length
    len$gb <- matched_gb_len
  }

  # get read counts on each position within pause peak (from kmin to kmax)
  bw_dfs$Xk <- map(rnap_grng,
                  ~ .x[(start(.x) >= 990000 + kmin) & (start(.x) <= 990000 + kmax), ]$score)

  ## Initial model: Poisson based MLEs ##
  # use read count within gene body to pre-estimate chi hat
  bw_dfs$chi <- bw_dfs$rc_gb / len$gb

  # TODO??
  # take care of single pause site or variable pause sites (pause peak)
  if (simpol@k_sd == 0) {
    bw_dfs$beta_org <- bw_dfs$chi / map_dbl(bw_dfs$Xk, k)
  } else {
    bw_dfs$beta_org <- bw_dfs$chi / (bw_dfs$rc_tss / len$tss)
  }
    
  bw_dfs$beta_max_rc <- bw_dfs$chi / map_dbl(bw_dfs$Xk, max)

  ## Adapted model: allows uncertainty in the pause site and steric hindrance ##
  # initialize beta using sum of read counts within pause peak
  bw_dfs$Xk_sum <- sapply(bw_dfs$Xk, sum)
  bw_dfs$beta_int <- bw_dfs$chi / bw_dfs$Xk_sum

  # initialize fk with some reasonable values based on heuristic
  fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
  fk_int <- fk_int / sum(fk_int)

  # estimate rates using EM
  em_ls <- list()
  # wrap the EM function in case there is an error
  main_EM <- possibly(main_EM, otherwise = NA)

  if (steric_hindrance) {
    f <- calculate_f(s = spacing, k = k)
    phi_int <- 0.5
    zeta <- 2000 # elongation rate
    # lambda used for scaling in EM, different from the one used to match coverage
    lambda1 <- 0.0505 * zeta ^ 2
  }


  for (i in 1:NROW(bw_dfs)) {
    rc <- bw_dfs[i, ]

    if (!steric_hindrance) {

      em_ls[[i]] <- main_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                            fk_int = fk_int, beta_int = rc$beta_int[[1]], chi_hat = rc$chi,
                            max_itr = 500, tor = 1e-4)

    } else {

      em_ls[[i]] <- main_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax, f1 = f[["f1"]], f2 = f[["f2"]],
                            fk_int = fk_int, beta_int = rc$beta_int[[1]], chi_hat = rc$chi,
                            phi_int = phi_int, lambda = lambda1, zeta = zeta,
                            max_itr = 500, tor = 1e-4)
    }
  }
  

  # get rate estimates and posterior distribution of pause sites
  bw_dfs$beta_adp <- map_dbl(em_ls, "beta", .default = NA)
  bw_dfs$Yk <- map(em_ls, "Yk", .default = NA)
  bw_dfs$fk <- map(em_ls, "fk", .default = NA)
  bw_dfs$fk_mean <- map_dbl(em_ls, "fk_mean", .default = NA)
  bw_dfs$fk_var <- map_dbl(em_ls, "fk_var", .default = NA)
  # calculate Yk / Xk
  # bw_dfs$proportion_Yk <- sapply(bw_dfs$Yk, sum) / sapply(bw_dfs$Xk, sum)
  bw_dfs$flag <- map_chr(em_ls, "flag", .default = NA)

  if (steric_hindrance) bw_dfs$phi <- map_dbl(em_ls, "phi", .default = NA)

  # add number of RNAP before pause site to output if it exists
  if (count_rnap) bw_dfs$rnap_n <- rnap_n_ls

  # Return simulation transcription rates object
  return(methods::new(Class = "simulation_transcription_rates",
                      simpol = simpol,
                      steric_hindrance = steric_hindrance,
                      trial = bw_dfs$trial,
                      chi = bw_dfs$chi,
                      beta_org = bw_dfs$beta_org,
                      beta_adp = bw_dfs$beta_adp,
                      phi = bw_dfs$phi,
                      fk = bw_dfs$fk,
                      fk_mean = bw_dfs$fk_mean,
                      fk_var = bw_dfs$fk_var,
                      flag = bw_dfs$flag,
                      rnap_n = bw_dfs$rnap_n,
                      ))
}
