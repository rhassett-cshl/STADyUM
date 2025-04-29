#' simulation_transcription_rates_valid
#'
#' Checks that a \code{simulation_transcription_rates} object is valid
#' @param object a \code{simulation_transcription_rates} object
#'
#' @return TRUE if valid, else errors
#' @keywords internal
simulation_transcription_rates_valid <- function(object) {
    errors <- character()

    # Check if required slots are present and have correct types
    if (!is(simpol(object), "simulate_polymerase")) {
        errors <- c(errors, "simpol must be a simulate_polymerase object")
    }
    if (!is.logical(steric_hindrance(object))) {
        errors <- c(errors, "steric_hindrance must be logical")
    }
    if (!is.numeric(trial(object))) {
        errors <- c(errors, "trial must be numeric")
    }
    if (!is.numeric(chi(object))) {
        errors <- c(errors, "chi must be numeric")
    }
    if (!is.numeric(beta_org(object))) {
        errors <- c(errors, "beta_org must be numeric")
    }
    if (!is.numeric(beta_adp(object))) {
        errors <- c(errors, "beta_adp must be numeric")
    }
    if (!is.numeric(phi(object))) {
        errors <- c(errors, "phi must be numeric")
    }
    if (!is.list(fk(object))) {
        errors <- c(errors, "fk must be a list")
    }
    if (!is.numeric(fk_mean(object))) {
        errors <- c(errors, "fk_mean must be numeric")
    }
    if (!is.numeric(fk_var(object))) {
        errors <- c(errors, "fk_var must be numeric")
    }
    if (!is.character(flag(object))) {
        errors <- c(errors, "flag must be character")
    }
    if (!is.list(rnap_n(object))) {
        errors <- c(errors, "rnap_n must be a list")
    }

    if (length(errors) == 0) TRUE else errors
}

#' Class simulation_transcription_rates
#'
#' Class \code{simulation_transcription_rates} estimates the transcription
#' rates, such as initiation, pause-release rates and landing pad occupancy,
#' from simulated data, such as nascent RNA sequencing read counts and genomic
#' coordinates simulated from polymerase movement over a specified time period
#' calculated by simulator, and contructs an object that holds these rates.
#'
#' @name simulation_transcription_rates-class
#' @rdname simulation_transcription_rates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @importClassesFrom data.table data.table
#' @exportClass simulation_transcription_rates
methods::setClass("simulation_transcription_rates",
    slots = c(
        simpol = "simulate_polymerase",
        steric_hindrance = "logical",
        trial = "numeric",
        chi = "numeric",
        beta_org = "numeric",
        beta_adp = "numeric",
        phi = "numeric",
        fk = "list",
        fk_mean = "numeric",
        fk_var = "numeric",
        flag = "character",
        rnap_n = "list"
    ),
    prototype = list(
        simpol = NULL,
        steric_hindrance = FALSE,
        trial = numeric(0),
        chi = numeric(0),
        beta_org = numeric(0),
        beta_adp = numeric(0),
        phi = 0,
        fk = list(),
        fk_mean = numeric(0),
        fk_var = numeric(0),
        flag = character(0),
        rnap_n = list()
    ),
    validity = simulation_transcription_rates_valid
)

#' estimate_simulation_transcription_rates
#'
#' Estimates the transcription rates, such as initiation, pause-release rates
#' and landing pad occupancy, from simulated data, such as nascent RNA
#' sequencing read counts and genomic coordinates, and contructs an object
#' that holds these rates
#'
#' @param simpol a \code{\link{simulate_polymerase-class}} object
#' @param steric_hindrance a logical value to determine whether to infer
#' landing-pad occupancy or not.
#' Defaults to FALSE.
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @return an \code{\link{_transcription_rates-class}} object
#'
#' @export
estimate_simulation_transcription_rates <-
    function(simpol, steric_hindrance = FALSE) {
        sample_cell <- 5000 # number of cells sampled each time
        sample_n <- 50 # number of times to sample
        matched_len <- 2e4 # gene length to be matched

        kmin <- 1
        kmax <- 200
        matched_gb_len <- matched_len - kmax
        count_rnap <- FALSE

        spacing <- spacing(simpol) + height(simpol)
        k <- k(simpol)

        # set a start coordinate for the simulated gene
        start_point <- 0.99 * 1e6
        # set lambda according to the read coverage of PRO-seq in the control
        # samples from Dukler et al. 2017
        lambda <- 102.1

        # last step of position matrix
        rnap_pos <- position_matrix(simpol)
        total_cell <- NCOL(rnap_pos)
        gene_len <- NROW(rnap_pos) - 1

        # count number of RNAP before the pause site
        # if (count_rnap) {
        # get probability vector
        # prob <- readRDS(str_replace(rds_in, ".RDS", "_prob_init.RDS"))
        # alpha <- as.double(gsub(".*a([0-9].*)b.*", "\\1", sel_sample))
        # beta <- as.double(gsub(".*b([0-9].*)g.*", "\\1", sel_sample))
        # zeta <- as.double(gsub(".*z([0-9].*)zsd.*", "\\1", sel_sample))

        # calculate time slice first then get corresponding probability for beta
        # beta_prob <- prob[1, 1] / alpha * beta
        # get pause position for every cell
        # idx <- which(prob == beta_prob, arr.ind = TRUE)
        # idx <- idx[idx[, 1] != 1, ]
        # }

        #### initiation and pause release rate estimates ####
        # generate regions for read counting
        gn_rng <- GRanges(
            seqnames = rep("chr1", 3),
            IRanges(
                start = c(1, kmax + 1, 1),
                end = c(kmax, gene_len, spacing)
            )
        )

        gn_rng <- IRanges::shift(gn_rng, shift = start_point)

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

        for (i in seq_len(sample_n)) {
            sel_cells <- sample(seq_len(total_cell),
                size = sample_cell,
                replace = TRUE
            )
            res_pos <- rnap_pos[, sel_cells]
            # get rid of position 1, which is always 1
            res_pos <- res_pos[-1, ]
            if (count_rnap) {
                # get pause sites
                pause_site <- idx[sel_cells, 1] - 1
                # generate data mask
                res_shape <- dim(res_pos)
                after_pause_len <- res_shape[1] - pause_site
                mask_mx <- map2(
                    pause_site,
                    after_pause_len,
                    function(x, y) c(rep(TRUE, x), rep(FALSE, y))
                )
                mask_mx <- Matrix::Matrix(
                    unlist(mask_mx),
                    nrow = res_shape[1],
                    ncol = res_shape[2]
                )
                # calculate rnap number before pause site for every cell
                rnap_n_ls[[i]] <- colSums(res_pos * mask_mx)
            }
            # combine rnap positions across all cells
            res_all <- rowSums(res_pos)
            # generate bigwig for positive strand
            rnap_grng[[i]] <- GRanges(
                seqnames = "chr1",
                IRanges(
                    start = (1 + start_point):(gene_len + start_point),
                    width = 1
                ),
                score = res_all,
                strand = "+",
                seqlengths = c("chr1" = gene_len * 10) + start_point
            )

            rm(res_pos, res_all)
        }

        #' @keywords internal
        summarise_bw <- function(bw, grng) {
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

        bw_dfs <- tibble(trial = seq_len(sample_n))
        bw_dfs$rc_region <- map(rnap_grng, ~ summarise_bw(.x, gn_rng))

        bw_dfs$rc_tss <- map_dbl(bw_dfs$rc_region, "tss")
        bw_dfs$rc_gb <- map_dbl(bw_dfs$rc_region, "gb")
        bw_dfs$rc_landing <- map_dbl(bw_dfs$rc_region, "landing")

        #### empirical way to calculate steric hindrance at pause site ####
        bw_dfs <- bw_dfs %>%
            mutate(
                R = (rc_tss + rc_gb) / sample_cell,
                R_pause = rc_tss / sample_cell,
                rnap_prop = rc_landing / sample_cell
            )

        # whether to match the simulated number of RNAPs to read coverage in
        # experimental data or not
        # here match RNAP number within kmin to kmax, RNAP in gene body will be
        # taken care of afterwards
        if (!is.null(lambda)) {
            rnap_grng <- map(rnap_grng, function(grng) {
                grng$score[kmin:kmax] <- rpois(
                    length(kmin:kmax),
                    grng$score[kmin:kmax] / sample_cell * lambda
                )
                # first 20bp get removed because they are usually not seen in
                # sequencing
                grng$score[seq_len(20)] <- 0
                return(grng)
            })
            bw_dfs$rc_region <- map(rnap_grng, ~ summarise_bw(.x, gn_rng))
            bw_dfs$rc_tss <- map_dbl(bw_dfs$rc_region, "tss")
        }

        # match RNAP number within gene bodies to desired read coverage
        if (!is.null(lambda)) {
            pois_mean <- (lambda * bw_dfs$rc_gb / sample_cell) *
                (matched_gb_len / len$gb)
            bw_dfs$rc_gb <- rpois(length(pois_mean), pois_mean)
            # assign matched gene body length as gene body length
            len$gb <- matched_gb_len
        }

        # get read counts on each position within pause peak (from kmin to kmax)
        bw_dfs$Xk <- map(
            rnap_grng,
            ~ .x[(start(.x) >= 990000 + kmin) &
                (start(.x) <= 990000 + kmax), ]$score
        )

        ## Initial model: Poisson based MLEs ##
        # use read count within gene body to pre-estimate chi hat
        bw_dfs$chi <- bw_dfs$rc_gb / len$gb

        # TODO??
        # take care of single pause site or variable pause sites (pause peak)
        if (ksd(simpol) == 0) {
            bw_dfs$beta_org <- bw_dfs$chi / map_dbl(bw_dfs$Xk, k)
        } else {
            bw_dfs$beta_org <- bw_dfs$chi / (bw_dfs$rc_tss / len$tss)
        }

        bw_dfs$beta_max_rc <- bw_dfs$chi / map_dbl(bw_dfs$Xk, max)

        ## Adapted model: allows uncertainty in the pause site and steric
        #  hindrance
        # initialize beta using sum of read counts within pause peak
        bw_dfs$Xk_sum <- vapply(bw_dfs$Xk, sum, numeric(1))

        # Handle cases where Xk_sum is 0 or chi is 0/Inf
        valid_indices <- which(bw_dfs$Xk_sum > 0 &
            bw_dfs$chi > 0 &
            !is.infinite(bw_dfs$chi))

        if (length(valid_indices) == 0) {
            stop("No valid data points found - check if there are any RNA
    polymerases in the pause region or gene body")
        }

        bw_dfs$beta_int <- rep(NA, nrow(bw_dfs))
        bw_dfs$beta_int[valid_indices] <- bw_dfs$chi[valid_indices] /
            bw_dfs$Xk_sum[valid_indices]

        # initialize fk with some reasonable values based on heuristic
        fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
        fk_int <- fk_int / sum(fk_int)

        # estimate rates using EM
        em_ls <- list()
        # wrap the EM function in case there is an error
        # main_EM <- possibly(main_EM, otherwise = NA)

        if (steric_hindrance) {
            f <- calculate_f(s = spacing, k = k)
            phi_int <- 0.5
            zeta <- 2000 # elongation rate
            # lambda used for scaling in EM, different from one used to match
            # coverage
            lambda1 <- 0.0505 * zeta^2
        }

        message("Starting EM algorithm...")

        for (i in seq_len(NROW(bw_dfs))) {
            rc <- bw_dfs[i, ]

            if (!steric_hindrance) {
                # rc$beta_int[[1]] is NaN

                em_ls[[i]] <- pause_escape_EM(
                    Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                    fk_int = fk_int, beta_int = rc$beta_int[[1]],
                    chi_hat = rc$chi, max_itr = 500, tor = 1e-4
                )
            } else {
                em_ls[[i]] <- steric_hindrance_EM(
                    Xk = rc$Xk[[1]], kmin = kmin,
                    kmax = kmax, f1 = f[["f1"]], f2 = f[["f2"]],
                    fk_int = fk_int, beta_int = rc$beta_int[[1]],
                    chi_hat = rc$chi, phi_int = phi_int,
                    lambda = lambda1, zeta = zeta, max_itr = 500,
                    tor = 1e-4
                )
            }
        }


        # get rate estimates and posterior distribution of pause sites
        bw_dfs$beta_adp <- map_dbl(em_ls, "beta", .default = NA)
        bw_dfs$Yk <- map(em_ls, "Yk", .default = NA)
        bw_dfs$fk <- map(em_ls, "fk", .default = NA)
        bw_dfs$fk_mean <- map_dbl(em_ls, "fk_mean", .default = NA)
        bw_dfs$fk_var <- map_dbl(em_ls, "fk_var", .default = NA)
        # calculate Yk / Xk
        # bw_dfs$proportion_Yk <- vapply(bw_dfs$Yk, sum, numeric(1)) / vapply
        # (bw_dfs$Xk, sum, numeric(1))
        bw_dfs$flag <- map_chr(em_ls, "flag", .default = NA)

        if (steric_hindrance) bw_dfs$phi <- map_dbl(em_ls, "phi", .default = NA)

        # add number of RNAP before pause site to output if it exists
        if (count_rnap) bw_dfs$rnap_n <- rnap_n_ls

        # calculate initiation rate
        bw_dfs$init_rate <- bw_dfs$R / (bw_dfs$R_pause + bw_dfs$R)

        # calculate pause release rate
        bw_dfs$pause_release_rate <- bw_dfs$R_pause / (bw_dfs$R_pause +
            bw_dfs$R)

        # calculate landing pad occupancy
        if (steric_hindrance) {
            bw_dfs$landing_pad_occupancy <- bw_dfs$rnap_prop
        }

        # calculate mean and variance of initiation rate
        init_rate_mean <- mean(bw_dfs$init_rate)
        init_rate_var <- var(bw_dfs$init_rate)

        # calculate mean and variance of pause release rate
        pause_release_rate_mean <- mean(bw_dfs$pause_release_rate)
        pause_release_rate_var <- var(bw_dfs$pause_release_rate)

        # calculate mean and variance of landing pad occupancy
        if (steric_hindrance) {
            landing_pad_occupancy_mean <- mean(bw_dfs$landing_pad_occupancy)
            landing_pad_occupancy_var <- var(bw_dfs$landing_pad_occupancy)
        }

        # create a list to hold the results
        results <- list(
            init_rate = init_rate_mean,
            init_rate_var = init_rate_var,
            pause_release_rate = pause_release_rate_mean,
            pause_release_rate_var = pause_release_rate_var
        )

        if (steric_hindrance) {
            results$landing_pad_occupancy <- landing_pad_occupancy_mean
            results$landing_pad_occupancy_var <- landing_pad_occupancy_var
        }

        # create and return the simulation_transcription_rates object
        new("simulation_transcription_rates",
            simpol = simpol,
            steric_hindrance = steric_hindrance,
            trial = sample_n,
            chi = init_rate_mean,
            beta_org = pause_release_rate_mean,
            beta_adp = if (steric_hindrance) landing_pad_occupancy_mean else 0,
            phi = if (steric_hindrance) landing_pad_occupancy_mean else 0,
            fk = results,
            fk_mean = c(init_rate_mean, pause_release_rate_mean),
            fk_var = c(init_rate_var, pause_release_rate_var),
            flag = if (steric_hindrance) {
                "with_steric_hindrance"
            } else {
                "no_steric_hindrance"
            },
            rnap_n = if (count_rnap) rnap_n_ls else list()
        )
    }

# Accessor methods
#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_simpol", function(object) standardGeneric("get_simpol"))
setMethod(
    "get_simpol", "simulation_transcription_rates",
    function(object) slot(object, "simpol")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric(
    "get_steric_hindrance",
    function(object) standardGeneric("get_steric_hindrance")
)
setMethod(
    "get_steric_hindrance", "simulation_transcription_rates",
    function(object) slot(object, "steric_hindrance")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_trial", function(object) standardGeneric("get_trial"))
setMethod(
    "get_trial", "simulation_transcription_rates",
    function(object) slot(object, "trial")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_chi", function(object) standardGeneric("get_chi"))
setMethod(
    "get_chi", "simulation_transcription_rates",
    function(object) slot(object, "chi")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_beta_org", function(object) standardGeneric("get_beta_org"))
setMethod(
    "get_beta_org", "simulation_transcription_rates",
    function(object) slot(object, "beta_org")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_beta_adp", function(object) standardGeneric("get_beta_adp"))
setMethod(
    "get_beta_adp", "simulation_transcription_rates",
    function(object) slot(object, "beta_adp")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_phi", function(object) standardGeneric("get_phi"))
setMethod(
    "get_phi", "simulation_transcription_rates",
    function(object) slot(object, "phi")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_fk", function(object) standardGeneric("get_fk"))
setMethod(
    "get_fk", "simulation_transcription_rates",
    function(object) slot(object, "fk")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_fk_mean", function(object) standardGeneric("get_fk_mean"))
setMethod(
    "get_fk_mean", "simulation_transcription_rates",
    function(object) slot(object, "fk_mean")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_fk_var", function(object) standardGeneric("get_fk_var"))
setMethod(
    "get_fk_var", "simulation_transcription_rates",
    function(object) slot(object, "fk_var")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_flag", function(object) standardGeneric("get_flag"))
setMethod(
    "get_flag", "simulation_transcription_rates",
    function(object) slot(object, "flag")
)

#' @rdname simulation_transcription_rates-class
#' @export
setGeneric("get_rnap_n", function(object) standardGeneric("get_rnap_n"))
setMethod(
    "get_rnap_n", "simulation_transcription_rates",
    function(object) slot(object, "rnap_n")
)

# Plotting methods
#' Plot transcription rates
#' @param object A simulation_transcription_rates object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the transcription rates
#' @export
setGeneric("plot_transcription_rates", function(
    object, file = NULL,
    width = 8, height = 6) {
    standardGeneric("plot_transcription_rates")
})
setMethod(
    "plot_transcription_rates", "simulation_transcription_rates",
    function(object, file = NULL, width = 8, height = 6) {
        # Create data frame for plotting
        df <- data.frame(
            trial = get_trial(object),
            chi = get_chi(object),
            beta_org = get_beta_org(object),
            beta_adp = get_beta_adp(object),
            phi = get_phi(object),
            fk_mean = get_fk_mean(object),
            fk_var = get_fk_var(object)
        )

        # Create plot
        p <- ggplot(df, aes(x = trial)) +
            geom_line(aes(y = chi, color = "chi")) +
            geom_line(aes(y = beta_org, color = "beta_org")) +
            geom_line(aes(y = beta_adp, color = "beta_adp")) +
            geom_line(aes(y = phi, color = "phi")) +
            theme_minimal() +
            labs(
                title = "Transcription Rates Across Trials",
                x = "Trial",
                y = "Rate",
                color = "Rate Type"
            )

        if (!is.null(file)) {
            ggsave(file, p, width = width, height = height)
        }

        return(p)
    }
)

#' Plot pause site distribution
#' @param object A simulation_transcription_rates object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the pause site distribution
#' @export
setGeneric("plot_pause_site_distribution", function(
    object, file = NULL,
    width = 8, height = 6) {
    standardGeneric("plot_pause_site_distribution")
})
setMethod(
    "plot_pause_site_distribution", "simulation_transcription_rates",
    function(object, file = NULL, width = 8, height = 6) {
        # Create data frame for plotting
        df <- data.frame(
            trial = get_trial(object),
            fk_mean = get_fk_mean(object),
            fk_var = get_fk_var(object)
        )

        # Create plot
        p <- ggplot(df, aes(x = fk_mean, y = fk_var)) +
            geom_point() +
            theme_minimal() +
            labs(
                title = "Pause Site Distribution",
                x = "Mean Position",
                y = "Variance"
            )

        if (!is.null(file)) {
            ggsave(file, p, width = width, height = height)
        }

        return(p)
    }
)

#' Plot RNAP counts
#' @param object A simulation_transcription_rates object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the RNAP counts
#' @export
setGeneric("plot_rnap_counts", function(
    object, file = NULL, width = 8,
    height = 6) {
    standardGeneric("plot_rnap_counts")
})
setMethod(
    "plot_rnap_counts", "simulation_transcription_rates",
    function(object, file = NULL, width = 8, height = 6) {
        rnap_n <- get_rnap_n(object)
        if (length(rnap_n) == 0) {
            stop("No RNAP counts available")
        }

        # Create data frame for plotting
        df <- data.frame(
            trial = rep(get_trial(object), each = length(rnap_n[[1]])),
            rnap_count = unlist(rnap_n)
        )

        # Create plot
        p <- ggplot(df, aes(x = rnap_count)) +
            geom_histogram(bins = 30) +
            theme_minimal() +
            labs(
                title = "RNAP Count Distribution",
                x = "RNAP Count",
                y = "Frequency"
            )

        if (!is.null(file)) {
            ggsave(file, p, width = width, height = height)
        }

        return(p)
    }
)

#' @rdname simulation_transcription_rates-class
#' @export
setMethod("show", "simulation_transcription_rates", function(object) {
    # Create a data frame for display
    df <- data.frame(
        Parameter = c(
            "Number of Trials", "Initiation Rate",
            "Pause Release Rate", "Steric Hindrance"
        ),
        Value = c(
            get_trial(object),
            round(get_chi(object), 4),
            round(get_beta_org(object), 4),
            get_steric_hindrance(object)
        )
    )

    # Print the object information
    cat("A simulation_transcription_rates object:\n\n")
    print(df, row.names = FALSE)
})

#' @examples
#' # Create a simulation_transcription_rates object
#' simpol <- new("simulate_polymerase",
#'     k = 50, ksd = 10, k_min = 30, k_max = 70,
#'     gene_len = 1000, alpha = 0.1, beta = 0.2, zeta = 1000,
#'     zeta_sd = 100, zeta_min = 800, zeta_max = 1200,
#'     cell_num = 1000, pol_size = 35, add_space = 15,
#'     time = 10, steps_to_record = 100
#' )
#'
#' # Estimate transcription rates
#' rates <- estimate_simulation_transcription_rates(simpol)
#'
#' # Show the object
#' show(rates)
