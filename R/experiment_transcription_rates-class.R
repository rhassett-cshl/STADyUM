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
    errors <- c(
        validate_counts(object),
        validate_bigwig_files(object),
        validate_regions(object),
        validate_gene_name_column(object),
        validate_steric_hindrance(object),
        validate_rates(object)
    )

    if (length(errors) == 0) TRUE else errors
}

# Helper validation functions
validate_counts <- function(object) {
    errors <- character()
    if (!is.data.frame(counts(object))) {
        errors <- c(errors, "counts must be a data.frame")
    }

    required_cols <- c("gene_id", "tss", "gb", "landing")
    missing_cols <- setdiff(required_cols, colnames(counts(object)))
    if (length(missing_cols) > 0) {
        errors <- c(errors, paste(
            "Missing required columns in counts:",
            paste(missing_cols, collapse = ", ")
        ))
    }
    return(errors)
}

validate_bigwig_files <- function(object) {
    errors <- character()
    if (!file.exists(bigwig_plus(object))) {
        errors <- c(errors, "bigwig_plus file does not exist")
    }
    if (!file.exists(bigwig_minus(object))) {
        errors <- c(errors, "bigwig_minus file does not exist")
    }
    return(errors)
}

validate_regions <- function(object) {
    errors <- character()
    if (!inherits(pause_regions(object), "GRanges")) {
        errors <- c(errors, "pause_regions must be a GRanges object")
    }
    if (!inherits(gene_body_regions(object), "GRanges")) {
        errors <- c(errors, "gene_body_regions must be a GRanges object")
    }
    return(errors)
}

validate_gene_name_column <- function(object) {
    errors <- character()
    if (!is.character(gene_name_column(object)) ||
        length(gene_name_column(object)) != 1) {
        errors <- c(errors, "gene_name_column must be a single character
        string")
    }
    return(errors)
}

validate_steric_hindrance <- function(object) {
    errors <- character()
    if (!is.logical(steric_hindrance(object)) ||
        length(steric_hindrance(object)) != 1) {
        errors <- c(errors, "steric_hindrance must be a single logical value")
    }

    if (steric_hindrance(object)) {
        if (!is.numeric(omega_scale(object)) ||
            length(omega_scale(object)) != 1 ||
            omega_scale(object) <= 0) {
            errors <- c(errors, "omega_scale must be a single positive numeric
            value")
        }
    }
    return(errors)
}

validate_rates <- function(object) {
    errors <- character()
    if (!inherits(rates(object), "tbl_df")) {
        errors <- c(errors, "rates must be a tibble")
    }

    required_cols <- c("gene_id", "initiation_rate", "pause_release_rate")
    if (steric_hindrance(object)) {
        required_cols <- c(required_cols, "landing_pad_occupancy")
    }
    missing_cols <- setdiff(required_cols, colnames(rates(object)))
    if (length(missing_cols) > 0) {
        errors <- c(errors, paste(
            "Missing required columns in rates:",
            paste(missing_cols, collapse = ", ")
        ))
    }
    return(errors)
}

#' Class experiment_transcription_rates
#'
#' Class \code{experiment_transcription_rates} has read counts, pause and gene
#' body genomic region coordinates, steric hindrance and omega scale factor
#' parameters used to estimate the transcription rates
#'
#' @slot counts a \code{data.frame} with five columns gene_id,
#' summarized_pause_counts, pause_length, summarized_gb_counts, gb_length
#' @slot bigwig_plus a path to bigwig for plus strand
#' @slot bigwig_minus a path to bigwig for minus strand
#' @slot pause_regions a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the pause region coordinates
#' @slot gene_body_regions a \code{\link[GenomicRanges]{GRanges-class}}
#' that holds all the gene body region coordinates
#' @slot gene_name_column a string for the gene name column in the GRanges
#' @slot steric_hindrance a logical value representing whether landing-pad
#' occupancy was inferred when estimating the rates
#' @slot omega_scale a numeric for the scale factor used to calculate omega
#' @slot rates a \code{tbl_df} containing estimated transcription rates such as
#' chi estimates, beta_org estimates from the initial model, beta_adp estimates
#' from the model with varying pause sites, fk_mean giving the mean position
#' of pause sites, fk_var for variance of pause sites, phi estimates for
#' landing-pad occupancy, omega_zeta for the effective initiation rate,
#' beta_zeta for the pause-escape rate, alpha_zeta for the potential initiation
#' rate, and likelihoods
#'
#' @name experiment_transcription_rates-class
#' @rdname experiment_transcription_rates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @exportClass experiment_transcription_rates
methods::setClass("experiment_transcription_rates",
    slots = c(
        counts = "data.frame",
        bigwig_plus = "character",
        bigwig_minus = "character",
        pause_regions = "GRanges",
        gene_body_regions = "GRanges",
        gene_name_column = "character",
        steric_hindrance = "logical",
        omega_scale = "ANY",
        rates = "tbl_df"
    ),
    validity = experiment_transcription_rates_valid
)

#' @keywords internal
input_validation_checks <- function(bigwig_plus, bigwig_minus, pause_regions,
    gene_body_regions, gene_name_column, steric_hindrance, omega_scale) {
    if (!file.exists(bigwig_plus) || !file.exists(bigwig_minus)) {
        stop("bigwig_plus or bigwig_minus file does not exist")
    }
    if (!file.access(bigwig_plus, 4) == 0 ||
        !file.access(bigwig_minus, 4) == 0) {
        stop("bigwig_plus or bigwig_minus file is not readable")
    }
    if (length(pause_regions) == 0 || length(gene_body_regions) == 0) {
        stop("pause_regions or gene_body_regions is empty")
    }
    if (!gene_name_column %in%
        colnames(S4Vectors::elementMetadata(pause_regions)) ||
        !gene_name_column %in%
            colnames(S4Vectors::elementMetadata(gene_body_regions))) {
        stop(sprintf("pause_regions or gene_body_regions does not have a column
        matching %s", gene_name_column))
    }
    if (!is.character(
        S4Vectors::elementMetadata(pause_regions)[, gene_name_column]
    ) || !is.character(
        S4Vectors::elementMetadata(gene_body_regions)[, gene_name_column]
    )) {
        stop(sprintf("gene name column %s must be of class 'character' in
        pause_regions and of gene_body_regions object", gene_name_column))
    }
    duplicated_pause_region_gene_names <-
        any(duplicated(
            S4Vectors::elementMetadata(pause_regions)[,gene_name_column
        ]))
    if (duplicated_pause_region_gene_names) {
        stop("One or more gene names are
        duplicated in pause region, gene names must be unique")
    }
    duplicated_gene_body_region_gene_names <-
        any(duplicated(
            S4Vectors::elementMetadata(gene_body_regions)[,gene_name_column
        ]))
    if (duplicated_gene_body_region_gene_names) {
        stop("One or more gene names are duplicated in gene body region, gene
        names must be unique")
    }
    if (steric_hindrance && (is.null(omega_scale) || !is.numeric(omega_scale) ||
        omega_scale <= 0)) {
        stop("For steric hindrance case, omega_scale parameter must be set to
        numeric greater than 0")
    }
}

#' @keywords internal
prepare_read_count_table <- function(bigwig_plus, bigwig_minus, pause_regions,
                                    gene_body_regions, kmax) {
    rc_cutoff <- 20 # read count cut-off for both gene body and pause peak
    bwp1_p3 <- import.bw(bigwig_plus); bwm1_p3 <- import.bw(bigwig_minus)

    if (sum(bwp1_p3$score) == 0 || sum(bwm1_p3$score) == 0) {
        stop("No reads found in plus or minus strand bigwig file")
    }

    bwp1_p3 <- process_bw(bw = bwp1_p3, strand = "+")
    bwm1_p3 <- process_bw(bw = bwm1_p3, strand = "-")
    bw1_p3 <- c(bwp1_p3, bwm1_p3)
    rm(bwp1_p3, bwm1_p3)

    # make sure pause region is the same as kmax used in EM
    pause_regions <- promoters(pause_regions, upstream = 0, downstream = kmax)

    rc1_pause <- summarise_bw(bw1_p3, pause_regions, "summarized_pause_counts")
    rc1_pause$pause_length <- kmax

    rc1_gb <- summarise_bw(bw1_p3, gene_body_regions, "summarized_gb_counts")
    rc1_gb$gb_length <- width(gene_body_regions)[match(
            rc1_gb$gene_id,gene_body_regions$gene_id)]

    rc1 <- Reduce(
        function(x, y) merge(x, y, by = "gene_id", all = TRUE),
        list(rc1_pause, rc1_gb)
    )

    # clean up some genes with missing values in tss length or gene body length
    rc1 <- rc1[!(is.na(rc1$pause_length) | is.na(rc1$gb_length)), ]
    rc1 <- rc1[(rc1$summarized_pause_counts > rc_cutoff) &
        (rc1$summarized_gb_counts > rc_cutoff), ]

    return(list(rc1 = rc1, bw1_p3 = bw1_p3))
}

#' @keywords internal
prepare_em_data <- function(rc1, bw1_p3, pause_regions, kmin, kmax, 
                            steric_hindrance, omega_scale, zeta) {
    # Create base DataFrame
    em_rate <- DataFrame(
        gene_id = rc1$gene_id, 
        s = rc1$summarized_gb_counts, 
        N = rc1$gb_length
    )
    
    # Calculate chi and get read counts
    em_rate$chi <- em_rate$s / em_rate$N
    Xk <- GenomicRanges::coverage(bw1_p3, weight = "score")[pause_regions]
    
    # Process read counts for each region
    Xk_list <- lapply(seq_along(pause_regions), function(i) {
        region <- pause_regions[i]
        counts <- as.numeric(Xk[[seqnames(region)]][ranges(region)])
        names(counts) <- start(region):end(region)
        counts
    })
    names(Xk_list) <- pause_regions$gene_id
    em_rate$Xk <- Xk_list[em_rate$gene_id]
    
    # Initialize beta
    em_rate$Xk_sum <- vapply(em_rate$Xk, sum, numeric(1))
    em_rate$beta_int <- em_rate$chi / em_rate$Xk_sum
    
    # Handle steric hindrance
    if (steric_hindrance) {
        em_rate$omega_zeta <- em_rate$chi * omega_scale
        em_rate$omega <- em_rate$omega_zeta / zeta
    }
    
    return(em_rate)
}

#' @keywords internal
run_em_algorithm <- function(em_rate, kmin, kmax, fk_int, steric_hindrance, 
                            zeta, lambda = NULL) {
    em_ls <- list()
    for (i in seq_len(NROW(em_rate))) {
        rc <- em_rate[i, ]
        if (!steric_hindrance) {
            em_ls[[i]] <- pause_escape_EM(
                Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                fk_int = fk_int, beta_int = rc$beta_int[[1]],
                chi_hat = rc$chi, max_itr = 500, tor = 1e-4
            )
        } else {
            em_ls[[i]] <- steric_hindrance_EM(
                Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax, 
                f1 = 0.517, f2 = 0.024, fk_int = fk_int, 
                beta_int = rc$beta_int[[1]], phi_int = 0.5, 
                chi_hat = rc$chi, lambda = lambda, zeta = zeta,
                max_itr = 500, tor = 1e-4
            )
        }
    }
    names(em_ls) <- em_rate$gene_id
    return(em_ls)
}

#' @keywords internal
process_em_results <- function(em_rate, em_ls, steric_hindrance, zeta) {
    # Extract results
    em_rate$beta_adp <- map_dbl(em_ls, "beta", .default = NA)
    em_rate$Yk <- map(em_ls, "Yk", .default = NA)
    em_rate$fk <- map(em_ls, "fk", .default = NA)
    em_rate$fk_mean <- map_dbl(em_ls, "fk_mean", .default = NA)
    em_rate$fk_var <- map_dbl(em_ls, "fk_var", .default = NA)
    
    # Calculate proportions and likelihood
    em_rate$t <- vapply(em_rate$Yk, sum, numeric(1))
    em_rate$proportion_Yk <- em_rate$t / vapply(em_rate$Xk, sum, numeric(1))
    em_rate$likelihood <- map_dbl(em_ls, ~ .x$likelihoods[[length(
        x$likelihoods)]])
    
    # Convert to tibble and handle steric hindrance
    em_rate <- as_tibble(em_rate)
    if (steric_hindrance) {
        em_rate$phi <- map_dbl(em_ls, "phi", .default = NA)
        em_rate <- em_rate %>%
            mutate(
                beta_zeta = beta_adp * zeta,
                alpha_zeta = omega_zeta / (1 - phi)
            )
    }
    
    return(em_rate)
}

#' @keywords internal
estimate_em_rates <- function(rc1, bw1_p3, pause_regions, kmin, kmax, fk_int,
                            steric_hindrance, omega_scale, zeta) {
    # Prepare data
    em_rate <- prepare_em_data(rc1, bw1_p3, pause_regions, kmin, kmax,
                            steric_hindrance, omega_scale, zeta)
    
    # Run EM algorithm
    lambda <- if (steric_hindrance) zeta^2 / omega_scale else NULL
    em_ls <- run_em_algorithm(em_rate, kmin, kmax, fk_int, steric_hindrance,
                            zeta, lambda)
    
    # Process results
    em_rate <- process_em_results(em_rate, em_ls, steric_hindrance, zeta)
    
    return(em_rate)
}

#' @keywords internal
prepare_rate_table <- function(em_rate, analytical_rate_tbl, steric_hindrance) {
    em_rate <- em_rate %>% left_join(analytical_rate_tbl, by = "gene_id")

    if (!steric_hindrance) {
        em_rate <- em_rate %>%
            select(gene_id, chi, beta_org, beta_adp, fk_mean, fk_var)
    } else {
        em_rate <- em_rate %>%
            select(
                gene_id, chi, beta_org, beta_adp, fk_mean, fk_var, phi,
                omega_zeta
            ) %>%
            mutate(
                beta_zeta = beta_adp * zeta, alpha_zeta = omega_zeta / (1 - phi)
            )
    }

    return(em_rate)
}

#' estimate_experiment_transcription_rates
#'
#' Estimates the transcription rates, such as initiation, pause-release rates
#' and landing pad occupancy, from experimental data, such as nascent RNA
#' sequencing read counts and genomic coordinates, and contructs an object
#' that holds these rates
#'
#' @param bigwig_plus the path to a bigwig file from the plus strand recording
#' PRO-seq read counts
#' @param bigwig_minus the path to a bigwig file from the minus strand recording
#' PRO-seq read counts
#' @param pause_regions a \link[GenomicRanges]{GRanges-class} object that must
#' contain a gene_id
#' @param gene_body_regions a \link[GenomicRanges]{GRanges-class} object that
#' must contain a gene_id
#' @param gene_name_column a string that indicates which column in the GRanges
#' gene names
#' @param steric_hindrance a logical value to determine whether to infer
#' landing-pad occupancy or not. Defaults to FALSE.
#' @param omega_scale a numeric value for scaling omega. Defaults to NULL.
#'
#' @return an \code{\link{experiment_transcription_rates-class}} object
#'
#' @export
estimate_experiment_transcription_rates <- function(bigwig_plus, bigwig_minus,
pause_regions, gene_body_regions, gene_name_column, steric_hindrance = FALSE,
omega_scale = NULL) {

    input_validation_checks(
        bigwig_plus, bigwig_minus, pause_regions,
        gene_body_regions, gene_name_column, steric_hindrance, omega_scale
    )

    # Force copy underlying GRanges obj to prevent modify in place side effects
    pause_regions <- GenomicRanges::makeGRangesFromDataFrame(
        data.table::copy(data.table::as.data.table(pause_regions)),
        keep.extra.columns = TRUE
    )
    gene_body_regions <- GenomicRanges::makeGRangesFromDataFrame(
        data.table::copy(data.table::as.data.table(gene_body_regions)),
        keep.extra.columns = TRUE
    )

    k <- 50; kmin <- 1; kmax <- 200; rnap_size <- 50; zeta <- 2000

    processed_data <- prepare_read_count_table(bigwig_plus, bigwig_minus,
    pause_regions, gene_body_regions, kmax)
    rc1 <- processed_data$rc1; bw1_p3 <- processed_data$bw1_p3

    message("estimating rates...")

    #### Initial model: Poisson-based Maximum Likelihood Estimation ####
    analytical_rate_tbl <- tibble::tibble(gene_id = rc1$gene_id, beta_org =
    (rc1$summarized_gb_counts / rc1$gb_length) / (rc1$summarized_pause_counts /
    rc1$pause_length))

    fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
    fk_int <- fk_int / sum(fk_int)

    em_rate <- estimate_em_rates(rc1, bw1_p3, pause_regions, kmin, kmax, fk_int,
    steric_hindrance, omega_scale, zeta)
    em_rate <- prepare_rate_table(em_rate, analytical_rate_tbl,
    steric_hindrance)

    return(methods::new(
        Class = "experiment_transcription_rates",
        counts = as.data.frame(rc1), bigwig_plus = bigwig_plus,
        bigwig_minus = bigwig_minus, pause_regions = pause_regions,
        gene_body_regions = gene_body_regions, gene_name_column =
        gene_name_column, steric_hindrance = steric_hindrance, omega_scale =
        omega_scale, rates = em_rate
    ))
}

#' Show method for experiment_transcription_rates objects
#' @param object An experiment_transcription_rates object
#' @return NULL (invisibly)
#' @export
#' @examples
#' # Create an experiment_transcription_rates object
#' exp_rates <- estimate_experiment_transcription_rates(
#'     bigwig_plus = "path/to/plus.bw",
#'     bigwig_minus = "path/to/minus.bw",
#'     pause_regions = GRanges("chr1:1-1000"),
#'     gene_body_regions = GRanges("chr1:1-2000"),
#'     gene_name_column = "gene_id"
#' )
#'
#' # Show the object
#' show(exp_rates)
methods::setMethod("show",
    signature = "experiment_transcription_rates",
    function(object) {
        cat("An experiment_transcription_rates object with:\n")
        cat("  -", length(unique(counts(object)$gene_id)), "genes\n")
        cat("  -", nrow(rates(object)), "rate estimates\n")
        cat("  - Steric hindrance:", steric_hindrance(object), "\n")
        if (steric_hindrance(object)) {
            cat("  - Omega scale:", omega_scale(object), "\n")
        }
    }
)


#' Get transcription rates from experiment_transcription_rates object
#'
#' @rdname experiment_transcription_rates-class
#' @description Retrieves the transcription rates from an
#' experiment_transcription_rates object
#' @param object An experiment_transcription_rates object
#' @return A tibble containing the transcription rates
#' @export
setGeneric("getRates", function(object) standardGeneric("getRates"))
setMethod("getRates", "experiment_transcription_rates", function(object) {
    slot(object, "rates")
})

#' Get read counts
#'
#' @param object An experiment_transcription_rates object
#' @return A data.frame containing the read counts
#' @export
setGeneric("getCounts", function(object) standardGeneric("getCounts"))
setMethod("getCounts", "experiment_transcription_rates", function(object) {
    counts(object)
})

#' Get genomic regions
#'
#' @param object An experiment_transcription_rates object
#' @param type Either "pause" or "gene_body" to specify which regions to return
#' @return A GRanges object containing the specified regions
#' @export
setGeneric("getRegions", function(object, type) standardGeneric("getRegions"))
setMethod(
    "getRegions", "experiment_transcription_rates",
    function(object, type) {
        if (!type %in% c("pause", "gene_body")) {
            stop("type must be either 'pause' or 'gene_body'")
        }
        if (type == "pause") {
            pause_regions(object)
        } else {
            gene_body_regions(object)
        }
    }
)

#' Export rates to CSV
#'
#' @param object An experiment_transcription_rates object
#' @param file Path to output CSV file
#' @return Outputs a CSV file with the rates
#' @export
setGeneric("exportRatesToCSV", function(object, file) {
    standardGeneric("exportRatesToCSV")
})
setMethod(
    "exportRatesToCSV", "experiment_transcription_rates",
    function(object, file) {
        write.csv(rates(object), file = file, row.names = FALSE)
    }
)

#' @keywords internal
create_scatter_plot <- function(data, rate_type) {
    ggplot2::ggplot(data, ggplot2::aes(
        x = .data$beta_org,
        y = .data[[rate_type]]
    )) +
        ggplot2::geom_point(color = "#1E88E5", alpha = 0.7, size = 2) +
        ggplot2::labs(x = "Original Beta", y = rate_type) +
        apply_common_theme()
}

#' @keywords internal
create_histogram_plot <- function(data, rate_type) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[rate_type]])) +
        ggplot2::geom_histogram(
            bins = 30, fill = "#1E88E5",
            color = "white", alpha = 0.7
        ) +
        ggplot2::labs(x = rate_type, y = "Count") +
        apply_common_theme()
}

#' @keywords internal
create_density_plot <- function(data, rate_type) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[rate_type]])) +
        ggplot2::geom_density(
            fill = "#1E88E5", color = "#0D47A1",
            alpha = 0.7
        ) +
        ggplot2::labs(x = rate_type, y = "Density") +
        apply_common_theme()
}

#' @keywords internal
apply_common_theme <- function() {
    ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_line(color = "gray90"),
            panel.grid.minor = ggplot2::element_line(color = "gray95"),
            axis.text = ggplot2::element_text(color = "black", size = 12),
            axis.title = ggplot2::element_text(color = "black", size = 14)
        )
}

#' Plot transcription rates
#'
#' @param object An experiment_transcription_rates object
#' @param type Type of plot to create ("scatter", "histogram", or "density")
#' @param rate_type Which rate to plot ("beta_org", "beta_adp", "chi", etc.)
#' @param file Optional path to save the plot. If provided, the plot will be
#' saved to this location.
#' @param width Width of the saved plot in inches. Default is 8.
#' @param height Height of the saved plot in inches. Default is 6.
#' @param dpi Resolution of the saved plot. Default is 300.
#' @param ... Additional arguments passed to the plotting function
#' @return A ggplot object
#' @export
setGeneric("plotRates", function(
    object, type = "scatter", rate_type = "beta_adp", file = NULL, width = 8,
    height = 6, dpi = 300, ...) {
    standardGeneric("plotRates")
})

setMethod("plotRates", "experiment_transcription_rates", function(
    object, type = "scatter", rate_type = "beta_adp", file = NULL, width = 8,
    height = 6, dpi = 300, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    if (!rate_type %in% colnames(rates(object))) {
        stop(sprintf("rate_type %s not found in rates data", rate_type))
    }

    data <- rates(object)

    p <- switch(type,
        scatter = create_scatter_plot(data, rate_type),
        histogram = create_histogram_plot(data, rate_type),
        density = create_density_plot(data, rate_type),
        stop("Invalid plot type. Choose from 'scatter', 'histogram', or
        'density'")
    )

    if (!is.null(file)) {
        ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi)
    }

    return(p)
})

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("rates", function(object) standardGeneric("rates"))

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("counts", function(object) standardGeneric("counts"))

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("pause_regions", function(object) standardGeneric("pause_regions"))

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("gene_body_regions", function(object) {
    standardGeneric("gene_body_regions")
})

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("gene_name_column", function(object) {
    standardGeneric("gene_name_column")
})

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("steric_hindrance", function(object) {
    standardGeneric("steric_hindrance")
})

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("omega_scale", function(object) standardGeneric("omega_scale"))

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("bigwig_plus", function(object) standardGeneric("bigwig_plus"))

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("bigwig_minus", function(object) standardGeneric("bigwig_minus"))

#' @rdname experiment_transcription_rates-class
#' @export
setMethod("rates", "experiment_transcription_rates", function(object) {
    slot(object, "rates")
})

#' @rdname experiment_transcription_rates-class
#' @export
setMethod("counts", "experiment_transcription_rates", function(object) {
    slot(object, "counts")
})

#' @rdname experiment_transcription_rates-class
#' @export
setMethod("pause_regions", "experiment_transcription_rates", function(object) {
    slot(object, "pause_regions")
})

#' @rdname experiment_transcription_rates-class
#' @export
setMethod(
    "gene_body_regions", "experiment_transcription_rates",
    function(object) {
        slot(object, "gene_body_regions")
    }
)

#' @rdname experiment_transcription_rates-class
#' @export
setMethod(
    "gene_name_column", "experiment_transcription_rates",
    function(object) {
        slot(object, "gene_name_column")
    }
)

#' @rdname experiment_transcription_rates-class
#' @export
setMethod(
    "steric_hindrance", "experiment_transcription_rates",
    function(object) {
        slot(object, "steric_hindrance")
    }
)

#' @rdname experiment_transcription_rates-class
#' @export
setMethod("omega_scale", "experiment_transcription_rates", function(object) {
    slot(object, "omega_scale")
})

#' @rdname experiment_transcription_rates-class
#' @export
setMethod("bigwig_plus", "experiment_transcription_rates", function(object) {
    slot(object, "bigwig_plus")
})

#' @rdname experiment_transcription_rates-class
#' @export
setMethod("bigwig_minus", "experiment_transcription_rates", function(object) {
    slot(object, "bigwig_minus")
})

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("save_rates", function(object, file) standardGeneric("save_rates"))

#' @rdname experiment_transcription_rates-class
#' @export
setMethod(
    "save_rates", "experiment_transcription_rates",
    function(object, file) {
        write.csv(rates(object), file = file, row.names = FALSE)
    }
)

#' @rdname experiment_transcription_rates-class
#' @export
setGeneric("plot_rate", function(object, rate_type) {
    standardGeneric("plot_rate")
})

#' @rdname experiment_transcription_rates-class
#' @export
setMethod(
    "plot_rate", "experiment_transcription_rates",
    function(object, rate_type) {
        if (!rate_type %in% colnames(rates(object))) {
            stop("rate_type must be a column name in rates")
        }
        data <- rates(object)
        ggplot(data, aes_string(x = rate_type)) +
            geom_histogram(bins = 30) +
            theme_minimal() +
            labs(
                title = paste("Distribution of", rate_type),
                x = rate_type,
                y = "Count"
            )
    }
)
