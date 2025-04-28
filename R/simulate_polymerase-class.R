# Validation function for vector/matrix dimensions
simulate_polymerase_valid <- function(object) {
    errors <- character()

    # Check vector lengths
    # if (length(object@pause_sites) != object@n) {
    # errors <- c(errors, "pause_sites length must match n")
    # }
    # if (length(object@probability_vector) != object@gene_len + 1) {
    # errors <- c(errors, "probability_vector length must match gene_len + 1")
    # }
    # if (length(object@combined_cells_data) != object@gene_len + 1) {
    # errors <- c(errors, "combined_cells_data length must match gene_len + 1")
    # }
    if (ncol(position_matrix(object)) != cell_num(object) || nrow(position_matrix(object)) != gene_len(object) + 1) {
        errors <- c(errors, "position_matrix dimensions do not match cell_num and gene_len")
    }

    if (length(errors) == 0) TRUE else errors
}

#' Class simulate_polymerase
#'
#' Class \code{simulate_polymerase} tracks movement of RNAPs along the DNA templates
#' of a large number of cells. It accepts several key user-specified parameters,
#' including the initiation rate, pause-escape rate, a constant or variable elongation
#' rate, the mean and variance of pause sites across cells, as well as the center-to-center
#' spacing constraint between RNAPs, the number of cells being simulated, the gene length,
#' and the total time of transcription. The simulator simply allows each RNAP to move forward
#' or not, in time slices of 1e-4 minutes, according to the specified position-specific rate parameters.
#' It assumes that at most one movement of each RNAP can occur per time slice. The simulator monitors for
#' collisions between adjacent RNAPs, prohibiting one RNAP to advance if it is at the boundary
#' of the allowable distance from the next. After running for specified time, it outputs
#' RNAPs for the last specified number of steps.
#'
#' @slot k an integer value for the mean of pause sites across cells.
#' @slot ksd a numeric value for the standard deviation of pause sites across cells.
#' @slot k_min an integer value for the upper bound of pause sites allowed.
#' @slot k_max an integer value for the lower bound of pause sites allowed.
#' @slot gene_len an integer value for the length of the gene.
#' @slot alpha a numeric value for the initiation rate.
#' @slot beta a numeric value for the pause release rate.
#' @slot zeta a numeric value for the mean elongation rate across sites.
#' @slot zeta_sd a numeric value for the standard deviation of pause sites across sites.
#' @slot zeta_min a numeric value for the minimum elongation rate.
#' @slot zeta_max a numeric value for the maximum elongation rate.
#' @slot cell_num an integer value for the number of cells to simulate.
#' @slot pol_size an integer value for the polymerase II size.
#' @slot add_space an integer value for the additional space in addition to RNAP size.
#' @slot time a numeric value for the time to simulate.
#' @slot steps_to_record an integer value for the number of steps to record in position matrix.
#' @slot pause_sites a numeric vector of pause sites
#' @slot probability_vector a numeric vector
#' @slot combined_cells_data an integer vector
#' @slot position_matrix a matrix of position of polymerase
#' @slot read_counts a numeric vector for read counts

#' @name simulate_polymerase-class
#' @rdname simulate_polymerase-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs
#' @importFrom reshape2 melt
#' @exportClass simulate_polymerase
methods::setClass("simulate_polymerase",
    slots = c(
        k = "integer", ksd = "numeric", k_min = "integer", k_max = "integer",
        gene_len = "integer", alpha = "numeric", beta = "numeric", zeta = "numeric",
        zeta_sd = "numeric", zeta_min = "numeric", zeta_max = "numeric", cell_num = "integer",
        pol_size = "integer", add_space = "integer", time = "numeric", delta_t = "numeric",
        steps_to_record = "integer", pause_sites = "numeric", probability_vector = "numeric",
        combined_cells_data = "integer", position_matrix = "matrix", read_counts = "numeric",
        avg_read_density = "numeric"
    ),
    validity = simulate_polymerase_valid
)

# Accessor methods
#' @rdname simulate_polymerase-class
#' @export
setGeneric("pause_sites", function(object) standardGeneric("pause_sites"))
setMethod("pause_sites", "simulate_polymerase", function(object) {
    slot(object, "pause_sites")
})

#' @rdname simulate_polymerase-class
#' @export
setGeneric("probability_vector", function(object) standardGeneric("probability_vector"))
setMethod("probability_vector", "simulate_polymerase", function(object) {
    slot(object, "probability_vector")
})

#' @rdname simulate_polymerase-class
#' @export
setGeneric("combined_cells_data", function(object) standardGeneric("combined_cells_data"))
setMethod("combined_cells_data", "simulate_polymerase", function(object) {
    slot(object, "combined_cells_data")
})

#' @rdname simulate_polymerase-class
#' @export
setGeneric("position_matrix", function(object) standardGeneric("position_matrix"))
setMethod("position_matrix", "simulate_polymerase", function(object) {
    slot(object, "position_matrix")
})

#' Get all simulation parameters
#' @param object A simulate_polymerase object
#' @return A list containing all simulation parameters
#' @export
setGeneric("get_parameters", function(object) standardGeneric("get_parameters"))
setMethod("get_parameters", "simulate_polymerase", function(object) {
    list(
        k = object@k,
        ksd = object@ksd,
        k_min = object@k_min,
        k_max = object@k_max,
        gene_len = object@gene_len,
        alpha = object@alpha,
        beta = object@beta,
        zeta = object@zeta,
        zeta_sd = object@zeta_sd,
        zeta_min = object@zeta_min,
        zeta_max = object@zeta_max,
        cell_num = object@cell_num,
        pol_size = object@pol_size,
        add_space = object@add_space,
        time = object@time,
        steps_to_record = object@steps_to_record
    )
})

#' Get pause sites as a data frame
#' @param object A simulate_polymerase object
#' @return A data frame with cell numbers and their pause sites
#' @export
setGeneric("get_pause_sites_df", function(object) standardGeneric("get_pause_sites_df"))
setMethod("get_pause_sites_df", "simulate_polymerase", function(object) {
    data.frame(
        cell = seq_len(object@cell_num),
        pause_site = object@pause_sites
    )
})

#' Get probability vector as a data frame
#' @param object A simulate_polymerase object
#' @return A data frame with positions and their transition probabilities
#' @export
setGeneric("get_probability_df", function(object) standardGeneric("get_probability_df"))
setMethod("get_probability_df", "simulate_polymerase", function(object) {
    data.frame(
        position = 0:object@gene_len,
        probability = object@probability_vector
    )
})

#' Get combined cells data as a data frame
#' @param object A simulate_polymerase object
#' @return A data frame with positions and their polymerase counts
#' @export
setGeneric("get_polymerase_counts_df", function(object) standardGeneric("get_polymerase_counts_df"))
setMethod("get_polymerase_counts_df", "simulate_polymerase", function(object) {
    data.frame(
        position = 0:object@gene_len,
        count = object@combined_cells_data
    )
})

#' Get position matrix as a tidy data frame
#' @param object A simulate_polymerase object
#' @return A data frame with cell, position, and polymerase presence
#' @export
setGeneric("get_position_df", function(object) standardGeneric("get_position_df"))
setMethod("get_position_df", "simulate_polymerase", function(object) {
    df <- melt(object@position_matrix)
    colnames(df) <- c("cell", "position", "polymerase_present")
    df
})

# Plotting methods
#' Plot polymerase distribution
#' @param object A simulate_polymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the distribution of polymerases across the gene
#' @export
setGeneric("plot_polymerase_distribution", function(object, file = NULL, width = 8, height = 6) standardGeneric("plot_polymerase_distribution"))
setMethod("plot_polymerase_distribution", "simulate_polymerase", function(object, file = NULL, width = 8, height = 6) {
    df <- data.frame(
        position = 0:object@gene_len,
        count = object@combined_cells_data
    )

    p <- ggplot(df, aes(x = position, y = count)) +
        geom_line() +
        theme_minimal() +
        labs(
            title = "Polymerase Distribution Across Gene",
            x = "Position",
            y = "Number of Polymerases"
        )

    if (!is.null(file)) {
        ggsave(file, p, width = width, height = height)
    }

    return(p)
})

#' Plot pause site distribution
#' @param object A simulate_polymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the distribution of pause sites
#' @export
setGeneric("plot_pause_sites", function(object, file = NULL, width = 8, height = 6) standardGeneric("plot_pause_sites"))
setMethod("plot_pause_sites", "simulate_polymerase", function(object, file = NULL, width = 8, height = 6) {
    df <- data.frame(
        cell = seq_len(object@cell_num),
        pause_site = object@pause_sites
    )

    p <- ggplot(df, aes(x = pause_site)) +
        geom_histogram(bins = 30) +
        theme_minimal() +
        labs(
            title = "Distribution of Pause Sites",
            x = "Pause Site Position",
            y = "Count"
        )

    if (!is.null(file)) {
        ggsave(file, p, width = width, height = height)
    }

    return(p)
})

#' Plot transition probabilities
#' @param object A simulate_polymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the transition probabilities across the gene
#' @export
setGeneric("plot_transition_probabilities", function(object, file = NULL, width = 8, height = 6) standardGeneric("plot_transition_probabilities"))
setMethod("plot_transition_probabilities", "simulate_polymerase", function(object, file = NULL, width = 8, height = 6) {
    df <- data.frame(
        position = 0:object@gene_len,
        probability = object@probability_vector
    )

    p <- ggplot(df, aes(x = position, y = probability)) +
        geom_line() +
        theme_minimal() +
        labs(
            title = "Transition Probabilities Across Gene",
            x = "Position",
            y = "Transition Probability"
        )

    if (!is.null(file)) {
        ggsave(file, p, width = width, height = height)
    }

    return(p)
})

#' Plot position matrix heatmap
#' @param object A simulate_polymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the position matrix as a heatmap
#' @export
setGeneric("plot_position_matrix", function(object, file = NULL, width = 8, height = 6) standardGeneric("plot_position_matrix"))
setMethod("plot_position_matrix", "simulate_polymerase", function(object, file = NULL, width = 8, height = 6) {
    df <- melt(object@position_matrix)
    colnames(df) <- c("Cell", "Position", "Value")

    p <- ggplot(df, aes(x = Position, y = Cell, fill = Value)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "blue") +
        theme_minimal() +
        labs(
            title = "Polymerase Position Matrix",
            x = "Position",
            y = "Cell"
        )

    if (!is.null(file)) {
        ggsave(file, p, width = width, height = height)
    }

    return(p)
})

#' Save all data frames to CSV files
#' @param object A simulate_polymerase object
#' @param dir Directory to save the files (default: "results")
#' @export
setGeneric("save_data_frames", function(object, dir = "results") standardGeneric("save_data_frames"))
setMethod("save_data_frames", "simulate_polymerase", function(object, dir = "results") {
    # Create directory if it doesn't exist
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }

    # Save each data frame
    write.csv(get_pause_sites_df(object), file.path(dir, "pause_sites.csv"), row.names = FALSE)
    write.csv(get_probability_df(object), file.path(dir, "transition_probabilities.csv"), row.names = FALSE)
    write.csv(get_polymerase_counts_df(object), file.path(dir, "polymerase_counts.csv"), row.names = FALSE)
    write.csv(get_position_df(object), file.path(dir, "position_matrix.csv"), row.names = FALSE)

    # Save parameters
    write.csv(as.data.frame(get_parameters(object)), file.path(dir, "parameters.csv"), row.names = TRUE)
})

#' Save all plots to files
#' @param object A simulate_polymerase object
#' @param dir Directory to save the plots (default: "results")
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @export
setGeneric("save_plots", function(object, dir = "results", width = 8, height = 6) standardGeneric("save_plots"))
setMethod("save_plots", "simulate_polymerase", function(object, dir = "results", width = 8, height = 6) {
    # Create directory if it doesn't exist
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }

    # Save each plot
    plot_polymerase_distribution(object, file.path(dir, "polymerase_distribution.pdf"), width, height)
    plot_pause_sites(object, file.path(dir, "pause_sites_distribution.pdf"), width, height)
    plot_transition_probabilities(object, file.path(dir, "transition_probabilities.pdf"), width, height)
    plot_position_matrix(object, file.path(dir, "position_matrix.pdf"), width, height)
})


#' @name simulate_polymerase
#' @rdname simulate_polymerase-class
#' @export
simulate_polymerase <- function(k, ksd, k_min, k_max, gene_len,
                                alpha, beta, zeta, zeta_sd, zeta_min,
                                zeta_max, cell_num, pol_size, add_space,
                                time, steps_to_record) {
    # Validate parameters
    errors <- character()

    # Check parameter ranges
    if (k_min <= k_max) {
        errors <- c(errors, "k_min must be less than k_max")
    }
    if (k > k_min || k < k_max) {
        errors <- c(errors, "k must be between k_min and k_max")
    }
    if (zeta_min <= zeta_max) {
        errors <- c(errors, "zeta_min must be less than zeta_max")
    }
    if (zeta > zeta_min || zeta < zeta_max) {
        errors <- c(errors, "zeta must be between zeta_min and zeta_max")
    }
    if (ksd <= 0) {
        errors <- c(errors, "ksd must be positive")
    }
    if (zeta_sd <= 0) {
        errors <- c(errors, "zeta_sd must be positive")
    }
    if (gene_len <= 0) {
        errors <- c(errors, "gene_len must be positive")
    }
    if (cell_num <= 0) {
        errors <- c(errors, "cell_num must be positive")
    }
    if (pol_size <= 0) {
        errors <- c(errors, "pol_size must be positive")
    }
    if (add_space < 0) {
        errors <- c(errors, "add_space must be non-negative")
    }
    if (time <= 0) {
        errors <- c(errors, "time must be positive")
    }
    if (steps_to_record < 0) {
        errors <- c(errors, "steps_to_record must be non-negative")
    }

    if (length(errors) > 0) {
        stop(sprintf("%s", paste(errors, collapse = "\n")))
    }

    # Call the C++ function
    result <- simulate_polymerase_cpp(
        k, k_min, k_max, ksd, gene_len,
        alpha, beta, zeta, zeta_sd,
        zeta_max, zeta_min, cell_num,
        pol_size, add_space, time, steps_to_record
    )

    # Create and return a simulate_polymerase object
    obj <- new("simulate_polymerase",
        probability_vector = result$probability_vector,
        pause_sites = result$pause_sites,
        combined_cells_data = result$combined_cells_data,
        position_matrix = result$position_matrix,
        k = as.integer(k),
        k_min = as.integer(k_min),
        k_max = as.integer(k_max),
        ksd = ksd,
        gene_len = as.integer(gene_len),
        alpha = alpha,
        beta = beta,
        zeta = zeta,
        zeta_sd = zeta_sd,
        zeta_max = zeta_max,
        zeta_min = zeta_min,
        cell_num = as.integer(cell_num),
        pol_size = as.integer(pol_size),
        add_space = as.integer(add_space),
        time = time,
        steps_to_record = as.integer(steps_to_record),
        read_counts = result$read_counts
    )

    # Validate the object
    validObject(obj)

    return(obj)
}

# TODO add another sample read counts method that takes a different gene length

#' Sample read counts from a simulate_polymerase object
#' @param object A simulate_polymerase object
#' @param read_density A numeric value for the read density within gene body in _Dukler et al._ (2017) for genes with median expression (i.e., 0.0489).
#' @return The read count per nucleotide value
#' @export
setGeneric("sample_read_counts_per_nucleotide", function(object, read_density = 0.0489) standardGeneric("sample_read_counts_per_nucleotide"))
setMethod("sample_read_counts_per_nucleotide", "simulate_polymerase", function(object, read_density = 0.0489) {
    cell_num <- object@cell_num
    k_max <- object@k_max
    total_rnap <- object@combined_cells_data

    N <- length(total_rnap)
    L <- N - k_max

    # To match our simulated read counts to reality, we need to compute a scaling factor lambda.
    # One way of doing it is computing the read density based on real experiments.
    # For example, we have computed the read density within gene body in _Dukler et al._
    # (2017) for genes with median expression (i.e., 0.0489).

    # If we assume the read counts following a Poisson distribution, we can then sample
    # the read counts with mean equals to the RNAP frequency multiplied by lambda.

    # TODO: handle case if lambda is INF because sum is 0
    lambda <- read_density / (sum(total_rnap[(k_max + 1):N]) / (L * cell_num))

    rc <- rpois(N, total_rnap / cell_num * lambda)

    gb_rc <- sum(rc[(k_max + 1):N])

    # read count per nucleotide
    rc_per_nt <- gb_rc / L

    # Update the read_counts field in the object
    object@read_counts <- rc_per_nt

    # Return the read count per nucleotide value
    return(rc_per_nt)
})

#' Sample read counts from a simulate_polymerase object
#' @param object A simulate_polymerase object
#' @param read_density A numeric value for the read density within gene body in _Dukler et al._ (2017) for genes with median expression (i.e., 0.0489).
#' @return The read count per nucleotide value
#' @export
setGeneric("sample_gene_body_avg_read_density", function(object, read_density = 0.0489) standardGeneric("sample_gene_body_avg_read_density"))
setMethod("sample_gene_body_avg_read_density", "simulate_polymerase", function(object, read_density = 0.0489) {
    cell_num <- object@cell_num
    k_max <- object@k_max
    total_rnap <- object@combined_cells_data

    N <- length(total_rnap)
    L <- N - k_max

    sim_avg_read_density <- sum(total_rnap[(k_max + 1):N]) / L

    object@avg_read_density <- sim_avg_read_density

    return(sim_avg_read_density)
})

#' @rdname simulate_polymerase-class
#' @export
setMethod("show", "simulate_polymerase", function(object) {
    cat("A simulate_polymerase object with:\n")
    cat("  - k =", k(object), "\n")
    cat("  - ksd =", ksd(object), "\n")
    cat("  - k_min =", k_min(object), "\n")
    cat("  - k_max =", k_max(object), "\n")
    cat("  - gene_len =", gene_len(object), "\n")
    cat("  - alpha =", alpha(object), "\n")
    cat("  - beta =", beta(object), "\n")
    cat("  - zeta =", zeta(object), "\n")
    cat("  - zeta_sd =", zeta_sd(object), "\n")
    cat("  - zeta_min =", zeta_min(object), "\n")
    cat("  - zeta_max =", zeta_max(object), "\n")
    cat("  - cell_num =", cell_num(object), "\n")
    cat("  - pol_size =", pol_size(object), "\n")
    cat("  - add_space =", add_space(object), "\n")
    cat("  - time =", time(object), "\n")
    cat("  - steps_to_record =", steps_to_record(object), "\n")
})

#' @rdname simulate_polymerase-class
#' @export
setGeneric("plot_probability", function(object) standardGeneric("plot_probability"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("plot_combined_cells", function(object) standardGeneric("plot_combined_cells"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("plot_avg_read_density", function(object) standardGeneric("plot_avg_read_density"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("simulate_read_counts", function(object) standardGeneric("simulate_read_counts"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("simulate_avg_read_density", function(object) standardGeneric("simulate_avg_read_density"))

#' @rdname simulate_polymerase-class
#' @export
setMethod("plot_probability", "simulate_polymerase", function(object) {
    df <- data.frame(
        position = 0:gene_len(object),
        probability = probability_vector(object)
    )
    ggplot(df, aes(x = position, y = probability)) +
        geom_line() +
        theme_minimal() +
        labs(
            title = "Probability Distribution",
            x = "Position",
            y = "Probability"
        )
})

#' @rdname simulate_polymerase-class
#' @export
setMethod("plot_combined_cells", "simulate_polymerase", function(object) {
    df <- data.frame(
        position = 0:gene_len(object),
        count = combined_cells_data(object)
    )
    ggplot(df, aes(x = position, y = count)) +
        geom_line() +
        theme_minimal() +
        labs(
            title = "Combined Cell Data",
            x = "Position",
            y = "Count"
        )
})

#' @rdname simulate_polymerase-class
#' @export
setMethod("plot_avg_read_density", "simulate_polymerase", function(object) {
    df <- data.frame(
        position = 0:gene_len(object),
        density = avg_read_density(object)
    )
    ggplot(df, aes(x = position, y = density)) +
        geom_line() +
        theme_minimal() +
        labs(
            title = "Average Read Density",
            x = "Position",
            y = "Density"
        )
})

#' @rdname simulate_polymerase-class
#' @export
setMethod("simulate_read_counts", "simulate_polymerase", function(object) {
    cell_num <- cell_num(object)
    k_max <- k_max(object)
    total_rnap <- combined_cells_data(object)
    
    # Calculate read counts per nucleotide
    rc_per_nt <- total_rnap / cell_num
    
    # Store read counts in object
    slot(object, "read_counts") <- rc_per_nt
    
    return(object)
})

#' @rdname simulate_polymerase-class
#' @export
setMethod("simulate_avg_read_density", "simulate_polymerase", function(object) {
    cell_num <- cell_num(object)
    k_max <- k_max(object)
    total_rnap <- combined_cells_data(object)
    
    # Calculate average read density
    sim_avg_read_density <- total_rnap / (cell_num * k_max)
    
    # Store average read density in object
    slot(object, "avg_read_density") <- sim_avg_read_density
    
    return(object)
})

#' @rdname simulate_polymerase-class
#' @export
setGeneric("k", function(object) standardGeneric("k"))
setMethod("k", "simulate_polymerase", function(object) slot(object, "k"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("ksd", function(object) standardGeneric("ksd"))
setMethod("ksd", "simulate_polymerase", function(object) slot(object, "ksd"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("k_min", function(object) standardGeneric("k_min"))
setMethod("k_min", "simulate_polymerase", function(object) slot(object, "k_min"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("k_max", function(object) standardGeneric("k_max"))
setMethod("k_max", "simulate_polymerase", function(object) slot(object, "k_max"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("gene_len", function(object) standardGeneric("gene_len"))
setMethod("gene_len", "simulate_polymerase", function(object) slot(object, "gene_len"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("alpha", function(object) standardGeneric("alpha"))
setMethod("alpha", "simulate_polymerase", function(object) slot(object, "alpha"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("beta", function(object) standardGeneric("beta"))
setMethod("beta", "simulate_polymerase", function(object) slot(object, "beta"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("zeta", function(object) standardGeneric("zeta"))
setMethod("zeta", "simulate_polymerase", function(object) slot(object, "zeta"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("zeta_sd", function(object) standardGeneric("zeta_sd"))
setMethod("zeta_sd", "simulate_polymerase", function(object) slot(object, "zeta_sd"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("zeta_min", function(object) standardGeneric("zeta_min"))
setMethod("zeta_min", "simulate_polymerase", function(object) slot(object, "zeta_min"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("zeta_max", function(object) standardGeneric("zeta_max"))
setMethod("zeta_max", "simulate_polymerase", function(object) slot(object, "zeta_max"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("cell_num", function(object) standardGeneric("cell_num"))
setMethod("cell_num", "simulate_polymerase", function(object) slot(object, "cell_num"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("pol_size", function(object) standardGeneric("pol_size"))
setMethod("pol_size", "simulate_polymerase", function(object) slot(object, "pol_size"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("add_space", function(object) standardGeneric("add_space"))
setMethod("add_space", "simulate_polymerase", function(object) slot(object, "add_space"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("time", function(object) standardGeneric("time"))
setMethod("time", "simulate_polymerase", function(object) slot(object, "time"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("steps_to_record", function(object) standardGeneric("steps_to_record"))
setMethod("steps_to_record", "simulate_polymerase", function(object) slot(object, "steps_to_record"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("read_counts", function(object) standardGeneric("read_counts"))
setMethod("read_counts", "simulate_polymerase", function(object) slot(object, "read_counts"))

#' @rdname simulate_polymerase-class
#' @export
setGeneric("avg_read_density", function(object) standardGeneric("avg_read_density"))
setMethod("avg_read_density", "simulate_polymerase", function(object) slot(object, "avg_read_density"))
