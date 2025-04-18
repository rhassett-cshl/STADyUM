# Validation function for vector/matrix dimensions
simulate_polymerase_valid <- function(object) {
    errors <- character()
    
    # Check vector lengths
    if (length(object@pause_sites) != object@n) {
        errors <- c(errors, "pause_sites length must match n")
    }
    if (length(object@probability_vector) != object@gene_len + 1) {
        errors <- c(errors, "probability_vector length must match gene_len + 1")
    }
    if (length(object@combined_cells_data) != object@gene_len + 1) {
        errors <- c(errors, "combined_cells_data length must match gene_len + 1")
    }
    if (nrow(object@position_matrix) != object@n || ncol(object@position_matrix) != object@gene_len + 1) {
        errors <- c(errors, "position_matrix dimensions must match n and gene_len + 1")
    }
    
    if (length(errors) == 0) TRUE else errors
}

#' Class simulate_polymerase
#'
#' Class \code{simulate_polymerase} has
#'
#' @slot k an integer value for the mean of pause sites across cells.
#' @slot k_min an integer value for the upper bound of pause sites allowed.
#' @slot k_max an integer value for the lower bound of pause sites allowed.
#' @slot ksd a numeric value for the standard deviation of pause sites across cells.
#' @slot gene_len an integer value for the length of the gene.
#' @slot alpha a numeric value for the initiation rate.
#' @slot beta a numeric value for the pause release rate.
#' @slot zeta a numeric value for the mean elongation rate across sites.
#' @slot zeta_sd a numeric value for the standard deviation of pause sites across sites.
#' @slot zeta_max a numeric value for the maximum elongation rate.
#' @slot zeta_min a numeric value for the minimum elongation rate.
#' @slot n an integer value for the number of cells to simulate.
#' @slot s an integer value for the polymerase II size.    
#' @slot h an integer value for the additional space in addition to RNAP size.
#' @slot time a numeric value for the time to simulate.
#' @slot delta_t a numeric value for the time step.
#' @slot csv_steps_to_record an integer value for the number of steps to record.
#' @slot pause_sites a numeric vector of pause sites
#' @slot probability_vector a numeric vector
#' @slot combined_cells_data an integer vector
#' @slot position_matrix a matrix of position of polymerase

#' @name simulate_polymerase-class
#' @rdname simulate_polymerase-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs
#' @importFrom reshape2 melt
#' @exportClass simulate_polymerase
methods::setClass("simulate_polymerase",
                  slots = c(k = "integer", k_min="integer", k_max="integer", ksd="numeric",
                            gene_len="integer", alpha="numeric", beta="numeric", zeta="numeric",
                            zeta_sd="numeric", zeta_max="numeric", zeta_min="numeric", n="integer",
                            s="integer", h="integer", time="numeric", delta_t="numeric",
                            csv_steps_to_record="integer", pause_sites="numeric", probability_vector="numeric", 
                            combined_cells_data="integer", position_matrix="matrix"),
                  validity = simulate_polymerase_valid
)

# Accessor methods
#' @rdname simulate_polymerase-class
#' @export
setGeneric("pause_sites", function(object) standardGeneric("pause_sites"))
setMethod("pause_sites", "simulate_polymerase", function(object) object@pause_sites)

#' @rdname simulate_polymerase-class
#' @export
setGeneric("probability_vector", function(object) standardGeneric("probability_vector"))
setMethod("probability_vector", "simulate_polymerase", function(object) object@probability_vector)

#' @rdname simulate_polymerase-class
#' @export
setGeneric("combined_cells_data", function(object) standardGeneric("combined_cells_data"))
setMethod("combined_cells_data", "simulate_polymerase", function(object) object@combined_cells_data)

#' @rdname simulate_polymerase-class
#' @export
setGeneric("position_matrix", function(object) standardGeneric("position_matrix"))
setMethod("position_matrix", "simulate_polymerase", function(object) object@position_matrix)

#' Get all simulation parameters
#' @param object A simulate_polymerase object
#' @return A list containing all simulation parameters
#' @export
setGeneric("get_parameters", function(object) standardGeneric("get_parameters"))
setMethod("get_parameters", "simulate_polymerase", function(object) {
    list(
        k = object@k,
        k_min = object@k_min,
        k_max = object@k_max,
        ksd = object@ksd,
        gene_len = object@gene_len,
        alpha = object@alpha,
        beta = object@beta,
        zeta = object@zeta,
        zeta_sd = object@zeta_sd,
        zeta_max = object@zeta_max,
        zeta_min = object@zeta_min,
        n = object@n,
        s = object@s,
        h = object@h,
        time = object@time,
        delta_t = object@delta_t,
        csv_steps_to_record = object@csv_steps_to_record
    )
})

#' Get pause sites as a data frame
#' @param object A simulate_polymerase object
#' @return A data frame with cell numbers and their pause sites
#' @export
setGeneric("get_pause_sites_df", function(object) standardGeneric("get_pause_sites_df"))
setMethod("get_pause_sites_df", "simulate_polymerase", function(object) {
    data.frame(
        cell = 1:object@n,
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
        cell = 1:object@n,
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
simulate_polymerase <- function(k, k_min, k_max, ksd, gene_len,
                              alpha, beta, zeta, zeta_sd,
                              zeta_max, zeta_min, total_cells,
                              s, h, time, delta_t, csv_steps_to_record) {
    # Validate parameters
    errors <- character()
    
    # Check parameter ranges
    if (k_min >= k_max) {
        errors <- c(errors, "k_min must be less than k_max")
    }
    if (k < k_min || k > k_max) {
        errors <- c(errors, "k must be between k_min and k_max")
    }
    if (zeta_min >= zeta_max) {
        errors <- c(errors, "zeta_min must be less than zeta_max")
    }
    if (zeta < zeta_min || zeta > zeta_max) {
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
    if (total_cells <= 0) {
        errors <- c(errors, "total_cells must be positive")
    }
    if (s <= 0) {
        errors <- c(errors, "s must be positive")
    }
    if (h < 0) {
        errors <- c(errors, "h must be non-negative")
    }
    if (time <= 0) {
        errors <- c(errors, "time must be positive")
    }
    if (delta_t <= 0) {
        errors <- c(errors, "delta_t must be positive")
    }
    if (csv_steps_to_record < 0) {
        errors <- c(errors, "csv_steps_to_record must be non-negative")
    }
    
    if (length(errors) > 0) {
        stop(paste(errors, collapse = "\n"))
    }

    # Call the C++ function
    result <- simulate_polymerase_cpp(k, k_min, k_max, ksd, gene_len,
                                    alpha, beta, zeta, zeta_sd,
                                    zeta_max, zeta_min, total_cells,
                                    s, h, time, delta_t, csv_steps_to_record)

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
        n = as.integer(total_cells),
        s = as.integer(s),
        h = as.integer(h),
        time = time,
        delta_t = delta_t,
        csv_steps_to_record = as.integer(csv_steps_to_record)
    )
    
    # Validate the object
    validObject(obj)
    
    return(obj)
}

# TODO add another sample read counts method that takes a different gene length

#' Sample read counts from a simulate_polymerase object
#' @param object A simulate_polymerase object
#' @param read_density A numeric value for the read density within gene body in _Dukler et al._ (2017) for genes with median expression (i.e., 0.0489).
#' @export
setGeneric("sample_read_counts", function(object, read_density=0.0489) standardGeneric("sample_read_counts"))
setMethod("sample_read_counts", "simulate_polymerase", function(object, read_density=0.0489) {
    cell_n <- object@n
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
    lambda <- read_density / (sum(total_rnap[(k_max + 1):N]) / (L * cell_n))

    set.seed(12345678)

    print(cell_n)
    print(total_rnap)
    print(N)
    print(lambda)
    
    rc <- rpois(N, total_rnap / cell_n * lambda)

    print(rc)

    gb_rc <- sum(rc[(k_max + 1):N])

    # read count per nucleotide
    rc_per_nt <- gb_rc / L

    return(rc_per_nt)

})