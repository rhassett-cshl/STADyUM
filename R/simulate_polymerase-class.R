#' @keywords internal
simulatePolymeraseValid <- function(object) {
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
    if (ncol(positionMatrix(object)) != cellNum(object) ||
        nrow(positionMatrix(object)) != geneLen(object) + 1) {
        errors <- c(errors, "positionMatrix dimensions do not match cellNum
        and geneLen")
    }

    if (length(errors) == 0) TRUE else errors
}

#' @keywords internal
validateSimulatePolymeraseParams <- function(
    k, ksd, kMin, kMax, geneLen,
    alpha, beta, zeta, zetaSd, zetaMin, zetaMax, cellNum, polSize,
    addSpace, time, stepsToRecord) {
    errors <- character()

    # Check parameter ranges
    if (kMin > kMax) {
        errors <- c(errors, "kMin must be less than kMax")
    }
    if (k < kMin || k > kMax) {
        errors <- c(errors, "k must be between kMin and kMax")
    }
    if (zetaMin > zetaMax) {
        errors <- c(errors, "zetaMin must be less than zetaMax")
    }
    if (zeta < zetaMin || zeta > zetaMax) {
        errors <- c(errors, "zeta must be between zetaMin and zetaMax")
    }
    if (ksd <= 0) {
        errors <- c(errors, "ksd must be positive")
    }
    if (zetaSd <= 0) {
        errors <- c(errors, "zetaSd must be positive")
    }
    if (geneLen <= 0) {
        errors <- c(errors, "geneLen must be positive")
    }
    if (cellNum <= 0) {
        errors <- c(errors, "cellNum must be positive")
    }
    if (polSize <= 0) {
        errors <- c(errors, "polSize must be positive")
    }
    if (addSpace < 0) {
        errors <- c(errors, "addSpace must be non-negative")
    }
    if (time <= 0) {
        errors <- c(errors, "time must be positive")
    }
    if (stepsToRecord < 0) {
        errors <- c(errors, "stepsToRecord must be non-negative")
    }

    if (length(errors) > 0) {
        stop(sprintf("%s", paste(errors, collapse = "\n")))
    }
}

#' Class simulatePolymerase
#'
#' Class \code{simulatePolymerase} tracks movement of RNAPs along the DNA
#' templates of a large number of cells. It accepts several key user-specified
#' parameters, including the initiation rate, pause-escape rate, a constant or
#' variable elongation rate, the mean and variance of pause sites across cells,
#' as well as the center-to-center spacing constraint between RNAPs, the number
#' of cells being simulated, the gene length, and the total time of
#' transcription. The simulator simply allows each RNAP to move forward or not,
#' in time slices of 1e-4 minutes, according to the specified position-specific
#' rate parameters. It assumes that at most one movement of each RNAP can occur
#' per time slice. The simulator monitors for collisions between adjacent
#' RNAPs, prohibiting one RNAP to advance if it is at the boundary of the
#' allowable distance from the next. After running for specified time, it
#' outputs the position of RNAPs for the last specified number of steps.
#'
#' @slot k an integer value for the mean of pause sites across cells.
#' @slot ksd a numeric value for the standard deviation of pause sites across
#' cells.
#' @slot kMin an integer value for the upper bound of pause sites allowed.
#' @slot kMax an integer value for the lower bound of pause sites allowed.
#' @slot geneLen an integer value for the length of the gene.
#' @slot alpha a numeric value for the initiation rate.
#' @slot beta a numeric value for the pause release rate.
#' @slot zeta a numeric value for the mean elongation rate across sites.
#' @slot zetaSd a numeric value for the standard deviation of pause sites
#' across sites.
#' @slot zetaMin a numeric value for the minimum elongation rate.
#' @slot zetaMax a numeric value for the maximum elongation rate.
#' @slot cellNum an integer value for the number of cells to simulate.
#' @slot polSize an integer value for the polymerase II size.
#' @slot addSpace an integer value for the additional space in addition to
#' RNAP size.
#' @slot time a numeric value for the time to simulate.
#' @slot stepsToRecord an integer value for the number of steps to record in
#' position matrix.
#' @slot pauseSites a numeric vector of pause sites
#' @slot probabilityVector a numeric vector
#' @slot combinedCellsData an integer vector
#' @slot positionMatrix a matrix of position of polymerase
#' @slot readCounts a numeric vector for read counts

#' @name simulatePolymerase-class
#' @rdname simulatePolymerase-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs
#' @importFrom reshape2 melt
#' @exportClass simulatePolymerase
methods::setClass("simulatePolymerase",
    slots = c(
        k = "integer", ksd = "numeric", kMin = "integer", kMax = "integer",
        geneLen = "integer", alpha = "numeric", beta = "numeric",
        zeta = "numeric", zetaSd = "numeric", zetaMin = "numeric",
        zetaMax = "numeric", cellNum = "integer", polSize = "integer",
        addSpace = "integer", time = "numeric", deltaT = "numeric",
        stepsToRecord = "integer", pauseSites = "numeric",
        probabilityVector = "numeric", combinedCellsData = "integer",
        positionMatrix = "matrix", readCounts = "numeric",
        avgReadDensity = "numeric"
    ),
    validity = simulatePolymeraseValid
)

# Accessor methods
#' @rdname simulatePolymerase-class
#' @export
setGeneric("pauseSites", function(object) standardGeneric("pauseSites"))
setMethod("pauseSites", "simulatePolymerase", function(object) {
    slot(object, "pauseSites")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("probabilityVector", function(object) {
    standardGeneric("probabilityVector")
})
setMethod("probabilityVector", "simulatePolymerase", function(object) {
    slot(object, "probabilityVector")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("combinedCellsData", function(object) {
    standardGeneric("combinedCellsData")
})
setMethod("combinedCellsData", "simulatePolymerase", function(object) {
    slot(object, "combinedCellsData")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("positionMatrix", function(object) {
    standardGeneric("positionMatrix")
})
setMethod("positionMatrix", "simulatePolymerase", function(object) {
    slot(object, "positionMatrix")
})

#' Get all simulation parameters
#' @param object A simulatePolymerase object
#' @return A list containing all simulation parameters
#' @export
setGeneric("getParameters", function(object) standardGeneric("getParameters"))
setMethod("getParameters", "simulatePolymerase", function(object) {
    list(
        k = object@k,
        ksd = object@ksd,
        kMin = object@kMin,
        kMax = object@kMax,
        geneLen = object@geneLen,
        alpha = object@alpha,
        beta = object@beta,
        zeta = object@zeta,
        zetaSd = object@zetaSd,
        zetaMin = object@zetaMin,
        zetaMax = object@zetaMax,
        cellNum = object@cellNum,
        polSize = object@polSize,
        addSpace = object@addSpace,
        time = object@time,
        stepsToRecord = object@stepsToRecord
    )
})

#' Get pause sites as a data frame
#' @param object A simulatePolymerase object
#' @return A data frame with cell numbers and their pause sites
#' @export
setGeneric("getPauseSitesDf", function(object) {
    standardGeneric("getPauseSitesDf")
})
setMethod("getPauseSitesDf", "simulatePolymerase", function(object) {
    data.frame(
        cell = seq_len(object@cellNum),
        pauseSite = object@pauseSites
    )
})

#' Get probability vector as a data frame
#' @param object A simulatePolymerase object
#' @return A data frame with positions and their transition probabilities
#' @export
setGeneric("getProbabilityDf", function(object) {
    standardGeneric("getProbabilityDf")
})
setMethod("getProbabilityDf", "simulatePolymerase", function(object) {
    data.frame(
        position = 0:object@geneLen,
        probability = object@probabilityVector
    )
})

#' Get combined cells data as a data frame
#' @param object A simulatePolymerase object
#' @return A data frame with positions and their polymerase counts
#' @export
setGeneric("getPolymeraseCountsDf", function(object) {
    standardGeneric("getPolymeraseCountsDf")
})
setMethod("getPolymeraseCountsDf", "simulatePolymerase", function(object) {
    data.frame(
        position = 0:object@geneLen,
        count = object@combinedCellsData
    )
})

# TODO: polymerase_present??

#' Get position matrix as a tidy data frame
#' @param object A simulatePolymerase object
#' @return A data frame with cell, position, and polymerase presence
#' @export
setGeneric("getPositionDf", function(object) {
    standardGeneric("getPositionDf")
})
setMethod("getPositionDf", "simulatePolymerase", function(object) {
    df <- reshape2::melt(object@positionMatrix)
    colnames(df) <- c("cell", "position", "polymerase_present")
    df
})

# Plotting methods
#' Plot polymerase distribution
#' @param object A simulatePolymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the distribution of polymerases across the
#' gene
#' @export
setGeneric("plotPolymeraseDistribution", function(
    object, file = NULL,
    width = 8, height = 6) {
    standardGeneric("plotPolymeraseDistribution")
})
setMethod(
    "plotPolymeraseDistribution", "simulatePolymerase",
    function(object, file = NULL, width = 8, height = 6) {
        df <- data.frame(
            position = 0:object@geneLen,
            count = object@combinedCellsData
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
    }
)

#' Plot pause site distribution
#' @param object A simulatePolymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the distribution of pause sites
#' @export
setGeneric("plotPauseSites", function(
    object, file = NULL, width = 8,
    height = 6) {
    standardGeneric("plotPauseSites")
})
setMethod("plotPauseSites", "simulatePolymerase", function(
    object,
    file = NULL, width = 8, height = 6) {
    df <- data.frame(
        cell = seq_len(object@cellNum),
        pauseSite = object@pauseSites
    )

    p <- ggplot(df, aes(x = pauseSite)) +
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
#' @param object A simulatePolymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the transition probabilities across the gene
#' @export
setGeneric("plotTransitionProbabilities", function(
    object, file = NULL,
    width = 8, height = 6) {
    standardGeneric("plotTransitionProbabilities")
})
setMethod(
    "plotTransitionProbabilities", "simulatePolymerase",
    function(object, file = NULL, width = 8, height = 6) {
        df <- data.frame(
            position = 0:object@geneLen,
            probability = object@probabilityVector
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
    }
)

#' Plot position matrix heatmap
#' @param object A simulatePolymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the position matrix as a heatmap
#' @export
setGeneric("plotPositionMatrix", function(
    object, file = NULL, width = 8,
    height = 6) {
    standardGeneric("plotPositionMatrix")
})
setMethod("plotPositionMatrix", "simulatePolymerase", function(
    object,
    file = NULL, width = 8, height = 6) {
    df <- reshape2::melt(object@positionMatrix)
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
    }
)

#' Save all data frames to CSV files
#' @param object A simulatePolymerase object
#' @param dir Directory to save the files (default: "results")
#' @return Outputs a CSV file with the data frames
#' @export
setGeneric("saveDataFrames", function(object, dir = "results") {
    standardGeneric("saveDataFrames")
})
setMethod("saveDataFrames", "simulatePolymerase", function(
    object,
    dir = "results") {
    # Create directory if it doesn't exist
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }

    # Save each data frame
    write.csv(getPauseSitesDf(object),
        file.path(dir, "pause_sites.csv"),
        row.names = FALSE
    )
    write.csv(getProbabilityDf(object), file.path(
        dir,
        "transition_probabilities.csv"
    ), row.names = FALSE)
    write.csv(getPolymeraseCountsDf(object), file.path(
        dir,
        "polymerase_counts.csv"
    ), row.names = FALSE)
    write.csv(getPositionDf(object), file.path(dir, "position_matrix.csv"),
        row.names = FALSE
    )

    # Save parameters
    write.csv(as.data.frame(getParameters(object)), file.path(
        dir,
        "parameters.csv"
    ), row.names = TRUE)
})

#' Save all plots to files
#' @param object A simulatePolymerase object
#' @param dir Directory to save the plots (default: "results")
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return Outputs directory where plots were saved to pdf files
#' @export
setGeneric("savePlots", function(
    object, dir = "results", width = 8,
    height = 6) {
    standardGeneric("savePlots")
})
setMethod("savePlots", "simulatePolymerase", function(
    object,
    dir = "results", width = 8, height = 6) {
    # Create directory if it doesn't exist
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }

    # Save each plot
    plotPolymeraseDistribution(object, file.path(
        dir,
        "polymerase_distribution.pdf"
    ), width, height)
    plotPauseSites(
        object, file.path(dir, "pause_sites_distribution.pdf"),
        width, height
    )
    plotTransitionProbabilities(object, file.path(
        dir,
        "transition_probabilities.pdf"
    ), width, height)
    plotPositionMatrix(
        object, file.path(dir, "position_matrix.pdf"), width,
        height
    )
})


#' @name simulatePolymerase
#' @rdname simulatePolymerase-class
#' @export
simulatePolymerase <- function(
    k, ksd, kMin, kMax, geneLen,
    alpha, beta, zeta, zetaSd, zetaMin, zetaMax, cellNum, polSize,
    addSpace, time, stepsToRecord) {
    # Validate parameters
    validateSimulatePolymeraseParams(
        k, ksd, kMin, kMax, geneLen,
        alpha, beta, zeta, zetaSd, zetaMin, zetaMax, cellNum, polSize,
        addSpace, time, stepsToRecord
    )

    # Call the C++ function
    result <- simulatePolymeraseCpp(
        k, kMin, kMax, ksd, geneLen,
        alpha, beta, zeta, zetaSd,
        zetaMax, zetaMin, cellNum,
        polSize, addSpace, time, stepsToRecord
    )

    # Create and return a simulatePolymerase object
    obj <- new("simulatePolymerase",
        probabilityVector = result$probabilityVector,
        pauseSites = result$pauseSites,
        combinedCellsData = result$combinedCellsData,
        positionMatrix = result$positionMatrix,
        k = as.integer(k),
        kMin = as.integer(kMin),
        kMax = as.integer(kMax),
        ksd = ksd,
        geneLen = as.integer(geneLen),
        alpha = alpha,
        beta = beta,
        zeta = zeta,
        zetaSd = zetaSd,
        zetaMax = zetaMax,
        zetaMin = zetaMin,
        cellNum = as.integer(cellNum),
        polSize = as.integer(polSize),
        addSpace = as.integer(addSpace),
        time = time,
        stepsToRecord = as.integer(stepsToRecord),
        readCounts = result$readCounts
    )

    # Validate the object
    validObject(obj)

    return(obj)
}

# TODO add another sample read counts method that takes a different gene length

#' Sample read counts from a simulatePolymerase object
#' @param object A simulatePolymerase object
#' @param readDensity A numeric value for the read density within gene body in
#' Dukler et al._ (2017) for genes with median expression (i.e., 0.0489).
#' @return The read count per nucleotide value
#' @export
setGeneric("sampleReadCountsPerNucleotide", function(
    object,
    readDensity = 0.0489) {
    standardGeneric("sampleReadCountsPerNucleotide")
})
setMethod("sampleReadCountsPerNucleotide", "simulatePolymerase", function(
    object, readDensity = 0.0489) {
    cellNum <- object@cellNum
    kMax <- object@kMax
    totalRnap <- object@combinedCellsData

    N <- length(totalRnap)
    L <- N - kMax

    # To match our simulated read counts to reality, we need to compute a
    # scaling factor lambda. One way of doing it is computing the read density
    # based on real experiments. For example, we have computed the read density
    # within gene body in _Dukler et al._ (2017) for genes with median
    # expression (i.e., 0.0489).

    # If we assume the read counts following a Poisson distribution, we can
    # then sample the read counts with mean equals to the RNAP frequency
    # multiplied by lambda.

    # TODO: handle case if lambda is INF because sum is 0
    lambda <- readDensity / (sum(totalRnap[(kMax + 1):N]) / (L * cellNum))

    rc <- rpois(N, totalRnap / cellNum * lambda)

    gbRc <- sum(rc[(kMax + 1):N])

    # read count per nucleotide
    rcPerNt <- gbRc / L

    # Update the readCounts field in the object
    object@readCounts <- rcPerNt

    # Return the read count per nucleotide value
    return(rcPerNt)
})

#' Sample read counts from a simulatePolymerase object
#' @param object A simulatePolymerase object
#' @param readDensity A numeric value for the read density within gene body in
#' _Dukler et al._ (2017) for genes with median expression (i.e., 0.0489).
#' @return The read count per nucleotide value
#' @export
setGeneric("sampleGeneBodyAvgReadDensity", function(
    object,
    readDensity = 0.0489) {
    standardGeneric("sampleGeneBodyAvgReadDensity")
})
setMethod(
    "sampleGeneBodyAvgReadDensity", "simulatePolymerase",
    function(object, readDensity = 0.0489) {
        cellNum <- object@cellNum
        kMax <- object@kMax
        totalRnap <- object@combinedCellsData

        N <- length(totalRnap)
        L <- N - kMax

        simAvgReadDensity <- sum(totalRnap[(kMax + 1):N]) / L

        object@avgReadDensity <- simAvgReadDensity

        return(simAvgReadDensity)
    }
)

#' @rdname simulatePolymerase-class
#' @export
setMethod("show", "simulatePolymerase", function(object) {
    cat("A simulatePolymerase object with:\n")
    cat("  - k =", k(object), "\n")
    cat("  - ksd =", ksd(object), "\n")
    cat("  - kMin =", kMin(object), "\n")
    cat("  - kMax =", kMax(object), "\n")
    cat("  - geneLen =", geneLen(object), "\n")
    cat("  - alpha =", alpha(object), "\n")
    cat("  - beta =", beta(object), "\n")
    cat("  - zeta =", zeta(object), "\n")
    cat("  - zetaSd =", zetaSd(object), "\n")
    cat("  - zetaMin =", zetaMin(object), "\n")
    cat("  - zetaMax =", zetaMax(object), "\n")
    cat("  - cellNum =", cellNum(object), "\n")
    cat("  - polSize =", polSize(object), "\n")
    cat("  - addSpace =", addSpace(object), "\n")
    cat("  - time =", time(object), "\n")
    cat("  - stepsToRecord =", stepsToRecord(object), "\n")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("plotProbability", function(object) {
    standardGeneric("plotProbability")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("plotCombinedCells", function(object) {
    standardGeneric("plotCombinedCells")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("plotAvgReadDensity", function(object) {
    standardGeneric("plotAvgReadDensity")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("simulateReadCounts", function(object) {
    standardGeneric("simulateReadCounts")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("simulateAvgReadDensity", function(object) {
    standardGeneric("simulateAvgReadDensity")
})

#' @rdname simulatePolymerase-class
#' @export
setMethod("plotProbability", "simulatePolymerase", function(object) {
    df <- data.frame(
        position = 0:geneLen(object),
        probability = probabilityVector(object)
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

#' @rdname simulatePolymerase-class
#' @export
setMethod("plotCombinedCells", "simulatePolymerase", function(object) {
    df <- data.frame(
        position = 0:geneLen(object),
        count = combinedCellsData(object)
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

#' @rdname simulatePolymerase-class
#' @export
setMethod("plotAvgReadDensity", "simulatePolymerase", function(object) {
    df <- data.frame(
        position = 0:geneLen(object),
        density = avgReadDensity(object)
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

#' @rdname simulatePolymerase-class
#' @export
setMethod("simulateReadCounts", "simulatePolymerase", function(object) {
    cellNum <- cellNum(object)
    kMax <- kMax(object)
    totalRnap <- combinedCellsData(object)

    # Calculate read counts per nucleotide
    rcPerNt <- totalRnap / cellNum

    # Store read counts in object
    slot(object, "readCounts") <- rcPerNt

    return(object)
})

#' @rdname simulatePolymerase-class
#' @export
setMethod("simulateAvgReadDensity", "simulatePolymerase", function(object) {
    cellNum <- cellNum(object)
    kMax <- kMax(object)
    totalRnap <- combinedCellsData(object)

    # Calculate average read density
    simAvgReadDensity <- totalRnap / (cellNum * kMax)

    # Store average read density in object
    slot(object, "avgReadDensity") <- simAvgReadDensity

    return(object)
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("k", function(object) standardGeneric("k"))
setMethod("k", "simulatePolymerase", function(object) slot(object, "k"))

#' @rdname simulatePolymerase-class
#' @export
setGeneric("ksd", function(object) standardGeneric("ksd"))
setMethod("ksd", "simulatePolymerase", function(object) slot(object, "ksd"))

#' @rdname simulatePolymerase-class
#' @export
setGeneric("kMin", function(object) standardGeneric("kMin"))
setMethod("kMin", "simulatePolymerase", function(object) {
    slot(object, "kMin")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("kMax", function(object) standardGeneric("kMax"))
setMethod("kMax", "simulatePolymerase", function(object) {
    slot(object, "kMax")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("geneLen", function(object) standardGeneric("geneLen"))
setMethod("geneLen", "simulatePolymerase", function(object) {
    slot(object, "geneLen")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("alpha", function(object) standardGeneric("alpha"))
setMethod("alpha", "simulatePolymerase", function(object) {
    slot(object, "alpha")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("beta", function(object) standardGeneric("beta"))
setMethod("beta", "simulatePolymerase", function(object) slot(object, "beta"))

#' @rdname simulatePolymerase-class
#' @export
setGeneric("zeta", function(object) standardGeneric("zeta"))
setMethod("zeta", "simulatePolymerase", function(object) slot(object, "zeta"))

#' @rdname simulatePolymerase-class
#' @export
setGeneric("zetaSd", function(object) standardGeneric("zetaSd"))
setMethod("zetaSd", "simulatePolymerase", function(object) {
    slot(object, "zetaSd")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("zetaMin", function(object) standardGeneric("zetaMin"))
setMethod("zetaMin", "simulatePolymerase", function(object) {
    slot(object, "zetaMin")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("zetaMax", function(object) standardGeneric("zetaMax"))
setMethod("zetaMax", "simulatePolymerase", function(object) {
    slot(object, "zetaMax")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("cellNum", function(object) standardGeneric("cellNum"))
setMethod("cellNum", "simulatePolymerase", function(object) {
    slot(object, "cellNum")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("polSize", function(object) standardGeneric("polSize"))
setMethod("polSize", "simulatePolymerase", function(object) {
    slot(object, "polSize")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("addSpace", function(object) standardGeneric("addSpace"))
setMethod("addSpace", "simulatePolymerase", function(object) {
    slot(object, "addSpace")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("time", function(object) standardGeneric("time"))
setMethod("time", "simulatePolymerase", function(object) slot(object, "time"))

#' @rdname simulatePolymerase-class
#' @export
setGeneric("stepsToRecord", function(object) {
    standardGeneric("stepsToRecord")
})
setMethod("stepsToRecord", "simulatePolymerase", function(object) {
    slot(object, "stepsToRecord")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("readCounts", function(object) standardGeneric("readCounts"))
setMethod("readCounts", "simulatePolymerase", function(object) {
    slot(object, "readCounts")
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("avgReadDensity", function(object) {
    standardGeneric("avgReadDensity")
})
setMethod("avgReadDensity", "simulatePolymerase", function(object) {
    slot(object, "avgReadDensity")
})
