#' @title Constructor for SimulatePolymerase object
#'
#' @description
#' Class \code{SimulatePolymerase} tracks the movement of RNAPs along the DNA
#' templates. It contains the parameters for the simulator as well as the
#' simulator results including the position of RNAPs for the last step, the
#' pause sites, a probability vector containing probability RNAPs move forward
#' or not, and a combined cells data vector containing the total number of
#' RNAPs at each site across all cells
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
#' @slot zetaVec a character value for the path to the zetaVec file.
#' @slot cellNum an integer value for the number of cells to simulate.
#' @slot polSize an integer value for the polymerase II size.
#' @slot addSpace an integer value for the additional space in addition to
#' RNAP size.
#' @slot time a numeric value for the time to simulate.
#' @slot timesToRecord a numeric vector of specific time points to record
#' position matrices for, or NULL to record no extra position matrices. Final
#' position matrix is always recorded. Default is NULL.
#' @slot pauseSites a numeric vector of pause sites
#' @slot siteProbabilities a numeric matrix representing the probability that
#' the polymerase move forward or not at each site for each cell (sites x cells)
#' @slot combinedCellsData an integer vector representing the total number of
#' RNAPs at each site across all cells
#' @slot positionMatrices a list of position matrices
#' @slot finalPositionMatrix a matrix representing the final position matrix
#' @slot readCounts a numeric vector for read counts per nucleotide
#'
#' @name SimulatePolymerase-class
#' @rdname SimulatePolymerase-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @importFrom methods slot new is slot<-
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs
#' @importFrom ggplot2 geom_tile scale_fill_gradient ggsave geom_histogram
#' @importFrom dplyr %>%
#' @importFrom stats prcomp sd
#' @importFrom utils head
#' @exportClass SimulatePolymerase
methods::setClass("SimulatePolymerase",
    slots = c(
        k = "integer", ksd = "numeric", kMin = "integer", kMax = "integer",
        geneLen = "integer", alpha = "numeric", beta = "numeric",
        zeta = "numeric", zetaSd = "numeric", zetaMin = "numeric",
        zetaMax = "numeric", zetaVec = "character", cellNum = "integer",
        polSize = "integer", addSpace = "integer", time = "numeric",
        timesToRecord = "ANY", pauseSites = "numeric",
        siteProbabilities = "matrix", combinedCellsData = "integer",
        positionMatrices = "list", finalPositionMatrix = "matrix",
        readCounts = "ANY"
    )
)

validateSimulatePolymeraseParams <- function(
    k, ksd, kMin, kMax, geneLen,
    alpha, beta, zeta, zetaSd, zetaMin, zetaMax, cellNum, polSize, addSpace) {
    errors <- character()
    if (!is.numeric(k) || k <= 0 || k %% 1 != 0 || k < kMin || k > kMax) {
        errors <- c(errors, "k must be a positive int between kMin and kMax")
    }
    if (!is.numeric(ksd) || ksd < 0 || ksd %% 1 != 0) {
        errors <- c(errors, "ksd must be a non-negative integer")
    }
    if (!is.numeric(kMin) || kMin <= 0 || kMin %% 1 != 0 || kMin >= kMax) {
        errors <- c(errors, "kMin must be a positive integer <= kMax")
    }
    if (!is.numeric(kMax) || kMax <= 0 || kMax %% 1 != 0) {
        errors <- c(errors, "kMax must be a positive integer")
    }
    if (!is.numeric(alpha) || alpha <= 0) {
        errors <- c(errors, "alpha must be a positive number")
    }
    if (!is.numeric(beta) || beta <= 0) {
        errors <- c(errors, "beta must be a positive number")
    }
    if (!is.numeric(zeta) || zeta <= 0 || zeta < zetaMin || zeta > zetaMax) {
        errors <- c(errors, "zeta must be positive and zetaMin<=zeta<=zetaMax")
    }
    if (!is.numeric(zetaSd) || zetaSd < 0) {
        errors <- c(errors, "zetaSd must be a non-negative number")
    }
    if (!is.numeric(zetaMin) || zetaMin <= 0 || zetaMin > zetaMax) {
        errors <- c(errors, "zetaMin must be positive and zetaMin<=zetaMax")
    }
    if (!is.numeric(zetaMax) || zetaMax <= 0) {
        errors <- c(errors, "zetaMax must be a positive number")
    }
    if (!is.numeric(geneLen) || geneLen <= 0 || geneLen %% 1 != 0 ||
        geneLen < kMax) {
        errors <- c(errors, "geneLen must be a positive int greater than kMax")
    }
    if (!is.numeric(cellNum) || cellNum <= 0 || cellNum %% 1 != 0) {
        errors <- c(errors, "cellNum must be a positive integer")
    }
    if (!is.numeric(polSize) || polSize <= 0 || polSize %% 1 != 0) {
        errors <- c(errors, "polSize must be a positive integer")
    }
    if (!is.numeric(addSpace) || addSpace < 0 || addSpace %% 1 != 0) {
        errors <- c(errors, "addSpace must be a non-negative integer")
    }
    if (length(errors) > 0) stop(sprintf("%s", paste(errors, collapse = "\n")))
}

validateTimesToRecord <- function(timesToRecord, time, geneLen, cellNum) {
    errors <- character()

    if (!is.numeric(time) || time < 1e-4) {
        errors <- c(errors, "time must be a positive number greater than 1e-4")
    }
    if (!is.null(timesToRecord)) {
        if (any(timesToRecord < 0) || any(timesToRecord > time)) {
            errors <- c(errors, "timesToRecord must be >= 0 and <= time")
        }
        numTimes <- length(timesToRecord)
        estMemMB <- (numTimes * (geneLen + 1) * cellNum * 4) / (1024 * 1024)
        if (estMemMB > 1024) {
            warning(sprintf("Large memory usage estimated: %.1f MB for position
            matrices. Reduce timesToRecord or cellNum", estMemMB))
        }
        if (estMemMB > 4096) {
            errors <- c(errors, sprintf("Memory usage too high: %.1f MB. Max
            allowed: 4GB", estMemMB))
        }
    }
    if (length(errors) > 0) stop(sprintf("%s", paste(errors, collapse = "\n")))
}

validateAndLoadZetaVec <- function(zetaVec, geneLen) {
    if (is.null(zetaVec)) {
        return(NULL)
    }

    # Check if file exists and is readable
    if (!file.exists(zetaVec)) {
        stop("zetaVec file does not exist")
    }
    if (file.access(zetaVec, 4) != 0) {
        stop("zetaVec file is not readable")
    }

    # Read first line to check column count
    first_line <- readLines(zetaVec, n = 1)
    if (length(strsplit(first_line, ",")[[1]]) != 1) {
        stop("zetaVec file must contain exactly one column")
    }

    # Read the full file
    zeta_data <- read.csv(zetaVec)

    # Convert to numeric vector and ensure correct length
    zeta_vec <- as.numeric(zeta_data[[1]])

    # Validate that all values are greater than 0 and non-negative
    if (any(is.na(zeta_vec)) || any(zeta_vec <= 0)) {
        stop("zetaVec file must contain only positive values")
    }

    # Check that the zeta_vec length matches geneLen
    if (length(zeta_vec) != geneLen) {
        stop(sprintf(
            "zetaVec file length (%d) does not match geneLen (%d)",
            length(zeta_vec), geneLen
        ))
    }

    return(zeta_vec)
}

#' @title Simulator for Tracking RNAP Movement
#'
#' @description
#' Runs the simulator that tracks the movement of RNAPs along the DNA
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
#' @rdname SimulatePolymerase-class
#' @param k an integer value for the mean of pause sites across cells.
#' @param ksd a numeric value for the standard deviation of pause sites across
#' cells.
#' @param kMin an integer value for the upper bound of pause sites allowed.
#' @param kMax an integer value for the lower bound of pause sites allowed.
#' @param geneLen an integer value for the length of the gene.
#' @param alpha a numeric value for the initiation rate.
#' @param beta a numeric value for the pause release rate.
#' @param zeta a numeric value for the mean elongation rate across sites.
#' @param zetaSd a numeric value for the standard deviation of pause sites
#' across sites.
#' @param zetaMin a numeric value for the minimum elongation rate.
#' @param zetaMax a numeric value for the maximum elongation rate.
#' @param zetaVec a character value for the path to the zetaVec file.
#' @param cellNum an integer value for the number of cells to simulate.
#' @param polSize an integer value for the polymerase II size.
#' @param addSpace an integer value for the additional space in addition to
#' RNAP size.
#' @param time a numeric value for the time to simulate.
#' @param timesToRecord a numeric vector of specific time points to record
#' position matrices for, or NULL to record no extra position matrices. Final
#' position matrix is always recorded. Default is NULL.
#' @return a \code{SimulatePolymerase} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500,
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 100, polSize = 33, 
#'     addSpace = 17, time = 1, timesToRecord = NULL
#' )
#' # Create a SimulatePolymerase object with specific time points recorded
#' sim2 <- simulatePolymerase(
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, zetaMax
#' = 2500,
#'     zetaVec = NULL, cellNum = 1000, polSize = 33, addSpace = 17, time = 1,
#'     timesToRecord = c(0.5, 1.0)
#' )
#' @export
simulatePolymerase <- function(
    k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
    alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, 
    zetaMax = 2500,
    zetaVec = NULL, cellNum = 1000, polSize = 33, addSpace = 17, time = 1,
    timesToRecord = NULL) {
    validateSimulatePolymeraseParams(
        k, ksd, kMin, kMax, geneLen, alpha, beta,
        zeta, zetaSd, zetaMin, zetaMax, cellNum, polSize, addSpace
    )
    validateTimesToRecord(timesToRecord, time, geneLen, cellNum)

    zeta_vec <- validateAndLoadZetaVec(zetaVec, geneLen)

    if (cellNum > 1000 && time > 200) {
        message("Note: The simulation may take a few minutes to run due to a
        large number of cells and long simulation time.")
    }

    result <- .Call("_STADyUM_simulate_polymerase_cpp",
        k, ksd, kMin, kMax, geneLen, alpha, beta, zeta, zetaSd,
        zetaMin, zetaMax, cellNum, polSize, addSpace, time, timesToRecord,
        zeta_vec,
        PACKAGE = "STADyUM"
    )

    obj <- new("SimulatePolymerase",
        siteProbabilities = result$probabilityVector,
        pauseSites = result$pauseSites,
        combinedCellsData = result$combinedCellsData,
        positionMatrices = result$positionMatrices,
        finalPositionMatrix = result$finalPositionMatrix,
        k = as.integer(k), kMin = as.integer(kMin), kMax = as.integer(kMax),
        ksd = ksd, geneLen = as.integer(geneLen), alpha = alpha, beta = beta,
        zeta = zeta, zetaSd = zetaSd, zetaMax = zetaMax, zetaMin = zetaMin,
        zetaVec = if (is.null(zetaVec)) "" else zetaVec,
        cellNum = as.integer(cellNum), polSize = as.integer(polSize),
        addSpace = as.integer(addSpace), time = time,
        timesToRecord = timesToRecord, readCounts = NULL
    )

    obj <- sampleReadCounts(obj)

    return(obj)
}

#' @note The `simulatePolymerase` function can take a few minutes to run for
#' large numbers of cells (> 1000) and long simulation times (> 200).

#' @rdname SimulatePolymerase-class
#' @title Sample Read Counts from SimulatePolymerase Object
#'
#' @description
#' Sample read counts from a SimulatePolymerase object. To match our simulated
#' read counts to reality, we need to compute a scaling factor lambda. One way
#' of doing it is computing the read density based on real experiments. For
#' example, we have computed the read density within gene body in
#' _Dukler et al_ (2017) for genes with median expression (i.e., 0.0489).
#' If we assume the read counts following a Poisson distribution, we can then
#' sample the read counts with mean equals to the RNAP frequency multiplied by
#' lambda.
#' @param object A SimulatePolymerase object
#' @param readDensity A numeric value for the read density within gene body in
#' _Dukler et al._ (2017) for genes with median expression (i.e., 0.0489).
#' @return The SimulatePolymerase object with read counts sampled given the
#' readDensity parameter used
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, 
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 100, polSize = 33, 
#'     addSpace = 17, time = 1, timesToRecord = NULL
#' )
#' # Sample read counts per nucleotide
#' sim <- sampleReadCounts(sim)
#' # Print the read counts per nucleotide
#' print(readCounts(sim))
#' @export
setGeneric("sampleReadCounts", function(
    object,
    readDensity = 0.0489) {
    standardGeneric("sampleReadCounts")
})
#' @rdname SimulatePolymerase-class
setMethod("sampleReadCounts", "SimulatePolymerase", function(
    object, readDensity = 0.0489) {
    cellNum <- slot(object, "cellNum")
    kMax <- slot(object, "kMax")
    totalRnap <- slot(object, "combinedCellsData")

    N <- length(totalRnap)
    L <- N - kMax

    # TODO: handle case if lambda is INF because sum is 0
    lambda <- readDensity / (sum(totalRnap[(kMax + 1):N]) / (L * cellNum))

    rc <- rpois(N, totalRnap / cellNum * lambda)

    slot(object, "readCounts") <- rc

    return(object)
})

#' @rdname SimulatePolymerase-class
#' @title Show Method for SimulatePolymerase Object
#'
#' @description
#' Show method for SimulatePolymerase object in human readable format
#' including summary statistics
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @export
setMethod("show", "SimulatePolymerase", function(object) {
    cat("A SimulatePolymerase object with:\n")
    cat("  - gene length =", slot(object, "geneLen"), "\n")
    cat("  - number of cells =", slot(object, "cellNum"), "\n")

    pauseSites <- slot(object, "pauseSites")
    if (length(pauseSites) > 0) {
        pauseMean <- mean(pauseSites)
        pauseSd <- sd(pauseSites)
        cat("  - pause sites mean =", round(pauseMean, 2), "\n")
        cat("  - pause sites std dev =", round(pauseSd, 2), "\n")
    }

    combinedData <- slot(object, "combinedCellsData")
    if (length(combinedData) > 0) {
        siteData <- data.frame(
            position = 0:(length(combinedData) - 1), count = combinedData
        )
        topSites <- siteData[order(siteData$count, decreasing = TRUE), ]
        top10 <- head(topSites, 10)
        cat("  - top 10 most occupied sites across all cells:\n")
        for (i in seq_len(nrow(top10))) {
            cat(sprintf(
                "    position %4d: %d polymerases\n",
                top10$position[i], top10$count[i]
            ))
        }
    }
    readCounts <- slot(object, "readCounts")
    kMax <- slot(object, "kMax")
    N <- length(combinedData)
    L <- N - kMax
    gbRc <- sum(readCounts[(kMax + 1):N])
    gbAvgRc <- gbRc / L
    pauseRc <- sum(readCounts[seq_len(kMax)])
    pauseAvgRc <- pauseRc / kMax
    cat("  - gene body average read counts =", round(gbAvgRc, 2), "\n")
    cat("  - pause region average read counts =", round(pauseAvgRc, 2), "\n")

    zeroReadsFraction <- sum(readCounts == 0) / length(readCounts)
    cat(
        "  - percent of nucleotides with zero reads =",
        round(zeroReadsFraction * 100, 1), "%\n"
    )
    cat("\nTo access the full simulation parameters, use: parameters(object)\n")
    cat("\nTo access the all sampled read counts, use: readCounts(object)\n")
})

# Accessor methods
#' @rdname SimulatePolymerase-class
#' @title Accessor for Pause Sites
#'
#' @description
#' Accessor for the pause sites numeric vector from a SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @export
setGeneric("pauseSites", function(object) standardGeneric("pauseSites"))
#' @rdname SimulatePolymerase-class
setMethod("pauseSites", "SimulatePolymerase", function(object) {
    slot(object, "pauseSites")
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Probability Matrix
#'
#' @description
#' Accessor for the probability numeric matrix from a SimulatePolymerase object.
#' The matrix has dimensions (sites x cells) where each element represents the
#' transition probability for a specific site and cell.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @export
setGeneric("siteProbabilities", function(object) {
    standardGeneric("siteProbabilities")
})
#' @rdname SimulatePolymerase-class
setMethod("siteProbabilities", "SimulatePolymerase", function(object) {
    slot(object, "siteProbabilities")
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Combined Cells Data
#'
#' @description
#' Accessor for the combined cells data numeric vector from a
#' SimulatePolymerase object.
#'
#' @export
setGeneric("combinedCellsData", function(object) {
    standardGeneric("combinedCellsData")
})
#' @rdname SimulatePolymerase-class
setMethod("combinedCellsData", "SimulatePolymerase", function(object) {
    slot(object, "combinedCellsData")
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Position Matrices
#'
#' @description
#' Accessor for the position matrices from a SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @export
setGeneric("positionMatrices", function(object) {
    standardGeneric("positionMatrices")
})
#' @rdname SimulatePolymerase-class
setMethod("positionMatrices", "SimulatePolymerase", function(object) {
    slot(object, "positionMatrices")
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Final Position Matrix
#'
#' @description
#' Accessor for the final position matrix from a SimulatePolymerase object.
#' This represents the last state of the simulation.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @export
setGeneric("finalPositionMatrix", function(object) {
    standardGeneric("finalPositionMatrix")
})
#' @rdname SimulatePolymerase-class
setMethod("finalPositionMatrix", "SimulatePolymerase", function(object) {
    slot(object, "finalPositionMatrix")
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Simulation Parameters
#'
#' @description
#' Accessor for the simulation parameters from a SimulatePolymerase object.
#'
#' @param object A SimulatePolymerase object
#' @return A list containing all simulation parameters including k, ksd, kMin,
#' kMax, geneLen, alpha, beta, zeta, zetaSd, zetaMin, zetaMax, zetaVec, cellNum,
#' polSize, addSpace, time, and stepsToRecord
#' @export
setGeneric("parameters", function(object) standardGeneric("parameters"))
#' @rdname SimulatePolymerase-class
setMethod("parameters", "SimulatePolymerase", function(object) {
    list(
        k = slot(object, "k"),
        ksd = slot(object, "ksd"),
        kMin = slot(object, "kMin"),
        kMax = slot(object, "kMax"),
        geneLen = slot(object, "geneLen"),
        alpha = slot(object, "alpha"),
        beta = slot(object, "beta"),
        zeta = slot(object, "zeta"),
        zetaSd = slot(object, "zetaSd"),
        zetaMin = slot(object, "zetaMin"),
        zetaMax = slot(object, "zetaMax"),
        zetaVec = slot(object, "zetaVec"),
        cellNum = slot(object, "cellNum"),
        polSize = slot(object, "polSize"),
        addSpace = slot(object, "addSpace"),
        time = slot(object, "time"),
        timesToRecord = if (is.null(slot(object, "timesToRecord"))) {
            NULL
        } else if (length(slot(object, "timesToRecord")) == 0) {
            NULL
        } else {
            slot(
                object,
                "timesToRecord"
            )
        }
    )
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Read Counts
#'
#' @description
#' Accessor for the read counts numeric vector sampled from a
#' SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase} object
#' @export
setGeneric("readCounts", function(object) standardGeneric("readCounts"))
#' @rdname SimulatePolymerase-class
setMethod("readCounts", "SimulatePolymerase", function(object) {
    slot(object, "readCounts")
})

#' @rdname SimulatePolymerase-class
#' @title Get Position Matrix for Specific Time Point
#'
#' @description
#' Get the position matrix for a specific time point from a SimulatePolymerase
#' object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @param timePoint a numeric value specifying the time point to retrieve
#' @return A 2D matrix of polymerase positions for the specified time point
#' 
#' @export
setGeneric("getPositionMatrixAtTime", function(object, timePoint) {
    standardGeneric("getPositionMatrixAtTime")
})
#' @rdname SimulatePolymerase-class
setMethod("getPositionMatrixAtTime", "SimulatePolymerase", function(
    object,
    timePoint) {
    matrices <- slot(object, "positionMatrices")
    if (length(matrices) == 0) {
        stop("No extra position matrices were recorded. Use
        finalPositionMatrix() to get the final state.")
    }

    # Extract available time points as numbers
    available_times <- as.numeric(gsub("t_", "", names(matrices)))

    # Find the closest time point within a small tolerance
    tolerance <- 1e-10 # Very small tolerance for floating point comparison
    closest_idx <- which(abs(available_times - timePoint) < tolerance)

    if (length(closest_idx) > 0) {
        # Use the first match if multiple
        time_name <- names(matrices)[closest_idx[1]]
        return(matrices[[time_name]])
    } else {
        stop(sprintf(
            "Time point %f not found. Available time points: %s",
            timePoint, paste(available_times, collapse = ", ")
        ))
    }
})

#' @rdname SimulatePolymerase-class
#' @title Get Available Time Points
#'
#' @description
#' Get all available time points from a SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @return A numeric vector of available time points
#' @export
setGeneric("getAvailableTimePoints", function(object) {
    standardGeneric("getAvailableTimePoints")
})
#' @rdname SimulatePolymerase-class
setMethod("getAvailableTimePoints", "SimulatePolymerase", function(object) {
    matrices <- slot(object, "positionMatrices")
    if (length(matrices) == 0) {
        return(numeric(0))
    }
    available_times <- as.numeric(gsub("t_", "", names(matrices)))
    return(sort(available_times))
})

#' @rdname SimulatePolymerase-class
#' @title Get Combined Cells Data for Specific Time Point
#'
#' @description
#' Get the combined cells data vector for a specific time point from a
#' SimulatePolymerase object. This function calculates the total number of
#' polymerases at each site across all cells for the specified time point,
#' similar to the final combinedCellsData but for intermediate time points.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @param timePoint a numeric value specifying the time point to retrieve.
#' If NULL, returns the final combinedCellsData.
#' @return An integer vector representing the total number of polymerases at
#' each site across all cells for the specified time point
#' @export
setGeneric("getCombinedCellsDataAtTime", function(object, timePoint = NULL) {
    standardGeneric("getCombinedCellsDataAtTime")
})
#' @rdname SimulatePolymerase-class
setMethod("getCombinedCellsDataAtTime", "SimulatePolymerase", function(
    object,
    timePoint = NULL) {
    if (is.null(timePoint)) {
        # Return the final combined cells data
        return(slot(object, "combinedCellsData"))
    } else {
        # Get position matrix for the specified time point
        pos_matrix <- getPositionMatrixAtTime(object, timePoint)

        # Calculate combined cells data by summing across cells for each site
        combined_data <- rowSums(pos_matrix)

        return(combined_data)
    }
})

## PLotting utilities

#' @rdname SimulatePolymerase-class
#' @title Plot Pause Site Distribution
#'
#' @description
#' Plot the distribution of pause sites across the gene as a histogram.
#'
#' @param object A SimulatePolymerase-class object
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
#' @rdname SimulatePolymerase-class
setMethod("plotPauseSites", "SimulatePolymerase", function(
    object,
    file = NULL, width = 8, height = 6) {
    df <- data.frame(
        cell = seq_len(slot(object, "cellNum")),
        pauseSite = slot(object, "pauseSites")
    )

    # Calculate mean and standard deviation
    pauseMean <- mean(df$pauseSite)
    pauseSd <- sd(df$pauseSite)

    p <- ggplot(df, aes(x = pauseSite)) +
        geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
        geom_vline(
            xintercept = pauseMean, color = "red",
            linetype = "dashed", size = 1
        ) +
        geom_vline(
            xintercept = pauseMean + pauseSd, color = "orange",
            linetype = "dotted", size = 0.8
        ) +
        geom_vline(
            xintercept = pauseMean - pauseSd, color = "orange",
            linetype = "dotted", size = 0.8
        ) +
        theme_minimal() +
        labs(
            title = "Distribution of Pause Sites",
            subtitle = sprintf("Mean: %.1f, SD: %.1f", pauseMean, pauseSd),
            x = "Pause Site Position",
            y = "Count"
        ) +
        theme(plot.title = element_text(hjust = 0.5))

    if (!is.null(file)) {
        ggsave(file, p, width = width, height = height)
    }

    return(p)
})

#' @rdname SimulatePolymerase-class
#' @title Plot Position Matrix Heatmap
#'
#' @description
#' Plot position matrices as a ggplot2 heatmap showing polymerase positions
#' across all cells at specific time points. Each cell in the heatmap
#' represents whether a polymerase is present (1) or absent (0) at a specific
#' site in a specific cell. By default, shows the final position matrix. Plots
#' for pause region which is from site 1 to  user-defined kMax.
#' @param object A SimulatePolymerase-class object
#' @param timePoint Optional time point to plot. If NULL, plots the final
#' position matrix.
#' @param maxCells Maximum number of cells to display (for performance with
#' large datasets). If NULL, shows all cells.
#' @param file Optional file path to save the plot to
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Plot resolution in DPI
#' @return A ggplot object showing heatmap of polymerase positions
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, 
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 100, polSize = 33, 
#'     addSpace = 17, time = 1, timesToRecord = c(0.5, 1.0))
#' # Plot final position matrix
#' plotPositionHeatmap(sim, file="position_heatmap.png")
#' # Plot position matrix at time 0.5
#' plotPositionHeatmap(sim, timePoint = 0.5, file="position_heatmap.png")
#' # Plot position matrix with maxCells = 100
#' plotPositionHeatmap(sim, timePoint = 0.5, maxCells = 100,
#' file="position_heatmap.png")
#' @export
setGeneric("plotPositionHeatmap", function(
    object, timePoint = NULL, maxCells = NULL, file = NULL, width = 10, 
    height = 8, dpi = 300) {
    standardGeneric("plotPositionHeatmap")
})
#' @rdname SimulatePolymerase-class
setMethod("plotPositionHeatmap", "SimulatePolymerase", function(
    object, timePoint = NULL, maxCells = NULL, file = NULL, width = 10, 
    height = 8, dpi = 300) {
    if (is.null(timePoint)) {
        matrix <- slot(object, "finalPositionMatrix")
        plotTitle <- "Final Polymerase Positions Heatmap"
    } else {
        matrix <- getPositionMatrixAtTime(object, timePoint)
        plotTitle <- sprintf("Polymerase Positions Heatmap at Time %.2f",
            timePoint
        )
    }
    if (nrow(matrix) == 0 || ncol(matrix) == 0) {
        stop("Position matrix is empty. Run the simulation first.")
    }
    
    # Limit cells for visualization if specified
    if (!is.null(maxCells) && ncol(matrix) > maxCells) {
        cellIndices <- seq(1, ncol(matrix), length.out = maxCells)
        matrix <- matrix[, cellIndices]
    }

    kMax <- slot(object, "kMax")
    matrix <- matrix[seq_len(kMax), ]
    df <- expand.grid(Site = seq_len(nrow(matrix)), 
                        Cell = seq_len(ncol(matrix)))
    df$Polymerase <- as.vector(matrix)

    p <- ggplot(df, aes(x = Cell, y = Site, fill = factor(Polymerase))) +
        geom_tile(color = "grey80") + 
        scale_fill_manual(values = c("0" = "white", "1" = "red"), 
        name = "Polymerase") +
        labs(title = plotTitle, x = "Cell", y = "Site") +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank()
        )

    if (!is.null(file)) {
        ggsave(file, p, width = width, height = height, dpi = dpi)
    } 

    return(p)
})

#' @rdname SimulatePolymerase-class
#' @title PCA Plot of Polymerase Position Matrix
#'
#' @description
#' Create a ggplot2 PCA plot of the polymerase position matrix, showing the
#' first two principal components for each cell. Useful for exploring
#' heterogeneity and clustering among cells. Useful for exploring the site
#' occupancy patterns for each cell.
#'
#' @param object A SimulatePolymerase object
#' @param timePoint Optional time point to plot. If NULL, plots the final
#' position matrix.
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Plot resolution in DPI
#' @return A ggplot object showing the PCA plot of cells
#' @export
setGeneric("plotPolymerasePCA", function(
    object, timePoint = NULL, file = NULL, width = 8, height = 6, dpi = 300) {
    standardGeneric("plotPolymerasePCA")
})
#' @rdname SimulatePolymerase-class
setMethod(
    "plotPolymerasePCA", "SimulatePolymerase",
    function(object, timePoint = NULL, file = NULL, width = 8, height = 6, 
    dpi = 300) {
        if (inherits(object, "SimulatePolymerase")) {
            if (is.null(timePoint)) {
                mat <- slot(object, "finalPositionMatrix")
                plotTitle <- "PCA of Polymerase Position Matrix (Final)"
            } else {
                mat <- getPositionMatrixAtTime(object, timePoint)
                plotTitle <- sprintf(
                    "PCA of Polymerase Position Matrix (Time %.2f)", timePoint
                )
            }
        } else {
            stop("object must be a SimulatePolymerase object")
        }

        if (ncol(mat) < 2) stop("Need at least 2 cells (columns) for PCA.")
        if (nrow(mat) < 2) stop("Need at least 2 sites (rows) for PCA.")

        # Transpose: cells as rows, sites as features
        pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
        varExplained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
        
        df <- data.frame(
            PC1 = pca$x[, 1], PC2 = pca$x[, 2], Cell = seq_len(ncol(mat)), 
            PolymeraseCount = colSums(mat)
        )

        p <- ggplot(df, aes(x = PC1, y = PC2, color = PolymeraseCount)) +
            geom_point(size = 3, alpha = 0.8) +
            scale_color_gradient(
                low = "lightblue", high = "darkred", name = "Polymerase\nCount"
            ) +
            labs(
                title = plotTitle,
                x = sprintf("PC1 (%.1f%% variance)", varExplained[1]),
                y = sprintf("PC2 (%.1f%% variance)", varExplained[2])
            ) +
            theme_minimal() +
            theme(plot.title = element_text(size = 14, face = "bold", 
                hjust = 0.5), legend.position = "right")

        if (!is.null(file)) {
            ggsave(file, p, width = width, height = height, dpi = dpi)
        }

        return(p)
    }
)

validatePlotRange <- function(start, end, data) {
    if (!is.null(start) || !is.null(end)) {
        plotStart <- if (!is.null(start)) start + 1 else 1
        plotEnd <- if (!is.null(end)) end else length(data)

        if (plotStart < 1 || plotStart > length(data)) {
            stop(sprintf("start must be between 0 and %d", 
            length(data) - 1))
        }
        if (plotEnd < plotStart || plotEnd > length(data)) {
            stop(sprintf("end must be between %d and %d", plotStart, 
            length(data)))
        }

        df <- data.frame(
            position = plotStart:plotEnd,
            count = data[plotStart:plotEnd]
        )
    } else {
        df <- data.frame(position = 2:length(data), 
        count = data[2:length(data)])
    }

    return(df)
}

#' @rdname SimulatePolymerase-class
#' @title Plot Combined Cells Data
#'
#' @description
#' Plot combined cells data as a ggplot2 plot. For SimulatePolymerase objects,
#' can plot data from specific time points by calculating combined cells data
#' from position matrices.
#'
#' @param object A SimulatePolymerase-class object
#' @param start Integer, starting position for plotting (default: NULL, uses
#' full range)
#' @param end Integer, ending position for plotting (default: NULL, uses full
#' range)
#' @param timePoint Optional numeric value specifying the time point being
#' plotted (used for automatic title generation and documentation)
#' @param file Optional file path to save the plot to
#' @return A ggplot object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, 
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 100, polSize = 33,
#'     addSpace = 17, time = 1, timesToRecord = c(0.5, 1.0))
#' # Plot final combined cells data
#' plotCombinedCells(sim, file="combined_cells.html")
#' # Plot specific time point data
#' plotCombinedCells(sim, timePoint = 0.5, file="combined_cells.html")
#' # Plot with custom range and title
#' plotCombinedCells(sim, timePoint = 0.5, start = 100, end = 500,
#' file="combined_cells.html")
#' @export
setGeneric("plotCombinedCells", function(
    object, start = NULL, end = NULL, timePoint = NULL, file = NULL) {
    standardGeneric("plotCombinedCells")
})
#' @rdname SimulatePolymerase-class
setMethod(
    "plotCombinedCells", "SimulatePolymerase",
    function(object, start = NULL, end = NULL, timePoint = NULL, file = NULL) {
        if (is.null(timePoint)) {
            data <- slot(object, "combinedCellsData")
            timePoint <- slot(object, "time") / 1000
        } else {
            pos_matrix <- getPositionMatrixAtTime(object, timePoint)
            data <- rowSums(pos_matrix)
        }
        if (length(data) == 0) stop("No data available for plotting.")

        pauseSites <- slot(object, "pauseSites")
        title <- sprintf("Polymerase Occupancy at Time %.4f Minutes", timePoint)
        df <- validatePlotRange(start, end, data)

        df <- subset(df, count > 0)
        p <- ggplot(df, aes(x=position, y=count)) +
        geom_segment( aes(x=position, xend=position, y=0, yend=count),
        color="lightgrey", size = 0.5) +
        geom_point( color="steelblue", size=3) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
        ) +
        labs(title = title, x = "Gene Site/Position", 
        y = "Number of Polymerases") +
        scale_y_continuous(expand = expansion(mult = c(0, .01))) +
        scale_x_continuous(expand = expansion(mult = c(.01, .05)))


        pmean <- mean(pauseSites)
        psd <- sd(pauseSites)
        if (pmean >= min(df$position) && pmean <= max(df$position)) {
            p <- p + geom_vline(
                xintercept = pmean, color = "red",
                linetype = "dashed", size = 1
            )
        }

        return(p)
    }
)