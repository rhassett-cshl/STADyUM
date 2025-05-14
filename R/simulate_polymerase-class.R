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
    if (ncol(positionMatrix(object)) != slot(object, "cellNum") ||
        nrow(positionMatrix(object)) != slot(object, "geneLen") + 1) {
        errors <- c(errors, "positionMatrix dimensions do not match cellNum
        and geneLen")
    }

    if (length(errors) == 0) TRUE else errors
}

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
    zeta_data <- read.csv(zetaVec, header = FALSE, 
                        colClasses = "numeric")
    
    # Convert to numeric vector and ensure correct length
    zeta_vec <- as.numeric(zeta_data[[1]])
    
    return(zeta_vec)
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
        zetaMax = "numeric", zetaVec="character", cellNum = "integer", 
        polSize = "integer", addSpace = "integer", time = "numeric", 
        deltaT = "numeric", stepsToRecord = "integer", pauseSites = "numeric",
        probabilityVector = "numeric", combinedCellsData = "integer",
        positionMatrix = "matrix", readCounts = "ANY",
        avgReadDensity = "ANY"
    ),
    validity = simulatePolymeraseValid
)


#' @name simulatePolymerase
#' @rdname simulatePolymerase-class
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
#' @param stepsToRecord an integer value for the number of steps to record in
#' position matrix.
#' @return a \code{simulatePolymerase} object
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Print the object
#' print(sim)
#' @export
simulatePolymerase <- function(
    k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
    alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
    zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
    stepsToRecord=1) {
    
    validateSimulatePolymeraseParams(
        k, ksd, kMin, kMax, geneLen, alpha, beta, zeta, zetaSd, zetaMin,
        zetaMax, cellNum, polSize, addSpace, time, stepsToRecord
    )
    
    zeta_vec <- validateAndLoadZetaVec(zetaVec, geneLen)
    
    result <- .Call("_STADyUM_simulate_polymerase_cpp",
        k, ksd, kMin, kMax, geneLen, alpha, beta, zeta, zetaSd,
        zetaMin, zetaMax, cellNum, polSize, addSpace, time, stepsToRecord,
        zeta_vec, PACKAGE = "STADyUM"
    )
    
    # Create and return a simulatePolymerase object
    obj <- new("simulatePolymerase",
        probabilityVector = result$probabilityVector,
        pauseSites = result$pauseSites,
        combinedCellsData = result$combinedCellsData,
        positionMatrix = result$positionMatrix,
        k = as.integer(k), kMin = as.integer(kMin), kMax = as.integer(kMax),
        ksd = ksd, geneLen = as.integer(geneLen), alpha = alpha, beta = beta,
        zeta = zeta, zetaSd = zetaSd, zetaMax = zetaMax, zetaMin = zetaMin,
        zetaVec = if (is.null(zetaVec)) "" else zetaVec,
        cellNum = as.integer(cellNum), polSize = as.integer(polSize),
        addSpace = as.integer(addSpace), time = time,
        stepsToRecord = as.integer(stepsToRecord), readCounts = NULL,
        avgReadDensity = NULL
    )
    
    validObject(obj)
    
    return(obj)
}

#' Sample read counts from a simulatePolymerase object. To match our simulated
#' read counts to reality, we need to compute a scaling factor lambda. One way
#' of doing it is computing the read density based on real experiments. For
#' example, we have computed the read density within gene body in 
#' _Dukler et al_ (2017) for genes with median expression (i.e., 0.0489). 
#' If we assume the read counts following a Poisson distribution, we can then
#' sample the read counts with mean equals to the RNAP frequency multiplied by
#' lambda.
#' @param object A simulatePolymerase object
#' @param readDensity A numeric value for the read density within gene body in
#' _Dukler et al._ (2017) for genes with median expression (i.e., 0.0489).
#' @return The read count per nucleotide value
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Sample read counts per nucleotide
#' readCounts <- sampleReadCountsPerNucleotide(sim)
#' # Print the read counts per nucleotide
#' print(readCounts)
#' @export
setGeneric("sampleReadCountsPerNucleotide", function(
    object,
    readDensity = 0.0489) {
    standardGeneric("sampleReadCountsPerNucleotide")
})
setMethod("sampleReadCountsPerNucleotide", "simulatePolymerase", function(
    object, readDensity = 0.0489) {
    cellNum <- slot(object, "cellNum")
    kMax <- slot(object, "kMax")
    totalRnap <- slot(object, "combinedCellsData")

    N <- length(totalRnap)
    L <- N - kMax

    # TODO: handle case if lambda is INF because sum is 0
    lambda <- readDensity / (sum(totalRnap[(kMax + 1):N]) / (L * cellNum))

    rc <- rpois(N, totalRnap / cellNum * lambda)

    gbRc <- sum(rc[(kMax + 1):N])

    ## read count per nucleotide
    rcPerNt <- gbRc / L

    slot(object, "readCounts") <- rcPerNt

    return(rcPerNt)
})

#' Sample read counts from a simulatePolymerase object
#' @param object A simulatePolymerase object
#' @param readDensity A numeric value for the read density within gene body in
#' _Dukler et al._ (2017) for genes with median expression (i.e., 0.0489).
#' @return The average read density within gene body
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1) 
#' # Sample average read density within gene body
#' avgReadDensity <- sampleGeneBodyAvgReadDensity(sim)
#' # Print the average read density within gene body
#' print(avgReadDensity)
#' @export
setGeneric("sampleGeneBodyAvgReadDensity", function(
    object,
    readDensity = 0.0489) {
    standardGeneric("sampleGeneBodyAvgReadDensity")
})
setMethod(
    "sampleGeneBodyAvgReadDensity", "simulatePolymerase",
    function(object, readDensity = 0.0489) {
        cellNum <- slot(object, "cellNum")
        kMax <- slot(object, "kMax")
        totalRnap <- slot(object, "combinedCellsData")

        N <- length(totalRnap)
        L <- N - kMax

        simAvgReadDensity <- sum(totalRnap[(kMax + 1):N]) / L

        slot(object, "avgReadDensity") <- simAvgReadDensity

        return(simAvgReadDensity)
    }
)


#' @rdname simulatePolymerase-class
#' @param object a \code{simulatePolymerase} object
#' @return a \code{simulatePolymerase} object
#' @export
setGeneric("simulateReadCounts", function(object) {
    standardGeneric("simulateReadCounts")
})
#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Simulate read counts
#' sim <- simulateReadCounts(sim)
#' # Print the read counts
#' print(sim)
#' @export
setMethod("simulateReadCounts", "simulatePolymerase", function(object) {
    cellNum <- slot(object, "cellNum")
    kMax <- slot(object, "kMax")
    totalRnap <- slot(object, "combinedCellsData")

    # Calculate read counts per nucleotide
    rcPerNt <- totalRnap / cellNum

    # Store read counts in object
    slot(object, "readCounts") <- rcPerNt

    return(object)
})

#' @rdname simulatePolymerase-class
#' @export
setGeneric("simulateAvgReadDensity", function(object) {
    standardGeneric("simulateAvgReadDensity")
})
#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Simulate average read density
#' sim <- simulateAvgReadDensity(sim)
#' # Print the average read density
#' print(sim)
#' @export
setMethod("simulateAvgReadDensity", "simulatePolymerase", function(object) {
    cellNum <- slot(object, "cellNum")
    kMax <- slot(object, "kMax")
    totalRnap <- slot(object, "combinedCellsData")

    # Calculate average read density
    simAvgReadDensity <- totalRnap / (cellNum * kMax)

    # Store average read density in object
    slot(object, "avgReadDensity") <- simAvgReadDensity

    return(object)
})

#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Show the object
#' show(sim)
#' @export
setMethod("show", "simulatePolymerase", function(object) {
    cat("A simulatePolymerase object with:\n")
    cat("  - k =", slot(object, "k"), "\n")
    cat("  - ksd =", slot(object, "ksd"), "\n")
    cat("  - kMin =", slot(object, "kMin"), "\n")
    cat("  - kMax =", slot(object, "kMax"), "\n")
    cat("  - geneLen =", slot(object, "geneLen"), "\n")
    cat("  - alpha =", slot(object, "alpha"), "\n")
    cat("  - beta =", slot(object, "beta"), "\n")
    cat("  - zeta =", slot(object, "zeta"), "\n")
    cat("  - zetaSd =", slot(object, "zetaSd"), "\n")
    cat("  - zetaMin =", slot(object, "zetaMin"), "\n")
    cat("  - zetaMax =", slot(object, "zetaMax"), "\n")
    cat("  - cellNum =", slot(object, "cellNum"), "\n")
    cat("  - polSize =", slot(object, "polSize"), "\n")
    cat("  - addSpace =", slot(object, "addSpace"), "\n")
    cat("  - time =", slot(object, "time"), "\n")
    cat("  - stepsToRecord =", slot(object, "stepsToRecord"), "\n")
})

# Plotting methods
#' Plot polymerase distribution
#' @param object A simulatePolymerase object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the distribution of polymerases across the
#' gene
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Plot polymerase distribution
#' plotPolymeraseDistribution(sim)
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
            position = 0:slot(object, "geneLen"),
            count = slot(object, "combinedCellsData")
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
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1) 
#' # Plot pause site distribution
#' plotPauseSites(sim)
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
        cell = seq_len(slot(object, "cellNum")),
        pauseSite = slot(object, "pauseSites")
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
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Plot transition probabilities
#' plotTransitionProbabilities(sim)
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
            position = 0:slot(object, "geneLen"),
            probability = slot(object, "probabilityVector")
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
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Plot position matrix
#' plotPositionMatrix(sim)
#' 
#' @export
setGeneric("plotPositionMatrix", function(
    object, file = NULL, width = 8,
    height = 6) {
    standardGeneric("plotPositionMatrix")
})
setMethod("plotPositionMatrix", "simulatePolymerase", function(
    object,
    file = NULL, width = 8, height = 6) {
    df <- reshape2::melt(slot(object, "positionMatrix"))
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
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Save data frames
#' saveDataFrames(sim)
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
    write.csv(pauseSites(object),
        file.path(dir, "pause_sites.csv"),
        row.names = FALSE
    )
    write.csv(probabilityVector(object), file.path(
        dir,
        "transition_probabilities.csv"
    ), row.names = FALSE)
    write.csv(combinedCellsData(object), file.path(
        dir,
        "polymerase_counts.csv"
    ), row.names = FALSE)
    write.csv(positionMatrix(object), file.path(dir, "position_matrix.csv"),
        row.names = FALSE
    )

    # Save parameters
    write.csv(as.data.frame(parameters(object)), file.path(
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
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Save plots
#' savePlots(sim)
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

#' @rdname simulatePolymerase-class
#' @export
setGeneric("plotAvgReadDensity", function(object) {
    standardGeneric("plotAvgReadDensity")
})
#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Plot average read density
#' plotAvgReadDensity(sim)
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

# Accessor methods
#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Get pause sites
#' pauseSites <- pauseSites(sim)
#' # Print the pause sites
#' print(pauseSites)
#' @export
setGeneric("pauseSites", function(object) standardGeneric("pauseSites"))
setMethod("pauseSites", "simulatePolymerase", function(object) {
    slot(object, "pauseSites")
})

#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Get probability vector
#' probabilityVector <- probabilityVector(sim)
#' # Print the probability vector
#' print(probabilityVector)
#' @export
setGeneric("probabilityVector", function(object) {
    standardGeneric("probabilityVector")
})
setMethod("probabilityVector", "simulatePolymerase", function(object) {
    slot(object, "probabilityVector")
})

#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Get combined cells data
#' combinedCellsData <- combinedCellsData(sim)
#' # Print the combined cells data
#' print(combinedCellsData)
#' @export
setGeneric("combinedCellsData", function(object) {
    standardGeneric("combinedCellsData")
})
setMethod("combinedCellsData", "simulatePolymerase", function(object) {
    slot(object, "combinedCellsData")
})

#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Get position matrix
#' positionMatrix <- positionMatrix(sim)
#' # Print the position matrix
#' print(positionMatrix)
#' @export
setGeneric("positionMatrix", function(object) {
    standardGeneric("positionMatrix")
})
setMethod("positionMatrix", "simulatePolymerase", function(object) {
    slot(object, "positionMatrix")
})

#' Get all simulation parameters
#' @param object A simulatePolymerase object
#' @return A list containing all simulation parameters including k, ksd, kMin,
#' kMax, geneLen, alpha, beta, zeta, zetaSd, zetaMin, zetaMax, zetaVec, cellNum,
#' polSize, addSpace, time, and stepsToRecord
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Get parameters
#' parameters <- parameters(sim)
#' # Print the parameters
#' print(parameters)
#' @export
setGeneric("parameters", function(object) standardGeneric("parameters"))

#' @rdname parameters
#' @export
setMethod("parameters", "simulatePolymerase", function(object) {
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
        stepsToRecord = slot(object, "stepsToRecord")
    )
})

#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Get read counts    
#' readCounts <- readCounts(sim)
#' # Print the read counts
#' print(readCounts)
#' @export
setGeneric("readCounts", function(object) standardGeneric("readCounts"))
setMethod("readCounts", "simulatePolymerase", function(object) {
    slot(object, "readCounts")
})

#' @rdname simulatePolymerase-class
#' @examples
#' # Create a simulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Get average read density
#' avgReadDensity <- avgReadDensity(sim)
#' # Print the average read density
#' print(avgReadDensity)
#' @export
setGeneric("avgReadDensity", function(object) {
    standardGeneric("avgReadDensity")
})
setMethod("avgReadDensity", "simulatePolymerase", function(object) {
    slot(object, "avgReadDensity")
})
