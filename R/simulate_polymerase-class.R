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
    zeta_data <- read_csv(zetaVec, 
                            col_names = FALSE,
                            col_types = "d",
                            show_col_types = FALSE)
    
    # Convert to numeric vector and ensure correct length
    zeta_vec <- as.numeric(zeta_data[[1]])
    
    return(zeta_vec)
}

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
#' @slot deltaT a numeric value for the time step size in the simulation.
#' @slot stepsToRecord an integer value for the number of steps to record in
#' position matrix.
#' @slot pauseSites a numeric vector of pause sites
#' @slot probabilityVector a numeric vector representing the probability that
#' the polymerase move forward or not at each site
#' @slot combinedCellsData an integer vector representing the total number of
#' RNAPs at each site across all cells
#' @slot positionMatrix a matrix of position of polymerase
#' @slot readCounts a numeric vector for read counts per nucleotide
#'
#' @name SimulatePolymerase-class
#' @rdname SimulatePolymerase-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @importFrom methods slot new is slot<- validObject
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs
#' @importFrom ggplot2 geom_tile scale_fill_gradient ggsave geom_histogram
#' @importFrom reshape2 melt
#' @importFrom readr read_csv
#' @exportClass SimulatePolymerase
methods::setClass("SimulatePolymerase",
    slots = c(
        k = "integer", ksd = "numeric", kMin = "integer", kMax = "integer",
        geneLen = "integer", alpha = "numeric", beta = "numeric",
        zeta = "numeric", zetaSd = "numeric", zetaMin = "numeric",
        zetaMax = "numeric", zetaVec="character", cellNum = "integer", 
        polSize = "integer", addSpace = "integer", time = "numeric", 
        deltaT = "numeric", stepsToRecord = "integer", pauseSites = "numeric",
        probabilityVector = "numeric", combinedCellsData = "integer",
        positionMatrix = "matrix", readCounts = "ANY"
    ),
    validity = simulatePolymeraseValid
)

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
#' @param stepsToRecord an integer value for the number of steps to record in
#' position matrix.
#' @return a \code{SimulatePolymerase} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
    
    # Create and return a SimulatePolymerase object
    obj <- new("SimulatePolymerase",
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
        stepsToRecord = as.integer(stepsToRecord), readCounts = NULL)

    sampleReadCountsPerNucleotide(obj)
    
    validObject(obj)
    
    return(obj)
}

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
#' @return The read count per nucleotide value
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("sampleReadCountsPerNucleotide", "SimulatePolymerase", function(
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

#' @rdname SimulatePolymerase-class
#' @title Show Method for SimulatePolymerase Object
#'
#' @description
#' Show method for SimulatePolymerase object in human readable format
#' including summary statistics
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Show the object
#' show(sim)
#' @export
setMethod("show", "SimulatePolymerase", function(object) {
    cat("A SimulatePolymerase object with:\n")
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

#' @rdname SimulatePolymerase-class
#' @title Plot Polymerase Distribution
#'
#' @description
#' Plot the distribution of polymerases across the gene.
#'
#' @param object A SimulatePolymerase-class object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the distribution of polymerases across the
#' gene
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
    "plotPolymeraseDistribution", "SimulatePolymerase",
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

#' @rdname SimulatePolymerase-class
#' @title Plot Pause Site Distribution
#'
#' @description
#' Plot the distribution of pause sites across the gene.
#'
#' @param object A SimulatePolymerase-class object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the distribution of pause sites
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("plotPauseSites", "SimulatePolymerase", function(
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

#' @rdname SimulatePolymerase-class
#' @title Plot Transition Probabilities
#'
#' @description
#' Plot the transition probabilities across the gene.
#'
#' @param object A SimulatePolymerase-class object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the transition probabilities across the gene
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
    "plotTransitionProbabilities", "SimulatePolymerase",
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

#' @rdname SimulatePolymerase-class
#' @title Plot Position Matrix Heatmap
#'
#' @description
#' Plot the position matrix as a heatmap.
#'
#' @param object A SimulatePolymerase-class object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the position matrix as a heatmap
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("plotPositionMatrix", "SimulatePolymerase", function(
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

#' @rdname SimulatePolymerase-class
#' @title Save all plots to files
#'
#' @description
#' Save all plots to files.
#'
#' @param object A SimulatePolymerase-class object
#' @param dir Directory to save the plots (default: "results")
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return Outputs directory where plots were saved to pdf files
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("savePlots", "SimulatePolymerase", function(
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

# Accessor methods
#' @rdname SimulatePolymerase-class
#' @title Accessor for Pause Sites
#'
#' @description
#' Accessor for the pause sites numeric vector from a SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("pauseSites", "SimulatePolymerase", function(object) {
    slot(object, "pauseSites")
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Probability Vector 
#'
#' @description
#' Accessor for the probability numeric vector from a SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("probabilityVector", "SimulatePolymerase", function(object) {
    slot(object, "probabilityVector")
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Combined Cells Data
#'
#' @description
#' Accessor for the combined cells data numeric vector from a
#' SimulatePolymerase object.
#'
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("combinedCellsData", "SimulatePolymerase", function(object) {
    slot(object, "combinedCellsData")
})

#' @rdname SimulatePolymerase-class
#' @title Accessor for Position Matrix
#'
#' @description
#' Accessor for the position matrix from a SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("positionMatrix", "SimulatePolymerase", function(object) {
    slot(object, "positionMatrix")
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
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
        stepsToRecord = slot(object, "stepsToRecord")
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
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
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
setMethod("readCounts", "SimulatePolymerase", function(object) {
    slot(object, "readCounts")
})