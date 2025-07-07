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
#' @slot timePointsToRecord a numeric vector of specific time points to record
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
#' @importFrom reshape2 melt
#' @importFrom readr read_csv
#' @importFrom plotly plot_ly layout add_segments add_annotations
#' @exportClass SimulatePolymerase
methods::setClass("SimulatePolymerase",
    slots = c(
        k = "integer", ksd = "numeric", kMin = "integer", kMax = "integer",
        geneLen = "integer", alpha = "numeric", beta = "numeric",
        zeta = "numeric", zetaSd = "numeric", zetaMin = "numeric",
        zetaMax = "numeric", zetaVec="character", cellNum = "integer", 
        polSize = "integer", addSpace = "integer", time = "numeric", 
        timePointsToRecord = "ANY", pauseSites = "numeric", 
        siteProbabilities = "matrix", combinedCellsData = "integer", positionMatrices = "list", finalPositionMatrix = "matrix", 
        readCounts = "ANY"
))

validateSimulatePolymeraseParams <- function(
    k, ksd, kMin, kMax, geneLen,
    alpha, beta, zeta, zetaSd, zetaMin, zetaMax, cellNum, polSize,
    addSpace, time, timePointsToRecord) {
    errors <- character()

    # Check parameter ranges
    if(!is.numeric(k) || k <= 0 || k %% 1 != 0 || k < kMin || k > kMax) {
        errors <- c(errors, "k must be a positive integer between kMin and kMax")
    }
    if(!is.numeric(ksd) || ksd < 0 || ksd %% 1 != 0) {
        errors <- c(errors, "ksd must be a non-negative integer")
    }
    if(!is.numeric(kMin) || kMin <= 0 || kMin %% 1 != 0) {
        errors <- c(errors, "kMin must be a positive integer")
    }
    if(!is.numeric(kMax) || kMax <= 0 || kMax %% 1 != 0) {
        errors <- c(errors, "kMax must be a positive integer")
    }
    if (kMin >= kMax) {
        errors <- c(errors, "kMin must be less than or equal to kMax")
    }
    if (k <= kMin || k >= kMax) {
        errors <- c(errors, "k must be between kMin and kMax (inclusive)")
    }
    if(!is.numeric(alpha) || alpha <= 0) {
        errors <- c(errors, "alpha must be a positive number")
    }
    if(!is.numeric(beta) || beta <= 0) {
        errors <- c(errors, "beta must be a positive number")
    }
    if(!is.numeric(zeta) || zeta <= 0 || zeta < zetaMin || zeta > zetaMax) {
        errors <- c(errors, "zeta must be a positive number between zetaMin and zetaMax")
    }
    if(!is.numeric(zetaSd) || zetaSd < 0) {
        errors <- c(errors, "zetaSd must be a non-negative number")
    }
    if(!is.numeric(zetaMin) || zetaMin <= 0 || zetaMin > zetaMax) {
        errors <- c(errors, "zetaMin must be a positive number less than or equal to zetaMax")
    }
    if(!is.numeric(zetaMax) || zetaMax <= 0) {
        errors <- c(errors, "zetaMax must be a positive number")
    }
    if(!is.numeric(geneLen) || geneLen <= 0 || geneLen %% 1 != 0 || geneLen < kMax) {
        errors <- c(errors, "geneLen must be a positive integer greater than kMax")
    }
    if(!is.numeric(cellNum) || cellNum <= 0 || cellNum %% 1 != 0) {
        errors <- c(errors, "cellNum must be a positive integer")
    }
    if(!is.numeric(polSize) || polSize <= 0 || polSize %% 1 != 0) {
        errors <- c(errors, "polSize must be a positive integer")
    }
    if(!is.numeric(addSpace) || addSpace < 0 || addSpace %% 1 != 0) {
        errors <- c(errors, "addSpace must be a non-negative integer")
    }
    if (!is.numeric(time) || time < 1e-4) {
        errors <- c(errors, "time must be a positive number greater than 1e-4")
    }
    
    # Memory usage checks for position matrices
    if (!is.null(timePointsToRecord)) {
        if (any(timePointsToRecord < 0)) {
            errors <- c(errors, "timePointsToRecord must contain non-negative values")
        }
        if (any(timePointsToRecord > time)) {
            errors <- c(errors, "timePointsToRecord must not exceed simulation time")
        }
        
        # Check memory usage for position matrices
        numTimePoints <- length(timePointsToRecord)
        totalSites <- geneLen + 1
        estimatedMemoryMB <- (numTimePoints * totalSites * cellNum * 4) / (1024 * 1024)
        
        # Warn if memory usage is high (> 1GB)
        if (estimatedMemoryMB > 1024) {
            warning(sprintf("Large memory usage estimated: %.1f MB for position matrices. Consider reducing timePointsToRecord or cellNum.", estimatedMemoryMB))
        }
        
        # Hard limit at 4GB to prevent system crashes
        if (estimatedMemoryMB > 4096) {
            errors <- c(errors, sprintf("Memory usage too high: %.1f MB. Reduce timePointsToRecord or cellNum. Maximum allowed: 4GB", estimatedMemoryMB))
        }
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
    
    # Validate that all values are greater than 0 and non-negative
    if (any(is.na(zeta_vec)) || any(zeta_vec <= 0)) {
        stop("zetaVec file must contain only positive values")
    }
    
    # Check that the zeta_vec length matches geneLen
    if (length(zeta_vec) != geneLen) {
        stop(sprintf("zetaVec file length (%d) does not match geneLen (%d)", 
                    length(zeta_vec), geneLen))
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
#' @param timePointsToRecord a numeric vector of specific time points to record
#' position matrices for, or NULL to record no extra position matrices. Final
#' position matrix is always recorded. Default is NULL.
#' @return a \code{SimulatePolymerase} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=NULL)
#' # Create a SimulatePolymerase object with specific time points recorded
#' sim2 <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=c(0.5, 1.0))
#' @export
simulatePolymerase <- function(
    k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
    alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
    zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
    timePointsToRecord=NULL) {
    
    validateSimulatePolymeraseParams(
        k, ksd, kMin, kMax, geneLen, alpha, beta, zeta, zetaSd, zetaMin,
        zetaMax, cellNum, polSize, addSpace, time, timePointsToRecord
    )
    
    zeta_vec <- validateAndLoadZetaVec(zetaVec, geneLen)
    
    result <- .Call("_STADyUM_simulate_polymerase_cpp",
        k, ksd, kMin, kMax, geneLen, alpha, beta, zeta, zetaSd,
        zetaMin, zetaMax, cellNum, polSize, addSpace, time, timePointsToRecord,
        zeta_vec, PACKAGE = "STADyUM"
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
        timePointsToRecord = timePointsToRecord, readCounts = NULL)

    obj <- sampleReadCounts(obj)
        
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
#' @return The SimulatePolymerase object with read counts sampled given the
#' readDensity parameter used
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=NULL)
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
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=c(0.5, 1.0))
#' # Show the object
#' show(sim)
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
            position = 0:(length(combinedData) - 1),
            count = combinedData
        )
        topSites <- siteData[order(siteData$count, decreasing = TRUE), ]
        top10 <- head(topSites, 10)
        
        cat("  - top 10 most occupied sites across all cells:\n")
        for (i in 1:nrow(top10)) {
            cat(sprintf("    position %4d: %d polymerases\n", 
                       top10$position[i], top10$count[i]))
        }
    }
    readCounts <- slot(object, "readCounts")
    kMax <- slot(object, "kMax")
    N <- length(combinedData)
    L <- N - kMax

    gbRc <- sum(readCounts[(kMax + 1):N])
    gbAvgRc <- gbRc / L
    pauseRc <- sum(readCounts[1:kMax])
    pauseAvgRc <- pauseRc / kMax

    cat("  - gene body average read counts =", round(gbAvgRc, 2), "\n")
    cat("  - pause region average read counts =", round(pauseAvgRc, 2), "\n")
    
    # Calculate fraction of nucleotides with zero reads to indicate sparsity
    zeroReadsFraction <- sum(readCounts == 0) / length(readCounts)
    cat("  - percent of nucleotides with zero reads =", 
        round(zeroReadsFraction * 100, 1), "%\n")

    cat("\nTo access the full simulation parameters, use: parameters(object)\n")
    cat("\nTo access the all sampled read counts, use: readCounts(object)\n")
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
#'     timePointsToRecord=c(0.5, 1.0))
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
#'     timePointsToRecord=NULL) 
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
  
  # Calculate mean and standard deviation
  pause_mean <- mean(df$pauseSite)
  pause_sd <- sd(df$pauseSite)
  
  p <- ggplot(df, aes(x = pauseSite)) +
    geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue", 
                   color = "black") +
    # Add vertical line for mean
    geom_vline(xintercept = pause_mean, color = "red", 
               linetype = "dashed", size = 1) +
    # Add vertical lines for mean ± 1 SD
    geom_vline(xintercept = pause_mean + pause_sd, color = "orange", 
               linetype = "dotted", size = 0.8) +
    geom_vline(xintercept = pause_mean - pause_sd, color = "orange", 
               linetype = "dotted", size = 0.8) +
    theme_minimal() +
    labs(
      title = "Distribution of Pause Sites",
      subtitle = sprintf("Mean: %.1f, SD: %.1f", pause_mean, pause_sd),
      x = "Pause Site Position",
      y = "Count"
    )
  
  if (!is.null(file)) {
    ggsave(file, p, width = width, height = height)
  }
  
  return(p)
})



#' @rdname SimulatePolymerase-class
#' @title Plot Combined Cells Data (Interactive Plotly)
#'
#' @description
#' Plot the combined cells data as an interactive plotly visualization, 
#' excluding the first site. This provides zoom, pan, and hover capabilities 
#' for exploring polymerase occupancy patterns. Red dashed line represents mean pause site position
#' across all cells calculated from pause site vector
#'
#' @param object A SimulatePolymerase-class object
#' @param start Integer, starting position for plotting (default: NULL, uses excludeFirstSite logic)
#' @param end Integer, ending position for plotting (default: NULL, uses full range)
#' @return A plotly object showing interactive polymerase occupancy
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=NULL)
#' # Plot interactive
#' plotCombinedCells(sim)
#' @export
setGeneric("plotCombinedCells", function(
    object, start = NULL, end = NULL) {
  standardGeneric("plotCombinedCells")
})
setMethod(
  "plotCombinedCells", "SimulatePolymerase",
  function(object, start = NULL, end = NULL) {
    data <- slot(object, "combinedCellsData")
    if (length(data) == 0) {
      stop("combinedCellsData is empty. Run the simulation first.")
    }
    
    if (!is.null(start) || !is.null(end)) {
      # Use explicit start/end if provided
      plot_start <- if (!is.null(start)) start + 1 else 2
      plot_end <- if (!is.null(end)) end else length(data)
      
      # Validate range
      if (plot_start < 1 || plot_start > length(data) - 1) {
        stop(sprintf("start must be between 1 and %d", length(data) - 1))
      }
      if (plot_end < plot_start || plot_end > length(data)) {
        stop(sprintf("end must be between %d and %d", plot_start, length(data)))
      }
      
      df <- data.frame(
        position = plot_start:plot_end,
        count = data[plot_start:plot_end]
      )
    }
    
    p <- ggplot(df, aes(x = position, y = count)) +
      geom_line(color = "steelblue", size = 1) +
      geom_point(color = "darkblue", size = 0.8, alpha = 0.7) +
      theme_minimal() +
      labs(
        title = "Polymerase Occupancy (Interactive)",
        x = "Position",
        y = "Number of Polymerases"
      ) +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")
      )
    
    pause_sites <- slot(object, "pauseSites")
    if (length(pause_sites) > 0) {
      # Calculate pause site statistics
      pause_mean <- mean(pause_sites)
      pause_sd <- sd(pause_sites)
      
      # Only add pause site annotations if they fall within the plotted range
      if (pause_mean >= min(df$position) && pause_mean <= max(df$position)) {
        p <- p + 
          # Add vertical line for mean pause site
          geom_vline(xintercept = pause_mean, 
                     color = "red", linetype = "dashed", 
                     size = 1) +
          # Add pause site label
          annotate("text", x = pause_mean, y = Inf,
                   label = sprintf("Pause Site\nMean: %.0f ± %.0f", 
                                   pause_mean, pause_sd),
                   vjust = 2, hjust = 0.5, color = "red",
                   fontface = "bold", size = 3)
      }
    }
    
    # Convert to plotly
    plotly_p <- ggplotly(p, tooltip = c("x", "y")) %>%
      layout(
        xaxis = list(title = "Position"),
        yaxis = list(title = "Number of Polymerases"),
        hovermode = "x unified"
      ) %>%
      config(displayModeBar = TRUE, 
             modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d"))
    
    return(plotly_p)
  }
)


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
#'     timePointsToRecord=NULL)
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
#' @title Accessor for Probability Matrix 
#'
#' @description
#' Accessor for the probability numeric matrix from a SimulatePolymerase object.
#' The matrix has dimensions (sites x cells) where each element represents the
#' transition probability for a specific site and cell.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=NULL)
#' # Get probability matrix
#' siteProbabilities <- siteProbabilities(sim)
#' # Print the probability matrix dimensions
#' print(dim(siteProbabilities))
#' @export
setGeneric("siteProbabilities", function(object) {
    standardGeneric("siteProbabilities")
})
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
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=NULL)
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
#' @title Accessor for Position Matrices
#'
#' @description
#' Accessor for the position matrices from a SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=c(0.5, 1.0))
#' # Get position matrices
#' positionMatrices <- positionMatrices(sim)
#' # Print the position matrices
#' print(positionMatrices)
#' @export
setGeneric("positionMatrices", function(object) {
    standardGeneric("positionMatrices")
})
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
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=c(0.5, 1.0))
#' # Get final position matrix
#' finalMatrix <- finalPositionMatrix(sim)
#' # Print the final position matrix
#' print(finalMatrix)
#' @export
setGeneric("finalPositionMatrix", function(object) {
    standardGeneric("finalPositionMatrix")
})
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
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=c(0.5, 1.0))
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
        timePointsToRecord = if (is.null(slot(object, "timePointsToRecord"))) NULL else if (length(slot(object, "timePointsToRecord")) == 0) NULL else slot(object, "timePointsToRecord")
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
#'     timePointsToRecord=NULL)
#' # Get read counts    
#' readCounts <- readCounts(sim)
#' # Print the read counts
#' print(readCounts)
#' @export
setGeneric("readCounts", function(object) standardGeneric("readCounts"))
setMethod("readCounts", "SimulatePolymerase", function(object) {
    slot(object, "readCounts")
})

#' @rdname SimulatePolymerase-class
#' @title Get Position Matrix for Specific Time Point
#'
#' @description
#' Get the position matrix for a specific time point from a SimulatePolymerase object.
#'
#' @param object a \code{SimulatePolymerase-class} object
#' @param timePoint a numeric value specifying the time point to retrieve
#' @return A 2D matrix of polymerase positions for the specified time point
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=c(0.5, 1.0))
#' # Get position matrix for time 0.5
#' posMatrix <- getPositionMatrixAtTime(sim, 0.5)
#' # Print the position matrix
#' print(posMatrix)
#' @export
setGeneric("getPositionMatrixAtTime", function(object, timePoint) {
    standardGeneric("getPositionMatrixAtTime")
})
setMethod("getPositionMatrixAtTime", "SimulatePolymerase", function(object, timePoint) {
    matrices <- slot(object, "positionMatrices")
    if (length(matrices) == 0) {
        stop("No extra position matrices were recorded. Use finalPositionMatrix() to get the final state.")
    }
    
    # Extract available time points as numbers
    available_times <- as.numeric(gsub("t_", "", names(matrices)))
    
    # Find the closest time point within a small tolerance
    tolerance <- 1e-10  # Very small tolerance for floating point comparison
    closest_idx <- which(abs(available_times - timePoint) < tolerance)
    
    if (length(closest_idx) > 0) {
        # Use the first match if multiple (shouldn't happen with proper tolerance)
        time_name <- names(matrices)[closest_idx[1]]
        return(matrices[[time_name]])
    } else {
        stop(sprintf("Time point %f not found. Available time points: %s", 
                    timePoint, paste(available_times, collapse=", ")))
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
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=NULL)
#' # Get available time points
#' availableTimes <- getAvailableTimePoints(sim)
#' # Print the available time points
#' print(availableTimes)
#' @export
setGeneric("getAvailableTimePoints", function(object) {
    standardGeneric("getAvailableTimePoints")
})
setMethod("getAvailableTimePoints", "SimulatePolymerase", function(object) {
    matrices <- slot(object, "positionMatrices")
    if (length(matrices) == 0) {
        return(numeric(0))
    }
    available_times <- as.numeric(gsub("t_", "", names(matrices)))
    return(sort(available_times))
})

#' @rdname SimulatePolymerase-class
#' @title Plot Position Matrix
#'
#' @description
#' Plot the position matrix as a heatmap showing polymerase positions across
#' sites and cells.
#'
#' @param object A SimulatePolymerase-class object
#' @param timePoint Optional time point to plot (if NULL, plots the first available time point)
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
#'     timePointsToRecord=c(0.5, 1.0))
#' # Plot position matrix
#' plotPositionMatrix(sim)
#' @export
setGeneric("plotPositionMatrix", function(
    object, timePoint = NULL, file = NULL, width = 8, height = 6) {
    standardGeneric("plotPositionMatrix")
})
setMethod(
    "plotPositionMatrix", "SimulatePolymerase",
    function(object, timePoint = NULL, file = NULL, width = 8, height = 6) {
        matrices <- slot(object, "positionMatrices")
        
        if (length(matrices) == 0) {
            stop("No extra position matrices were recorded. Use finalPositionMatrix() to access the final state.")
        }
        
        if (is.null(timePoint)) {
            # Use the first available time point
            timePoint <- as.numeric(gsub("t_", "", names(matrices)[1]))
        }
        
        pos_matrix <- getPositionMatrixAtTime(object, timePoint)
        
        # Convert to long format for plotting
        df <- reshape2::melt(pos_matrix)
        colnames(df) <- c("Site", "Cell", "Value")
        
        p <- ggplot(df, aes(x = Cell, y = Site, fill = factor(Value))) +
            geom_tile() +
            scale_fill_manual(values = c("0" = "white", "1" = "red")) +
            theme_minimal() +
            labs(
                title = sprintf("Polymerase Positions at Time %.2f", timePoint),
                x = "Cell",
                y = "Site",
                fill = "Polymerase Present"
            ) +
            theme(legend.position = "none")
        
        if (!is.null(file)) {
            ggsave(file, p, width = width, height = height)
        }
        
        return(p)
    }
)

#' @rdname SimulatePolymerase-class
#' @title Plot Polymerase Count per Cell
#'
#' @description
#' Plot the distribution of polymerase counts per cell from the final position matrix.
#' This shows how transcription activity is distributed across all cells using a histogram.
#'
#' @param object A SimulatePolymerase-class object
#' @param maxCells Maximum number of cells to display (for performance with 
#' large datasets). If NULL, shows all cells.
#' @param samplingMethod Method for sampling cells when maxCells is specified:
#' "random" (random sampling), "even" (evenly spaced), or "first" (first N cells)
#' @param bins Number of bins for the histogram (default: 30)
#' @param addStatistics Logical, whether to add summary statistics
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing histogram of polymerase counts per cell
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=NULL)
#' # Plot all cells
#' plotPolymeraseCountPerCell(sim)
#' # Plot sampled cells
#' plotPolymeraseCountPerCell(sim, maxCells=200, samplingMethod="even")
#' # Custom bins
#' plotPolymeraseCountPerCell(sim, bins=50)
#' @export
setGeneric("plotPolymeraseCountPerCell", function(
    object, maxCells = NULL, samplingMethod = "even", bins = 30, 
    addStatistics = TRUE, file = NULL, width = 8, height = 6) {
  standardGeneric("plotPolymeraseCountPerCell")
})
setMethod(
  "plotPolymeraseCountPerCell", "SimulatePolymerase",
  function(object, maxCells = NULL, samplingMethod = "even", bins = 30, 
           addStatistics = TRUE, file = NULL, width = 8, height = 6) {
    matrix <- slot(object, "finalPositionMatrix")
    
    if (nrow(matrix) == 0 || ncol(matrix) == 0) {
      stop("finalPositionMatrix is empty. Run the simulation first.")
    }
    
    # Calculate polymerase count per cell
    cell_counts <- colSums(matrix)
    total_cells <- length(cell_counts)
    
    # Sample cells if specified
    if (!is.null(maxCells) && total_cells > maxCells) {
      if (samplingMethod == "random") {
        cell_indices <- sample(1:total_cells, maxCells)
      } else if (samplingMethod == "even") {
        cell_indices <- seq(1, total_cells, length.out = maxCells)
      } else if (samplingMethod == "first") {
        cell_indices <- 1:maxCells
      } else {
        stop("samplingMethod must be 'random', 'even', or 'first'")
      }
      cell_counts <- cell_counts[cell_indices]
      cell_numbers <- cell_indices
    } else {
      cell_numbers <- 1:total_cells
    }
    
    # Create data frame for histogram
    df <- data.frame(count = cell_counts)
    
    # Create the histogram
    p <- ggplot(df, aes(x = count)) +
      geom_histogram(bins = bins, fill = "steelblue", 
                     color = "darkblue", alpha = 0.7) +
      theme_minimal() +
      labs(
        title = sprintf("Distribution of Polymerase Counts per Cell (%s)", 
                        if (!is.null(maxCells) && total_cells > maxCells) 
                          sprintf("Showing %d of %d cells", maxCells, total_cells) 
                        else sprintf("All %d cells", total_cells)),
        x = "Number of Polymerases per Cell",
        y = "Number of Cells"
      )
    
    # Add statistics if requested
    if (addStatistics) {
      mean_count <- mean(cell_counts)
      median_count <- median(cell_counts)
      sd_count <- sd(cell_counts)
      
      # Add vertical lines for mean and median
      p <- p + 
        geom_vline(xintercept = mean_count, color = "red", 
                   linetype = "dashed", size = 1) +
        geom_vline(xintercept = median_count, color = "orange", 
                   linetype = "dotted", size = 1) +
        annotate("text", x = mean_count, y = Inf, 
                 label = sprintf("Mean: %.1f", mean_count), 
                 vjust = 2, hjust = -0.1, color = "red", 
                 fontface = "bold") +
        annotate("text", x = median_count, y = Inf, 
                 label = sprintf("Median: %.1f", median_count), 
                 vjust = 4, hjust = -0.1, color = "orange", 
                 fontface = "bold") +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("SD: %.1f", sd_count), 
                 vjust = 6, hjust = 1.1, color = "black", 
                 fontface = "bold")
    }
    
    if (!is.null(file)) {
      ggsave(file, p, width = width, height = height)
    }
    
    return(p)
  }
)


#' @rdname SimulatePolymerase-class
#' @title Plot Position Matrix Heatmap (Interactive)
#'
#' @description
#' Plot position matrices as an interactive plotly heatmap showing 
#' polymerase positions across all cells at specific time points. Each cell in
#' the heatmap represents whether a polymerase is present (1) or absent (0) 
#' at a specific site in a specific cell. By default, shows the final position matrix.
#' 
#' The heatmap is interactive, allowing users to hover over cells to see the
#' exact polymerase positions. Users can also zoom in and out, pan, and click
#' on cells to get more detailed information.
#'
#' @param object A SimulatePolymerase-class object
#' @param timePoint Optional time point to plot. If NULL, plots the final
#' position matrix.
#' @param maxCells Maximum number of cells to display (for performance with 
#' large datasets). If NULL, shows all cells.
#' @param addPauseSites Logical, whether to add pause site annotations
#' @param addGeneBody Logical, whether to add gene body annotation
#' @return A plotly object showing interactive heatmap of polymerase positions
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- SimulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     timePointsToRecord=NULL)
#' # Plot final position heatmap
#' plotFinalPositionHeatmap(sim, maxCells=100)
#' # Plot specific time point
#' plotFinalPositionHeatmap(sim, timePoint=0.5, maxCells=100)
#' @export
setGeneric("plotFinalPositionHeatmap", function(
    object, timePoint = NULL, maxCells = NULL, addPauseSites = TRUE, addGeneBody = TRUE) {
  standardGeneric("plotFinalPositionHeatmap")
})
setMethod(
  "plotFinalPositionHeatmap", "SimulatePolymerase",
  function(object, timePoint = NULL, maxCells = NULL, addPauseSites = TRUE, addGeneBody = TRUE) {
    # Get the appropriate matrix based on timePoint
    if (is.null(timePoint)) {
      # Use final position matrix
      matrix <- slot(object, "finalPositionMatrix")
      plot_title <- "Final Polymerase Positions Heatmap"
    } else {
      # Use position matrix for specific time point
      matrix <- getPositionMatrixAtTime(object, timePoint)
      plot_title <- sprintf("Polymerase Positions Heatmap at Time %.2f", timePoint)
    }
    
    if (nrow(matrix) == 0 || ncol(matrix) == 0) {
      stop("Position matrix is empty. Run the simulation first.")
    }
    
    # Limit cells for visualization if specified
    if (!is.null(maxCells) && ncol(matrix) > maxCells) {
      # Sample cells evenly across the range
      cell_indices <- seq(1, ncol(matrix), length.out = maxCells)
      matrix <- matrix[, cell_indices]
    }
    
    # Create the plotly heatmap
    p <- plot_ly(
      z = matrix,
      type = "heatmap",
      colorscale = list(
        list(0, "white"),
        list(1, "red")
      ),
      showscale = TRUE,
      colorbar = list(
        title = "Polymerase",
        ticktext = c("Absent", "Present"),
        tickvals = c(0, 1)
      )
    ) %>%
      layout(
        title = list(
          text = plot_title,
          font = list(size = 16)
        ),
        xaxis = list(
          title = "Cell",
          showticklabels = FALSE,
          showgrid = FALSE
        ),
        yaxis = list(
          title = "Site",
          showgrid = FALSE
        ),
        margin = list(l = 60, r = 60, t = 80, b = 60)
      )
    
    # Add pause site annotation if requested
    if (addPauseSites) {
      pause_sites <- slot(object, "pauseSites")
      if (length(pause_sites) > 0) {
        pause_mean <- mean(pause_sites)
        pause_sd <- sd(pause_sites)
        
        # Add horizontal line for mean pause site
        p <- p %>% add_segments(
          x = 0, xend = ncol(matrix),
          y = pause_mean, yend = pause_mean,
          line = list(color = "blue", width = 2, dash = "dash"),
          showlegend = FALSE
        ) %>%
          add_annotations(
            x = ncol(matrix) * 0.5,
            y = pause_mean,
            text = sprintf("Pause Site<br>Mean: %.0f ± %.0f", pause_mean, pause_sd),
            showarrow = FALSE,
            font = list(color = "blue", size = 12),
            bgcolor = "rgba(255,255,255,0.8)",
            bordercolor = "blue",
            borderwidth = 1
          )
      }
    }
    
    # Add gene body annotation if requested
    if (addGeneBody) {
      gene_len <- slot(object, "geneLen")
      k_max <- slot(object, "kMax")
      
      # Add gene body region annotation
      p <- p %>% add_annotations(
        x = ncol(matrix) * 0.5,
        y = (k_max + 1 + gene_len) / 2,
        text = "Gene Body",
        showarrow = FALSE,
        font = list(color = "darkgreen", size = 12),
        bgcolor = "rgba(255,255,255,0.8)",
        bordercolor = "darkgreen",
        borderwidth = 1
      )
    }
    
    return(p)
  }
)