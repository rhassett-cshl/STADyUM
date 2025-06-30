#' @title Constructor for SimulationTranscriptionRates object
#'
#' @description
#' Class containing the simulated data, such as the nascent RNA sequencing reads
#' sampled from the simulated polymerase movement in the 
#' \code{SimulatePolymerase} object. It also contains the estimated average
#' read depths along the gene body and pause regions, given fixed or varied
#' pause sites, as well as the landing pad occupancy estimate. These can be
#' under a model with or without steric hindrance.
#' 
#' @slot simpol a \code{\linkS4class{SimulatePolymerase}} object
#' @slot stericHindrance a logical value to determine whether to infer
#' landing-pad occupancy or not. Defaults to FALSE.
#' @slot rates a \code{\link[tibble]{tbl_df}} containing the estimated rates
#' with columns:
#' \describe{
#'   \item{trial}{Numeric. Trial number}
#'   \item{chi}{Numeric. RNAP density along gene body}
#'   \item{betaOrg}{Numeric. RNAP density along pause region}
#'   \item{betaAdp}{Numeric. RNAP density along pause region from adapted model
#' which allows pause sites to vary across cells}
#'   \item{fkMean}{Numeric. Mean position of pause sites}
#'   \item{fkVar}{Numeric. Variance of pause site positions}
#'   \item{phi}{Numeric. Landing-pad occupancy (only if steric hindrance is
#'   enabled)}
#' }
#' @slot rnapN a list of \code{\link[GenomicRanges]{GRanges-class}} objects
#' that hold the RNAP positions
#'
#' @name SimulationTranscriptionRates-class
#' @rdname SimulationTranscriptionRates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @importClassesFrom data.table data.table
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom tibble tibble
#' @importFrom Matrix Matrix
#' @importFrom plyranges group_by_overlaps group_by summarise
#' @importFrom purrr map map_chr map_dbl map2
#' @importFrom dplyr mutate
#' @importFrom methods slot is slot<- validObject
#' @importFrom stats rpois optim rnorm var dnorm uniroot
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs ggsave
#' @exportClass SimulationTranscriptionRates
methods::setClass("SimulationTranscriptionRates",
    slots = c(
        simpol = "SimulatePolymerase",
        stericHindrance = "logical",
        rates = "tbl_df",
        rnapN = "list"
    ),
    contains = "TranscriptionRates"
)

summariseSimulationBw <- function(bw, grng, regionNames) {
    rc <- grng %>%
        plyranges::group_by_overlaps(bw) %>%
        plyranges::group_by(query) %>%
        plyranges::summarise(score = sum(score))
    if (!1 %in% rc$query) {
        rc <- rbind(DataFrame(list(query = 1, score = 0)), rc)
    }
    rc <- as.list(rc$score)
    names(rc) <- regionNames
    return(rc)
}

prepareSimulationParameters <- function(simpol) {
    sampleCell <- 5000
    sampleN <- 50
    matchedLen <- 2e4
    kmin <- 1
    kmax <- 200
    matchedGbLen <- matchedLen - kmax

    spacing <- slot(simpol, "polSize") + slot(simpol, "addSpace")
    k <- slot(simpol, "k")
    alpha <- slot(simpol, "alpha")
    beta <- slot(simpol, "beta")
    prob <- siteProbabilities(simpol)
    startPoint <- 0.99 * 1e6
    lambda <- 102.1  # from Dukler et al. 2017

    rnapPos <- finalPositionMatrix(simpol)
    totalCell <- NCOL(rnapPos)
    geneLen <- NROW(rnapPos) - 1

    return(list(
        sampleCell = sampleCell,
        sampleN = sampleN,
        matchedLen = matchedLen,
        kmin = kmin,
        kmax = kmax,
        matchedGbLen = matchedGbLen,
        spacing = spacing,
        k = k,
        startPoint = startPoint,
        lambda = lambda,
        rnapPos = rnapPos,
        totalCell = totalCell,
        geneLen = geneLen,
        alpha = alpha,
        beta = beta,
        prob = prob
    ))
}

createGenomicRegions <- function(params) {
    gnRng <- GRanges(seqnames = rep("chr1", 3),
        IRanges::IRanges(start = c(1, params$kmax + 1, 1), 
                end = c(params$kmax, params$geneLen, params$spacing))
    )

    gnRng <- IRanges::shift(gnRng, shift = params$startPoint)

    regionNames <- c("tss", "gb", "landing")
    names(gnRng) <- regionNames

    len <- as.list(width(gnRng))
    names(len) <- regionNames

    return(list(
        gnRng = gnRng,
        regionNames = regionNames,
        len = len
    ))
}

generateRnapPositions <- function(params, regions) {
    seeds <- seq(from = 2013, by = 1, length.out = params$sampleN)
    rnapGrng <- list()


    beta_prob <- params$prob[1, 1] / params$alpha * params$beta
    
    # Find pause sites for each cell by looking for the beta_prob value in each column
    pauseSitesPerCell <- numeric(params$totalCell)
    for (cell in 1:params$totalCell) {
        # Find where beta_prob occurs in this cell's probability vector
        cell_probs <- params$prob[, cell]
        pause_indices <- which(cell_probs == beta_prob)
        # Remove the first position (position 0) as it's not a pause site
        pause_indices <- pause_indices[pause_indices != 1]
        
        if (length(pause_indices) > 0) {
            # Take the first pause site found
            pauseSitesPerCell[cell] <- pause_indices[1] - 1
        } else {
            # If no pause site found, use a default value (e.g., k)
            pauseSitesPerCell[cell] <- params$k
        }
    }

    for (i in seq_len(params$sampleN)) {
        set.seed(seeds[i])
        selCells <- sample(seq_len(params$totalCell), 
            size = params$sampleCell, replace = TRUE)
        resPos <- params$rnapPos[, selCells]
        resPos <- resPos[-1, ]


        # Use the pre-calculated pause sites for the selected cells
        pauseSite <- pauseSitesPerCell[selCells]
        resShape <- dim(resPos)
        afterPauseLen <- resShape[1] - pauseSite

        
        maskMx <- map2(pauseSite, afterPauseLen,
            function(x, y) c(rep(TRUE, x), rep(FALSE, y))
        )
        maskMx <- matrix(unlist(maskMx), 
            nrow = resShape[1], ncol = resShape[2])
        

        # calculate rnap positions across all cells
        resAll <- rowSums(resPos)
        # generate bigwigs for positive strand
        rnapGrng[[i]] <- GRanges(seqnames = "chr1",
            IRanges::IRanges(start = (1 + params$startPoint):
                (params$geneLen + params$startPoint),
                width = 1), 
            score = resAll,
            strand = "+",
            seqlengths = c("chr1" = params$geneLen * 10) + 
                params$startPoint)

        rm(resPos, resAll)
    }

    return(rnapGrng)
}

calculateReadCounts <- function(rnapGrng, regions, params) {
    bwDfs <- tibble(
        trial = seq_len(params$sampleN),
        rcRegion = vector("list", params$sampleN),
        rcTss = numeric(params$sampleN),
        rcGb = numeric(params$sampleN),
        rcLanding = numeric(params$sampleN),
        R = numeric(params$sampleN),
        Rpause = numeric(params$sampleN),
        rnapProp = numeric(params$sampleN)
    )

    bwDfs$rcRegion <- map(rnapGrng, 
        ~ summariseSimulationBw(.x, regions$gnRng, regions$regionNames))

    bwDfs$rcTss <- map_dbl(bwDfs$rcRegion, "tss")
    bwDfs$rcGb <- map_dbl(bwDfs$rcRegion, "gb")
    bwDfs$rcLanding <- map_dbl(bwDfs$rcRegion, "landing")

    bwDfs$R <- (bwDfs$rcTss + bwDfs$rcGb) / params$sampleCell
    bwDfs$Rpause <- bwDfs$rcTss / params$sampleCell
    bwDfs$rnapProp <- bwDfs$rcLanding / params$sampleCell

    return(bwDfs)
}

adjustReadCoverage <- function(rnapGrng, regions, params, bwDfs) {
    if (!is.null(params$lambda)) {
        rnapGrng <- map(rnapGrng, function(grng) {
            grng$score[params$kmin:params$kmax] <- rpois(
                length(params$kmin:params$kmax),
                grng$score[params$kmin:params$kmax] / 
                    params$sampleCell * params$lambda
            )
            grng$score[seq_len(20)] <- 0
            return(grng)
        })
        
        bwDfs$rcRegion <- map(rnapGrng, 
            ~ summariseSimulationBw(.x, regions$gnRng, regions$regionNames))
        bwDfs$rcTss <- map_dbl(bwDfs$rcRegion, "tss")

        poisMean <- (params$lambda * bwDfs$rcGb / params$sampleCell) *
            (params$matchedGbLen / regions$len$gb)
        bwDfs$rcGb <- rpois(length(poisMean), poisMean)
        regions$len$gb <- params$matchedGbLen
    }

    return(list(bwDfs = bwDfs, rnapGrng = rnapGrng))
}

calculateInitialRates <- function(bwDfs, regions, params, simpol, rnapGrng) {
    bwDfs$Xk <- map(
        rnapGrng,
        ~ .x[(start(.x) >= 990000 + params$kmin) &
            (start(.x) <= 990000 + params$kmax), ]$score
    )

    bwDfs$chi <- bwDfs$rcGb / regions$len$gb

    if (parameters(simpol)$ksd == 0) {
        bwDfs$betaOrg <- bwDfs$chi / map_dbl(bwDfs$Xk, params$k)
    } else {
        bwDfs$betaOrg <- bwDfs$chi / (bwDfs$rcTss / regions$len$tss)
    }

    bwDfs$betaMaxRc <- bwDfs$chi / map_dbl(bwDfs$Xk, max)
    bwDfs$XkSum <- vapply(bwDfs$Xk, sum, numeric(1))

    validIndices <- which(bwDfs$XkSum > 0 &
        bwDfs$chi > 0 &
        !is.infinite(bwDfs$chi))

    if (length(validIndices) == 0) {
        stop("No RNAPs in the pause region or gene body")
    }

    bwDfs$betaInt <- rep(NA, nrow(bwDfs))
    bwDfs$betaInt[validIndices] <- bwDfs$chi[validIndices] /
        bwDfs$XkSum[validIndices]

    return(bwDfs)
}

runEmAlgorithm <- function(bwDfs, params, stericHindrance) {
    fkInt <- dnorm(params$kmin:params$kmax, mean = 50, sd = 100)
    fkInt <- fkInt / sum(fkInt)
    emLs <- list()

    if (stericHindrance) {
        f <- calculateF(s = params$spacing, k = params$k)
        phiInt <- 0.5
        zeta <- 2000
        lambda1 <- 0.0505 * zeta^2
    }

    message("Starting EM algorithm...")

    for (i in seq_len(NROW(bwDfs))) {
        rc <- bwDfs[i, ]
        if (!stericHindrance) {
            emLs[[i]] <- pauseEscapeEM(
                Xk = rc$Xk[[1]], 
                kmin = params$kmin, 
                kmax = params$kmax,
                fkInt = fkInt, 
                betaInt = rc$betaInt[[1]],
                chiHat = rc$chi, 
                maxItr = 500, 
                tor = 1e-4
            )
        } else {
            emLs[[i]] <- stericHindranceEM(
                Xk = rc$Xk[[1]], 
                kmin = params$kmin,
                kmax = params$kmax, 
                f1 = f[["f1"]], 
                f2 = f[["f2"]],
                fkInt = fkInt, 
                betaInt = rc$betaInt[[1]],
                chiHat = rc$chi, 
                phiInt = phiInt,
                lambda = lambda1, 
                zeta = zeta, 
                maxItr = 500,
                tor = 1e-4
            )
        }
    }

    return(emLs)
}

processEmResults <- function(bwDfs, emLs, stericHindrance) {
    bwDfs$betaAdp <- map_dbl(emLs, "beta", .default = NA)
    bwDfs$Yk <- map(emLs, "Yk", .default = NA)
    bwDfs$fk <- map(emLs, "fk", .default = NA)
    bwDfs$fkMean <- map_dbl(emLs, "fkMean", .default = NA)
    bwDfs$fkVar <- map_dbl(emLs, "fkVar", .default = NA)
    bwDfs$flag <- map_chr(emLs, "flag", .default = NA)
    
    if (stericHindrance) {
        bwDfs$phi <- map_dbl(emLs, "phi", .default = NA)
    }

    return(bwDfs)
}

#' @title Estimate transcription rates from simulation data
#'
#' @description
#' Estimates transcription rates from a SimulatePolymerase object using
#' Expectation Maximization likelihood formula.
#'
#' @param simpol a \code{\linkS4class{SimulatePolymerase}} object
#' @param stericHindrance a logical value to determine whether to infer
#' landing-pad occupancy or not. Defaults to FALSE.
#' @param ... Additional arguments (not used)
#' @return a \code{\linkS4class{SimulationTranscriptionRates}} object
#' 
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1)
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim)
#' # Print the estimated rates
#' print(estRates)
#'
#' @rdname SimulationTranscriptionRates-class
#' @export
setMethod("estimateTranscriptionRates", "SimulatePolymerase", 
function(x, stericHindrance=FALSE, ...) {
    
    simpol <- x  # x is the SimulatePolymerase object

    if (!is(simpol, "SimulatePolymerase")) {
        stop("simpol parameter must be a SimulatePolymerase object")
    }
    if (!is.logical(stericHindrance)) {
        stop("stericHindrance parameter must be a logical value")
    }
    
    # Prepare parameters and regions
    params <- prepareSimulationParameters(simpol)
    regions <- createGenomicRegions(params)

    # Generate RNAP positions
    rnapGrng <- generateRnapPositions(params, regions)

    # Calculate initial read counts
    bwDfs <- calculateReadCounts(rnapGrng, regions, params)

    # Adjust read coverage if needed
    adjustedData <- adjustReadCoverage(rnapGrng, regions, params, bwDfs)
    bwDfs <- adjustedData$bwDfs
    rnapGrng <- adjustedData$rnapGrng
    
    # Calculate initial rates
    bwDfs <- calculateInitialRates(bwDfs, regions, params, simpol, rnapGrng)
    
    # Run EM algorithm
    emLs <- runEmAlgorithm(bwDfs, params, stericHindrance)
    
    # Process EM results
    bwDfs <- processEmResults(bwDfs, emLs, stericHindrance)
    print(bwDfs)
    
    rates_tibble <- tibble(
        trial = bwDfs$trial,
        chi = bwDfs$chi,
        betaOrg = bwDfs$betaOrg,
        betaAdp = bwDfs$betaAdp,
        fkMean = bwDfs$fkMean,
        fkVar = bwDfs$fkVar,
        rcTss = bwDfs$rcTss,
        rcGb = bwDfs$rcGb,
        rcLanding = bwDfs$rcLanding,
        R = bwDfs$R,
        Rpause = bwDfs$Rpause,
        rnapProp = bwDfs$rnapProp
    )
    
    if (stericHindrance && "phi" %in% colnames(bwDfs)) {
        rates_tibble$phi <- bwDfs$phi
    }
    
    # Create and return the final object
    new("SimulationTranscriptionRates",
        simpol = simpol, 
        stericHindrance = stericHindrance,
        rates = rates_tibble
    )
})

#' @rdname SimulationTranscriptionRates-class
#' @title Show Method for SimulationTranscriptionRates Object
#'
#' @description
#' Show method for SimulationTranscriptionRates object in human readable format
#' including summary statistics
#'
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1)
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim)
#' # Show the object
#' show(estRates)
#' @export
setMethod("show", "SimulationTranscriptionRates", function(object) {
    cat("A SimulationTranscriptionRates object with:\n")
    cat("  - Steric hindrance:", stericHindrance(object), "\n")
    cat("  - number of trials:", nrow(rates(object)), "\n")

    ratesData <- rates(object)
    numericCols <- sapply(ratesData, is.numeric)
    excludeCols <- c("trial")
    if (!stericHindrance(object)) {
        excludeCols <- c(excludeCols, "phi")
    }
    validCols <- names(ratesData)[numericCols & !(names(ratesData) %in% excludeCols)]
    
    if (length(validCols) > 0) {
        cat("\nSummary statistics for rate estimates:\n")
        for (col in validCols) {
            values <- ratesData[[col]]
            values <- values[!is.na(values)]
            if (length(values) > 0) {
                meanVal <- mean(values, na.rm = TRUE)
                cat(sprintf("  - %s: %.4f (mean)\n", col, meanVal))
            }
        }
    }
    
})

#' @rdname SimulationTranscriptionRates-class
#' @title Plot Transcription Rates
#'
#' @description
#' Plot the transcription rates across trials.
#'
#' @param object A SimulationTranscriptionRates object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the transcription rates
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1)
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim)
#' # Plot transcription rates
#' plotTranscriptionRates(estRates)
#' @export
setGeneric("plotTranscriptionRates", function(
    object, file = NULL,
    width = 8, height = 6) {
    standardGeneric("plotTranscriptionRates")
})
setMethod(
    "plotTranscriptionRates", "SimulationTranscriptionRates",
    function(object, file = NULL, width = 8, height = 6) {
        # Get the rates tibble
        rates_df <- rates(object)
        
        # Create plot
        p <- ggplot(rates_df, aes(x = trial)) +
            geom_line(aes(y = chi, color = "chi")) +
            geom_line(aes(y = betaOrg, color = "betaOrg")) +
            geom_line(aes(y = betaAdp, color = "betaAdp")) +
            theme_minimal() +
            labs(
                title = "Transcription Rates Across Trials",
                x = "Trial",
                y = "Rate",
                color = "Rate Type"
            )

        # Add phi line if steric hindrance is enabled
        if (stericHindrance(object) && "phi" %in% colnames(rates_df)) {
            p <- p + geom_line(aes(y = phi, color = "phi"))
        }

        if (!is.null(file)) {
            ggsave(file, p, width = width, height = height)
        }

        return(p)
    }
)

#' @rdname SimulationTranscriptionRates-class
#' @title Plot Pause Site Distribution
#'
#' @description
#' Plot the pause site distribution across trials.
#'
#' @param object A SimulationTranscriptionRates object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the pause site distribution
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1)
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim)
#' # Plot pause site distribution
#' plotPauseSiteDistribution(estRates)
#' @export
setGeneric("plotPauseSiteDistribution", function(
    object, file = NULL,
    width = 8, height = 6) {
    standardGeneric("plotPauseSiteDistribution")
})
setMethod(
    "plotPauseSiteDistribution", "SimulationTranscriptionRates",
    function(object, file = NULL, width = 8, height = 6) {
        # Get the rates tibble
        rates_df <- rates(object)
        
        # Check if fkMean and fkVar columns exist
        if (!all(c("fkMean", "fkVar") %in% colnames(rates_df))) {
            stop("fkMean and fkVar columns not found in rates tibble")
        }

        # Create plot
        p <- ggplot(rates_df, aes(x = fkMean, y = fkVar)) +
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

# Accessor methods
#' @rdname SimulationTranscriptionRates-class
#' @title Accessor for SimulatePolymerase Object
#'
#' @description
#' Accessor for the SimulatePolymerase object from a
#' SimulationTranscriptionRates object.
#'
#' @param object a \code{SimulationTranscriptionRates} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1)
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim)
#' # Get simpol
#' simpol <- simpol(estRates)
#' # Print the simpol
#' print(simpol)
#' @export
setGeneric("simpol", function(object) standardGeneric("simpol"))
setMethod(
    "simpol", "SimulationTranscriptionRates",
    function(object) slot(object, "simpol")
)

#' @rdname SimulationTranscriptionRates-class
#' @title Accessor for Steric Hindrance
#'
#' @description
#' Accessor for the steric hindrance flag from a
#' SimulationTranscriptionRates object.
#'
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1)
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim)
#' # Get steric hindrance
#' stericHindrance <- stericHindrance(estRates)
#' # Print the steric hindrance
#' print(stericHindrance)
#' @export
setMethod("stericHindrance", "SimulationTranscriptionRates", function(object) {
    slot(object, "stericHindrance")
})

#' @rdname SimulationTranscriptionRates-class
#' @title Rates Tibble
#'
#' @description
#' Accessor for the tibble containing the estimated rates from a
#' SimulationTranscriptionRates object.
#'
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1)
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim)
#' # Get rates
#' rates <- rates(estRates)
#' # Print the rates
#' print(rates)
#' @export
setMethod("rates", "SimulationTranscriptionRates", function(object) {
    slot(object, "rates")
})

#' @examples
#' # Create a SimulationTranscriptionRates object
#' simpol <- simulatePolymerase(k = 50, ksd = 10, kMin = 30, kMax = 70,
#'     geneLen = 1000, alpha = 0.1, beta = 0.2, zeta = 1000,
#'     zetaSd = 100, zetaMin = 800, zetaMax = 1200,
#'     cellNum = 1000, polSize = 35, addSpace = 15,
#'     time = 10)
#'
#' # Estimate transcription rates
#' rates <- estimateTranscriptionRates(simpol)
#'
#' # Show the object
#' show(rates)
