simulationTranscriptionRatesValid <- function(object) {
    errors <- character()

    # Check if required slots are present and have correct types
    if (!is(slot(object, "simpol"), "simulatePolymerase")) {
        errors <- c(errors, "simpol must be a simulatePolymerase object")
    }
    if (!is.logical(slot(object, "stericHindrance"))) {
        errors <- c(errors, "stericHindrance must be logical")
    }
    if (!is.numeric(slot(object, "trial"))) {
        errors <- c(errors, "trial must be numeric")
    }
    if (!is.numeric(slot(object, "chi"))) {
        errors <- c(errors, "chi must be numeric")
    }
    if (!is.numeric(slot(object, "betaOrg"))) {
        errors <- c(errors, "betaOrg must be numeric")
    }
    if (!is.numeric(slot(object, "betaAdp"))) {
        errors <- c(errors, "betaAdp must be numeric")
    }
    if (!is.numeric(slot(object, "phi"))) {
        errors <- c(errors, "phi must be numeric")
    }
    if (!is.list(slot(object, "fk"))) {
        errors <- c(errors, "fk must be a list")
    }
    if (!is.numeric(slot(object, "fkMean"))) {
        errors <- c(errors, "fkMean must be numeric")
    }
    if (!is.numeric(slot(object, "fkVar"))) {
        errors <- c(errors, "fkVar must be numeric")
    }
    if (!is.character(slot(object, "flag"))) {
        errors <- c(errors, "flag must be character")
    }
    if (!is.list(slot(object, "rnapN"))) {
        errors <- c(errors, "rnapN must be a list")
    }

    if (length(errors) == 0) TRUE else errors
}

#' Class SimulationTranscriptionRates
#'
#' Class \code{SimulationTranscriptionRates} estimates the transcription
#' rates, such as initiation, pause-release rates and landing pad occupancy,
#' from simulated data, such as nascent RNA sequencing read counts and genomic
#' coordinates simulated from polymerase movement over a specified time period
#' calculated by simulator, and contructs an object that holds these rates.
#'
#' @name SimulationTranscriptionRates-class
#' @rdname SimulationTranscriptionRates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @importClassesFrom data.table data.table
#' @importFrom purrr map map_chr map_dbl
#' @exportClass SimulationTranscriptionRates
methods::setClass("SimulationTranscriptionRates",
    slots = c(
        simpol = "SimulatePolymerase",
        stericHindrance = "logical",
        trial = "numeric",
        chi = "numeric",
        betaOrg = "numeric",
        betaAdp = "numeric",
        phi = "numeric",
        fk = "list",
        fkMean = "numeric",
        fkVar = "numeric",
        flag = "character",
        rnapN = "list"
    ),
    prototype = list(
        simpol = NULL,
        stericHindrance = FALSE,
        trial = numeric(0),
        chi = numeric(0),
        betaOrg = numeric(0),
        betaAdp = numeric(0),
        phi = 0,
        fk = list(),
        fkMean = numeric(0),
        fkVar = numeric(0),
        flag = character(0),
        rnapN = list()
    ),
    validity = simulationTranscriptionRatesValid
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
    countRnap <- FALSE

    spacing <- slot(simpol, "polSize") + slot(simpol, "addSpace")
    k <- slot(simpol, "k")
    startPoint <- 0.99 * 1e6
    lambda <- 102.1  # from Dukler et al. 2017

    rnapPos <- positionMatrix(simpol)
    totalCell <- NCOL(rnapPos)
    geneLen <- NROW(rnapPos) - 1

    return(list(
        sampleCell = sampleCell,
        sampleN = sampleN,
        matchedLen = matchedLen,
        kmin = kmin,
        kmax = kmax,
        matchedGbLen = matchedGbLen,
        countRnap = countRnap,
        spacing = spacing,
        k = k,
        startPoint = startPoint,
        lambda = lambda,
        rnapPos = rnapPos,
        totalCell = totalCell,
        geneLen = geneLen
    ))
}

createGenomicRegions <- function(params) {
    gnRng <- GRanges(seqnames = rep("chr1", 3),
        IRanges(start = c(1, params$kmax + 1, 1), 
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
    rnapNls <- list()

    for (i in seq_len(params$sampleN)) {
        selCells <- sample(seq_len(params$totalCell), 
            size = params$sampleCell, replace = TRUE)
        resPos <- params$rnapPos[, selCells]
        resPos <- resPos[-1, ]

        if (params$countRnap) {
            pauseSite <- idx[selCells, 1] - 1
            resShape <- dim(resPos)
            afterPauseLen <- resShape[1] - pauseSite
            maskMx <- map2(pauseSite, afterPauseLen,
                function(x, y) c(rep(TRUE, x), rep(FALSE, y))
            )
            maskMx <- Matrix::Matrix(unlist(maskMx), 
                nrow = resShape[1], ncol = resShape[2])
            rnapNls[[i]] <- colSums(resPos * maskMx)
        }

        resAll <- rowSums(resPos)
        rnapGrng[[i]] <- GRanges(seqnames = "chr1",
            IRanges(start = (1 + params$startPoint):
                (params$geneLen + params$startPoint),
                width = 1), 
            score = resAll,
            strand = "+",
            seqlengths = c("chr1" = params$geneLen * 10) + 
                params$startPoint)

        rm(resPos, resAll)
    }

    return(list(rnapGrng = rnapGrng, rnapNls = rnapNls))
}

calculateReadCounts <- function(rnapData, regions, params) {
    bwDfs <- tibble(trial = seq_len(params$sampleN))
    bwDfs$rcRegion <- map(rnapData$rnapGrng, 
        ~ summariseSimulationBw(.x, regions$gnRng, regions$regionNames))

    bwDfs$rcTss <- map_dbl(bwDfs$rcRegion, "tss")
    bwDfs$rcGb <- map_dbl(bwDfs$rcRegion, "gb")
    bwDfs$rcLanding <- map_dbl(bwDfs$rcRegion, "landing")

    bwDfs <- bwDfs %>% mutate(
        R = (rcTss + rcGb) / params$sampleCell,
        Rpause = rcTss / params$sampleCell,
        rnapProp = rcLanding / params$sampleCell
    )

    return(bwDfs)
}

adjustReadCoverage <- function(rnapData, regions, params, bwDfs) {
    if (!is.null(params$lambda)) {
        rnapGrng <- map(rnapData$rnapGrng, function(grng) {
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

    if (slot(simpol, "ksd") == 0) {
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

    bwDfs$initRate <- bwDfs$R / (bwDfs$Rpause + bwDfs$R)
    bwDfs$pauseReleaseRate <- bwDfs$Rpause / (bwDfs$Rpause + bwDfs$R)

    if (stericHindrance) {
        bwDfs$landingPadOccupancy <- bwDfs$rnapProp
    }

    return(bwDfs)
}

calculateFinalRates <- function(bwDfs, stericHindrance) {
    initRateMean <- mean(bwDfs$initRate)
    initRateVar <- var(bwDfs$initRate)
    pauseReleaseRateMean <- mean(bwDfs$pauseReleaseRate)
    pauseReleaseRateVar <- var(bwDfs$pauseReleaseRate)

    results <- list(
        initRate = initRateMean,
        initRateVar = initRateVar,
        pauseReleaseRate = pauseReleaseRateMean,
        pauseReleaseRateVar = pauseReleaseRateVar
    )

    if (stericHindrance) {
        landingPadOccupancyMean <- mean(bwDfs$landingPadOccupancy)
        landingPadOccupancyVar <- var(bwDfs$landingPadOccupancy)
        results$landingPadOccupancy <- landingPadOccupancyMean
        results$landingPadOccupancyVar <- landingPadOccupancyVar
    }

    return(results)
}

#' estimateSimulationTranscriptionRates
#'
#' Estimates the transcription rates, such as initiation, pause-release rates
#' and landing pad occupancy, from simulated data, such as nascent RNA
#' sequencing read counts and genomic coordinates, and contructs an object
#' that holds these rates
#'
#' @param simpol a \code{\link{simulatePolymerase-class}} object
#' @param stericHindrance a logical value to determine whether to infer
#' landing-pad occupancy or not.
#' Defaults to FALSE.
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @return an \code{\link{SimulationTranscriptionRates-class}} object
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Print the estimated rates
#' print(estRates)
#'
#' @export
estimateSimulationTranscriptionRates <- function(simpol, stericHindrance=FALSE)
{
    # Prepare parameters and regions
    params <- prepareSimulationParameters(simpol)
    regions <- createGenomicRegions(params)

    # Generate RNAP positions
    rnapData <- generateRnapPositions(params, regions)

    # Calculate initial read counts
    bwDfs <- calculateReadCounts(rnapData, regions, params)

    # Adjust read coverage if needed
    adjustedData <- adjustReadCoverage(rnapData, regions, params, bwDfs)
    bwDfs <- adjustedData$bwDfs
    rnapGrng <- adjustedData$rnapGrng
    
    # Calculate initial rates
    bwDfs <- calculateInitialRates(bwDfs, regions, params, simpol, rnapGrng)
    
    # Run EM algorithm
    emLs <- runEmAlgorithm(bwDfs, params, stericHindrance)
    
    # Process EM results
    bwDfs <- processEmResults(bwDfs, emLs, stericHindrance)
    
    # Calculate final rates
    results <- calculateFinalRates(bwDfs, stericHindrance)
    
    # Create and return the final object
    new("SimulationTranscriptionRates",
        simpol = simpol, 
        stericHindrance = stericHindrance,
        trial = params$sampleN, 
        chi = results$initRate,
        betaOrg = results$pauseReleaseRate,
        betaAdp = if (stericHindrance) results$landingPadOccupancy else 0,
        phi = if (stericHindrance) results$landingPadOccupancy else 0,
        fk = results, 
        fkMean = c(results$initRate, results$pauseReleaseRate),
        fkVar = c(results$initRateVar, results$pauseReleaseRateVar),
        flag = if (stericHindrance) "withStericHindrance" else 
            "noStericHindrance",
        rnapN = if (params$countRnap) rnapData$rnapNls else list()
    )
}

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Show the object
#' show(estRates)
#' @export
setMethod("show", "SimulationTranscriptionRates", function(object) {
    # Create a data frame for display
    df <- data.frame(
        Parameter = c(
            "Number of Trials", "Initiation Rate",
            "Pause Release Rate", "Steric Hindrance"
        ),
        Value = c(
            trial(object),
            round(chi(object), 4),
            round(betaOrg(object), 4),
            stericHindrance(object)
        )
    )

    # Print the object information
    cat("A SimulationTranscriptionRates object:\n\n")
    print(df, row.names = FALSE)
})


# Plotting methods
#' Plot transcription rates
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
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
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
        # Create data frame for plotting
        df <- data.frame(
            trial = trial(object),
            chi = chi(object),
            betaOrg = betaOrg(object),
            betaAdp = betaAdp(object),
            phi = phi(object),
            fkMean = fkMean(object),
            fkVar = fkVar(object)
        )

        # Create plot
        p <- ggplot(df, aes(x = trial)) +
            geom_line(aes(y = chi, color = "chi")) +
            geom_line(aes(y = betaOrg, color = "betaOrg")) +
            geom_line(aes(y = betaAdp, color = "betaAdp")) +
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
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
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
        # Create data frame for plotting
        df <- data.frame(
            trial = trial(object),
            fkMean = fkMean(object),
            fkVar = fkVar(object)
        )

        # Create plot
        p <- ggplot(df, aes(x = fkMean, y = fkVar)) +
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
#' @param object A SimulationTranscriptionRates object
#' @param file Optional file path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return A ggplot object showing the RNAP counts
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Plot RNAP counts
#' plotRnapCounts(estRates)
#' @export
setGeneric("plotRnapCounts", function(
    object, file = NULL, width = 8,
    height = 6) {
    standardGeneric("plotRnapCounts")
})
setMethod(
    "plotRnapCounts", "SimulationTranscriptionRates",
    function(object, file = NULL, width = 8, height = 6) {
        rnapN <- rnapN(object)
        if (length(rnapN) == 0) {
            stop("No RNAP counts available")
        }

        # Create data frame for plotting
        df <- data.frame(
            trial = rep(trial(object), each = length(rnapN[[1]])),
            rnapCount = unlist(rnapN)
        )

        # Create plot
        p <- ggplot(df, aes(x = rnapCount)) +
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

# Accessor methods
#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
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
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get steric hindrance
#' stericHindrance <- stericHindrance(estRates)
#' # Print the steric hindrance
#' print(stericHindrance)
#' @export
setMethod("stericHindrance", "SimulationTranscriptionRates", function(object) {
    slot(object, "stericHindrance")
})

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get trial
#' trial <- trial(estRates)
#' # Print the trial
#' print(trial)
#' @export
setGeneric("trial", function(object) standardGeneric("trial"))
setMethod(
    "trial", "SimulationTranscriptionRates",
    function(object) slot(object, "trial")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get chi
#' chi <- chi(estRates)
#' # Print the chi
#' print(chi)
#' @export
setGeneric("chi", function(object) standardGeneric("chi"))
setMethod(
    "chi", "SimulationTranscriptionRates",
    function(object) slot(object, "chi")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get betaOrg
#' betaOrg <- betaOrg(estRates)
#' # Print the betaOrg
#' print(betaOrg)
#' @export
setGeneric("betaOrg", function(object) standardGeneric("betaOrg"))
setMethod(
    "betaOrg", "SimulationTranscriptionRates",
    function(object) slot(object, "betaOrg")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get betaAdp
#' betaAdp <- betaAdp(estRates)
#' # Print the betaAdp
#' print(betaAdp)
#' @export
setGeneric("betaAdp", function(object) standardGeneric("betaAdp"))
setMethod(
    "betaAdp", "SimulationTranscriptionRates",
    function(object) slot(object, "betaAdp")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get phi
#' phi <- phi(estRates)
#' # Print the phi
#' print(phi)
#' @export
setGeneric("phi", function(object) standardGeneric("phi"))
setMethod(
    "phi", "SimulationTranscriptionRates",
    function(object) slot(object, "phi")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get fk
#' fk <- fk(estRates)
#' # Print the fk
#' print(fk)
#' @export
setGeneric("fk", function(object) standardGeneric("fk"))
setMethod(
    "fk", "SimulationTranscriptionRates",
    function(object) slot(object, "fk")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get fkMean
#' fkMean <- fkMean(estRates)
#' # Print the fkMean
#' print(fkMean)
#' @export
setGeneric("fkMean", function(object) standardGeneric("fkMean"))
setMethod(
    "fkMean", "SimulationTranscriptionRates",
    function(object) slot(object, "fkMean")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get fkVar
#' fkVar <- fkVar(estRates)
#' # Print the fkVar
#' print(fkVar)
#' @export
setGeneric("fkVar", function(object) standardGeneric("fkVar"))
setMethod(
    "fkVar", "SimulationTranscriptionRates",
    function(object) slot(object, "fkVar")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get flag
#' flag <- flag(estRates)
#' # Print the flag
#' print(flag)
#' @export
setGeneric("flag", function(object) standardGeneric("flag"))
setMethod(
    "flag", "SimulationTranscriptionRates",
    function(object) slot(object, "flag")
)

#' @rdname SimulationTranscriptionRates-class
#' @examples
#' # Create a SimulatePolymerase object
#' sim <- simulatePolymerase(
#'     k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
#'     alpha=1, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
#'     zetaVec=NULL, cellNum=1000, polSize=33, addSpace=17, time=1, 
#'     stepsToRecord=1)
#' # Estimate transcription rates
#' estRates <- estimateSimulationTranscriptionRates(sim)
#' # Get rnapN
#' rnapN <- rnapN(estRates)
#' # Print the rnapN
#' print(rnapN)
#' @export
setGeneric("rnapN", function(object) standardGeneric("rnapN"))
setMethod(
    "rnapN", "SimulationTranscriptionRates",
    function(object) slot(object, "rnapN")
)

#' @examples
#' # Create a SimulationTranscriptionRates object
#' simpol <- simulatePolymerase(k = 50, ksd = 10, kMin = 30, kMax = 70,
#'     geneLen = 1000, alpha = 0.1, beta = 0.2, zeta = 1000,
#'     zetaSd = 100, zetaMin = 800, zetaMax = 1200,
#'     cellNum = 1000, polSize = 35, addSpace = 15,
#'     time = 10, stepsToRecord = 100
#' )
#'
#' # Estimate transcription rates
#' rates <- estimateSimulationTranscriptionRates(simpol)
#'
#' # Show the object
#' show(rates)
