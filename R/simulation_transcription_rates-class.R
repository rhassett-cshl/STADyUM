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
#' @slot name a character value for the name of the experiment
#' @slot stericHindrance a logical value to determine whether to infer
#' landing-pad occupancy or not. Defaults to FALSE.
#' @slot rates a \code{\link[tibble]{tbl_df}} containing the estimated rates
#' with columns:
#' \describe{
#' \item{trial}{Numeric. Trial number, each trial samples 5000 cells}
#' \item{chi}{Numeric. RNAP density along gene body given as estimate for the
#' gene body elongation rate \verb{[RNAPs/bp]}}
#' \item{betaOrg}{Numeric. Ratio of gene body RNAP density to pause region
#' RNAP density with fixed pause sites given as an estimate for the
#' pause-escape rate}
#' \item{betaAdp}{Numeric. Ratio of gene body RNAP density to pause region
#' RNAP density from adapted model which allows pause sites to vary across
#' cells given as an estimate for the pause-escape rate}
#' \item{phi}{Numeric. Landing-pad occupancy representing the fraction of time
#' in the simulation that the landing pad is occupied by RNA polymerase. The
#' landing pad, also known as the initiation site, is the region where the RNA
#' polymerase binds to the DNA and begins transcription. It is the polSize +
#' addSpace length (only applicable if steric hindrance is enabled)}
#' \item{fk}{list. Likelihood of pausing at each pause region position}
#' \item{fkMean}{Numeric. Mean position of pause sites}
#' \item{fkVar}{Numeric. Variance of pause site positions}
#' \item{totalTssRc}{Numeric. Total RNAP read counts in the TSS region (across
#' all cells)}
#' \item{totalGbRc}{Numeric. Total RNAP read counts in the gene body region
#' (across all cells)}
#' \item{totalLandingRc}{Numeric. Total RNAP read counts in the landing pad
#' region (across all cells)}
#' \item{avgRcPerCell}{Numeric. Average number of RNAPs per cell for gene
#' body and TSS regions}
#' \item{avgTssRcPerCell}{Numeric. Average number of RNAPs per cell for
#' pause regions (across all cells)}
#' \item{avgLandingRcPerCell}{Numeric. Average number of RNAPs per cell for
#' landing pad regions (across all cells)}
#' \item{actualPauseSiteCounts}{list. Simulated pause site counts for each
#' cell}
#' \item{expectedPauseSiteCounts}{list. Expected pause site counts for each
#' cell estimated from the EM algorithm}
#' \item{expectationMaximizationStatus}{character. Status of the expectation
#' maximization algorithm. "normal": converged normally, "single_site":
#' converged to single pause site, "max_iterations": reached maximum
#' iterations without convergence. Note max # of iterations is 500.}
#' \item{likelihood}{a numeric vector of the likelihood of the model}
#' }
#'
#' @name SimulationTranscriptionRates-class
#' @rdname SimulationTranscriptionRates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @importClassesFrom data.table data.table
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom tibble tibble
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
        name = "character",
        stericHindrance = "logical",
        rates = "tbl_df"
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
    ksd <- slot(simpol, "ksd")
    alpha <- slot(simpol, "alpha")
    beta <- slot(simpol, "beta")
    prob <- siteProbabilities(simpol)
    startPoint <- 0.99 * 1e6
    lambda <- 102.1 # from Dukler et al. 2017

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
        ksd = ksd,
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
    gnRng <- GRanges(
        seqnames = rep("chr1", 3),
        IRanges::IRanges(
            start = c(1, params$kmax + 1, 1),
            end = c(params$kmax, params$geneLen, params$spacing)
        )
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

    # Find pause sites for each cell from the beta_prob value in each column
    pauseSitesPerCell <- numeric(params$totalCell)
    for (cell in seq_len(params$totalCell)) {
        cell_probs <- params$prob[, cell]
        pause_indices <- which(cell_probs == beta_prob)
        pause_indices <- pause_indices[pause_indices != 1]

        if (length(pause_indices) > 0) {
            pauseSitesPerCell[cell] <- pause_indices[1] - 1
        } else {
            pauseSitesPerCell[cell] <- params$k
        }
    }

    for (i in seq_len(params$sampleN)) {
        selCells <- sample(seq_len(params$totalCell),
            size = params$sampleCell, replace = TRUE
        )
        resPos <- params$rnapPos[, selCells]
        resPos <- resPos[-1, ]
        # Use the pre-calculated pause sites for the selected cells
        pauseSite <- pauseSitesPerCell[selCells]
        resShape <- dim(resPos); afterPauseLen <- resShape[1] - pauseSite
        maskMx <- map2(
            pauseSite, afterPauseLen,
            function(x, y) c(rep(TRUE, x), rep(FALSE, y))
        )
        maskMx <- matrix(unlist(maskMx),
            nrow = resShape[1], ncol = resShape[2]
        )
        # calculate rnap positions across all cells
        resAll <- rowSums(resPos)
        rnapGrng[[i]] <- GRanges(
            seqnames = "chr1",
            IRanges::IRanges(start = (1 + params$startPoint):
            (params$geneLen + params$startPoint), width = 1),
            score = resAll, strand = "+",
            seqlengths = c("chr1" = params$geneLen * 10) + params$startPoint
        )
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

    bwDfs$rcRegion <- map(
        rnapGrng,
        ~ summariseSimulationBw(.x, regions$gnRng, regions$regionNames)
    )

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

        bwDfs$rcRegion <- map(
            rnapGrng,
            ~ summariseSimulationBw(.x, regions$gnRng, regions$regionNames)
        )
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

    if (params$ksd == 0) {
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
    bwDfs$likelihood <- map_dbl(
        emLs,
        ~ .x$likelihoods[[length(.x$likelihoods)]]
    )

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
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500,
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 1000, polSize = 33,
#'     addSpace = 17, time = 1)
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim, name="sim_beta1_k50")
#' # Print the estimated rates
#' print(estRates)
#'
#' @rdname SimulationTranscriptionRates-class
#' @export
setMethod(
    "estimateTranscriptionRates", "SimulatePolymerase",
    function(x, name, stericHindrance = FALSE, ...) {
        simpol <- x # x is the SimulatePolymerase object
        if (!is(simpol, "SimulatePolymerase")) {
            stop("simpol parameter must be a SimulatePolymerase object")
        }
        if (is.null(name) || !is.character(name) || length(name) < 1) {
            stop("name must be a string")
        }
        if (!is.logical(stericHindrance)) {
            stop("stericHindrance parameter must be a logical value")
        }
        params <- prepareSimulationParameters(simpol)
        regions <- createGenomicRegions(params)

        rnapGrng <- generateRnapPositions(params, regions)
        bwDfs <- calculateReadCounts(rnapGrng, regions, params)

        adjustedData <- adjustReadCoverage(rnapGrng, regions, params, bwDfs)
        bwDfs <- adjustedData$bwDfs
        rnapGrng <- adjustedData$rnapGrng

        bwDfs <- calculateInitialRates(bwDfs, regions, params, simpol, rnapGrng)
        emLs <- runEmAlgorithm(bwDfs, params, stericHindrance)
        bwDfs <- processEmResults(bwDfs, emLs, stericHindrance)

        rates_tibble <- tibble(
            trial = bwDfs$trial, chi = bwDfs$chi, betaOrg = bwDfs$betaOrg,
            betaAdp = bwDfs$betaAdp, fk = bwDfs$fk, fkMean = bwDfs$fkMean,
            fkVar = bwDfs$fkVar, totalTssRc = bwDfs$rcTss, 
            totalGbRc = bwDfs$rcGb, totalLandingRc = bwDfs$rcLanding,
            avgRcPerCell = bwDfs$R, avgTssRcPerCell = bwDfs$Rpause,
            avgLandingRcPerCell = bwDfs$rnapProp,
            actualPauseSiteCounts = bwDfs$Xk, 
            expectedPauseSiteCounts = bwDfs$Yk,
            expectationMaximizationStatus = bwDfs$flag,
            likelihood = bwDfs$likelihood
        )

        if (stericHindrance && "phi" %in% colnames(bwDfs)) {
            rates_tibble$phi <- bwDfs$phi
        }

        new("SimulationTranscriptionRates",
            simpol = simpol, name = name, stericHindrance = stericHindrance, 
            rates = rates_tibble
        )
    }
)

showEmStatus <- function(emStatus) {
}

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
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, 
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 1000, polSize = 33,
#'     addSpace = 17, time = 1
#' )
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim, name="sim_beta1_k50")
#' # Show the object
#' show(estRates)
#' @export
setMethod("show", "SimulationTranscriptionRates", function(object) {
    cat("A SimulationTranscriptionRates object with:\n")
    cat("  - Steric hindrance:", stericHindrance(object), "\n")
    cat("  - Number of Trials with 5000 cells:", nrow(rates(object)), "\n")
    cat("\nSummary statistics for rate estimates:\n")
    ratesData <- rates(object)

    chi_mean <- mean(ratesData$chi, na.rm = TRUE)
    cat(sprintf("  - chi (gene body RNAP density): %.2f RNAPs/bp\n", chi_mean))

    betaOrg_mean <- mean(ratesData$betaOrg, na.rm = TRUE)
    cat(sprintf("  - betaOrg (ratio of gene body RNAP density to pause region
        RNAP density, fixed sites): %.4f\n", betaOrg_mean))
    betaAdp_mean <- mean(ratesData$betaAdp, na.rm = TRUE)
    cat(sprintf("  - betaAdp (ratio of gene body RNAP density to pause region
        RNAP density, adapted model): %.4f\n", betaAdp_mean))

    rcTssTrialMean <- mean(ratesData$totalTssRc, na.rm = TRUE)
    cat(sprintf("  - TSS Read Counts Averaged Across Trials: ~ %.0f read
        counts\n", rcTssTrialMean))
    rcGbTrialMean <- mean(ratesData$totalGbRc, na.rm = TRUE)
    cat(sprintf("  - Gene Body Read Counts Averaged Across Trials: ~ %.0f read
        counts\n", rcGbTrialMean))
    rcLandingTrialMean <- mean(ratesData$totalLandingRc, na.rm = TRUE)
    cat(sprintf("  - Landing Pad Read Counts Averaged Across Trials: ~ %.0f
        read counts\n", rcLandingTrialMean))

    if (stericHindrance(object) && "phi" %in% colnames(ratesData)) {
        phi_mean <- mean(ratesData$phi, na.rm = TRUE)
        cat(sprintf("  - phi (landing pad occupancy): %.2f\n", phi_mean))
    }

    em_status <- ratesData$expectationMaximizationStatus
    total_trials <- length(em_status)
    max_iter_count <- sum(em_status == "max_iterations", na.rm = TRUE)
    normal_count <- sum(em_status == "normal", na.rm = TRUE)
    single_site_count <- sum(em_status == "single_site", na.rm = TRUE)
    cat("\nEM Algorithm Convergence:\n")
    cat(sprintf("  - Trials converged normally: %d/%d (%.1f%%)\n",
        normal_count, total_trials, (normal_count / total_trials) * 100
    ))
    cat(sprintf("  - Trials converged to single site: %d/%d (%.1f%%)\n",
        single_site_count, total_trials, 
        (single_site_count / total_trials) * 100
    ))
    cat(sprintf("  - Trials reached max iterations without convergence: %d/%d
    (%.1f%%)\n", max_iter_count, total_trials,
        (max_iter_count / total_trials) * 100
    ))
})

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
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, 
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 1000, polSize = 33,
#'     addSpace = 17, time = 1
#' )
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim, name="sim_beta1_k50")
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
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, 
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 1000, polSize = 33,
#'     addSpace = 17, time = 1
#' )
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim, name="sim_beta1_k50")
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
#'     k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 1950,
#'     alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500, 
#'     zetaMax = 2500, zetaVec = NULL, cellNum = 1000, polSize = 33,
#'     addSpace = 17, time = 1
#' )
#' # Estimate transcription rates
#' estRates <- estimateTranscriptionRates(sim, name="sim_beta1_k50")
#' # Get rates
#' rates <- rates(estRates)
#' # Print the rates
#' print(rates)
#' @export
setMethod("rates", "SimulationTranscriptionRates", function(object) {
    slot(object, "rates")
})