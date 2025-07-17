#' @title Constructor for ExperimentTranscriptionRates object
#'
#' @description
#' Class \code{ExperimentTranscriptionRates} has read counts, pause and gene
#' body genomic region coordinates, steric hindrance and omega scale factor
#' parameters used to estimate the transcription rates as well as the estimated
#' rates
#'
#' @slot counts a \code{data.frame} with five columns gene_id,
#' summarizedPauseCounts, pauseLength, summarizedGbCounts, gbLength
#' @slot bigwigPlus a path to bigwig for plus strand recording PRO-seq read
#' counts. This can be generated with the proseq2.0 pipeline.
#' @slot bigwigMinus a path to bigwig for minus strand recording PRO-seq read
#' counts. This can be generated with the proseq2.0 pipeline.
#' @slot pauseRegions a \code{\link[GenomicRanges]{GRanges-class}} object that
#' holds the pause regions coordinates
#' @slot geneBodyRegions a \code{\link[GenomicRanges]{GRanges-class}} object
#' that holds the gene body regions coordinates
#' @slot name a character value for the name of the experiment
#' @slot stericHindrance a logical value indicating whether to infer
#' landing-pad occupancy
#' @slot omegaScale a numeric value for scaling omega, or NULL if steric
#' hindrance is disabled
#' @slot rates a \code{\link[tibble]{tbl_df}} containing the estimated rates
#' with columns:
#' \describe{
#' \item{geneId}{Character. Gene identifier}
#' \item{chi}{Numeric. RNAP density along gene body given as estimate for the
#' gene body elongation rate [RNAPs/bp]}
#' \item{betaOrg}{Numeric. Ratio of gene body RNAP density to pause region
#' RNAP density with fixed pause sites given as an estimate for the
#' pause-escape rate}
#' \item{betaAdp}{Numeric. Ratio of gene body RNAP density to pause region
#' RNAP density from adapted model which allows pause sites to vary across
#' cells given as an estimate for the pause-escape rate}
#' \item{fkMean}{Numeric. Mean position of pause sites}
#' \item{fkVar}{Numeric. Variance of pause site positions}
#' \item{phi}{Numeric. Landing-pad occupancy (only if steric hindrance is
#' enabled)}
#' \item{betaZeta}{Numeric. Pause-escape rate (only if steric hindrance is
#' enabled)}
#' \item{alphaZeta}{Numeric. Potential initiation rate (only if steric
#' hindrance is enabled)}
#' \item{omegaZeta}{Numeric. Effective initiation rate (only if steric
#' hindrance is enabled)}
#' \item{totalGbRc}{Numeric. Total gene body read counts}
#' \item{gbLength}{Numeric. Gene body length}
#' \item{expectedPauseSiteCounts}{Numeric. Expected pause site counts}
#' \item{actualPauseSiteCounts}{Numeric. Observed pause site counts}
#' }
#'
#' @name ExperimentTranscriptionRates-class
#' @rdname ExperimentTranscriptionRates-class
#' @import GenomicRanges
#' @importClassesFrom tibble tbl_df
#' @importFrom dplyr mutate select left_join rename %>%
#' @importFrom stats dnorm uniroot density
#' @importFrom methods slot is slot<- validObject
#' @importFrom S4Vectors DataFrame splitAsList
#' @importFrom tibble as_tibble
#' @importFrom rtracklayer import.bw
#' @importFrom GenomeInfoDb seqnames
#' @exportClass ExperimentTranscriptionRates
methods::setClass("ExperimentTranscriptionRates",
    slots = c(
        counts = "data.frame",
        bigwigPlus = "character",
        bigwigMinus = "character",
        pauseRegions = "GRanges",
        geneBodyRegions = "GRanges",
        stericHindrance = "logical",
        omegaScale = "ANY",
        name = "character" 
    ),
    contains = "TranscriptionRates"
)

checkForOverlappingRegions <- function(pauseRegions, gbRegions) {
    pauseOverlaps <- findOverlaps(pauseRegions, drop.self = TRUE)
    gbOverlaps <- findOverlaps(gbRegions, drop.self = TRUE)
    overlapType <- switch(paste(length(pauseOverlaps) > 0, 
    length(gbOverlaps) > 0),
        "TRUE FALSE" = "pauseRegion",
        "FALSE TRUE" = "geneBodyRegion",
        "TRUE TRUE" = "both pauseRegion and geneBodyRegion",
        "none"
    )

    if (overlapType != "none") {
        stop(sprintf("Overlapping coordinates detected in %s. Handle multiple
        isoforms by selecting a single representative isoform per gene or
        providing non-overlapping regions", overlapType))
    }
}

inputValidationChecks <- function(
    bwPlus, bwMinus, pauseRegions,
    gbRegions, stericHindrance, omegaScale) {
    if (!file.exists(bwPlus) || !file.access(bwPlus, 4) == 0) {
        stop("bigwigPlus file does not exist or is not readable")
    }
    if (!file.exists(bwMinus) || !file.access(bwMinus, 4) == 0) {
        stop("bigwigMinus file does not exist or is not readable")
    }
    if (!methods::is(pauseRegions, "GRanges") || length(pauseRegions) == 0) {
        stop("pauseRegions must be GRanges object with at least one region")
    }
    if (!methods::is(gbRegions, "GRanges") || length(gbRegions) == 0) {
        stop("geneBodyRegions must be GRanges object with at least one region")
    }
    if (!is.logical(stericHindrance)) {
        stop("stericHindrance must be a single logical value")
    }
    if (!is.null(omegaScale) && !is.numeric(omegaScale)) {
        stop("omegaScale must be NULL or a numeric value")
    }
    if (!"gene_id" %in% colnames(S4Vectors::elementMetadata(pauseRegions)) ||
        !"gene_id" %in% colnames(S4Vectors::elementMetadata(gbRegions))) {
        stop("pauseRegions or geneBodyRegions does not have gene_id column")
    }

    duplicatedPauseRegionGeneNames <-
        any(duplicated(S4Vectors::elementMetadata(pauseRegions)[, "gene_id"]))
    if (duplicatedPauseRegionGeneNames) {
        stop("Gene names must be unique in pauseRegion")
    }
    duplicatedGbRegionGeneNames <-
        any(duplicated(S4Vectors::elementMetadata(gbRegions)[, "gene_id"]))
    if (duplicatedGbRegionGeneNames) {
        stop("Gene names must be unique in geneBodyRegion")
    }
    checkForOverlappingRegions(pauseRegions, gbRegions)
    if (stericHindrance && (is.null(omegaScale) || !is.numeric(omegaScale) ||
        omegaScale <= 0)) {
        stop("For steric hindrance, omegaScale must be set to a numeric > 0")
    }
}

prepareReadCountTable <- function(bwPlus, bwMinus, pauseRegions, gbRegions,   
                                    kmax) {
    pb <- progress::progress_bar$new(
        format = "Processing [:bar] :percent eta: :eta", total = 4
    )
    message("\nImporting bigwig files..."); pb$tick(0)
    bwp1P3 <- import.bw(bwPlus)
    bwm1P3 <- import.bw(bwMinus)
    if (sum(bwp1P3$score) == 0 || sum(bwm1P3$score) == 0) {
        stop("No reads found in plus or minus strand bigwig file")
    }
    bigwigChrs <- unique(c(
        as.character(seqnames(bwp1P3)),
        as.character(seqnames(bwm1P3))
    ))
    grangesChrs <- unique(c(
        as.character(seqnames(pauseRegions)),
        as.character(seqnames(gbRegions))
    ))
    pb$tick(); message("\nProcessing plus and minus strands bigwig...")
    bwp1P3 <- processBw(bw = bwp1P3, strand = "+")
    bwm1P3 <- processBw(bw = bwm1P3, strand = "-")
    bw1P3 <- c(bwp1P3, bwm1P3)
    rm(bwp1P3, bwm1P3)
    pauseRegions <- promoters(pauseRegions, upstream = 0, downstream = kmax)
    pb$tick(); message("\nSummarizing pause and gene body regions...")
    rc1Pause <- summariseBw(
        bw = bw1P3, grng = pauseRegions,
        colName = "summarizedPauseCounts"
    )
    rc1Gb <- summariseBw(
        bw = bw1P3, grng = gbRegions,
        colName = "summarizedGbCounts"
    )
    rc1Pause$pauseLength <- kmax
    rc1Gb$gbLength <- width(gbRegions)[match(rc1Gb$gene_id, gbRegions$gene_id)]
    pb$tick(); message("\nGenerating read counts table...")
    rc1 <- Reduce(
        function(x, y) merge(x, y, by = "gene_id", all = TRUE),
        list(rc1Pause, rc1Gb)
    )
    rc1 <- rc1[!(is.na(rc1$pauseLength) | is.na(rc1$gbLength)), ]
    rc1 <- rc1[(rc1$summarizedPauseCounts > 20) &
        (rc1$summarizedGbCounts > 20), ]
    pb$tick()
    return(list(rc1 = rc1, bw1P3 = bw1P3))
}

## For each region, extract per-base signal and store position + region ID
findValidPauseIndices <- function(pauseRegions, regionsChr, starts, ends,
bwCov, geneIds) {
    validIndices <- vapply(seq_along(pauseRegions), function(i) {
        chr <- regionsChr[i]; rStart <- starts[i]; rEnd <- ends[i]
        if (!chr %in% names(bwCov)) {
            warning(sprintf("Chrom %s not in coverage for %s", chr, geneIds[i]))
            return(FALSE)
        }
        if (rStart > rEnd || rStart < 1 || rEnd > length(bwCov[[chr]])) {
            warning(sprintf(
                "Invalid region coordinates for gene %s: %s:%d-%d",
                geneIds[i], chr, rStart, rEnd
            ))
            return(FALSE)
        }
        if (length(as.numeric(bwCov[[chr]][rStart:rEnd])) == 0) {
            warning(sprintf("No signal data extracted for gene %s", geneIds[i]))
            return(FALSE)
        }
        return(TRUE)
    }, logical(1))
    if (sum(validIndices) == 0) {
        stop("No valid signal data could be extracted")
    }
    return(validIndices)
}

prepareEmData <- function(rc1, bw1P3, pauseRegions, kmin, kmax,
                            stericHindrance, omegaScale, zeta) {
    emRate <- DataFrame(
        geneId = rc1$gene_id,
        totalGbRc = rc1$summarizedGbCounts, gbLength = rc1$gbLength
    )
    emRate$chi <- emRate$totalGbRc / emRate$gbLength
    bwCov <- coverage(bw1P3, weight = "score")
    geneIds <- pauseRegions$gene_id
    strands <- as.character(strand(pauseRegions))
    regionsChr <- as.character(seqnames(pauseRegions))
    starts <- start(pauseRegions); ends <- end(pauseRegions)
    validIndices <- findValidPauseIndices(pauseRegions, regionsChr, starts,
    ends, bwCov, geneIds)
    Xk <- do.call(rbind, lapply(which(validIndices), function(i) {
        chr <- regionsChr[i]
        rStart <- starts[i]
        rEnd <- ends[i]
        regionLen <- rEnd - rStart + 1
        signal <- as.numeric(bwCov[[chr]][rStart:rEnd])
        if (strands[i] == "-") signal <- rev(signal)
        data.frame(
            region = geneIds[i], position = seq_len(regionLen),
            signal = signal
        )
    }))
    Xk <- splitAsList(Xk$signal, Xk$region)
    emRate$Xk <- Xk[emRate$geneId]
    emRate$XkSum <- vapply(emRate$Xk, sum, numeric(1))
    emRate$betaInt <- emRate$chi / emRate$XkSum
    if (stericHindrance) {
        emRate$omegaZeta <- emRate$chi * omegaScale
        emRate$omega <- emRate$omegaZeta / zeta
    }
    return(emRate)
}

experimentRunEmAlgorithm <- function(emRate, kmin, kmax, fkInt,
                                        stericHindrance, zeta, lambda = NULL) {
    emLs <- list()
    for (i in seq_len(NROW(emRate))) {
        rc <- emRate[i, ]
        if (!stericHindrance) {
            emLs[[i]] <- pauseEscapeEM(
                Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                fkInt = fkInt, betaInt = rc$betaInt[[1]],
                chiHat = rc$chi, maxItr = 500, tor = 1e-4
            )
        } else {
            emLs[[i]] <- stericHindranceEM(
                Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                f1 = 0.517, f2 = 0.024, fkInt = fkInt,
                betaInt = rc$betaInt[[1]], phiInt = 0.5,
                chiHat = rc$chi, lambda = lambda, zeta = zeta,
                maxItr = 500, tor = 1e-4
            )
        }
    }
    names(emLs) <- emRate$geneId
    return(emLs)
}

experimentProcessEmResults <- function(emRate, emLs, stericHindrance, zeta) {
    emRate$betaAdp <- map_dbl(emLs, "beta", .default = NA)
    emRate$Yk <- map(emLs, "Yk", .default = NA)
    emRate$fk <- map(emLs, "fk", .default = NA)
    emRate$fkMean <- map_dbl(emLs, "fkMean", .default = NA)
    emRate$fkVar <- map_dbl(emLs, "fkVar", .default = NA)

    emRate$t <- vapply(emRate$Yk, sum, numeric(1))
    emRate$proportionYk <- emRate$t / vapply(emRate$Xk, sum, numeric(1))
    emRate$likelihood <- map_dbl(
        emLs,
        ~ .x$likelihoods[[length(.x$likelihoods)]]
    )

    ## Convert to tibble and handle steric hindrance
    emRate <- as_tibble(emRate)
    if (stericHindrance) {
        emRate$phi <- map_dbl(emLs, "phi", .default = NA)
        emRate <- emRate %>%
            mutate(
                betaZeta = betaAdp * zeta,
                alphaZeta = omegaZeta / (1 - phi)
            )
    }

    return(emRate)
}

estimateEmRates <- function(rc1, bw1P3, pauseRegions, kmin, kmax, fkInt,
                            stericHindrance, omegaScale, zeta) {
    emRate <- prepareEmData(
        rc1, bw1P3, pauseRegions, kmin, kmax,
        stericHindrance, omegaScale, zeta
    )
    lambda <- if (stericHindrance) zeta^2 / omegaScale else NULL
    emLs <- experimentRunEmAlgorithm(
        emRate, kmin, kmax, fkInt,
        stericHindrance, zeta, lambda
    )
    emRate <- experimentProcessEmResults(emRate, emLs, stericHindrance, zeta)
    return(emRate)
}

prepareRateTable <- function(emRate, analyticalRateTbl, stericHindrance) {
    emRate <- emRate %>% left_join(analyticalRateTbl, by = "geneId")

    if (!stericHindrance) {
        emRate <- emRate %>%
            select(
                geneId, chi, betaOrg, betaAdp, fkMean, fkVar, totalGbRc,
                gbLength, Yk, Xk, likelihood
            )%>%
            rename(expectedPauseSiteCounts = Yk, actualPauseSiteCounts = Xk)
    } else {
        emRate <- emRate %>%
            select(
                geneId, chi, betaOrg, betaAdp, fkMean, fkVar, phi,
                omegaZeta, betaZeta, alphaZeta, totalGbRc, gbLength, Yk, Xk,
                likelihood
            )%>%
            rename(expectedPauseSiteCounts = Yk, actualPauseSiteCounts = Xk)
    }

    return(emRate)
}

#' @title Estimate transcription rates from real PRO-seq data
#'
#' @description
#' Estimates the transcription rates using Expectation Maximization likelihood
#' formula. Estimates initiation, pause-release rates, average read depth along
#' gene body and pause regions. If steric hindrance is enabled, also estimates
#' the landing pad occupancy. Estimated from experimental data, such as nascent
#' RNA sequencing read counts and genomic coordinates, and contructs an object
#' that holds these rates. Considers models where pause sites are fixed or
#' varied. Considers models where steric hindrance is enabled or disabled.
#'
#' @param x The path to a bigwig file from the plus strand recording PRO-seq
#' read counts
#' @param bigwigMinus the path to a bigwig file from the minus strand recording
#' PRO-seq read counts
#' @param pauseRegions a \link[GenomicRanges]{GRanges-class} object that must
#' contain a gene_id
#' @param geneBodyRegions a \link[GenomicRanges]{GRanges-class} object that
#' must contain a gene_id
#' @param stericHindrance a logical value to determine whether to infer
#' landing-pad occupancy or not. Defaults to FALSE.
#' @param omegaScale a numeric value for scaling omega. Defaults to NULL.
#' @param ... Additional arguments (not used)
#'
#' @return an \code{\link{ExperimentTranscriptionRates-class}} object
#'
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#'     stericHindrance = TRUE,
#'     omegaScale = 1000,
#' )
#'
#' @rdname ExperimentTranscriptionRates-class
#' @export
setMethod(
    "estimateTranscriptionRates", "character",
    function(
        x, bigwigMinus, pauseRegions, geneBodyRegions, name, 
        stericHindrance = FALSE, omegaScale = NULL) { 
        bigwigPlus <- x # x is the first bigwig file path
        inputValidationChecks(
            bigwigPlus, bigwigMinus, pauseRegions,
            geneBodyRegions, stericHindrance, omegaScale
        )
        ## Force copy underlying GRanges obj to prevent side effects
        pauseRegions <- GenomicRanges::makeGRangesFromDataFrame(
            data.table::copy(data.table::as.data.table(pauseRegions)),
            keep.extra.columns = TRUE
        )
        geneBodyRegions <- GenomicRanges::makeGRangesFromDataFrame(
            data.table::copy(data.table::as.data.table(geneBodyRegions)),
            keep.extra.columns = TRUE
        )

        kmin <- 1; kmax <- 200; rnapSize <- 50; zeta <- 2000
        processedData <- prepareReadCountTable(
            bigwigPlus, bigwigMinus,
            pauseRegions, geneBodyRegions, kmax
        )
        rc1 <- processedData$rc1; bw1P3 <- processedData$bw1P3

        message("estimating rates...")
        ## Initial model: Poisson-based Maximum Likelihood Estimation
        analyticalRateTbl <- tibble::tibble(
            geneId = rc1$gene_id, betaOrg =
                (rc1$summarizedGbCounts / rc1$gbLength) /
                (rc1$summarizedPauseCounts /
                    rc1$pauseLength)
        )
        fkInt <- dnorm(kmin:kmax, mean = 50, sd = 100)
        fkInt <- fkInt / sum(fkInt)
        emRate <- estimateEmRates(
            rc1, bw1P3, pauseRegions, kmin, kmax, fkInt,
            stericHindrance, omegaScale, zeta
        )
        emRate <- prepareRateTable(emRate, analyticalRateTbl, stericHindrance)

        return(methods::new(
            Class = "ExperimentTranscriptionRates",
            counts = as.data.frame(rc1), bigwigPlus = bigwigPlus,
            bigwigMinus = bigwigMinus, pauseRegions = pauseRegions,
            geneBodyRegions = geneBodyRegions, 
            stericHindrance = stericHindrance, omegaScale = omegaScale, 
            rates = emRate, name = name  
        ))
    }
)

#' @title Show method for ExperimentTranscriptionRates objects
#'
#' @description
#' Enhanced show method for ExperimentTranscriptionRates objects that displays
#' summary statistics and provides guidance on data access
#'
#' @param object An ExperimentTranscriptionRates object
#' @return NULL (invisibly)
#' @export
methods::setMethod("show", "ExperimentTranscriptionRates", function(object) {
    cat("An ExperimentTranscriptionRates object with:\n")
    cat("  -", length(unique(counts(object)$gene_id)), "genes\n")
    cat("  -", nrow(rates(object)), "rate estimates\n")
    cat("  - Steric hindrance:", stericHindrance(object), "\n")
    if (stericHindrance(object)) {
        cat("  - Omega scale:", omegaScale(object), "\n")
    }

    cat("\nSummary statistics for rate estimates across all genes/features:\n")
    ratesData <- rates(object)

    chi_mean <- mean(ratesData$chi, na.rm = TRUE)
    cat(sprintf("  - chi (gene body RNAP density): %.2f RNAPs/bp\n", chi_mean))

    betaOrg_mean <- mean(ratesData$betaOrg, na.rm = TRUE)
    cat(sprintf("  - betaOrg (ratio of gene body RNAP density to pause region
    RNAP density, fixed sites): %.4f\n", betaOrg_mean))

    betaAdp_mean <- mean(ratesData$betaAdp, na.rm = TRUE)
    cat(sprintf("  - betaAdp (ratio of gene body RNAP density to pause region
    RNAP density, adapted model): %.4f\n", betaAdp_mean))

    fkMean_mean <- mean(ratesData$fkMean, na.rm = TRUE)
    cat(sprintf(
        "  - fkMean (mean position of pause sites): ~ %.0f bp\n",
        fkMean_mean
    ))

    fkVar_mean <- mean(ratesData$fkVar, na.rm = TRUE)
    cat(sprintf(
        "  - fkVar (variance of pause site positions): %.2f bp^2\n",
        fkVar_mean
    ))

    cat("\nTo access the full rates data, use: rates(object)\n")
})
## Accessors

#' @rdname ExperimentTranscriptionRates-class
#'
#' @title Accessor for Estimated Rates
#'
#' @description
#' Accessor for the estimated rates
#'
#' @param object an \code{ExperimentTranscriptionRates} object
#' @return a \code{DataFrame} with the following columns:
#' \item{geneId}{a character vector of gene IDs}
#' \item{chi}{a numeric vector of RNAP density along gene body}
#' \item{betaOrg}{a numeric vector }
#' \item{betaAdp}{a numeric vector }
#' \item{fkMean}{a numeric vector of the mean position of pause sites}
#' \item{fkVar}{a numeric vector of the variance of pause site positions}
#' \item{phi}{a numeric vector of the landing-pad occupancy}
#' \item{betaZeta}{a numeric vector of the pause-escape rate}
#' \item{alphaZeta}{a numeric vector of the potential initiation rate}
#' \item{omega}{a numeric vector of the effective initiation rate}
#' \item{omegaZeta}{a numeric vector of the effective initiation rate}
#' \item{expectedPauseSiteCounts}{a numeric vector of the expected pause site
#' counts}
#' \item{actualPauseSiteCounts}{a numeric vector of the actual pause site
#' counts that are observed}
#' \item{likelihood}{a numeric vector of the likelihood of the model}
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#'
#' # Get the rates from the object
#' rates(expRates)
#' @export
setMethod("rates", "ExperimentTranscriptionRates", function(object) {
    slot(object, "rates")
})

#' @rdname ExperimentTranscriptionRates-class
#' @title Accessor for Read Counts
#' @description
#' Accessor for the read counts
#'
#' @param object an \code{ExperimentTranscriptionRates} object
#' @return a \code{DataFrame} with the following columns:
#' \item{geneId}{a character vector of gene IDs}
#' \item{summarizedGbCounts}{a numeric vector of the sum of gene body read
#' counts}
#' \item{gbLength}{a numeric vector of the gene body length}
#' \item{summarizedPauseCounts}{a numeric vector of the sum of pause region
#' read counts}
#' @export
setGeneric("counts", function(object) standardGeneric("counts"))
setMethod("counts", "ExperimentTranscriptionRates", function(object) {
    slot(object, "counts")
})

#' @rdname ExperimentTranscriptionRates-class
#' @title Accessor for Pause Regions Coordinates
#' @description
#' Accessor for the pause region coords \link[GenomicRanges]{GRanges-class}
#'
#' @param object an \code{ExperimentTranscriptionRates} object
#' @return a \link[GenomicRanges]{GRanges-class} object with the pause regions
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' pauseRegions(expRates)
#' @export
setGeneric("pauseRegions", function(object) standardGeneric("pauseRegions"))
setMethod("pauseRegions", "ExperimentTranscriptionRates", function(object) {
    slot(object, "pauseRegions")
})

#' @rdname ExperimentTranscriptionRates-class
#' @title Accessor for Gene Body Regions Coordinates
#' @description
#' Accessor for the gene body regions coords \link[GenomicRanges]{GRanges-class}
#'
#' @param object an \code{ExperimentTranscriptionRates} object
#' @return a \link[GenomicRanges]{GRanges-class} object with the gene body
#' regions
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' geneBodyRegions(expRates)
#' @export
setGeneric("geneBodyRegions", function(object) {
    standardGeneric("geneBodyRegions")
})
setMethod(
    "geneBodyRegions", "ExperimentTranscriptionRates",
    function(object) {
        slot(object, "geneBodyRegions")
    }
)

#' @rdname ExperimentTranscriptionRates-class
#' @title Accessor for Omega Scale
#' @description
#' Accessor for the omega scale used to scale the omega value based on prior
#' studies.
#'
#' @param object an \code{ExperimentTranscriptionRates} object
#' @return a numeric value for the scaling factor for omega
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' omegaScale(expRates)
#' @export
setGeneric("omegaScale", function(object) standardGeneric("omegaScale"))
setMethod("omegaScale", "ExperimentTranscriptionRates", function(object) {
    slot(object, "omegaScale")
})

#' @rdname ExperimentTranscriptionRates-class
#' @title Accessor for Steric Hindrance
#' @description
#' Accessor for the steric hindrance flag. If TRUE, the landing-pad
#' occupancy is inferred in the rates held in this object.
#'
#' @param object an \code{ExperimentTranscriptionRates} object
#' @return a logical value to determine whether to infer landing-pad occupancy
#' or not
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expRates <- estimateTranscriptionRates(
#'     "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' stericHindrance(expRates)
#' @export
setMethod("stericHindrance", "ExperimentTranscriptionRates", function(object) {
    slot(object, "stericHindrance")
})