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
#' @slot stericHindrance a logical value indicating whether to infer
#' landing-pad occupancy
#' @slot omegaScale a numeric value for scaling omega, or NULL if steric
#' hindrance is disabled
#' @slot rates a \code{\link[tibble]{tbl_df}} containing the estimated rates
#' with columns:
#' \describe{
#'   \item{geneId}{Character. Gene identifier}
#'   \item{chi}{Numeric. RNAP density along gene body}
#'   \item{betaOrg}{Numeric. RNAP density along pause region}
#'   \item{betaAdp}{Numeric. RNAP density along pause region from adapted model
#' which allows pause sites to vary across cells}
#'   \item{fkMean}{Numeric. Mean position of pause sites}
#'   \item{fkVar}{Numeric. Variance of pause site positions}
#'   \item{phi}{Numeric. Landing-pad occupancy (only if steric hindrance is
#'   enabled)}
#'   \item{betaZeta}{Numeric. Pause-escape rate (only if steric hindrance is
#'   enabled)}
#'   \item{alphaZeta}{Numeric. Potential initiation rate (only if steric
#'   hindrance is enabled)}
#'   \item{omegaZeta}{Numeric. Effective initiation rate (only if steric
#'   hindrance is enabled)}
#' }
#'
#' @name ExperimentTranscriptionRates-class
#' @rdname ExperimentTranscriptionRates-class
#' @import GenomicRanges
#' @importClassesFrom tibble tbl_df
#' @importFrom dplyr mutate select left_join
#' @importFrom stats dnorm uniroot density
#' @importFrom magrittr %>%
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
        rates = "tbl_df"
    ))

inputValidationChecks <- function(bigwigPlus, bigwigMinus, pauseRegions,
    geneBodyRegions, stericHindrance, omegaScale) {
    if (!file.exists(bigwigPlus) || !file.access(bigwigPlus, 4) == 0) {
        stop("bigwigPlus file does not exist or is not readable")
    }
    if (!file.exists(bigwigMinus) || !file.access(bigwigMinus, 4) == 0) {
        stop("bigwigMinus file does not exist or is not readable")
    }
    if(!methods::is(pauseRegions, "GRanges") || length(pauseRegions) == 0) {
        stop("pauseRegions must be a GRanges object with at least one region")
    }
    if(!methods::is(geneBodyRegions, "GRanges") || length(geneBodyRegions) == 0) {
        stop("geneBodyRegions must be a GRanges object with at least one region")
    }
    if (!is.logical(stericHindrance)) {
        stop("stericHindrance must be a single logical value")
    }
    if (!is.null(omegaScale) && !is.numeric(omegaScale)) {
        stop("omegaScale must be NULL or a numeric value")
    }
    if (!"gene_id" %in%
        colnames(S4Vectors::elementMetadata(pauseRegions)) ||
        !"gene_id" %in%
            colnames(S4Vectors::elementMetadata(geneBodyRegions))) {
        stop("pauseRegions or geneBodyRegions does not have gene_id column")
    }

    duplicatedPauseRegionGeneNames <-
        any(duplicated(
            S4Vectors::elementMetadata(pauseRegions)[,"gene_id"]))
    if (duplicatedPauseRegionGeneNames) {
        stop("One or more gene names are
        duplicated in pause region, gene names must be unique")
    }
    duplicatedGeneBodyRegionGeneNames <-
        any(duplicated(
            S4Vectors::elementMetadata(geneBodyRegions)[,"gene_id"]))
    if (duplicatedGeneBodyRegionGeneNames) {
        stop("One or more gene names are duplicated in gene body region, gene
        names must be unique")
    }
    
    # Check for overlapping coordinates in GRanges objects
    # This could indicate multiple isoforms that need to be handled separately
    pauseOverlaps <- findOverlaps(pauseRegions, drop.self = TRUE)
    if (length(pauseOverlaps) > 0) {
        stop("Overlapping coordinates detected in pauseRegions. This may indicate multiple isoforms. Please handle multiple isoforms by selecting a single representative isoform per gene or by providing non-overlapping regions.")
    }
    
    geneBodyOverlaps <- findOverlaps(geneBodyRegions, drop.self = TRUE)
    if (length(geneBodyOverlaps) > 0) {
        stop("Overlapping coordinates detected in geneBodyRegions. This may 
        indicate multiple isoforms. Please handle multiple isoforms by selecting a single representative isoform per gene or by providing non-overlapping regions.")
    }
    
    if (stericHindrance && (is.null(omegaScale) || !is.numeric(omegaScale) ||
        omegaScale <= 0)) {
        stop("For steric hindrance case, omegaScale parameter must be set to
        numeric greater than 0")
    }
}

prepareReadCountTable <- function(bigwigPlus, bigwigMinus, pauseRegions,
                                    geneBodyRegions, kmax) {
    rcCutoff <- 20
    pb <- progress::progress_bar$new(
        format = "Processing [:bar] :percent eta: :eta",
        total = 4,  
    )
    
    message("\nImporting bigwig files...")
    pb$tick(0)
    bwp1P3 <- import.bw(bigwigPlus)
    bwm1P3 <- import.bw(bigwigMinus)
    if (sum(bwp1P3$score) == 0 || sum(bwm1P3$score) == 0) {
        stop("No reads found in plus or minus strand bigwig file")
    }
    # Check chromosome name compatibility between bigwig files and GRanges
    bigwigChrs <- unique(c(
        as.character(seqnames(bwp1P3)),
        as.character(seqnames(bwm1P3))
    ))
    grangesChrs <- unique(c(
        as.character(seqnames(pauseRegions)),
        as.character(seqnames(geneBodyRegions))
    ))
    missingChrs <- setdiff(grangesChrs, bigwigChrs)
    if (length(missingChrs) > 0) {
        stop("The following chromosomes in GRanges objects are not found in 
        bigwig files: ", paste(missingChrs, collapse = ", "))
    }
    
    pb$tick()

    message("\nProcessing plus and minus strands bigwig...")
    bwp1P3 <- processBw(bw = bwp1P3, strand = "+")    
    bwm1P3 <- processBw(bw = bwm1P3, strand = "-")
    bw1P3 <- c(bwp1P3, bwm1P3)
    rm(bwp1P3, bwm1P3)
    pb$tick()  

    pauseRegions <- promoters(pauseRegions, upstream = 0, downstream = kmax)

    message("\nSummarizing pause and gene body regions...")
    rc1Pause <- summariseBw(bw = bw1P3, grng = pauseRegions, 
                            colName = "summarizedPauseCounts")    
    rc1Gb <- summariseBw(bw = bw1P3, grng = geneBodyRegions, 
                        colName = "summarizedGbCounts")
    rc1Pause$pauseLength <- kmax
    rc1Gb$gbLength <- width(geneBodyRegions)[match(
            rc1Gb$gene_id, geneBodyRegions$gene_id)]
    pb$tick()


    message("\nGenerating read counts table...")
    rc1 <- Reduce(
        function(x, y) merge(x, y, by="gene_id", all=TRUE),
        list(rc1Pause, rc1Gb)
    )
    rc1 <- rc1[!(is.na(rc1$pauseLength) | is.na(rc1$gbLength)), ]
    rc1 <- rc1[(rc1$summarizedPauseCounts > rcCutoff) &
        (rc1$summarizedGbCounts > rcCutoff), ]
    pb$tick()

    return(list(rc1 = rc1, bw1P3 = bw1P3))
}

prepareEmData <- function(rc1, bw1P3, pauseRegions, kmin, kmax, 
                            stericHindrance, omegaScale, zeta) {
    emRate <- DataFrame(
        geneId = rc1$gene_id, 
        s = rc1$summarizedGbCounts, 
        N = rc1$gbLength
    )
    
    emRate$chi <- emRate$s / emRate$N

    bwCov <- coverage(bw1P3, weight = "score")

    geneIds <- pauseRegions$gene_id
    strands <- as.character(strand(pauseRegions))
    regionsChr <- as.character(seqnames(pauseRegions))
    starts <- start(pauseRegions)
    ends <- end(pauseRegions)

    ## For each region, extract per-base signal and store position + region ID
    # First, identify valid regions
    validIndices <- sapply(seq_along(pauseRegions), function(i) {
        chr <- regionsChr[i]
        regionStart <- starts[i]
        regionEnd <- ends[i]
        
        # Check if chromosome exists in coverage data
        if (!chr %in% names(bwCov)) {
            warning(sprintf("Chromosome %s not found in coverage data for gene %s", 
                          chr, geneIds[i]))
            return(FALSE)
        }
        
        # Check if region coordinates are valid
        if (regionStart > regionEnd || regionStart < 1 || 
            regionEnd > length(bwCov[[chr]])) {
            warning(sprintf("Invalid region coordinates for gene %s: %s:%d-%d", 
                          geneIds[i], chr, regionStart, regionEnd))
            return(FALSE)
        }
        
        signal <- as.numeric(bwCov[[chr]][regionStart:regionEnd])
        
        # Check if signal extraction was successful
        if (length(signal) == 0) {
            warning(sprintf("No signal data extracted for gene %s", geneIds[i]))
            return(FALSE)
        }
        
        return(TRUE)
    })
    
    # Check if we have any valid data
    if (sum(validIndices) == 0) {
        stop("No valid signal data could be extracted from any regions")
    }
    
    # Process only valid regions
    Xk <- do.call(rbind, lapply(which(validIndices), function(i) {
        chr <- regionsChr[i]
        regionStart <- starts[i]
        regionEnd <- ends[i]
        regionLen <- regionEnd - regionStart + 1
        signal <- as.numeric(bwCov[[chr]][regionStart:regionEnd])
        
        if (strands[i] == "-") signal <- rev(signal)

        data.frame(
            region = geneIds[i],
            position = seq_len(regionLen),
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
    emRate <- prepareEmData(rc1, bw1P3, pauseRegions, kmin, kmax,
                            stericHindrance, omegaScale, zeta)
    lambda <- if (stericHindrance) zeta^2 / omegaScale else NULL
    emLs <- experimentRunEmAlgorithm(emRate, kmin, kmax, fkInt,
                                    stericHindrance, zeta, lambda)
    emRate <- experimentProcessEmResults(emRate, emLs, stericHindrance, zeta)
    return(emRate)
}

prepareRateTable <- function(emRate, analyticalRateTbl, stericHindrance) {
    emRate <- emRate %>% left_join(analyticalRateTbl, by = "geneId")

    if (!stericHindrance) {
        emRate <- emRate %>%
            select(geneId, chi, betaOrg, betaAdp, fkMean, fkVar)
    } else {
        emRate <- emRate %>%
            select(
                geneId, chi, betaOrg, betaAdp, fkMean, fkVar, phi,
                omegaZeta, betaZeta, alphaZeta
            )
    }

    return(emRate)
}

#' @title Generic function for estimating transcription rates
#' 
#' @description
#' Generic function that estimates transcription rates from either simulation
#' data (SimulatePolymerase object) or experimental data (bigwig files and
#' genomic regions).
#'
#' @param x The input data (either a SimulatePolymerase object or bigwig files)
#' @param ... Additional arguments passed to the specific methods
#' @return An object containing estimated transcription rates
#' @export
setGeneric("estimateTranscriptionRates", function(x, ...) {
    standardGeneric("estimateTranscriptionRates")
})

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
#' expRates <- estimateTranscriptionRates(
#'     "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",    
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#' )
#'
#' @rdname ExperimentTranscriptionRates-class
#' @export
setMethod("estimateTranscriptionRates", "character", 
function(x, bigwigMinus, pauseRegions, geneBodyRegions, stericHindrance=FALSE, omegaScale=NULL, ...) {
    
    bigwigPlus <- x  # x is the first bigwig file path
    
    inputValidationChecks(
        bigwigPlus, bigwigMinus, pauseRegions,
        geneBodyRegions, stericHindrance, omegaScale
    )

    ## Force copy underlying GRanges obj to prevent modify in place side effects
    pauseRegions <- GenomicRanges::makeGRangesFromDataFrame(
        data.table::copy(data.table::as.data.table(pauseRegions)),
        keep.extra.columns = TRUE
    )
    geneBodyRegions <- GenomicRanges::makeGRangesFromDataFrame(
        data.table::copy(data.table::as.data.table(geneBodyRegions)),
        keep.extra.columns = TRUE
    )

    kmin <- 1; kmax <- 200; rnapSize <- 50; zeta <- 2000

    processedData <- prepareReadCountTable(bigwigPlus, bigwigMinus,
    pauseRegions, geneBodyRegions, kmax)
    rc1 <- processedData$rc1; bw1P3 <- processedData$bw1P3

    message("estimating rates...")

    ## Initial model: Poisson-based Maximum Likelihood Estimation 
    analyticalRateTbl <- tibble::tibble(geneId = rc1$gene_id, betaOrg =
    (rc1$summarizedGbCounts / rc1$gbLength) / (rc1$summarizedPauseCounts /
    rc1$pauseLength))

    fkInt <- dnorm(kmin:kmax, mean = 50, sd = 100)
    fkInt <- fkInt / sum(fkInt)

    emRate <- estimateEmRates(rc1, bw1P3, pauseRegions, kmin, kmax, fkInt,
    stericHindrance, omegaScale, zeta)
    emRate <- prepareRateTable(emRate, analyticalRateTbl, stericHindrance)

    return(methods::new(
        Class = "ExperimentTranscriptionRates",
        counts = as.data.frame(rc1), bigwigPlus = bigwigPlus,
        bigwigMinus = bigwigMinus, pauseRegions = pauseRegions,
        geneBodyRegions = geneBodyRegions, stericHindrance = stericHindrance,
        omegaScale = omegaScale, rates = emRate
    ))
})

#' @title Show method for ExperimentTranscriptionRates
#' @description
#' Custom display of \code{\link{ExperimentTranscriptionRates-class}} objects
#' showing summary statistics
#' 
#' @param object An ExperimentTranscriptionRates object
#' @return NULL (invisibly)
#' @rdname ExperimentTranscriptionRates-class
#' @export
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' expRates <- estimateTranscriptionRates(
#'     "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#' )
#'
#' # Show the object
#' show(expRates)
methods::setMethod("show",
    signature = "ExperimentTranscriptionRates",
    function(object) {
        cat("An ExperimentTranscriptionRates object with:\n")
        cat("  -", length(unique(counts(object)$gene_id)), "genes\n")
        cat("  -", nrow(rates(object)), "rate estimates\n")
        cat("  - Steric hindrance:", stericHindrance(object), "\n")
        if (stericHindrance(object)) {
            cat("  - Omega scale:", omegaScale(object), "\n")
        }
    }
)

createScatterPlot <- function(data, rateType) {
    ggplot2::ggplot(data, ggplot2::aes(
        x = .data$betaOrg,
        y = .data[[rateType]]
    )) +
        ggplot2::geom_point(color = "#1E88E5", alpha = 0.7, size = 2) +
        ggplot2::labs(x = "Original Beta", y = rateType) +
        applyCommonTheme()
}

createHistogramPlot <- function(data, rateType) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[rateType]])) +
        ggplot2::geom_histogram(
            bins = 30, fill = "#1E88E5",
            color = "white", alpha = 0.7
        ) +
        ggplot2::labs(x = rateType, y = "Count") +
        applyCommonTheme()
}

createDensityPlot <- function(data, rateType) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[rateType]])) +
        ggplot2::geom_density(
            fill = "#1E88E5", color = "#0D47A1",
            alpha = 0.7
        ) +
        ggplot2::labs(x = rateType, y = "Density") +
        applyCommonTheme()
}

applyCommonTheme <- function() {
    ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_line(color = "gray90"),
            panel.grid.minor = ggplot2::element_line(color = "gray95"),
            axis.text = ggplot2::element_text(color = "black", size = 12),
            axis.title = ggplot2::element_text(color = "black", size = 14)
        )
}

#' @title Plot transcription rates
#'
#' @description
#' Plots the transcription rates using ggplot2.
#'
#' @param object An ExperimentTranscriptionRates object
#' @param type Type of plot to create ("scatter", "histogram", or "density").
#' Defaults to "scatter"
#' @param rateType Which rate to plot ("betaOrg", "betaAdp", "chi", etc.).
#' Defaults to "betaAdp"
#' @param file Optional path to save the plot. If provided, the plot will be
#' saved to this location.
#' @param width Width of the saved plot in inches. Default is 8.
#' @param height Height of the saved plot in inches. Default is 6.
#' @param dpi Resolution of the saved plot. Default is 300.
#' @param ... Additional arguments passed to the plotting function
#' @return A ggplot object
#' @rdname ExperimentTranscriptionRates-class
#' 
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' expRates <- estimateTranscriptionRates(
#'     "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#' )
#'
#' # Plot rates as a histogram
#' plotRates(expRates, type = "histogram")

#' @export
setGeneric("plotRates", function(
    object, type = "scatter", rateType = "betaAdp", file = NULL, width = 8,
    height = 6, dpi = 300, ...) {
    standardGeneric("plotRates")
})

setMethod("plotRates", "ExperimentTranscriptionRates", function(
    object, type = "scatter", rateType = "betaAdp", file = NULL, width = 8,
    height = 6, dpi = 300, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    if (!rateType %in% colnames(rates(object))) {
        stop(sprintf("rateType %s not found in rates data", rateType))
    }

    data <- rates(object)

    p <- switch(type,
        scatter = createScatterPlot(data, rateType),
        histogram = createHistogramPlot(data, rateType),
        density = createDensityPlot(data, rateType),
        stop("Invalid plot type. Choose from 'scatter', 'histogram', or
        'density'")
    )

    if (!is.null(file)) {
        ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi)
    }

    return(p)
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
#' @export
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' expRates <- estimateTranscriptionRates(
#'     "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#' )
#' 
#' # Get the rates from the object
#' rates(expRates)
setGeneric("rates", function(object) standardGeneric("rates"))
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
#' expRates <- estimateTranscriptionRates(
#'     "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
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
#' expRates <- estimateTranscriptionRates(
#'     "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
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
#' expRates <- estimateTranscriptionRates(
#'     "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
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
#' expRates <- estimateTranscriptionRates(
#'     "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#' )    
#' stericHindrance(expRates)
#' @export
setMethod("stericHindrance", "ExperimentTranscriptionRates", function(object) {
    slot(object, "stericHindrance")
})

