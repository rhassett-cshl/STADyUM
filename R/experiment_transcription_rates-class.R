#' @importFrom dplyr mutate select left_join

experimentTranscriptionRatesValid <- function(object) {
    errors <- c(
        validateCounts(object),
        validateBigwigFiles(object),
        validateRegions(object),
        validateGeneNameColumn(object),
        validateStericHindrance(object),
        validateRates(object)
    )

    if (length(errors) == 0) TRUE else errors
}

validateCounts <- function(object) {
    errors <- character()
    if (!is.data.frame(counts(object))) {
        errors <- c(errors, "counts must be a data.frame")
    }

    return(errors)
}

validateBigwigFiles <- function(object) {
    errors <- character()
    if (!file.exists(bigwigPlus(object))) {
        errors <- c(errors, "bigwigPlus file does not exist")
    }
    if (!file.exists(bigwigMinus(object))) {
        errors <- c(errors, "bigwigMinus file does not exist")
    }
    return(errors)
}

validateRegions <- function(object) {
    errors <- character()
    if (!inherits(pauseRegions(object), "GRanges")) {
        errors <- c(errors, "pauseRegions must be a GRanges object")
    }
    if (!inherits(geneBodyRegions(object), "GRanges")) {
        errors <- c(errors, "geneBodyRegions must be a GRanges object")
    }
    return(errors)
}

validateGeneNameColumn <- function(object) {
    errors <- character()
    if (!is.character(geneNameColumn(object)) ||
        length(geneNameColumn(object)) != 1) {
        errors <- c(errors, "geneNameColumn must be a single character
        string")
    }
    return(errors)
}

validateStericHindrance <- function(object) {
    errors <- character()
    if (!is.logical(stericHindrance(object)) ||
        length(stericHindrance(object)) != 1) {
        errors <- c(errors, "stericHindrance must be a single logical value")
    }

    if (stericHindrance(object)) {
        if (!is.numeric(omegaScale(object)) ||
            length(omegaScale(object)) != 1 ||
            omegaScale(object) <= 0) {
            errors <- c(errors, "omegaScale must be a single positive numeric
            value")
        }
    }
    return(errors)
}

validateRates <- function(object) {
    errors <- character()
    if (!inherits(rates(object), "tbl_df")) {
        errors <- c(errors, "rates must be a tibble")
    }

    #required_cols <- c("gene_id", "initiation_rate", "pause_release_rate")
    #if (steric_hindrance(object)) {
    #    required_cols <- c(required_cols, "landing_pad_occupancy")
    #}
    #missing_cols <- setdiff(required_cols, colnames(rates(object)))
    #if (length(missing_cols) > 0) {
    #    errors <- c(errors, paste(
    #        "Missing required columns in rates:",
    #        paste(missing_cols, collapse = ", ")
    #    ))
    #}
    return(errors)
}

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
#' @slot pauseRegions a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the pause region coordinates for every gene to handle read counting in
#' this region
#' @slot geneBodyRegions a \code{\link[GenomicRanges]{GRanges-class}}
#' that holds all the gene body region coordinates for every gene to handle
#' read counting in this region
#' @slot geneNameColumn a string for the gene name column in the GRanges
#' @slot stericHindrance a logical value representing whether landing-pad
#' occupancy was inferred when estimating the rates. Landing pad occupancy
#' represents the probability that the landing pad required for a new
#' initiation event is already occupied by an RNAP If TRUE, the omegaScale
#' slot must be set to a numeric value greater than 0.
#' @slot omegaScale a numeric for the scale factor used to calculate omega.
#' Omega represents the effective initiation rate, the initiation rate after a
#' portion of initiation events are blocked by RNAPs in the pause region due to
#' steric hindrance. Scale factors calibrate omega based on prior knowledge
#' 
#' @slot data A \code{tbl_df} containing the estimated transcription rates:
#' \describe{
#'   \item{\eqn{\chi}}{Numeric. Maximum likelihood estimate of the average read
#' depth in the gene body region.}
#'   \item{\eqn{\beta_{org}}}{Numeric. Maximum likelihood estimate of average
#' read depth in pause region without varying pause sites}
#'   \item{\eqn{\beta_{adp}}}{Numeric. Maximum likelihood estimate of average
#' read depth in pause region with varying pause sites}
#'   \item{fkMean}{Numeric. Mean position of pause sites}
#'   \item{fkVar}{Numeric. Variance of pause sites}
#'   \item{\eqn{\phi}}{Numeric. Landing-pad occupancy estimates representing
#' probability of RNAP occupying the landing pad required for a new initiation
#' event}
#'   \item{\eqn{\omega_{\zeta}}}{Numeric. Effective initiation rate, considering
#' steric hindrance}
#'   \item{\eqn{\beta_{\zeta}}}{Numeric. Pause-escape rate}
#'   \item{\eqn{\alpha_{\zeta}}}{Numeric. Potential initiation rate}
#' }
#'
#' @name ExperimentTranscriptionRates-class
#' @rdname ExperimentTranscriptionRates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @exportClass ExperimentTranscriptionRates
methods::setClass("ExperimentTranscriptionRates",
    slots = c(
        counts = "data.frame",
        bigwigPlus = "character",
        bigwigMinus = "character",
        pauseRegions = "GRanges",
        geneBodyRegions = "GRanges",
        geneNameColumn = "character",
        stericHindrance = "logical",
        omegaScale = "ANY",
        rates = "tbl_df"
    ),  
    validity = experimentTranscriptionRatesValid
)

inputValidationChecks <- function(bigwigPlus, bigwigMinus, pauseRegions,
    geneBodyRegions, geneNameColumn, stericHindrance, omegaScale) {
    if (!file.exists(bigwigPlus) || !file.exists(bigwigMinus)) {
        stop("bigwigPlus or bigwigMinus file does not exist")
    }
    if (!file.access(bigwigPlus, 4) == 0 ||
        !file.access(bigwigMinus, 4) == 0) {
        stop("bigwigPlus or bigwigMinus file is not readable")
    }
    if (length(pauseRegions) == 0 || length(geneBodyRegions) == 0) {
        stop("pauseRegions or geneBodyRegions is empty")
    }
    if (!geneNameColumn %in%
        colnames(S4Vectors::elementMetadata(pauseRegions)) ||
        !geneNameColumn %in%
            colnames(S4Vectors::elementMetadata(geneBodyRegions))) {
        stop(sprintf("pauseRegions or geneBodyRegions does not have a column
        matching %s", geneNameColumn))
    }
    if (!is.character(
        S4Vectors::elementMetadata(pauseRegions)[, geneNameColumn]
    ) || !is.character(
        S4Vectors::elementMetadata(geneBodyRegions)[, geneNameColumn]
    )) {
        stop(sprintf("gene name column %s must be of class 'character' in
        pauseRegions and of geneBodyRegions object", geneNameColumn))
    }
    duplicatedPauseRegionGeneNames <-
        any(duplicated(
            S4Vectors::elementMetadata(pauseRegions)[,geneNameColumn
        ]))
    if (duplicatedPauseRegionGeneNames) {
        stop("One or more gene names are
        duplicated in pause region, gene names must be unique")
    }
    duplicatedGeneBodyRegionGeneNames <-
        any(duplicated(
            S4Vectors::elementMetadata(geneBodyRegions)[,geneNameColumn
        ]))
    if (duplicatedGeneBodyRegionGeneNames) {
        stop("One or more gene names are duplicated in gene body region, gene
        names must be unique")
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
            rc1Gb$geneId, geneBodyRegions$geneId)]
    pb$tick()


    message("\nGenerating read counts table...")
    rc1 <- Reduce(
        function(x, y) merge(x, y, by=geneNameColumn, all=TRUE),
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
        geneId = rc1$geneId, 
        s = rc1$summarizedGbCounts, 
        N = rc1$gbLength
    )
    
    emRate$chi <- emRate$s / emRate$N

    bwCov <- coverage(bw1P3, weight = "score")

    geneIds <- pauseRegions$geneId
    strands <- as.character(strand(pauseRegions))
    regionsChr <- as.character(seqnames(pauseRegions))
    starts <- start(pauseRegions)
    ends <- end(pauseRegions)

    ## For each region, extract per-base signal and store position + region ID
    Xk <- do.call(rbind, lapply(seq_along(pauseRegions), function(i) {
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
    names(Xk) <- pauseRegions$geneId
    
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
            emLs[[i]] <- pauseEscapeEm(
                Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                fkInt = fkInt, betaInt = rc$betaInt[[1]],
                chiHat = rc$chi, maxItr = 500, tor = 1e-4
            )
        } else {
            emLs[[i]] <- stericHindranceEm(
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

#' estimateExperimentTranscriptionRates
#' 
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
#' @param bigwigPlus the path to a bigwig file from the plus strand recording
#' PRO-seq read counts
#' @param bigwigMinus the path to a bigwig file from the minus strand recording
#' PRO-seq read counts
#' @param pauseRegions a \link[GenomicRanges]{GRanges-class} object that must
#' contain a geneId
#' @param geneBodyRegions a \link[GenomicRanges]{GRanges-class} object that
#' must contain a geneId
#' @param geneNameColumn a string that indicates which column in the GRanges
#' represents gene names information. Defaults to "gene_id"
#' @param stericHindrance a logical value to determine whether to infer
#' landing-pad occupancy or not. Defaults to FALSE.
#' @param omegaScale a numeric value for scaling omega. Defaults to NULL.
#'
#' @return an \code{\link{ExperimentTranscriptionRates-class}} object
#' 
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",    
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
#' )
#'
#' @rdname ExperimentTranscriptionRates-class
#' @export
estimateExperimentTranscriptionRates <- function(bigwigPlus, bigwigMinus,
pauseRegions, geneBodyRegions, geneNameColumn="gene_id", stericHindrance=FALSE,
omegaScale=NULL) {

    inputValidationChecks(
        bigwigPlus, bigwigMinus, pauseRegions,
        geneBodyRegions, geneNameColumn, stericHindrance, omegaScale
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
    analyticalRateTbl <- tibble::tibble(geneId = rc1$geneId, betaOrg =
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
        geneBodyRegions = geneBodyRegions, geneNameColumn =
        geneNameColumn, stericHindrance = stericHindrance, omegaScale =
        omegaScale, rates = emRate
    ))
}

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
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
#' )
#'
#' # Show the object
#' show(expRates)
methods::setMethod("show",
    signature = "ExperimentTranscriptionRates",
    function(object) {
        cat("An ExperimentTranscriptionRates object with:\n")
        cat("  -", length(unique(counts(object)$geneId)), "genes\n")
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
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
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

#' @rdname ExperimentTranscriptionRates-class
#' @export
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
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
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
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
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
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
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
#' @title Accessor for Gene Name Column
#' @description
#' Accessor for the gene name column in the GRanges object that contains gene
#' names
#'
#' @param object an \code{ExperimentTranscriptionRates} object
#' @return a string that indicates which column in the GRanges represents gene
#' names information
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
#' )    
#' geneNameColumn(expRates)
#' @export
setGeneric("geneNameColumn", function(object) {
    standardGeneric("geneNameColumn")
})
setMethod(
    "geneNameColumn", "ExperimentTranscriptionRates",
    function(object) {
        slot(object, "geneNameColumn")
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
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
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
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "gene_id"
#' )    
#' stericHindrance(expRates)
#' @export
setMethod("stericHindrance", "ExperimentTranscriptionRates", function(object) {
    slot(object, "stericHindrance")
})

