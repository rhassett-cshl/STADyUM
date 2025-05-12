#' @importFrom rtracklayer import.bw
#' @import GenomicRanges
#' @import tibble
#' @importFrom methods new
#' @importFrom S4Vectors elementMetadata DataFrame splitAsList
#' @importFrom dplyr mutate select %>% left_join
#' @importFrom purrr map map_dbl
#' @importFrom stats dnorm
#' @importFrom utils write.csv
#' @import ggplot2
#' @import progress
#'
#' @title Check validity of experimentTranscriptionRates object
#' @description Validates an experimentTranscriptionRates object
#' @param object A \code{experimentTranscriptionRates} object
#' @return TRUE if valid, else errors
#' @keywords internal
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

# Helper validation functions
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

#' Class experimentTranscriptionRates
#'
#' Class \code{experimentTranscriptionRates} has read counts, pause and gene
#' body genomic region coordinates, steric hindrance and omega scale factor
#' parameters used to estimate the transcription rates
#'
#' @slot counts a \code{data.frame} with five columns gene_id,
#' summarizedPauseCounts, pauseLength, summarizedGbCounts, gbLength
#' @slot bigwigPlus a path to bigwig for plus strand
#' @slot bigwigMinus a path to bigwig for minus strand
#' @slot pauseRegions a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the pause region coordinates
#' @slot geneBodyRegions a \code{\link[GenomicRanges]{GRanges-class}}
#' that holds all the gene body region coordinates
#' @slot geneNameColumn a string for the gene name column in the GRanges
#' @slot stericHindrance a logical value representing whether landing-pad
#' occupancy was inferred when estimating the rates
#' @slot omegaScale a numeric for the scale factor used to calculate omega
#' @slot rates a \code{tbl_df} containing estimated transcription rates such as
#' chi estimates, betaOrg estimates from the initial model, betaAdp estimates
#' from the model with varying pause sites, fkMean giving the mean position
#' of pause sites, fkVar for variance of pause sites, phi estimates for
#' landing-pad occupancy, omegaZeta for the effective initiation rate,
#' betaZeta for the pause-escape rate, alphaZeta for the potential initiation
#' rate, and likelihoods
#'
#' @name experimentTranscriptionRates-class
#' @rdname experimentTranscriptionRates-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom tibble tbl_df
#' @exportClass experimentTranscriptionRates
methods::setClass("experimentTranscriptionRates",
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

#' @keywords internal
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

#' @keywords internal
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
    rc1Pause <- summariseBw(bw = bw1P3, regions = pauseRegions, 
                            "summarizedPauseCounts")    
    rc1Gb <- summariseBw(bw = bw1P3, regions = geneBodyRegions, 
                        "summarizedGbCounts")
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

#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
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
#' Estimates the transcription rates, such as initiation, pause-release rates
#' and landing pad occupancy, from experimental data, such as nascent RNA
#' sequencing read counts and genomic coordinates, and contructs an object
#' that holds these rates
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
#' @return an \code{\link{experimentTranscriptionRates-class}} object
#'
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
        Class = "experimentTranscriptionRates",
        counts = as.data.frame(rc1), bigwigPlus = bigwigPlus,
        bigwigMinus = bigwigMinus, pauseRegions = pauseRegions,
        geneBodyRegions = geneBodyRegions, geneNameColumn =
        geneNameColumn, stericHindrance = stericHindrance, omegaScale =
        omegaScale, rates = emRate
    ))
}

#' Show method for experimentTranscriptionRates objects
#' @param object An experimentTranscriptionRates object
#' @return NULL (invisibly)
#' @export
#' @examples
#' # Create an experimentTranscriptionRates object
#' expRates <- estimateExperimentTranscriptionRates(
#'     bigwigPlus = "path/to/plus.bw",
#'     bigwigMinus = "path/to/minus.bw",
#'     pauseRegions = GRanges("chr1:1-1000"),
#'     geneBodyRegions = GRanges("chr1:1-2000"),
#'     geneNameColumn = "geneId"
#' )
#'
#' # Show the object
#' show(expRates)
methods::setMethod("show",
    signature = "experimentTranscriptionRates",
    function(object) {
        cat("An experimentTranscriptionRates object with:\n")
        cat("  -", length(unique(counts(object)$geneId)), "genes\n")
        cat("  -", nrow(rates(object)), "rate estimates\n")
        cat("  - Steric hindrance:", stericHindrance(object), "\n")
        if (stericHindrance(object)) {
            cat("  - Omega scale:", omegaScale(object), "\n")
        }
    }
)

#' Export rates to CSV
#'
#' @param object An experimentTranscriptionRates object
#' @param file Path to output CSV file. Defaults to "experiment_rates.csv"
#' @return Outputs a CSV file with the rates
#' @export
setGeneric("exportRatesToCSV", function(object, file="experiment_rates.csv") {
    standardGeneric("exportRatesToCSV")
})
setMethod(
    "exportRatesToCSV", "experimentTranscriptionRates",
    function(object, file) {
        write.csv(rates(object), file = file, row.names = FALSE)
    }
)

#' @keywords internal
createScatterPlot <- function(data, rateType) {
    ggplot2::ggplot(data, ggplot2::aes(
        x = .data$betaOrg,
        y = .data[[rateType]]
    )) +
        ggplot2::geom_point(color = "#1E88E5", alpha = 0.7, size = 2) +
        ggplot2::labs(x = "Original Beta", y = rateType) +
        applyCommonTheme()
}

#' @keywords internal
createHistogramPlot <- function(data, rateType) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[rateType]])) +
        ggplot2::geom_histogram(
            bins = 30, fill = "#1E88E5",
            color = "white", alpha = 0.7
        ) +
        ggplot2::labs(x = rateType, y = "Count") +
        applyCommonTheme()
}

#' @keywords internal
createDensityPlot <- function(data, rateType) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[rateType]])) +
        ggplot2::geom_density(
            fill = "#1E88E5", color = "#0D47A1",
            alpha = 0.7
        ) +
        ggplot2::labs(x = rateType, y = "Density") +
        applyCommonTheme()
}

#' @keywords internal
applyCommonTheme <- function() {
    ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_line(color = "gray90"),
            panel.grid.minor = ggplot2::element_line(color = "gray95"),
            axis.text = ggplot2::element_text(color = "black", size = 12),
            axis.title = ggplot2::element_text(color = "black", size = 14)
        )
}

#' Plot transcription rates
#'
#' @param object An experimentTranscriptionRates object
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
#' @export
setGeneric("plotRates", function(
    object, type = "scatter", rateType = "betaAdp", file = NULL, width = 8,
    height = 6, dpi = 300, ...) {
    standardGeneric("plotRates")
})

#' @rdname experimentTranscriptionRates-class
#' @export
setMethod("plotRates", "experimentTranscriptionRates", function(
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

#' @rdname experimentTranscriptionRates-class
#' @export
setGeneric("rates", function(object) standardGeneric("rates"))
#' @rdname experimentTranscriptionRates-class
#' @export
setMethod("rates", "experimentTranscriptionRates", function(object) {
    slot(object, "rates")
})

#' @rdname experimentTranscriptionRates-class
#' @export
setGeneric("counts", function(object) standardGeneric("counts"))
#' @rdname experimentTranscriptionRates-class
#' @export
setMethod("counts", "experimentTranscriptionRates", function(object) {
    slot(object, "counts")
})

#' @rdname experimentTranscriptionRates-class
#' @export
setGeneric("pauseRegions", function(object) standardGeneric("pauseRegions"))
#' @rdname experimentTranscriptionRates-class
#' @export
setMethod("pauseRegions", "experimentTranscriptionRates", function(object) {
    slot(object, "pauseRegions")
})

#' @rdname experimentTranscriptionRates-class
#' @export
setGeneric("geneBodyRegions", function(object) {
    standardGeneric("geneBodyRegions")
})
#' @rdname experimentTranscriptionRates-class
#' @export
setMethod(
    "geneBodyRegions", "experimentTranscriptionRates",
    function(object) {
        slot(object, "geneBodyRegions")
    }
)

#' @rdname experimentTranscriptionRates-class
#' @export
setGeneric("geneNameColumn", function(object) {
    standardGeneric("geneNameColumn")
})
#' @rdname experimentTranscriptionRates-class
#' @export
setMethod(
    "geneNameColumn", "experimentTranscriptionRates",
    function(object) {
        slot(object, "geneNameColumn")
    }
)

#' @rdname experimentTranscriptionRates-class
#' @export
setGeneric("stericHindrance", function(object) {
    standardGeneric("stericHindrance")
})
#' @rdname experimentTranscriptionRates-class
#' @export
setMethod(
    "stericHindrance", "experimentTranscriptionRates",
    function(object) {
        slot(object, "stericHindrance")
    }
)

#' @rdname experimentTranscriptionRates-class
#' @export
setGeneric("omegaScale", function(object) standardGeneric("omegaScale"))

#' @rdname experimentTranscriptionRates-class
#' @export
setMethod("omegaScale", "experimentTranscriptionRates", function(object) {
    slot(object, "omegaScale")
})