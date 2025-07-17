#' @importFrom dplyr mutate bind_rows bind_cols
#' @importFrom tibble tibble
#' @importFrom purrr map2 pmap map_dbl
#' @importFrom stats pchisq p.adjust
#' @importFrom methods slot is slot<- validObject
#' @title Constructor for TranscriptionRatesLRT object
#'
#' @description
#' Constructs results of likelihood ratio test comparing transcription data
#' estimated from two different sets of experimental read counts.
#'
#' @name TranscriptionRatesLRT-class
#' @rdname TranscriptionRatesLRT-class
#' @exportClass TranscriptionRatesLRT
methods::setClass("TranscriptionRatesLRT",
    slots = c(
        expData1 = "ExperimentTranscriptionRates",
        expData2 = "ExperimentTranscriptionRates",
        spikeInScalingFactor = "character",
        omegaTbl = "tbl_df",
        betaTbl = "tbl_df"
    )
)

computeOmegaLRT <- function(lambda1, lambda2, rc1, rc2) {
    tao1 <- lambda1 / (lambda1 + lambda2)
    tao2 <- 1 - tao1

    chi1 <- rc1$chi
    chi2 <- rc2$chi * lambda1 / lambda2
    p <- numeric(length(rc1$geneId)) # is this correct?

    omegaTbl <-
        tibble(
            geneId = rc1$geneId, chi1 = chi1, chi2 = chi2,
            lfc = log2(chi2 / chi1)
        )

    omegaTbl <- omegaTbl %>%
        bind_cols(bind_rows(map2(rc1$totalGbRc, rc2$totalGbRc, omegaLRT,
            tao1 = tao1, tao2 = tao2
        )))

    omegaTbl <- omegaTbl %>% mutate(padj = p.adjust(p, method = "BH"))
    return(omegaTbl)
}

computeBetaLRTParams <- function(rc1, rc2, kmin, kmax) {
    fkInt <- dnorm(kmin:kmax, mean = 50, sd = 100)
    fkInt <- fkInt / sum(fkInt)
    s1 <- rc1$totalGbRc
    s2 <- rc2$totalGbRc
    t1H1 <- map_dbl(rc1$expectedPauseSiteCounts, sum)
    t2H1 <- map_dbl(rc2$expectedPauseSiteCounts, sum)
    Xk1 <- rc1$actualPauseSiteCounts
    Xk2 <- rc2$actualPauseSiteCounts
    M <- rc1$gbLength
    chiHat <- (s1 + s2) / M
    betaInt <- chiHat / (map_dbl(rc1$actualPauseSiteCounts, sum) 
    + map_dbl(rc2$actualPauseSiteCounts, sum))
    chiHat1 <- rc1$chi
    chiHat2 <- rc2$chi
    chiHat2 <- rc2$chi

    list(
        fkInt = fkInt, s1 = s1, s2 = s2, t1H1 = t1H1, t2H1 = t2H1,
        Xk1 = Xk1, Xk2 = Xk2, M = M, chiHat = chiHat, betaInt = betaInt,
        chiHat1 = chiHat1, chiHat2 = chiHat2
    )
}

runEMH0BetaLRT <- function(params, kmin, kmax, maxItr, tor) {
    emRes <- pmap(
        list(
            params$Xk1, params$Xk2, params$betaInt, params$chiHat,
            params$chiHat1, params$chiHat2
        ),
        function(x, y, z, k, m, n) {
            tryCatch(
                mainExpectationMaximizationH0(
                    params$fkInt,
                    Xk1 = x, Xk2 = y, kmin, kmax,
                    betaInt = z, chiHat = k, chiHat1 = m, chiHat2 = n,
                    maxItr = maxItr, tor = tor
                ),
                error = function(err) {
                    list(
                        "beta" = NA, "Yk1" = NA, "Yk2" = NA, "likelihoods" =
                            list(NA)
                    )
                }
            )
        }
    )

    list(
        emRes = emRes,
        h0Likelihood = map_dbl(
            emRes,
            ~ .x$likelihoods[[length(.x$likelihoods)]]
        )
    )
}

runEMH1BetaLRT <- function(params, h0Results, kmin, kmax, maxItr, tor) {
    idx <- h0Results$tStats < 0; h0Beta <- map_dbl(h0Results$emRes, "beta")
    h0Fk1 <- map(h0Results$emRes, "fk1"); h0Fk2 <- map(h0Results$emRes, "fk2")
    emHc <- pmap(
        list(h0Fk1[idx], params$Xk1[idx], h0Beta[idx], params$chiHat1[idx]),
        function(x, y, z, k) {
            tryCatch(
                pauseEscapeEM(
                    fkInt = x, Xk = y, kmin = kmin, kmax = kmax, betaInt = z,
                    chiHat = k, maxItr = maxItr, tor = tor
                ),
                error = function(err) {
                    list(
                        "beta" = NA, "Yk" = NA, "fkMean" = NA, "fkVar" = NA,
                        "likelihoods" = list(NA)
                    )
                }
            )
        }
    )
    emHt <- pmap(
        list(h0Fk2[idx], params$Xk2[idx], h0Beta[idx], params$chiHat2[idx]),
        function(x, y, z, k) {
            tryCatch(
                pauseEscapeEM(
                    fkInt = x, Xk = y, kmin = kmin, kmax = kmax, betaInt = z,
                    chiHat = k, maxItr = maxItr, tor = tor
                ),
                error = function(err) {
                    list(
                        "beta" = NA, "Yk" = NA, "fkMean" = NA, "fkVar" = NA,
                        "likelihoods" = list(NA)
                    )
                }
            )
        }
    )
    list(
        emHc = emHc, emHt = emHt,
        h1Likelihood1 = map_dbl(
            emHc,
            ~ .x$likelihoods[[length(.x$likelihoods)]]
        ),
        h1Likelihood2 = map_dbl(
            emHt,
            ~ .x$likelihoods[[length(.x$likelihoods)]]
        )
    )
}

constructBetaLRTTable <- function(rc1, rc2, h0Results, h1Results) {
    beta1 <- rc1$betaAdp
    beta2 <- rc2$betaAdp
    tStats <- rc1$likelihood + rc2$likelihood - h0Results$h0Likelihood
    p <- numeric(length(rc1$geneId)) # is this correct?

    betaTbl <- tibble(
        geneId = rc1$geneId, beta1 = beta1, beta2 = beta2,
        lfc = log2(beta2 / beta1), fkMean1 = rc1$fkMean, fkMean2 = rc2$fkMean,
        fkVar1 = rc1$fkVar, fkVar2 = rc2$fkVar, tStats = tStats
    )

    idx <- betaTbl$tStats < 0

    betaTblIdx <- tibble(
        geneId = names(h1Results$emHc),
        beta1 = map_dbl(h1Results$emHc, "beta"),
        beta2 = map_dbl(h1Results$emHt, "beta"),
        lfc = log2(beta2 / beta1),
        fkMean1 = map_dbl(h1Results$emHc, "fkMean"),
        fkMean2 = map_dbl(h1Results$emHt, "fkMean"),
        fkVar1 = map_dbl(h1Results$emHc, "fkVar"),
        fkVar2 = map_dbl(h1Results$emHt, "fkVar"),
        tStats = h1Results$h1Likelihood1 + h1Results$h1Likelihood2 -
            h0Results$h0Likelihood[idx]
    )

    betaTbl <- bind_rows(betaTbl[!idx, ], betaTblIdx)

    betaTbl <- betaTbl %>%
        mutate(
            p = pchisq(2 * tStats,
                df = 1, ncp = 0, lower.tail = FALSE,
                log.p = FALSE
            )
        ) %>%
        mutate(padj = p.adjust(p, method = "BH"))

    return(betaTbl)
}

computeBetaLRT <- function(rc1, rc2, kmin, kmax) {
    maxItr <- 500
    tor <- 1e-6

    params <- computeBetaLRTParams(rc1, rc2, kmin, kmax)
    h0Results <- runEMH0BetaLRT(params, kmin, kmax, maxItr, tor)
    h1Results <- runEMH1BetaLRT(params, h0Results, kmin, kmax, maxItr, tor)
    betaTbl <- constructBetaLRTTable(rc1, rc2, h0Results, h1Results)

    return(betaTbl)
}

#' @rdname TranscriptionRatesLRT-class
#' @title Likelihood Ratio Test
#'
#' @description
#' Likelihood ratio test comparing aspects of transcriptional
#' dynamics of the same TU under different conditions using Likelihood Ratio
#' Test Statistics. Uses read counts and rate estimates estimated from
#' estimateTranscriptionRates. The method also requires
#' scaling factors to determine changes in χ estimates. They can be the numbers
#' of total mapped reads or spike-in reads from the samples. Likelihood ratio
#' test computes the log 2 fold change in χ estimates between conditions, the
#' log 2 fold change in beta estimates between conditions, the t-statistics for
#' the likelihood ratio tests, and the adjusted p-values based on the "BH"
#' method.
#' @param expData1 an \code{\linkS4class{ExperimentTranscriptionRates}} object
#' @param expData2 an \code{\linkS4class{ExperimentTranscriptionRates}} object
#' @param spikeInScalingFactor path to a csv file containing scale factors
#' based on total or spike-in reads
#'
#' Note: Gene body length assumed to be the same between conditions
#'
#' @return a \code{\linkS4class{TranscriptionRatesLRT}} object
#'
#' @examples
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expData1 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_treated",
#' )
#' spikeInScalingFactor <- "inst/extdata/spikein_scaling_factor.csv"
#' lrts <- likelihoodRatioTest(expData1, expData2, spikeInScalingFactor)
#' # Print the likelihood ratio test object
#' print(lrts)
#' @export
likelihoodRatioTest <- function(expData1, expData2, spikeInScalingFactor) {
    if (!is(expData1, "ExperimentTranscriptionRates")) {
        stop("expData1 must be an ExperimentTranscriptionRates object")
    }
    if (!is(expData2, "ExperimentTranscriptionRates")) {
        stop("expData2 must be an ExperimentTranscriptionRates object")
    }
    k <- 50; kmin <- 1; kmax <- 200; rnapSize <- 50; zeta <- 2000; sigP <- 0.05
    lfc1 <- 0; lfc2 <- 0; maxItr <- 500; tor <- 1e-6
    rc1 <- rates(expData1); rc2 <- rates(expData2)

    ## Get union set of genes being analyzed
    gnUnion <- intersect(rc1$geneId, rc2$geneId)
    rc1 <- rc1[match(gnUnion, rc1$geneId), ]
    rc2 <- rc2[match(gnUnion, rc2$geneId), ]

    ## Poisson-based Likelihood Ratio Tests
    ## Use # of spike-in or total # of mappable reads as scaling factor
    scaleTbl <- read.csv(spikeInScalingFactor)

    required_cols <- c("control_1", "control_2", "treated_1", "treated_2")
    missing_cols <- setdiff(required_cols, colnames(scaleTbl))
    if (length(missing_cols) > 0) {
        stop(
            "scaleTbl is missing required columns: ",
            paste(missing_cols, collapse = ", "),
            "\nExpected columns: ", paste(required_cols, collapse = ", ")
        )
    }

    ## Cancel out M and zeta since they are the same between conditions
    lambda1 <- scaleTbl$control_1 + ifelse(is.na(scaleTbl$control_2), 0,
        scaleTbl$control_2
    )
    lambda2 <- scaleTbl$treated_1 + ifelse(is.na(scaleTbl$treated_2), 0,
        scaleTbl$treated_2
    )

    omegaTbl <- computeOmegaLRT(lambda1, lambda2, rc1, rc2)
    betaTbl <- computeBetaLRT(rc1, rc2, kmin, kmax)

    return(new("TranscriptionRatesLRT",
        expData1 = expData1,
        expData2 = expData2,
        spikeInScalingFactor = spikeInScalingFactor,
        omegaTbl = omegaTbl,
        betaTbl = betaTbl
    ))
}

#' @rdname TranscriptionRatesLRT-class
#' @title Accessor for ExperimentTranscriptionRates Object
#'
#' @description
#' Accessor for the first ExperimentTranscriptionRates object from a
#' TranscriptionRatesLRT object.
#'
#' @param object a \code{\linkS4class{TranscriptionRatesLRT}} object
#'
#' @examples
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expData1 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_treated",
#' )
#' spikeInScalingFactor <- "inst/extdata/spikein_scaling_factor.csv"
#' lrts <- likelihoodRatioTest(expData1, expData2, spikeInScalingFactor)
#' expData1(lrts)
#' @export
setGeneric("expData1", function(object) standardGeneric("expData1"))
setMethod("expData1", "TranscriptionRatesLRT", function(object) {
    slot(object, "expData1")
})

#' @rdname TranscriptionRatesLRT-class
#' @title Accessor for ExperimentTranscriptionRates Object
#'
#' @description
#' Accessor for the second ExperimentTranscriptionRates object from a
#' TranscriptionRatesLRT object.
#'
#' @param object a \code{\linkS4class{TranscriptionRatesLRT}} object
#'
#' @examples
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expData1 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_treated",
#' )
#' spikeInScalingFactor <- "inst/extdata/spikein_scaling_factor.csv"
#' lrts <- likelihoodRatioTest(expData1, expData2, spikeInScalingFactor)
#' expData2(lrts)
#' @export
setGeneric("expData2", function(object) standardGeneric("expData2"))
setMethod("expData2", "TranscriptionRatesLRT", function(object) {
    slot(object, "expData2")
})

#' @rdname TranscriptionRatesLRT-class
#' @title Accessor for Spike-In Scaling Factor
#'
#' @description
#' Accessor for the spike-in scaling factor from a
#' TranscriptionRatesLRT object.
#'
#' @param object a \code{\linkS4class{TranscriptionRatesLRT}} object
#'
#' @examples
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expData1 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_treated",
#' )
#' spikeInScalingFactor <- "inst/extdata/spikein_scaling_factor.csv"
#' lrts <- likelihoodRatioTest(expData1, expData2, spikeInScalingFactor)
#' spikeInScalingFactor(lrts)
#' @export
setGeneric("spikeInScalingFactor", function(object) {
    standardGeneric("spikeInScalingFactor")
})
setMethod(
    "spikeInScalingFactor", "TranscriptionRatesLRT",
    function(object) slot(object, "spikeInScalingFactor")
)

#' @rdname TranscriptionRatesLRT-class
#' @title Accessor for Omega Table
#'
#' @description
#' Accessor for the omega table from a TranscriptionRatesLRT object.
#'
#' @param object a \code{\linkS4class{TranscriptionRatesLRT}} object
#'
#' @return tbl_df 
#'
#' @examples
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expData1 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_treated",
#' )
#' spikeInScalingFactor <- "inst/extdata/spikein_scaling_factor.csv"
#' lrts <- likelihoodRatioTest(expData1, expData2, spikeInScalingFactor)
#' omegaTbl(lrts)
#' @export
setGeneric("omegaTbl", function(object) {
    standardGeneric("omegaTbl")
})
setMethod(
    "omegaTbl", "TranscriptionRatesLRT",
    function(object) slot(object, "omegaTbl")
)


#' @rdname TranscriptionRatesLRT-class
#' @title Accessor for Beta Table
#'
#' @description
#' Accessor for the beta table from a TranscriptionRatesLRT object.
#'
#' @param object a \code{\linkS4class{TranscriptionRatesLRT}} object
#'
#' @return tbl_df 
#'
#' @examples
#' load("inst/extdata/granges_for_read_counting_chr21_subset.RData")
#' expData1 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-control-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_control",
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
#'     name = "K562_treated",
#' )
#' spikeInScalingFactor <- "inst/extdata/spikein_scaling_factor.csv"
#' lrts <- likelihoodRatioTest(expData1, expData2, spikeInScalingFactor)
#' betaTbl(lrts)
#' @export
setGeneric("betaTbl", function(object) {
    standardGeneric("betaTbl")
})
setMethod(
    "betaTbl", "TranscriptionRatesLRT",
    function(object) slot(object, "betaTbl")
)

# Plotting Utilities

#' @title Plot pause site contour map comparison between two conditions
#'
#' @description
#' Plot a contour map with mean pause site position on the x-axis and pause site
#' standard deviation on the y-axis. The plot is a comparison between two
#' conditions.
#'
#' @param object an \code{\link{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
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
#'     stericHindrance = TRUE,
#'     omegaScale = 1000,
#'     name = "K562_control",
#' )
#' plotPauseSiteContourMap(expRates, file="pause_sites_contour_map.png")
#'
#' @rdname TranscriptionRatesLRT-class
#' @export
setGeneric("plotPauseSiteContourMapTwoConditions", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotPauseSiteContourMapTwoConditions")
})

setMethod(
    "plotPauseSiteContourMapTwoConditions", "TranscriptionRatesLRT",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {

        betaTbl <- betaTbl(object)
        name1 <- slot(expData1(object), "name")
        name2 <- slot(expData2(object), "name")

        # Calculate standard deviation from variance
        betaTbl$fkSd1 <- sqrt(betaTbl$fkVar1)
        betaTbl$fkSd2 <- sqrt(betaTbl$fkVar2)

        # Determine the limits for the plot with a buffer
        buffer <- 0.05
        x_min <- min(betaTbl$fkMean1, betaTbl$fkMean2, na.rm = TRUE)
        x_max <- max(betaTbl$fkMean1, betaTbl$fkMean2, na.rm = TRUE)
        y_min <- min(betaTbl$fkSd1, betaTbl$fkSd2, na.rm = TRUE)
        y_max <- max(betaTbl$fkSd1, betaTbl$fkSd2, na.rm = TRUE)

        x_range <- x_max - x_min
        y_range <- y_max - y_min

        p <- ggplot(betaTbl) +
            geom_density_2d(aes(x = fkMean1, y = fkSd1, color = name1)) +
            geom_density_2d(aes(x = fkMean2, y = fkSd2, color = name2)) +
            geom_point(aes(x = fkMean1, y = fkSd1, color = name1), alpha = 0.6,
            size = 1.5) +
            geom_point(aes(x = fkMean2, y = fkSd2, color = name2), alpha = 0.6,
            size = 1.5) +
            labs(
                x = "Mean Pause Site Position (bp)",
                y = "Pause Site Standard Deviation (bp)",
                title = "Pause Site Mean vs Standard Deviation Distributions",
                color = "Condition"
            ) +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5)) +
            xlim(x_min - buffer * x_range, x_max + buffer * x_range) +
            ylim(0, y_max + buffer * y_range)  # Set lower limit to 0

        if (!is.null(file)) {
            ggsave(file, p,
                width = width, height = height,
                dpi = dpi
            )
        }
        return(p)
    }
)

#' @title Plot Beta Violin Plot
#'
#' @description
#' Plot a violin plot comparing the beta values between two conditions.
#'
#' @param object an \code{\link{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
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
#'     stericHindrance = TRUE,
#'     omegaScale = 1000,
#'     name = "K562_control",
#' )
#' BetaViolinPlot(expRates, file="boxplot.png")
#'
#' @rdname TranscriptionRatesLRT-class
#' @export
setGeneric("BetaViolinPlot", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("BetaViolinPlot")
})

setMethod(
    "BetaViolinPlot", "TranscriptionRatesLRT",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {

        betaTbl <- betaTbl(object)
        name1 <- slot(expData1(object), "name")
        name2 <- slot(expData2(object), "name")

        p <- betaTbl %>%
            dplyr::select(beta1, beta2) %>%
            tidyr::pivot_longer(cols = c(beta1, beta2),
                                names_to = "Condition",
                                values_to = "Beta",
                                names_prefix = "beta") %>%
            dplyr::mutate(Condition = 
            dplyr::recode(Condition, "1" = name1, "2" = name2)) %>%
            ggpubr::ggviolin(x = "Condition", y = "Beta", fill = "Condition",
                            palette = c("#00AFBB", "#E7B800"),
                            add = "boxplot", add.params = list(fill = "white"))
                            +
            labs(title = "Beta Values Comparison", x = "Condition", 
            y = "Beta Value") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5))

        if (!is.null(file)) {
            ggsave(file, p,
                width = width, height = height,
                dpi = dpi
            )
        }
        return(p)
    }
)

#' @title Plot Chi Violin Plot
#'
#' @description
#' Plot a violin plot comparing the chi values between two conditions.
#'
#' @param object an \code{\link{TranscriptionRates}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
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
#'     stericHindrance = TRUE,
#'     omegaScale = 1000,
#'     name = "K562_control",
#' )
#' ChiViolinPlot(expRates, file="boxplot.png")
#'
#' @rdname TranscriptionRatesLRT-class
#' @export
setGeneric("ChiViolinPlot", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("ChiViolinPlot")
})

setMethod(
    "ChiViolinPlot", "TranscriptionRatesLRT",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {

        omegaTbl <- omegaTbl(object)
        name1 <- slot(expData1(object), "name")
        name2 <- slot(expData2(object), "name")

        p <- omegaTbl %>%
            dplyr::select(chi1, chi2) %>%
            tidyr::pivot_longer(cols = c(chi1, chi2),
                                names_to = "Condition",
                                values_to = "Chi",
                                names_prefix = "chi") %>%
            dplyr::mutate(Condition = 
            dplyr::recode(Condition, "1" = name1, "2" = name2)) %>%
            ggpubr::ggviolin(x = "Condition", y = "Chi", fill = "Condition",
                            palette = c("#00AFBB", "#E7B800"),
                            add = "boxplot", add.params = list(fill = "white"))
                            +
            labs(title = "Chi Values Comparison", x = "Condition", 
            y = "Chi Value") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5))

        if (!is.null(file)) {
            ggsave(file, p,
                width = width, height = height,
                dpi = dpi
            )
        }
        return(p)
    }
)