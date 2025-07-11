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
    t1H1 <- map_dbl(rc1$Yk, sum)
    t2H1 <- map_dbl(rc2$Yk, sum)
    Xk1 <- rc1$Xk
    Xk2 <- rc2$Xk
    M <- rc1$gbLength
    chiHat <- (s1 + s2) / M
    betaInt <- chiHat / (map_dbl(rc1$Xk, sum) + map_dbl(rc2$Xk, sum))
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
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
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
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
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
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
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
#' )
#' expData2 <- estimateTranscriptionRates(
#'     bigwigPlus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_plus_chr21_subset.bw",
#'     bigwigMinus = 
#'      "inst/extdata/PROseq-K562-vihervaara-treated-SE_minus_chr21_subset.bw",
#'     pauseRegions = bw_pause_21_subset,
#'     geneBodyRegions = bw_gene_body_21_subset,
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
