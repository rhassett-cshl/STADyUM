#' @importFrom dplyr mutate bind_rows bind_cols case_when
#' @importFrom tibble tibble
#' @importFrom purrr map2 pmap map_dbl
#' @importFrom stats pchisq p.adjust
#' @importFrom methods slot is slot<- validObject
#' @importFrom utils read.csv
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
        transcriptionRates1 = "TranscriptionRates",
        transcriptionRates2 = "TranscriptionRates",
        spikeInFile = "ANY",
        chiTbl = "tbl_df",
        betaTbl = "tbl_df",
        fkTbl  = "tbl_df"
    )
)

computeChiLRT <- function(lambda1, lambda2, rc1, rc2, isExperiment) {
    tao1 <- lambda1 / (lambda1 + lambda2)
    tao2 <- 1 - tao1

    chi1 <- rc1$chi
    chi2 <- rc2$chi * lambda1 / lambda2
    if(isExperiment) {
        p <- numeric(length(rc1$geneId)) 
    } else {
        p <- numeric(length(rc1$trial))
    }

    if(isExperiment) {
        chiTbl <- tibble(
            geneId = rc1$geneId, chi1 = chi1, chi2 = chi2,
            lfc = log2(chi2 / chi1)
        )
    } else {
        chiTbl <- tibble(
            trial = rc1$trial, chi1 = chi1, chi2 = chi2,
            lfc = log2(chi2 / chi1)
        )
    }

    chiTbl <- chiTbl %>%
        bind_cols(bind_rows(map2(rc1$totalGbRc, rc2$totalGbRc, chiLRT,
        tao1 = tao1, tao2 = tao2)))

    chiTbl <- chiTbl %>% mutate(padj = p.adjust(p, method = "BH"))
    return(chiTbl)
}

computeBetaLRTParams <- function(rc1, rc2, scaleFactor, kmin, kmax, gbLength) {
    fkInt <- dnorm(kmin:kmax, mean = 50, sd = 100)
    fkInt <- fkInt / sum(fkInt)
    s1 <- rc1$totalGbRc
    s2 <- rc2$totalGbRc
    t1H1 <- map_dbl(rc1$expectedPauseSiteCounts, sum)
    t2H1 <- map_dbl(rc2$expectedPauseSiteCounts, sum)
    Xk1 <- rc1$actualPauseSiteCounts
    Xk2 <- rc2$actualPauseSiteCounts
    M <- gbLength
    chiHat <- (s1 + s2*scaleFactor) / M
    betaInt <- chiHat / (map_dbl(rc1$actualPauseSiteCounts, sum) 
    + map_dbl(rc2$actualPauseSiteCounts, sum))
    chiHat1 <- rc1$chi
    chiHat2 <- rc2$chi
    rc1Likelihood <- rc1$likelihood
    rc2Likelihood <- rc2$likelihood

    list(
        fkInt = fkInt, s1 = s1, s2 = s2, t1H1 = t1H1, t2H1 = t2H1,
        Xk1 = Xk1, Xk2 = Xk2, M = M, chiHat = chiHat, betaInt = betaInt,
        chiHat1 = chiHat1, chiHat2 = chiHat2,
        rc1Likelihood = rc1Likelihood, rc2Likelihood = rc2Likelihood
    )
}

computeFkLRTParams <- function(rc1, rc2, scaleFactor, kmin, kmax, gbLength) {
    fkInt <- dnorm(kmin:kmax, mean = 50, sd = 100)
    fkInt <- fkInt / sum(fkInt)
    s1 <- rc1$totalGbRc
    s2 <- rc2$totalGbRc
    t1H1 <- map_dbl(rc1$expectedPauseSiteCounts, sum)
    t2H1 <- map_dbl(rc2$expectedPauseSiteCounts, sum)
    Xk1 <- rc1$actualPauseSiteCounts
    Xk2 <- rc2$actualPauseSiteCounts
    M <- gbLength
    chiHat <- (s1 + s2*scaleFactor) / M
    chiHat1 <- rc1$chi
    chiHat2 <- rc2$chi

    betaInt1 <- rc1$betaAdp
    betaInt2 <- rc2$betaAdp
    rc1Likelihood <- rc1$likelihood
    rc2Likelihood <- rc2$likelihood

    list(
        fkInt = fkInt, s1 = s1, s2 = s2, t1H1 = t1H1, t2H1 = t2H1,
        Xk1 = Xk1, Xk2 = Xk2, M = M, chiHat = chiHat, betaInt1 = betaInt1, betaInt2 = betaInt2,
        chiHat1 = chiHat1, chiHat2 = chiHat2,
        rc1Likelihood = rc1Likelihood, rc2Likelihood = rc2Likelihood
    )
}

runEMH0BetaLRT <- function(params, kmin, kmax, scaleFactor, maxItr, tor) {
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
                    betaInt = z, chiHat = k, chiHat1 = m, chiHat2 = n, scaleFactor = scaleFactor, maxItr = maxItr, tor = tor
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

runEMH0FkLRT <- function(params, kmin, kmax, scaleFactor, maxItr, tor) {
    emRes <- pmap(
        list(
            params$Xk1, params$Xk2, params$betaInt1, params$betaInt2, params$chiHat,
            params$chiHat1, params$chiHat2
        ),
        function(x, y, z1, z2, k, m, n) {
            tryCatch(
                mainExpectationMaximizationFkH0(
                    params$fkInt,
                    Xk1 = x, Xk2 = y, kmin, kmax,
                    betaInt1 = z1, betaInt2 = z2, chiHat = k, chiHat1 = m, chiHat2 = n, scaleFactor = scaleFactor, maxItr = maxItr, tor = tor
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

runEMH1BetaLRT <- function(params, h0Results, kmin, kmax, maxItr, tor, scaleFactor) {
    tStats <- params$rc1Likelihood + params$rc2Likelihood * scaleFactor - h0Results$h0Likelihood
    idx <- tStats < 0; h0Beta <- map_dbl(h0Results$emRes, "beta")
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

runEMH1FkLRT <- function(params, h0Results, kmin, kmax, maxItr, tor, scaleFactor) {
    tStats <- params$rc1Likelihood + params$rc2Likelihood * scaleFactor - h0Results$h0Likelihood
    h0Beta1 <- map_dbl(h0Results$emRes, "beta1");  h0Beta2 <- map_dbl(h0Results$emRes, "beta2")
    idx <- tStats < 0; h0Fk <- map(h0Results$emRes, "fk")
    
    emHc <- pmap(
        list(h0Fk[idx], params$Xk1[idx], h0Beta1[idx], params$chiHat1[idx]),
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
        list(h0Fk[idx], params$Xk2[idx], h0Beta2[idx], params$chiHat2[idx]),
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


constructBetaLRTTable <- function(rc1, rc2, scaleFactor, h0Results, h1Results,
isExperiment)
{
    beta1 <- rc1$betaAdp; beta2 <- rc2$betaAdp
    tStats <- rc1$likelihood + rc2$likelihood * scaleFactor - h0Results$h0Likelihood
    if(isExperiment) p <- numeric(length(rc1$geneId))
    else p <- numeric(length(rc1$trial))

    if(isExperiment) {
        betaTbl <- tibble(geneId = rc1$geneId, beta1 = beta1, beta2 = beta2,
            lfc = log2(beta2 / beta1), fkMean1 = rc1$fkMean, 
            fkMean2 = rc2$fkMean, fkVar1 = rc1$fkVar, fkVar2 = rc2$fkVar,
            tStats = tStats)
    } else {
        betaTbl <- tibble(trial = rc1$trial, beta1 = beta1, beta2 = beta2,
            lfc = log2(beta2 / beta1), fkMean1 = rc1$fkMean, 
            fkMean2 = rc2$fkMean, fkVar1 = rc1$fkVar, fkVar2 = rc2$fkVar,
            tStats = tStats)
    }

    idx <- betaTbl$tStats < 0
    betaTblIdx <- tibble(
        geneId = rc1$geneId[idx],
        beta1 = map_dbl(h1Results$emHc, "beta"),
        beta2 = map_dbl(h1Results$emHt, "beta"),
        lfc = log2(beta2 / beta1),
        fkMean1 = map_dbl(h1Results$emHc, "fkMean"),
        fkMean2 = map_dbl(h1Results$emHt, "fkMean"),
        fkVar1 = map_dbl(h1Results$emHc, "fkVar"),
        fkVar2 = map_dbl(h1Results$emHt, "fkVar"),
        tStats = h1Results$h1Likelihood1 + h1Results$h1Likelihood2 * scaleFactor - h0Results$h0Likelihood[idx]
    )

    betaTbl <- bind_rows(betaTbl[!idx, ], betaTblIdx)
    betaTbl <- betaTbl %>%
        mutate(p = pchisq(2 * tStats, df = 1, ncp = 0, lower.tail = FALSE,
                log.p = FALSE)
        ) %>%
        mutate(padj = p.adjust(p, method = "BH"))

    return(betaTbl)
}

constructFkLRTTable <- function(rc1, rc2, scaleFactor, h0Results, h1Results,
isExperiment)
{
    tStats <- rc1$likelihood + rc2$likelihood * scaleFactor - h0Results$h0Likelihood
    if(isExperiment) p <- numeric(length(rc1$geneId))
    else p <- numeric(length(rc1$trial))

    if(isExperiment) {
        fkTbl <- tibble(geneId = rc1$geneId, 
            fkMean1 = rc1$fkMean, 
            fkMean2 = rc2$fkMean, fkVar1 = rc1$fkVar, fkVar2 = rc2$fkVar,
            tStats = tStats)
    } else {
        fkTbl <- tibble(trial = rc1$trial,
            fkMean1 = rc1$fkMean, 
            fkMean2 = rc2$fkMean, fkVar1 = rc1$fkVar, fkVar2 = rc2$fkVar,
            tStats = tStats)
    }

    idx <- fkTbl$tStats < 0
    fkTblIdx <- tibble(
        geneId = rc1$geneId[idx],
        fkMean1 = map_dbl(h1Results$emHc, "fkMean"),
        fkMean2 = map_dbl(h1Results$emHt, "fkMean"),
        fkVar1 = map_dbl(h1Results$emHc, "fkVar"),
        fkVar2 = map_dbl(h1Results$emHt, "fkVar"),
        tStats = h1Results$h1Likelihood1 + h1Results$h1Likelihood2 * scaleFactor - h0Results$h0Likelihood[idx]
    )
    

    fkTbl <- bind_rows(fkTbl[!idx, ], fkTblIdx)
    fkTbl <- fkTbl %>%
        mutate(p = pchisq(2 * tStats, df = 2, ncp = 0, lower.tail = FALSE,
                log.p = FALSE)
        ) %>%
        mutate(padj = p.adjust(p, method = "BH")
               )
    
    return(fkTbl)
}



computeBetaLRT <- function(rc1, rc2, scaleFactor, kmin, kmax, gbLength, isExperiment) {
    maxItr <- 500
    tor <- 1e-6

    params <- computeBetaLRTParams(rc1, rc2, scaleFactor, kmin, kmax, gbLength)
    h0Results <- runEMH0BetaLRT(params, kmin, kmax, scaleFactor, maxItr, tor)
    h1Results <- runEMH1BetaLRT(params, h0Results, kmin, kmax, maxItr, tor, scaleFactor)
    betaTbl <- constructBetaLRTTable(rc1, rc2, scaleFactor, h0Results, h1Results, isExperiment)

    return(betaTbl)
}

computeFkLRT <- function(rc1, rc2, scaleFactor, kmin, kmax, gbLength, isExperiment) {
    maxItr <- 500
    tor <- 1e-6

    params <- computeFkLRTParams(rc1, rc2, scaleFactor, kmin, kmax, gbLength)
    h0Results <- runEMH0FkLRT(params, kmin, kmax, scaleFactor, maxItr, tor)
    h1Results <- runEMH1FkLRT(params, h0Results, kmin, kmax, maxItr, tor, scaleFactor)
    FkTbl <- constructFkLRTTable(rc1, rc2, scaleFactor, h0Results, h1Results, isExperiment)

    return(FkTbl)
}

#' @rdname TranscriptionRatesLRT-class
#' @title Likelihood Ratio Test
#'
#' @description
#' Likelihood ratio test comparing aspects of transcriptional
#' dynamics of the same TU under different conditions using Likelihood Ratio
#' Test Statistics. Uses read counts and rate estimates estimated from
#' estimateTranscriptionRates. The method also requires
#' scaling factors to determine changes in \eqn{\chi} estimates. They can be the numbers
#' of total mapped reads or spike-in reads from the samples. Likelihood ratio
#' test computes the log 2 fold change in \eqn{\chi} estimates between conditions, the
#' log 2 fold change in beta estimates between conditions, the t-statistics for
#' the likelihood ratio tests, and the adjusted p-values based on the "BH"
#' method.
#' @param transcriptionRates1 an \code{\linkS4class{TranscriptionRates}} object
#' @param transcriptionRates2 an \code{\linkS4class{TranscriptionRates}} object
#' @param scaleFactor a numeric value to scale the beta estimates. Defaults to 1
#' @param spikeInFile path to a csv file containing scale factors
#' based on total or spike-in reads or NULL if not provided. Defaults to NULL.
#'
#'
#' @return a \code{\linkS4class{TranscriptionRatesLRT}} object
#'
#' @examples
#' load(system.file("extdata", "granges_for_read_counting_DLD1_chr21.RData",
#' package = "STADyUM"))
#' transcriptionRates1 <- estimateTranscriptionRates(system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
#' bigwigMinus = system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw", package = "STADyUM"),
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control"
#' )
#' transcriptionRates2 <- estimateTranscriptionRates(system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
#' bigwigMinus = system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw", package = "STADyUM"),
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Treated"
#' )
#' spikeInFile <- system.file("extdata", "spikein_scaling_factor.csv", 
#' package = "STADyUM")
#' lrts <- likelihoodRatioTest(transcriptionRates1, transcriptionRates2,
#' scaleFactor=1, spikeInFile)
#' # Print the likelihood ratio test object
#' print(lrts)
#' @export
likelihoodRatioTest <- function(transcriptionRates1, transcriptionRates2, scaleFactor=1, spikeInFile = NULL) {
    if (!is(transcriptionRates1, "TranscriptionRates") ||
        !is(transcriptionRates2, "TranscriptionRates")) {
        stop("transcriptionRates1 and transcriptionRates2 must be an
            TranscriptionRates object")
    }
    k <- 50; rnapSize <- 50; zeta <- 2000; sigP <- 0.05; lfc1 <- 0; lfc2 <- 0;
    maxItr <- 500; tor <- 1e-6; rc1 <- rates(transcriptionRates1); 
    rc2 <- rates(transcriptionRates2)
    kmin <- 1; kmax <- length(rc1$actualPauseSiteCounts[[1]])
    isExperiment <- is(transcriptionRates1, "ExperimentTranscriptionRates")
    if(isExperiment) {
        gnUnion <- intersect(rc1$geneId, rc2$geneId)
        rc1 <- rc1[match(gnUnion, rc1$geneId), ]
        rc2 <- rc2[match(gnUnion, rc2$geneId), ]
        gbLength <- rc1$gbLength
    }
    else {
        simpol1_params <- parameters(simpol(transcriptionRates1))
        simpol2_params <- parameters(simpol(transcriptionRates2))
        gbLength <- simpol1_params$geneLen - simpol1_params$kMax
        if(gbLength != simpol2_params$geneLen - simpol2_params$kMax) {
            stop("Gene body length is not the same between conditions")
        }
    }
    if(!is.null(spikeInFile)) {
        scaleTbl <- read.csv(spikeInFile)
        required_cols <- c("control_1", "control_2", "treated_1", "treated_2")
        missing_cols <- setdiff(required_cols, colnames(scaleTbl))
        if (length(missing_cols) > 0) {
            stop("scaleTbl is missing required columns: ",
                paste(missing_cols, collapse = ", "),
                "\nExpected columns: ", paste(required_cols, collapse = ", "))
        }
    }
    else {
        scaleTbl <- tibble(control_1 = 1, control_2 = 1, treated_1 = 1,
        treated_2 = 1)
    }
    lambda1 <- scaleTbl$control_1 + scaleTbl$control_2
    lambda2 <- scaleTbl$treated_1 + scaleTbl$treated_2
    chiTbl <- computeChiLRT(lambda1, lambda2, rc1, rc2, isExperiment)
    betaTbl <- computeBetaLRT(rc1, rc2, scaleFactor, kmin, kmax, gbLength, isExperiment)
    fkTbl <- computeFkLRT(rc1, rc2, scaleFactor, kmin, kmax, gbLength, isExperiment)
    return(new("TranscriptionRatesLRT",
        transcriptionRates1 = transcriptionRates1,
        transcriptionRates2 = transcriptionRates2, chiTbl = chiTbl,
        spikeInFile = spikeInFile, betaTbl = betaTbl, fkTbl = fkTbl
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
#' @export
setGeneric("transcriptionRates1", 
function(object) standardGeneric("transcriptionRates1"))
#' @rdname TranscriptionRatesLRT-class
setMethod("transcriptionRates1", "TranscriptionRatesLRT", function(object) {
    slot(object, "transcriptionRates1")
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
#' @export
setGeneric("transcriptionRates2", 
function(object) standardGeneric("transcriptionRates2"))
#' @rdname TranscriptionRatesLRT-class
setMethod("transcriptionRates2", "TranscriptionRatesLRT", function(object) {
    slot(object, "transcriptionRates2")
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
#' @export
setGeneric("spikeInFile", function(object) {
    standardGeneric("spikeInFile")
})
#' @rdname TranscriptionRatesLRT-class
setMethod(
    "spikeInFile", "TranscriptionRatesLRT",
    function(object) slot(object, "spikeInFile")
)

#' @rdname TranscriptionRatesLRT-class
#' @title Accessor for Chi Table
#'
#' @description
#' Accessor for the chi table from a TranscriptionRatesLRT object.
#'
#' @param object a \code{\linkS4class{TranscriptionRatesLRT}} object
#'
#' @return tbl_df 
#'
#' @export
setGeneric("chiTbl", function(object) {
    standardGeneric("chiTbl")
})
#' @rdname TranscriptionRatesLRT-class
setMethod(
    "chiTbl", "TranscriptionRatesLRT",
    function(object) slot(object, "chiTbl")
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
#' @export
setGeneric("betaTbl", function(object) {
    standardGeneric("betaTbl")
})
#' @rdname TranscriptionRatesLRT-class
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
#' @param object an \code{\linkS4class{TranscriptionRatesLRT}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @examples
#' # Create an ExperimentTranscriptionRates object
#' load(system.file("extdata", "granges_for_read_counting_DLD1_chr21.RData",
#' package = "STADyUM"))
#' transcriptionRates1 <- estimateTranscriptionRates(system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
#' bigwigMinus = system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw", package = "STADyUM"),
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Control"
#' )
#' transcriptionRates2 <- estimateTranscriptionRates(system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
#' bigwigMinus = system.file("extdata",
#' "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw", package = "STADyUM"),
#'     pauseRegions = bw_pause_filtered,
#'     geneBodyRegions = bw_gb_filtered,
#'     name = "Treated"
#' )
#' spikeInFile <- system.file("extdata", "spikein_scaling_factor.csv", 
#' package = "STADyUM")
#' lrts <- likelihoodRatioTest(transcriptionRates1, transcriptionRates2,
#' scaleFactor=1, spikeInFile)
#' plotPauseSiteContourMapTwoConditions(lrts,
#' file="pause_sites_contour_map.png")
#'
#' @rdname TranscriptionRatesLRT-class
#' @export
setGeneric("plotPauseSiteContourMapTwoConditions", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotPauseSiteContourMapTwoConditions")
})
#' @rdname TranscriptionRatesLRT-class
setMethod(
    "plotPauseSiteContourMapTwoConditions", "TranscriptionRatesLRT",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {

        betaTbl <- betaTbl(object)
        name1 <- slot(transcriptionRates1(object), "name")
        name2 <- slot(transcriptionRates2(object), "name")

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
            geom_density_2d(aes(x = .data$fkMean1, y = .data$fkSd1, color = name1)) +
            geom_density_2d(aes(x = .data$fkMean2, y = .data$fkSd2, color = name2)) +
            geom_point(aes(x = .data$fkMean1, y = .data$fkSd1, color = name1), alpha = 0.6,
            size = 1.5) +
            geom_point(aes(x = .data$fkMean2, y = .data$fkSd2, color = name2), alpha = 0.6,
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

#' @rdname TranscriptionRatesLRT-class
#' @title BetaViolinPlot
#'
#' @description
#' Plot a violin plot comparing the beta values between two conditions.
#'
#' @param object an \code{\linkS4class{TranscriptionRatesLRT}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @export
setGeneric("BetaViolinPlot", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("BetaViolinPlot")
})
#' @rdname TranscriptionRatesLRT-class
setMethod(
    "BetaViolinPlot", "TranscriptionRatesLRT",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {

        betaTbl <- betaTbl(object)
        name1 <- slot(transcriptionRates1(object), "name")
        name2 <- slot(transcriptionRates2(object), "name")

        p <- betaTbl %>%
            dplyr::select(.data$beta1, .data$beta2) %>%
            tidyr::pivot_longer(cols = c(.data$beta1, .data$beta2),
                                names_to = "Condition",
                                values_to = "Beta",
                                names_prefix = "beta") %>%
            dplyr::mutate(Condition = 
            dplyr::recode(.data$Condition, "1" = name1, "2" = name2),
                ScaledBeta = .data$Beta * 2000,  # Scale by zeta
                LogScaledBeta = log(.data$ScaledBeta)) %>%  # Log of scaled beta
            ggplot(aes(x = .data$Condition, y = .data$LogScaledBeta, fill = .data$Condition)) +
            geom_violin(trim = FALSE) +
            geom_boxplot(width = 0.1, fill = "white") +
            scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
            labs(title = "Log of Scaled Beta Values Comparison", 
            x = "Condition", 
            y = expression(log(beta * zeta))) +  
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
#' @param object an \code{\linkS4class{TranscriptionRatesLRT}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @rdname TranscriptionRatesLRT-class
#' @export
setGeneric("ChiViolinPlot", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("ChiViolinPlot")
})
#' @rdname TranscriptionRatesLRT-class
setMethod(
    "ChiViolinPlot", "TranscriptionRatesLRT",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {

        chiTbl <- chiTbl(object)
        name1 <- slot(transcriptionRates1(object), "name")
        name2 <- slot(transcriptionRates2(object), "name")

        p <- chiTbl %>%
            dplyr::select(.data$chi1, .data$chi2) %>%
            tidyr::pivot_longer(cols = c(.data$chi1, .data$chi2),
                                names_to = "Condition",
                                values_to = "Chi",
                                names_prefix = "chi") %>%
            dplyr::mutate(Condition = 
            dplyr::recode(.data$Condition, "1" = name1, "2" = name2),
                ScaledChi = .data$Chi * 2000,  # Scale by zeta
                LogScaledChi = log(.data$ScaledChi)) %>%  # Log of scaled chi
            ggplot(aes(x = .data$Condition, y = .data$LogScaledChi, fill = .data$Condition)) +
            geom_violin(trim = FALSE) +
            geom_boxplot(width = 0.1, fill = "white") +
            scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
            labs(title = "Log of Scaled Chi Values Comparison", 
            x = "Condition", 
            y = expression(log(chi * zeta))) +  
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

#' @title MA plot for log fold change treated/control vs log mean beta*zeta
#'
#' @description
#' Plot a MA plot for log fold change treated/control vs log mean beta*zeta.
#'
#' @param object an \code{\linkS4class{TranscriptionRatesLRT}} object
#' @param file the path to a file to save the plot to
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @param dpi the resolution of the plot in dpi
#'
#' @return an \code{\link{ggplot2}} object
#'
#' @rdname TranscriptionRatesLRT-class
#' @export
setGeneric("plotLfcMa", function(
    object, file = NULL, width = 8,
    height = 6, dpi = 300) {
    standardGeneric("plotLfcMa")
})
#' @rdname TranscriptionRatesLRT-class
setMethod(
    "plotLfcMa", "TranscriptionRatesLRT",
    function(object, file = NULL, width = 8,
            height = 6, dpi = 300) {

        betaTbl <- betaTbl(object)
        name1 <- slot(transcriptionRates1(object), "name")
        name2 <- slot(transcriptionRates2(object), "name")
        zeta <- 2000

        betaTbl <- betaTbl %>%
        mutate(beta = (.data$beta1 + .data$beta2) / 2,
                logbeta_zeta = log2(.data$beta * zeta),
                category =
                case_when(
                    (padj < 0.05) & (lfc > log2(1.2)) ~ "Up",
                    (padj < 0.05) & (lfc < log2(0.8)) ~ "Down",
                    TRUE ~ "Others"
                ),
                category = factor(.data$category, levels = c("Up", "Down", "Others")))

            p <- betaTbl %>%
            ggplot(aes(x = .data$logbeta_zeta, 
            y = .data$lfc, color = .data$category)) +
            geom_point(alpha = 0.5, size = 0.5) +
            scale_color_manual(values=c("#E41A1C", "#377EB8", "gray")) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
            ylim(-6, 6) +
            xlim(-10, 10) +
            labs(y = bquote(log[2]*"("*beta[.(name1)]/beta[.(name2)]*")"),
                x = expression(log[2]*bar(beta*zeta)),
                color = "Category") +
            theme_minimal()


        if (!is.null(file)) {
            ggsave(file, p,
                width = width, height = height,
                dpi = dpi
            )
        }
        return(p)
    }
)
