#' Class transcriptionRatesLrts for comparing aspects of transcriptional
#' dynamics of the same TU under different conditions using Likelihood Ratio
#' Test Statistics
#'
#' Class \code{transcriptionRatesLrts}
#'
#' @name transcriptionRatesLrts-class
#' @rdname transcriptionRatesLrts-class
#' @exportClass transcriptionRatesLrts
methods::setClass("transcriptionRatesLrts",
    slots = c(
        expData1 = "experimentTranscriptionRates",
        expData2 = "experimentTranscriptionRates",
        sc = "character"
    ),
    # validity =
)


#' likelihoodRatioTest
#'
#' DESCRIPTION
#' @param expData1 an \code{\link{experimentTranscriptionRates-class}} object
#' @param expData2 an \code{\link{experimentTranscriptionRates-class}} object
#' @param sc path to a csv file containing scale factors based on total or
#' spike-in reads
#'
#' @return a \code{\link{likelihoodRatioTest-class}} object
#'
#' @export
likelihoodRatioTest <- function(expData1, expData2, sc) {
    helperTcIn <- "helper_function_em_two_condition.R"
    helperPrIn <- "helper_function_em_pause_escape.R"

    omegaOut <- file.path(resultDir, "omega.csv")
    betaOut <- file.path(resultDir, "beta.csv")

    theme_set(cowplot::theme_cowplot())

    # set up parameters
    k <- 50
    kmin <- 1
    kmax <- 200 # also used as k on the poisson case

    rnapSize <- 50
    zeta <- 2000

    # criteria for significance, p.adj < 0.05 and 2 fold differences
    sigP <- 0.05
    lfc1 <- 0
    lfc2 <- 0

    # read in summary statistics from the varied pause site model expData1 and
    # expData2
    # union set of genes being analyzed
    gnUnion <- intersect(expData1$geneId, expData2$geneId)
    expData1 <- expData1[match(gnUnion, expData1$geneId), ]
    expData2 <- expData2[match(gnUnion, expData2$geneId), ]

    #### Poisson-based Likelihood Ratio Tests ####
    # read in number of spike-in or total number of mappable reads
    # use them as scaling factor
    scaleTbl <- readCsv(spikeIn, showColTypes = FALSE)
    # subset the right table
    scaleTbl <- scaleTbl[strDetect(expData2, scaleTbl$sample), ]
    #
    # based on formula (27) and (28), cancel out M and zeta since they are the
    # same between conditions
    lambda1 <- scaleTbl$control1 + ifelse(is.na(scaleTbl$control2), 0,
        scaleTbl$control2
    )
    lambda2 <- scaleTbl$treated1 + ifelse(is.na(scaleTbl$treated2), 0,
        scaleTbl$treated2
    )

    ## LRT for omega ##
    tao1 <- lambda1 / (lambda1 + lambda2)
    tao2 <- 1 - tao1

    omegaTbl <-
        tibble(
            geneId = expData1$geneId,
            chi1 = expData1$chi,
            chi2 = expData2$chi * lambda1 / lambda2,
            lfc = log2(chi2 / chi1)
        )

    omegaTbl <- omegaTbl %>%
        bind_cols(bind_rows(map2(expData1$s, expData2$s, omegaLrt,
            tao1 = tao1, tao2 = tao2
        )))

    omegaTbl <- omegaTbl %>%
        mutate(padj = p.adjust(p, method = "BH"))

    # LRT for beta #
    # need to jointly do EM one more time for H0, which assume betas are the
    # same between conditions, initialize fk with some reasonable values based
    # on heuristic
    fkInt <- dnorm(kmin:kmax, mean = 50, sd = 100)
    fkInt <- fkInt / sum(fkInt)
    # try uniform distribution which gives similar results
    # fk_int <- rep(1 / kmax, kmax)

    # collect and construct stats
    s1 <- expData1$s
    s2 <- expData2$s
    t1H1 <- mapDbl(expData1$Yk, sum)
    t2H1 <- mapDbl(expData2$Yk, sum)
    Xk1 <- expData1$Xk
    Xk2 <- expData2$Xk
    # gene body length, assumed to be same between conditions
    M <- expData1$N

    # some values inherit from EM for H1, could be further integrated into EM
    # here
    chiHat <- (s1 + s2) / M
    # beta_int <- chi_hat / (t1_h1 + t2_h1)
    betaInt <- chiHat / (mapDbl(expData1$Xk, sum)
    + mapDbl(expData2$Xk, sum))

    # chiHat for control and test sets
    chiHat1 <- expData1$chi
    chiHat2 <- expData2$chi
    # max iterations and tolerance for EM
    maxItr <- 500
    tor <- 1e-6
    # run EM for multiple combinations of parameters
    emRes <- pmap(
        list(Xk1, Xk2, betaInt, chiHat, chiHat1, chiHat2),
        function(x, y, z, k, m, n) {
            tryCatch(
                mainEMH0(fkInt,
                    Xk1 = x, Xk2 = y, kmin, kmax, betaInt = z,
                    chiHat = k, chiHat1 = m, chiHat2 = n,
                    maxItr, tor
                ),
                error = function(err) {
                    # handling the error, one of the cases is when
                    # there is no read counts in the pause region
                    list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
                }
            )
        }
    )

    h0Likelihood <- map_dbl(emRes, ~ .x$likelihoods[[length(.x$likelihoods)]])

    # use eq (25) instead of (31) to compute T stats
    betaTbl <-
        tibble(
            geneId = expData1$geneId,
            beta1 = expData1$betaAdp,
            beta2 = expData2$betaAdp,
            lfc = log2(beta2 / beta1),
            fkMean1 = expData1$fkMean,
            fkMean2 = expData2$fkMean,
            fkVar1 = expData1$fkVar,
            fkVar2 = expData2$fkVar,
            tStats = expData1$likelihood + expData2$likelihood -
                h0Likelihood
        )

    # some genes with negative T stats, fix them
    idx <- betaTbl$tStats < 0

    # use parameter estimates from h0 as initial values for EM in h1
    h0Beta <- map_dbl(emRes, "beta")
    h0Fk1 <- map(emRes, "fk1")
    h0Fk2 <- map(emRes, "fk2")

    emHc <- pmap(
        list(h0Fk1[idx], Xk1[idx], h0Beta[idx], chiHat1[idx]),
        function(x, y, z, k) {
            tryCatch(
                mainEM(
                    fkInt = x, Xk = y, kmin = kmin, kmax = kmax,
                    betaInt = z, chiHat = k,
                    maxItr = maxItr, tor = tor
                ),
                error = function(err) {
                    list(
                        "beta" = NA, "Yk" = NA,
                        "fkMean" = NA, "fkVar" = NA, "likelihoods" = NA
                    )
                }
            )
        }
    )

    emHt <- pmap(
        list(h0Fk2[idx], Xk2[idx], h0Beta[idx], chiHat2[idx]),
        function(x, y, z, k) {
            tryCatch(
                mainEM(
                    fkInt = x, Xk = y, kmin = kmin, kmax = kmax,
                    betaInt = z, chiHat = k,
                    maxItr = maxItr, tor = tor
                ),
                error = function(err) {
                    list(
                        "beta" = NA, "Yk" = NA,
                        "fkMean" = NA, "fkVar" = NA, "likelihoods" = NA
                    )
                }
            )
        }
    )

    h1Likelihood1 <- map_dbl(emHc, ~ .x$likelihoods[[length(.x$likelihoods)]])
    h1Likelihood2 <- map_dbl(emHt, ~ .x$likelihoods[[length(.x$likelihoods)]])

    betaTblIdx <-
        tibble(
            geneId = names(emHc),
            beta1 = map_dbl(emHc, "beta"),
            beta2 = map_dbl(emHt, "beta"),
            lfc = log2(beta2 / beta1),
            fkMean1 = map_dbl(emHc, "fkMean"),
            fkMean2 = map_dbl(emHt, "fkMean"),
            fkVar1 = map_dbl(emHc, "fkVar"),
            fkVar2 = map_dbl(emHt, "fkVar"),
            tStats = h1Likelihood1 + h1Likelihood2 - h0Likelihood[idx]
        )

    betaTbl <- bind_rows(betaTbl[!idx, ], betaTblIdx)

    betaTbl <- betaTbl %>%
        mutate(p = pchisq(2 * tStats,
            df = 1, ncp = 0, lower.tail = FALSE,
            log.p = FALSE
        ))

    betaTbl <- betaTbl %>% mutate(padj = p.adjust(p, method = "BH"))

    p <- omegaTbl %>%
        select(contains("chi")) %>%
        pivot_longer(cols = contains("chi")) %>%
        mutate(value = log2(value)) %>%
        ggboxplot(
            x = "name", y = "value", fill = "name",
            palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            outlier.shape = NA, notch = TRUE
        ) +
        stat_compare_means() +
        scale_x_discrete(labels = c("chi1" = "Control", "chi2" = "Treated")) +
        labs(x = "", y = expression(log[2] * chi)) +
        # coord_cartesian(ylim = c(scale_tbl$chi_ymin, scale_tbl$chi_ymax)) +
        theme(legend.position = "none")

    ggsave(file.path(resultDir, "chi_distribution.png"),
        plot = p,
        width = 4, height = 5
    )

    p <- betaTbl %>%
        select(contains("beta")) %>%
        pivot_longer(cols = contains("beta")) %>%
        mutate(value = log2(value * zeta)) %>%
        ggboxplot(
            x = "name", y = "value", fill = "name",
            palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            outlier.shape = NA, notch = TRUE
        ) +
        stat_compare_means(label.x.npc = "left", label.y = scaleTbl$betaYmax) +
        scale_x_discrete(labels = c("beta1" = "Control", "beta2" = "Treated")) +
        labs(x = "", y = expression(log[2] * beta * zeta)) +
        coord_cartesian(ylim = c(scaleTbl$betaYmin, scaleTbl$betaYmax)) +
        theme(legend.position = "none")

    ggsave(file.path(resultDir, "beta_distribution.png"),
        plot = p,
        width = 4, height = 5
    )

    # mean vs. LFC
    omegaTbl <- omegaTbl %>%
        mutate(
            chi = (chi1 + chi2) / 2,
            logchi = log2(chi),
            category =
                case_when(
                    (padj < 0.05) & (lfc > lfc1) ~ "Up",
                    (padj < 0.05) & (lfc < lfc2) ~ "Down",
                    TRUE ~ "Others"
                ),
            category = factor(category, levels = c("Up", "Down", "Others"))
        )

    p <- omegaTbl %>%
        ggplot(aes(x = logchi, y = lfc, color = category)) +
        geom_point(alpha = 0.5, size = 0.5) +
        scale_color_manual(values = c("#E41A1C", "#377EB8", "gray")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
        ylim(-6, 6) +
        labs(
            y = expression(log[2] * "FC(Treated/Control)"),
            x = expression(log[2] * "mean(" * chi * ")"),
            color = "Category"
        ) +
        cowplot::theme_cowplot()

    ggsave(file.path(resultDir, "chi_mean_vs_lfc.png"),
        plot = p,
        width = 6, height = 3
    )


    betaTbl <- betaTbl %>%
        mutate(
            beta = (beta1 + beta2) / 2,
            logbetaZeta = log2(beta * zeta),
            category =
                case_when(
                    (padj < 0.05) & (lfc > lfc1) ~ "Up",
                    (padj < 0.05) & (lfc < lfc2) ~ "Down",
                    TRUE ~ "Others"
                ),
            category = factor(category, levels = c("Up", "Down", "Others"))
        )

    p <- betaTbl %>%
        ggplot(aes(x = logbetaZeta, y = lfc, color = category)) +
        geom_point(alpha = 0.5, size = 0.5) +
        scale_color_manual(values = c("#E41A1C", "#377EB8", "gray")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
        ylim(-6, 6) +
        xlim(-10, 10) +
        labs(
            y = expression(log[2] * "FC(Treated/Control)"),
            x = expression(log[2] * "mean(" * beta * zeta * ")"),
            color = "Category"
        ) +
        cowplot::theme_cowplot()

    ggsave(file.path(resultDir, "beta_mean_vs_lfc.png"),
        plot = p,
        width = 6, height = 3
    )


    # Changes of initiation and pause release rates
    lfcTbl <- omegaTbl %>%
        select(geneId, lfc, category) %>%
        inner_join(betaTbl %>% select(geneId, lfc, category),
            by = "geneId", suffix = c("_chi", "_beta")
        )
    #
    # Number of genes with differential rates
    lfcSummary <- lfcTbl %>%
        select(contains("category")) %>%
        pivot_longer(everything()) %>%
        group_by(name, value) %>%
        summarise(count = n()) %>%
        mutate(name = str_remove(name, "category_"))

    p <- lfcSummary %>%
        ggplot(aes(x = name, y = count, fill = value)) +
        geom_col(position = "dodge") +
        scale_x_discrete(labels = c(
            "alpha" = expression(alpha),
            "beta" = expression(beta),
            "chi" = expression(chi)
        )) +
        geom_text(aes(label = count),
            position = position_dodge(width = 0.9),
            vjust = -0.25
        ) +
        labs(x = "", y = "Number of Genes", fill = "Category")

    ggsave(file.path(resultDir, "gene_number_with_differential_rates.png"),
        plot = p,
        width = 6, height = 4
    )

    # mean and variance of pause sites in different beta categories
    betaTbl <- betaTbl %>%
        mutate(
            fkStd1 = fkVar1^0.5,
            fkStd2 = fkVar2^0.5
        )

    # plot mean and std of pause sites side by side
    p1 <- betaTbl %>%
        ggplot(aes(x = fkMean1, y = fkStd1)) +
        geom_pointdensity() +
        scale_color_viridis() +
        geom_vline(
            xintercept = c(50, 100), linetype = "dashed",
            color = "gray"
        ) +
        coord_cartesian(xlim = c(0, 200), ylim = c(0, 70)) +
        labs(x = "Mean of k", y = "SD of k")

    p2 <- betaTbl %>%
        ggplot(aes(x = fkMean2, y = fkStd2)) +
        geom_pointdensity() +
        scale_color_viridis() +
        geom_vline(
            xintercept = c(50, 100), linetype = "dashed",
            color = "gray"
        ) +
        coord_cartesian(xlim = c(0, 200), ylim = c(0, 70)) +
        labs(x = "Mean of k", y = "SD of k")

    p <- cowplot::plot_grid(p1, p2)

    ggsave(
        filename = file.path(resultDir, "mean_vs_std_of_k_pointdensity.png"),
        plot = p,
        width = 14, height = 5
    )

    # output tables
    write_csv(omegaTbl %>% select(geneId, chi1, chi2, lfc, tStats, padj),
        file = omegaOut
    )
    write_csv(betaTbl %>%
        select(
            geneId, beta1, beta2, lfc, tStats, padj,
            fkMean1, fkMean2, fkVar1, fkVar2
        ), file = betaOut)
}

#' @rdname transcriptionRatesLrts-class
#' @export
setGeneric("expData1", function(object) standardGeneric("expData1"))
setMethod("expData1", "transcriptionRatesLrts", function(object) {
    slot(object, "expData1")
})

#' @rdname transcriptionRatesLrts-class
#' @export
setGeneric("expData2", function(object) standardGeneric("expData2"))
setMethod("expData2", "transcriptionRatesLrts", function(object) {
    slot(object, "expData2")
})

#' @rdname transcriptionRatesLrts-class
#' @export
setGeneric("sc", function(object) standardGeneric("sc"))
setMethod("sc", "transcriptionRatesLrts", function(object) slot(object, "sc"))

#' @rdname transcriptionRatesLrts-class
#' @export
setGeneric("likelihoodRatioTest", function(expData1, expData2, sc) {
    standardGeneric("likelihoodRatioTest")
})
