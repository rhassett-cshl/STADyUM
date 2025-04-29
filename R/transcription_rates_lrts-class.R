#' Class transcription_rates_lrts for comparing aspects of transcriptional
#' dynamics of the same TU under different conditions using Likelihood Ratio
#' Test Statistics
#'
#' Class \code{transcription_rates_lrts}
#'
#' @name transcription_rates_lrts-class
#' @rdname transcription_rates_lrts-class
#' @exportClass transcription_rates_lrts
methods::setClass("transcription_rates_lrts",
    slots = c(
        exp_data1 = "experiment_transcription_rates",
        exp_data2 = "experiment_transcription_rates",
        sc = "character"
    ),
    # validity =
)


#' likelihood_ratio_test
#'
#' DESCRIPTION
#' @param exp_data1 an \code{\link{experiment_transcription_rates-class}} object
#' @param exp_data2 an \code{\link{experiment_transcription_rates-class}} object
#' @param sc path to a csv file containing scale factors based on total or
#' spike-in reads
#'
#' @return a \code{\link{likelihood_ratio_test-class}} object
#'
#' @export
likelihood_ratio_test <- function(exp_data1, exp_data2, sc) {
    helper_tc_in <- "helper_function_em_two_condition.R"
    helper_pr_in <- "helper_function_em_pause_escape.R"

    omega_out <- file.path(result_dir, "omega.csv")
    beta_out <- file.path(result_dir, "beta.csv")

    theme_set(cowplot::theme_cowplot())

    # set up parameters
    k <- 50
    kmin <- 1
    kmax <- 200 # also used as k on the poisson case

    rnap_size <- 50
    zeta <- 2000

    # criteria for significance, p.adj < 0.05 and 2 fold differences
    sig_p <- 0.05
    lfc1 <- 0
    lfc2 <- 0

    # read in summary statistics from the varied pause site model exp_data1 and
    # exp_data2
    # union set of genes being analyzed
    gn_union <- intersect(exp_data1$gene_id, exp_data2$gene_id)
    exp_data1 <- exp_data1[match(gn_union, exp_data1$gene_id), ]
    exp_data2 <- exp_data2[match(gn_union, exp_data2$gene_id), ]

    #### Poisson-based Likelihood Ratio Tests ####
    # read in number of spike-in or total number of mappable reads
    # use them as scaling factor
    scale_tbl <- read_csv(spike_in, show_col_types = FALSE)
    # subset the right table
    scale_tbl <- scale_tbl[str_detect(exp_data2, scale_tbl$sample), ]
    #
    # based on formula (27) and (28), cancel out M and zeta since they are the
    # same between conditions
    lambda1 <- scale_tbl$control_1 + ifelse(is.na(scale_tbl$control_2), 0,
        scale_tbl$control_2
    )
    lambda2 <- scale_tbl$treated_1 + ifelse(is.na(scale_tbl$treated_2), 0,
        scale_tbl$treated_2
    )

    ## LRT for omega ##
    tao1 <- lambda1 / (lambda1 + lambda2)
    tao2 <- 1 - tao1

    omega_tbl <-
        tibble(
            gene_id = exp_data1$gene_id,
            chi1 = exp_data1$chi,
            chi2 = exp_data2$chi * lambda1 / lambda2,
            lfc = log2(chi2 / chi1)
        )

    omega_tbl <- omega_tbl %>%
        bind_cols(bind_rows(map2(exp_data1$s, exp_data2$s, omega_lrt,
            tao1 = tao1, tao2 = tao2
        )))

    omega_tbl <- omega_tbl %>%
        mutate(padj = p.adjust(p, method = "BH"))

    ## LRT for beta ##
    # need to jointly do EM one more time for H0, which assume betas are the
    # same between conditions, initialize fk with some reasonable values based
    # on heuristic
    fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
    fk_int <- fk_int / sum(fk_int)
    # try uniform distribution which gives similar results
    # fk_int <- rep(1 / kmax, kmax)

    # collect and construct stats
    s1 <- exp_data1$s
    s2 <- exp_data2$s
    t1_h1 <- map_dbl(exp_data1$Yk, sum)
    t2_h1 <- map_dbl(exp_data2$Yk, sum)
    Xk1 <- exp_data1$Xk
    Xk2 <- exp_data2$Xk
    # gene body length, assumed to be same between conditions
    M <- exp_data1$N

    # some values inherit from EM for H1, could be further integrated into EM
    # here
    chi_hat <- (s1 + s2) / M
    # beta_int <- chi_hat / (t1_h1 + t2_h1)
    beta_int <- chi_hat / (map_dbl(exp_data1$Xk, sum)
    + map_dbl(exp_data2$Xk, sum))

    # chi_hat for control and test sets
    chi_hat1 <- exp_data1$chi
    chi_hat2 <- exp_data2$chi
    # max iterations and tolerance for EM
    max_itr <- 500
    tor <- 1e-6
    # run EM for multiple combinations of parameters
    em_res <- pmap(
        list(Xk1, Xk2, beta_int, chi_hat, chi_hat1, chi_hat2),
        function(x, y, z, k, m, n) {
            tryCatch(
                main_EM_h0(fk_int,
                    Xk1 = x, Xk2 = y, kmin, kmax, beta_int = z,
                    chi_hat = k, chi_hat1 = m, chi_hat2 = n,
                    max_itr, tor
                ),
                error = function(err) {
                    # handling the error, one of the cases is when
                    # there is no read counts in the pause region
                    list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
                }
            )
        }
    )

    h0_likelihood <- map_dbl(em_res, ~ .x$likelihoods[[length(.x$likelihoods)]])

    beta_tbl <-
        tibble(
            gene_id = exp_data1$gene_id,
            beta1 = exp_data1$beta_adp,
            beta2 = exp_data2$beta_adp,
            lfc = log2(beta2 / beta1),
            fk_mean1 = exp_data1$fk_mean,
            fk_mean2 = exp_data2$fk_mean,
            fk_var1 = exp_data1$fk_var,
            fk_var2 = exp_data2$fk_var,
            # use eq (25) instead of (31) to compute T stats
            t_stats = exp_data1$likelihood + exp_data2$likelihood -
                h0_likelihood
        )

    # some genes with negative T stats, fix them
    idx <- beta_tbl$t_stats < 0

    # use parameter estimates from h0 as initial values for EM in h1
    h0_beta <- map_dbl(em_res, "beta")
    h0_fk1 <- map(em_res, "fk1")
    h0_fk2 <- map(em_res, "fk2")

    em_hc <- pmap(
        list(h0_fk1[idx], Xk1[idx], h0_beta[idx], chi_hat1[idx]),
        function(x, y, z, k) {
            tryCatch(
                main_EM(
                    fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
                    beta_int = z, chi_hat = k,
                    max_itr = max_itr, tor = tor
                ),
                error = function(err) {
                    list(
                        "beta" = NA, "Yk" = NA,
                        "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA
                    )
                }
            )
        }
    )

    em_ht <- pmap(
        list(h0_fk2[idx], Xk2[idx], h0_beta[idx], chi_hat2[idx]),
        function(x, y, z, k) {
            tryCatch(
                main_EM(
                    fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
                    beta_int = z, chi_hat = k,
                    max_itr = max_itr, tor = tor
                ),
                error = function(err) {
                    list(
                        "beta" = NA, "Yk" = NA,
                        "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA
                    )
                }
            )
        }
    )

    h1_likelihood1 <- map_dbl(em_hc, ~ .x$likelihoods[[length(.x$likelihoods)]])
    h1_likelihood2 <- map_dbl(em_ht, ~ .x$likelihoods[[length(.x$likelihoods)]])

    beta_tbl_idx <-
        tibble(
            gene_id = names(em_hc),
            beta1 = map_dbl(em_hc, "beta"),
            beta2 = map_dbl(em_ht, "beta"),
            lfc = log2(beta2 / beta1),
            fk_mean1 = map_dbl(em_hc, "fk_mean"),
            fk_mean2 = map_dbl(em_ht, "fk_mean"),
            fk_var1 = map_dbl(em_hc, "fk_var"),
            fk_var2 = map_dbl(em_ht, "fk_var"),
            t_stats = h1_likelihood1 + h1_likelihood2 - h0_likelihood[idx]
        )

    beta_tbl <- bind_rows(beta_tbl[!idx, ], beta_tbl_idx)

    beta_tbl <- beta_tbl %>%
        mutate(p = pchisq(2 * t_stats,
            df = 1, ncp = 0, lower.tail = FALSE,
            log.p = FALSE
        ))

    beta_tbl <- beta_tbl %>% mutate(padj = p.adjust(p, method = "BH"))

    p <- omega_tbl %>%
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

    ggsave(file.path(result_dir, "chi_distribution.png"),
        plot = p,
        width = 4, height = 5
    )

    p <- beta_tbl %>%
        select(contains("beta")) %>%
        pivot_longer(cols = contains("beta")) %>%
        mutate(value = log2(value * zeta)) %>%
        ggboxplot(
            x = "name", y = "value", fill = "name",
            palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            outlier.shape = NA, notch = TRUE
        ) +
        stat_compare_means(label.x.npc = "left", label.y = scale_tbl$beta_ymax)
    +
        scale_x_discrete(labels = c("beta1" = "Control", "beta2" = "Treated")) +
        labs(x = "", y = expression(log[2] * beta * zeta)) +
        coord_cartesian(ylim = c(scale_tbl$beta_ymin, scale_tbl$beta_ymax)) +
        theme(legend.position = "none")

    ggsave(file.path(result_dir, "beta_distribution.png"),
        plot = p,
        width = 4, height = 5
    )

    # mean vs. LFC
    omega_tbl <- omega_tbl %>%
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

    p <- omega_tbl %>%
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

    ggsave(file.path(result_dir, "chi_mean_vs_lfc.png"),
        plot = p,
        width = 6, height = 3
    )


    beta_tbl <- beta_tbl %>%
        mutate(
            beta = (beta1 + beta2) / 2,
            logbeta_zeta = log2(beta * zeta),
            category =
                case_when(
                    (padj < 0.05) & (lfc > lfc1) ~ "Up",
                    (padj < 0.05) & (lfc < lfc2) ~ "Down",
                    TRUE ~ "Others"
                ),
            category = factor(category, levels = c("Up", "Down", "Others"))
        )

    p <- beta_tbl %>%
        ggplot(aes(x = logbeta_zeta, y = lfc, color = category)) +
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

    ggsave(file.path(result_dir, "beta_mean_vs_lfc.png"),
        plot = p,
        width = 6, height = 3
    )


    # Changes of initiation and pause release rates
    lfc_tbl <- omega_tbl %>%
        select(gene_id, lfc, category) %>%
        inner_join(beta_tbl %>% select(gene_id, lfc, category),
            by = "gene_id", suffix = c("_chi", "_beta")
        )
    #
    # Number of genes with differential rates
    lfc_summary <- lfc_tbl %>%
        select(contains("category")) %>%
        pivot_longer(everything()) %>%
        group_by(name, value) %>%
        summarise(count = n()) %>%
        mutate(name = str_remove(name, "category_"))

    p <- lfc_summary %>%
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

    ggsave(file.path(result_dir, "gene_number_with_differential_rates.png"),
        plot = p,
        width = 6, height = 4
    )

    # mean and variance of pause sites in different beta categories
    beta_tbl <- beta_tbl %>%
        mutate(
            fk_std1 = fk_var1^0.5,
            fk_std2 = fk_var2^0.5
        )

    # plot mean and std of pause sites side by side
    p1 <- beta_tbl %>%
        ggplot(aes(x = fk_mean1, y = fk_std1)) +
        geom_pointdensity() +
        scale_color_viridis() +
        geom_vline(
            xintercept = c(50, 100), linetype = "dashed",
            color = "gray"
        ) +
        coord_cartesian(xlim = c(0, 200), ylim = c(0, 70)) +
        labs(x = "Mean of k", y = "SD of k")

    p2 <- beta_tbl %>%
        ggplot(aes(x = fk_mean2, y = fk_std2)) +
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
        filename = file.path(result_dir, "mean_vs_std_of_k_pointdensity.png"),
        plot = p,
        width = 14, height = 5
    )

    # output tables
    write_csv(omega_tbl %>% select(gene_id, chi1, chi2, lfc, t_stats, padj),
        file = omega_out
    )
    write_csv(beta_tbl %>%
        select(
            gene_id, beta1, beta2, lfc, t_stats, padj,
            fk_mean1, fk_mean2, fk_var1, fk_var2
        ), file = beta_out)
}

#' @rdname transcription_rates_lrts-class
#' @export
setGeneric("exp_data1", function(object) standardGeneric("exp_data1"))
setMethod("exp_data1", "transcription_rates_lrts", function(object) {
    slot(object, "exp_data1")
})

#' @rdname transcription_rates_lrts-class
#' @export
setGeneric("exp_data2", function(object) standardGeneric("exp_data2"))
setMethod("exp_data2", "transcription_rates_lrts", function(object) {
    slot(object, "exp_data2")
})

#' @rdname transcription_rates_lrts-class
#' @export
setGeneric("sc", function(object) standardGeneric("sc"))
setMethod("sc", "transcription_rates_lrts", function(object) slot(object, "sc"))

#' @rdname transcription_rates_lrts-class
#' @export
setGeneric("likelihood_ratio_test", function(exp_data1, exp_data2, sc) {
    standardGeneric("likelihood_ratio_test")
})
