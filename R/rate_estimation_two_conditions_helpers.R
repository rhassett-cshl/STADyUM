getFk <- function(Xk, Yk, fk, kmin, kmax) {
    t <- sum(Yk)
    u <- sum(Yk * seq(kmin, kmax))
    v <- sum(Yk * seq(kmin, kmax)^2)

    w <- sum(fk / (1 - fk) * (Xk - Yk), na.rm = TRUE)
    z <- sum(fk / (1 - fk) * (Xk - Yk) * seq(kmin, kmax), na.rm = TRUE)
    r <- sum(fk / (1 - fk) * (Xk - Yk) * seq(kmin, kmax)^2, na.rm = TRUE)

    fkMean <- (u - z) / (t - w)
    fkVar <- (v - r) / (t - w) - fkMean^2

    # avoid negative values
    if (fkVar < 1e-10) {
        fk[seq_len(length(fk))] <- 0
        # sometimes it looks like an integer but actually it's not
        fk[round(fkMean)] <- 1
        return(list("t" = t, "fk" = fk))
    }

    fk <- dnorm(kmin:kmax, mean = fkMean, sd = fkVar^0.5)
    fk <- fk / sum(fk)

    return(list("t" = t, "fk" = fk))
}

getMaximizationH0 <- function(chiHat, Xk1, Xk2, Yk1, Yk2, fk1, fk2, kmin, kmax)
{
    fk1Ls <- getFk(Xk1, Yk1, fk1, kmin, kmax)
    fk2Ls <- getFk(Xk2, Yk2, fk2, kmin, kmax)

    beta <- chiHat / (fk1Ls$t + fk2Ls$t)

    return(list("beta" = beta, "fk1" = fk1Ls$fk, "fk2" = fk2Ls$fk))
}

# EM function to estimate parameters when assuming beta is the same between
# conditions
mainExpectationMaximizationH0 <- function(fkInt, Xk1, Xk2, kmin, kmax, betaInt,
                                            chiHat, chiHat1, chiHat2, 
                                            maxItr = 100, tor = 1e-3) {
    betas <- list(); likelihoods <- list(); yks <- list(); flag <- "normal"
    for (i in seq_len(maxItr)) {
        if (i == 1) {
            Yk1 <- getExpectation(fkInt, Xk1, betaInt)
            Yk2 <- getExpectation(fkInt, Xk2, betaInt)
            hats <- getMaximizationH0(
                chiHat, Xk1, Xk2, Yk1, Yk2,
                fkInt, fkInt, kmin, kmax
            )
        }
        if (i != 1) {
            Yk1 <- getExpectation(hats$fk1, Xk1, hats$beta)
            Yk2 <- getExpectation(hats$fk2, Xk2, hats$beta)
            hats <- getMaximizationH0(
                chiHat, Xk1, Xk2, Yk1, Yk2,
                hats$fk1, hats$fk2, kmin, kmax
            )
        }
        likelihoods[[i]] <-
            getLikelihood(
                beta = hats$beta, chi = chiHat1, Xk = Xk1, Yk = Yk1,
                fk = hats$fk1
            ) +
            getLikelihood(
                beta = hats$beta, chi = chiHat2, Xk = Xk2, Yk = Yk2,
                fk = hats$fk2
            )
        betas[[i]] <- hats$beta
        if (any(hats$fk1 == 1) & any(hats$fk2 == 1)) {
            hats$beta <- chiHat / (Xk1[which(hats$fk1 == 1)] +
                Xk2[which(hats$fk2 == 1)])
            flag <- "single_site"
            break
        }
        if (i > 1) {
            diff <- likelihoods[[i]] - likelihoods[[i - 1]]
            if (diff <= tor) break
        }
    }
    if (i == maxItr) flag <- "max_iteration"
    return(list(
        "beta" = hats$beta, "Yk1" = Yk1, "Yk2" = Yk2,
        "fk1" = hats$fk1, "fk2" = hats$fk2, "betas" = betas,
        "likelihoods" = likelihoods, "yks" = yks, "flag" = flag
    ))
}

#### Functions for LRTs ####
# formulas are based on the unified model preprint v5
# formula (27), calculate t stats for omega
omegaLRT <- function(s1, s2, tao1, tao2) {
    ## compute T statistic and p values
    tStats <- s1 * log(s1 / (tao1 * (s1 + s2))) + s2 * log(s2 / (tao2 * (s1 +
        s2)))
    p <- pchisq(2 * tStats, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    return(c("tStats" = tStats, "p" = p))
}
