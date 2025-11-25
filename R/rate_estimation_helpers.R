processBw <- function(bw, strand) {
    strand(bw) <- strand
    bw$score <- abs(bw$score)
    bw <- bw[bw$score > 0]
    bw <- GenomeInfoDb::keepStandardChromosomes(bw, pruning.mode = "coarse")

    ## Ensure each range is a single basepair
    bw <- GenomicRanges::GRanges(
        seqnames = seqnames(bw),
        ranges = IRanges::IRanges(
            start = start(bw),
            end = start(bw),
            names = names(bw)
        ),
        strand = strand(bw),
        score = bw$score
    )

    return(bw)
}

summariseBw <- function(bw, grng, colName) {
    hits <- findOverlaps(bw, grng, ignore.strand = FALSE)
    
    if (length(hits) == 0) {
        rc <- DataFrame(gene_id = character(0), score = numeric(0))
    } else {
        qh <- queryHits(hits)
        sh <- subjectHits(hits)
        
        weighted <- mcols(bw)$score[qh] * width(bw)[qh]
        sums <- tapply(weighted, sh, sum)
        
        gene_ids <- mcols(grng)$gene_id[as.integer(names(sums))]
        
        rc <- DataFrame(
            gene_id = gene_ids,
            score = as.numeric(sums)
        )
    }
    
    colnames(rc) <- c("gene_id", colName)
    return(rc)
}

getExpectation <- function(fk, Xk, beta) {
    Yk <- Xk / (1 - beta + beta / fk)
    return(Yk)
}

getLikelihood <- function(beta, chi, Xk, Yk, fk) {
    ## part of the original likelihood function associated with beta, Xk and Yk
    ## used as criteria to terminate EM
    t <- sum(Yk)

    idx1 <- fk != 0
    idx2 <- 1 - fk != 0
    likelihood <- -t * log(beta) - chi / beta +
        sum(Yk[idx1] * log(fk[idx1])) + sum((Xk - Yk)[idx2] *
            log(1 - fk[idx2]))
    return(likelihood)
}

######## functions for adapted model #########
# model allows pause sites to vary across cells
# EM doesn't include phi estimates
# functions for EM based on Gaussian distributed k
pauseEscapeMaximization <- function(chiHat, Xk, Yk, fk, kmin, kmax) {
    t <- sum(Yk)
    u <- sum(Yk * seq(kmin, kmax))
    v <- sum(Yk * seq(kmin, kmax)^2)

    w <- sum(fk / (1 - fk) * (Xk - Yk))
    z <- sum(fk / (1 - fk) * (Xk - Yk) * seq(kmin, kmax))
    r <- sum(fk / (1 - fk) * (Xk - Yk) * seq(kmin, kmax)^2)

    fkMean <- (u - z) / (t - w)
    fkVar <- (v - r) / (t - w) - fkMean^2

    if (fkVar < 1e-10) {
        fk[seq_len(length(fk))] <- 0
        fk[round(fkMean)] <- 1
    } else {
        fk <- dnorm(kmin:kmax, mean = fkMean, sd = fkVar^0.5)
        fk <- fk / sum(fk)
    }

    beta <- chiHat / t

    return(list(
        "beta" = beta, "fk" = fk, "fkMean" = fkMean,
        "fkVar" = fkVar
    ))
}

pauseEscapeEM <- function(
    fkInt, Xk, kmin, kmax, betaInt, chiHat, maxItr = 100, tor = 1e-3) {
    betas <- list(); likelihoods <- list(); flag <- "normal"

    for (i in seq_len(maxItr)) {
        if (i == 1) {
            Yk <- getExpectation(fkInt, Xk, betaInt)
            hats <- pauseEscapeMaximization(
                chiHat, Xk, Yk, fkInt, kmin,
                kmax
            )
            beta <- betaInt
        }
        if (i != 1) {
            Yk <- getExpectation(hats$fk, Xk, hats$beta)
            hats <- pauseEscapeMaximization(
                chiHat, Xk, Yk, hats$fk, kmin,
                kmax
            )
        }
        likelihoods[[i]] <-
            getLikelihood(
                beta = hats$beta, chi = chiHat, Xk = Xk, Yk = Yk,
                fk = hats$fk
            )

        betas[[i]] <- hats$beta
        if (any(hats$fk == 1)) {
            hats$beta <- chiHat / Xk[which(hats$fk == 1)]
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
        "beta" = hats$beta, "Yk" = Yk, "fk" = hats$fk,
        "fkMean" = hats$fkMean, "fkVar" = hats$fkVar,
        "betas" = betas, "likelihoods" = likelihoods, "flag" = flag
    ))
}

######## functions for steric hindrance #########
# model allows both varied pause sites and steric hindrance, EM contains phi
# estimations
# functions for EM based on Gaussian distributed k
multRnapPhi <-
    function(alpha, beta, f1, f2) {
        return(
            (1 - f1 - f2) * alpha / (alpha + beta) +
                f1 * alpha^2 / (alpha^2 + beta^2 + alpha * beta) +
                f2 * alpha^3 / (beta^2 * alpha + beta^3 + alpha^2 * beta +
                    alpha^3)
        )
    }

phiPolynom <- function(phi, beta, omega, f1, f2) {
    alpha <- omega / (1 - phi)

    return(multRnapPhi(alpha, beta, f1, f2) - phi)
}

# find phi corresponding to omega and beta by solving polynomial that results
# from substituting alpha = omega/(1-phi) into the equation that defines phi in
# terms of alpha and beta
multRnapPhiOmega <- function(omega, beta, f1, f2) {
    ## set bounds for solution
    lb <- 1e-6
    ub <- 1 - 1e-6
    epsilon <- 1e-3

    ## make sure opposite signs at bounds; if not treat as an edge case
    phi1 <- phiPolynom(lb, beta, omega, f1, f2)
    phi2 <- phiPolynom(ub, beta, omega, f1, f2)

    if ((phi1 > 0 & phi2 > 0) | (phi1 < 0 & phi2 < 0)) {
        ## in this case phi is almost certainly close to 0 or 1
        ## pick the closer case
        if (abs(phi1) < epsilon) {
            return(epsilon)
        } else if (abs(phi2) < epsilon) {
            return(1 - epsilon)
        } else {
            simpleError("Edge case fail")
        }
    }

    try(phiRoot <- uniroot(phiPolynom, c(lb, ub), beta, omega, f1, f2),
        silent = FALSE
    )

    return(phiRoot$root)
}

# version of above that uses log parameterization of beta
betaEcllOmegaLog <- function(args, omega, chi, t, f1, f2) {
    beta <- exp(args[1])

    ## fix Beta prior with a=b=2
    a <- 2
    b <- 2

    phi <- multRnapPhiOmega(omega, beta, f1, f2)

    if (phi <= 0 | phi >= 1) {
        retval <- -Inf
    } else {
        retval <- -t * log(beta) - chi / beta + (a - 1) * log(phi) + (b - 1) *
            log(1 - phi)
    }
    return(-retval)
}

betaMStepOmega <- function(chi, t, f1, f2, oldphi, oldbeta, lambda, zeta) {
    omega <- chi * zeta / lambda
    ret <- list("par" = NA_integer_)

    try(ret <-
        optim(c(log(oldbeta)), betaEcllOmegaLog,
            gr = NULL, omega = omega, chi = chi,
            t = t, f1 = f1, f2 = f2, method = "L-BFGS-B", lower = log(omega)
        ), silent = FALSE)

    beta <- exp(ret$par)
    phi <- multRnapPhiOmega(omega, beta, f1, f2)

    return(list("beta" = beta, "phi" = phi))
}

stericHindranceMaximization <- function(
    chiHat, Xk, Yk, fk, kmin, kmax, f1,
    f2, phi, beta, lambda, zeta) {
    t <- sum(Yk)
    u <- sum(Yk * seq(kmin, kmax))
    v <- sum(Yk * seq(kmin, kmax)^2)

    w <- sum(fk / (1 - fk) * (Xk - Yk))
    z <- sum(fk / (1 - fk) * (Xk - Yk) * seq(kmin, kmax))
    r <- sum(fk / (1 - fk) * (Xk - Yk) * seq(kmin, kmax)^2)

    fkMean <- (u - z) / (t - w)
    fkVar <- (v - r) / (t - w) - fkMean^2

    if (fkVar < 1e-10) {
        fk[seq_len(length(fk))] <- 0
        fk[round(fkMean)] <- 1
    } else {
        fk <- dnorm(kmin:kmax, mean = fkMean, sd = fkVar^0.5)
        fk <- fk / sum(fk)
    }

    param <- betaMStepOmega(
        chi = chiHat, t = t, f1 = f1, f2 = f2, oldphi = phi,
        oldbeta = beta, lambda = lambda, zeta = zeta
    )

    return(list(
        "beta" = param$beta, "phi" = param$phi, "fk" = fk,
        "fkMean" = fkMean, "fkVar" = fkVar
    ))
}

calculateF <- function(s, k) {
    x <- round(rnorm(1e7, mean = k, sd = 25))
    x <- x[x >= 17 & x <= 200]
    f <- mean(x > s)
    f1 <- mean((x > s) & (x <= 2 * s))
    f2 <- mean(x > 2 * s)
    return(c("f" = f, "f1" = f1, "f2" = f2))
}

stericHindranceEM <- function(
    Xk, kmin, kmax, f1, f2, fkInt, betaInt, phiInt, chiHat, lambda, zeta,
    maxItr = 100, tor = 1e-3) {
    betas <- list(); likelihoods <- list(); flag <- "normal"

    for (i in seq_len(maxItr)) {
        if (i == 1) {
            Yk <- getExpectation(fkInt, Xk, betaInt)
            hats <- stericHindranceMaximization(
                chiHat, Xk, Yk, fkInt, kmin, kmax, f1, f2,
                phiInt, betaInt, lambda, zeta
            )
            beta <- betaInt
        }
        if (i != 1) {
            Yk <- getExpectation(hats$fk, Xk, hats$beta)
            hats <- stericHindranceMaximization(
                chiHat, Xk, Yk, hats$fk, kmin, kmax, f1, f2,
                hats$phi, hats$beta, lambda, zeta
            )
        }
        likelihoods[[i]] <-
            getLikelihood(
                beta = hats$beta, chi = chiHat, Xk = Xk, Yk = Yk,
                fk = hats$fk
            )

        betas[[i]] <- hats$beta
        if (any(hats$fk == 1)) {
            hats$beta <- chiHat / Xk[which(hats$fk == 1)]
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
        "beta" = hats$beta, "Yk" = Yk, "fk" = hats$fk, "fkMean" =
            hats$fkMean, "fkVar" = hats$fkVar, "betas" = betas,
        "likelihoods" = likelihoods, "phi" = hats$phi, "flag" = flag
    ))
}
