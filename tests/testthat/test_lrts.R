test_that("Transcription Rates Likelihood Ratio Test works as intented", {
    simpol1 <- simulatePolymerase(k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
    alpha=0.5, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500,
    cellNum=1000, polSize=33, addSpace=17, time=40, timesToRecord=NULL)

    simpol2 <- simulatePolymerase(k=100, ksd=25, kMin=75, kMax=200, geneLen=1950, alpha=0.5, beta=10, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500, cellNum=1000, polSize=33, addSpace=17, time=40, timesToRecord=NULL)

    simRates1 <- estimateTranscriptionRates(simpol1, name="low_k_low_beta")

    simRates2 <- estimateTranscriptionRates(simpol2, name="high_k_high_beta")

    lrts <- likelihoodRatioTest(simRates1, simRates2)

    betaTbl <- betaTbl(lrts)

    avgBeta1 <- mean(betaTbl$beta1, na.rm = TRUE)
    avgBeta2 <- mean(betaTbl$beta2, na.rm = TRUE)
    expect_gt(avgBeta2, avgBeta1)  

    avgFkMean1 <- mean(betaTbl$fkMean1, na.rm = TRUE)
    avgFkMean2 <- mean(betaTbl$fkMean2, na.rm = TRUE)
    expect_gt(avgFkMean2, avgFkMean1)  
})