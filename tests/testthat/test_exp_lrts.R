test_that("Experimental Data Transcription Rates Likelihood Ratio Test works as intented", {

    load(system.file("extdata", "granges_for_read_counting_DLD1_chr21.RData", 
    package = "STADyUM"))

    controlRates <- estimateTranscriptionRates(system.file("extdata","PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus_chr21.bw", package = "STADyUM"),
    system.file("extdata", "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus_chr21.bw",
    package = "STADyUM"), bw_pause_filtered, bw_gb_filtered, "Control", stericHindrance=TRUE, omegaScale=12.3768278981277)

    treatedRates <- estimateTranscriptionRates(system.file("extdata",
    "PROseq-DLD1-aoi-NELFC_Auxin-SE_plus_chr21.bw", package = "STADyUM"),
    system.file("extdata", "PROseq-DLD1-aoi-NELFC_Auxin-SE_minus_chr21.bw", 
    package = "STADyUM"), bw_pause_filtered, bw_gb_filtered, "Auxin Treated", stericHindrance=TRUE, omegaScale=11.0318379571379)

    lrts <- likelihoodRatioTest(controlRates, treatedRates, spikeInFile=system.file("extdata", "spikein_scaling_factor.csv",  package = "STADyUM"))

    betaTbl <- betaTbl(lrts)

    avgBeta1 <- mean(betaTbl$beta1, na.rm = TRUE)
    avgBeta2 <- mean(betaTbl$beta2, na.rm = TRUE)
    expect_gt(avgBeta1, avgBeta2)  

    avgFkMean1 <- mean(betaTbl$fkMean1, na.rm = TRUE)
    avgFkMean2 <- mean(betaTbl$fkMean2, na.rm = TRUE)
    expect_gt(avgFkMean2, avgFkMean1)  

    plotMeanPauseDistrib <- plotMeanPauseDistrib(controlRates)
    expect_s3_class(plotMeanPauseDistrib, "ggplot")
    plotExpAct <- plotExpectedVsActualPauseSiteCounts(controlRates)
    expect_s3_class(plotExpAct, "ggplot")
    chiPlot <- plotChiDistrib(controlRates)
    expect_s3_class(chiPlot, "ggplot")

    expShowTest <- capture.output(show(controlRates))
    expect_true(any(grepl("- Omega scale:", expShowTest)))

    chiViolinPlot <- ChiViolinPlot(lrts)
    expect_s3_class(chiViolinPlot, "ggplot")
})