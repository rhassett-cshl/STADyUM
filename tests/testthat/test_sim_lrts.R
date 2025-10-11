test_that("Simulated Transcription Rates Likelihood Ratio Test works as intented", {

    simpol1 <- simulatePolymerase(k=50, ksd=25, kMin=17, kMax=200, geneLen=1950,
    alpha=0.5, beta=1, zeta=2000, zetaSd=1000, zetaMin=1500, zetaMax=2500, cellNum=100, polSize=33, addSpace=17, time=30, timesToRecord=c(10,20,30))

    simpol2 <- simulatePolymerase(k=100, ksd=25, kMin=75, kMax=200,
    geneLen=1950, alpha=0.5, beta=10, zeta=2000, zetaSd=1000, zetaMin=1500,
    zetaMax=2500, cellNum=100, polSize=33, addSpace=17, time=30,
    timesToRecord=c(10,20,30))

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

    pauseSitesPlot <- plotPauseSites(simpol1)
    expect_s3_class(pauseSitesPlot, "ggplot")
    positionHeatmap <- plotPositionHeatmap(simpol1, timePoint=10)
    expect_s3_class(positionHeatmap, "ggplot")
    cellsPlot <- plotCombinedCells(simpol1, timePoint=10)
    expect_s3_class(cellsPlot, "ggplot")
    pcaPlot <- plotPolymerasePCA(simpol1, timePoint=10)
    expect_s3_class(pcaPlot, "ggplot")

    sampleReadCounts(simpol1)
    readCounts <- readCounts(simpol1)
    expect_true(all(readCounts >= 0))
    expect_gt(mean(readCounts), 0)
    expect_lt(mean(readCounts), 1)

    betaPlot <- BetaViolinPlot(lrts)
    expect_s3_class(betaPlot, "ggplot")
    pauseContourMap <- plotPauseSiteContourMapTwoConditions(lrts)
    expect_s3_class(pauseContourMap, "ggplot")
    maPlot <- plotLfcMa(lrts)
    expect_s3_class(maPlot, "ggplot")

    timePoints <- getAvailableTimePoints(simpol1)
    expect_equal(timePoints, c(10,20))
    combinedCellsData10 <- getCombinedCellsDataAtTime(simpol1, timePoint=10)
    expect_true(all(combinedCellsData10 >= 0))
    expect_true(combinedCellsData10[1] == 100)

    betaVsChiPlot <- plotBetaVsChi(simRates1)
    expect_s3_class(betaVsChiPlot, "ggplot")
    pauseSiteContourMap <- plotPauseSiteContourMap(simRates1)
    expect_s3_class(pauseSiteContourMap, "ggplot")

    simpolShowTest <- capture.output(show(simpol1))
    expect_true(any(grepl("top 10 most occupied sites across all cells", simpolShowTest)))
    expect_true(parameters(simpol1)$k == 50)

    simRatesShowTest <- capture.output(show(simRates1))
    expect_true(any(grepl("chi", simRatesShowTest)))
    expect_true(any(grepl("betaOrg", simRatesShowTest)))
    expect_true(any(grepl("betaAdp", simRatesShowTest)))
})

test_that("validateSimulatePolymeraseParams accepts good params and rejects bad", {
  expect_silent(validateSimulatePolymeraseParams(
    k = 50, ksd = 25, kMin = 17, kMax = 200, geneLen = 500,
    alpha = 1, beta = 1, zeta = 2000, zetaSd = 1000, zetaMin = 1500,
    zetaMax = 2500, cellNum = 10, polSize = 5, addSpace = 0
  ))
  expect_error(validateSimulatePolymeraseParams(k = -1, ksd = 1, kMin = 1, kMax = 2,
                                               geneLen = 10, alpha = 1, beta = 1,
                                               zeta = 100, zetaSd = 1, zetaMin = 50,
                                               zetaMax = 200, cellNum = 1, polSize = 1,
                                               addSpace = 0))
})

test_that("validateTimesToRecord checks ranges and warns on large memory", {
  expect_silent(validateTimesToRecord(timesToRecord = c(0.1, 0.5),
                                      time = 1, geneLen = 100, cellNum = 10))
  expect_error(validateTimesToRecord(timesToRecord = c(-1), time = 1, geneLen = 10, cellNum = 1))
})

test_that("validateAndLoadZetaVec reads a single-column CSV and errors on bad input", {
  tf <- tempfile(fileext = ".csv")
  write.csv(data.frame(z = rep(2, 50)), tf, row.names = FALSE)
  zz <- validateAndLoadZetaVec(tf, geneLen = 50)
  expect_type(zz, "double")
  expect_length(zz, 50)
  unlink(tf)
  # wrong length
  tf2 <- tempfile(fileext = ".csv")
  write.csv(data.frame(z = rep(2, 10)), tf2, row.names = FALSE)
  expect_error(validateAndLoadZetaVec(tf2, geneLen = 20))
  unlink(tf2)
})