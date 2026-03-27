#' @title Constructor for EpigenomicFeatures object
#'
#' @description
#' Class \code{EpigenomicFeatures}...
#'
#' @slot bigwigPlus a path to bigwig for plus strand QC recording PRO-seq read
#' counts. This can be generated with the proseq2.0 pipeline.
#' @slot bigwigMinus a path to bigwig for minus strand QC recording PRO-seq read
#' counts. This can be generated with the proseq2.0 pipeline.
#' @slot geneBodyRegions a \code{\link[GenomicRanges]{GRanges-class}} object
#' that holds the gene body regions coordinates used for gene coordinate 
#' binning for feature smoothing
#' @slot histoneInputFolder a path to folder containing the histone bigwig
#' input files
#' @slot chromosomes a character vector of the chromosomes to process (default: # chr22). NULL means all chromosomes.
#' @slot repeatsBedfile a path to a bed file containing the repeat coordinates #' for low complexity regions
#' @slot ctcfBigwig a path to a bigwig file containing the CTCF signal
#' @slot wgbsBedfile a path to a bed file containing the WGBS signal
#' @slot spliceJunctionsBedfile a path to a bed file containing the splice 
#' junction coordinates and strand information
#' @slot stemLoopsBedfile a path to a bed file containing the stem-loop 
#' coordinates and strand information
#'
#' @name EpigenomicFeatures-class
#' @rdname EpigenomicFeatures-class
#' @importClassesFrom tibble tbl_df
#' @importFrom dplyr mutate select left_join rename %>% relocate filter
#' @importFrom methods slot is slot<- validObject
#' @importFrom GenomicRanges width resize 
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors DataFrame splitAsList queryHits subjectHits mcols
#' @importFrom tibble as_tibble
#' @importFrom rtracklayer import.bw
#' @importFrom dplyr mutate select left_join rename %>% relocate filter
#' @importFrom plyranges as_granges find_overlaps_directed tile_ranges
#' @importFrom plyranges anchor_5p anchor_3p shift_upstream shift_downstream
#' @importFrom S4Vectors DataFrame splitAsList queryHits subjectHits mcols
#' @importFrom tibble as_tibble
#' @importFrom rtracklayer import.bw
#' @importFrom readr read_tsv
#' @importFrom data.table :=
#' @importFrom GenomicRanges width
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @exportClass EpigenomicFeatures
methods::setClass("EpigenomicFeatures",
    slots = c(
        bigwigPlus = "character",
        bigwigMinus = "character",
        geneBodyRegions = "GRanges",
        histoneInputFolder = "character",
        ctcfBigwig = "character",
        wgbsBedfile = "character",
        spliceJunctionsBedfile = "character",
        stemLoopsBedfile = "character",
        repeatsBedfile = "character",
        chromosomes = "character"
    ),
    contains = "VIRTUAL"
)

## Validate input files

validateEpigenomicInputFiles <- function(
    bigwigPlus,
    bigwigMinus,
    histoneInputFolder,
    ctcfBigwig,
    wgbsBedfile,
    spliceJunctionsBedfile,
    repeatsBedfile
) {
    if (!dir.exists(histoneInputFolder)) {
        stop("histoneInputFolder does not exist: ", histoneInputFolder)
    }

    histoneBwFiles <- list.files(
        path = histoneInputFolder,
        pattern = "\\.(bw|bigWig)$",
        ignore.case = TRUE,
        full.names = TRUE
    )
    if (length(histoneBwFiles) == 0) {
        stop(
            "No histone bigWig files found in histoneInputFolder: ",
            histoneInputFolder
        )
    }

    requiredFiles <- c(
        bigwigPlus = bigwigPlus,
        bigwigMinus = bigwigMinus,
        ctcfBigwig = ctcfBigwig,
        wgbsBedfile = wgbsBedfile,
        spliceJunctionsBedfile = spliceJunctionsBedfile,
        repeatsBedfile = repeatsBedfile
    )
    missingFiles <- names(requiredFiles)[!file.exists(requiredFiles)]
    if (length(missingFiles) > 0) {
        missingDetails <- paste(
            paste0(missingFiles, ": ", requiredFiles[missingFiles]),
            collapse = "\n"
        )
        stop("Missing required input file(s):\n", missingDetails)
    }
}

## Helper functions for preprocessing epigenomic features

buildChromosomeTargets <- function(chromosomes) {
    chromosomes <- as.character(chromosomes)
    chromosomes <- chromosomes[!is.na(chromosomes) & nzchar(chromosomes)]
    if (length(chromosomes) == 0) {
        return(character())
    }

    bareChromosomes <- sub("^chr", "", chromosomes, ignore.case = TRUE)
    unique(c(
        chromosomes,
        bareChromosomes,
        paste0("chr", bareChromosomes)
    ))
}

filterByChromosomes <- function(grng, chromosomes) {
    chromTargets <- buildChromosomeTargets(chromosomes)
    if (length(chromTargets) == 0) {
        return(grng)
    }

    grng[as.character(GenomicRanges::seqnames(grng)) %in% chromTargets]
}

# Preprocess histone features from ENCODE by extracting bigwig files from the folder, standardizing the chromosome style, and normalizing the signal by RPM, and filtering by chromosomes.
preprocessHistoneFeatures <- function(histoneInputFolder, chromosomes) {
    # List all histone bigWig files in the directory.
    bwFiles <- list.files(
        path = histoneInputFolder,
        pattern = "\\.(bw|bigWig)$",
        ignore.case = TRUE,
        full.names = TRUE
    )

    # Build one normalized GRanges per histone.
    histoneData <- lapply(bwFiles, function(file) {
        histoneName <- sub(
            "\\.(bw|bigWig)$",
            "",
            basename(file),
            ignore.case = TRUE
        )
        data <- rtracklayer::import.bw(file)
        data$score <- data$score * 1e6 / sum(data$score)
        GenomeInfoDb::seqlevelsStyle(data) <- "NCBI"
        chromosomes <- as.character(chromosomes)
        if (length(chromosomes) > 0) {
            data <- data[as.character(seqnames(data)) %in% chromosomes]
        }
        list(name = histoneName, data = data)
    })
    histones <- lapply(histoneData, function(x) x$data)
    names(histones) <- vapply(histoneData, function(x) x$name, character(1))
    return(histones)

}

# Preprocess CTCF feature from ENCODE by extracting bigwig file, standardizing the chromosome style, and normalizing the signal by RPM, and filtering by chromosomes
preprocessCTCFFeature <- function(ctcfBigwig, chromosomes) {
    data <- rtracklayer::import.bw(ctcfBigwig)
    data$score <- data$score * 1e6 / sum(data$score)
        
    ctcfGrng = plyranges::as_granges(data)
    GenomeInfoDb::seqlevelsStyle(ctcfGrng) <- 'NCBI'

    ctcfGrng <- filterByChromosomes(ctcfGrng, chromosomes)

    return(ctcfGrng)
}

# Preprocess WGBS file from ENCODE by reading in the bed file, standardizing the chromosome style, and filtering by chromosomes, calculate position of methylated C with coverage greater than 5, and returns GRanges object including coverage and methylation score
preprocessWGBSFeature <- function(wgbsBedfile, chromosomes) {
    wgbs = readr::read_tsv(wgbsBedfile, col_names = FALSE)
    ## set a cut for coverage
    cvrgCut <- 5
    chromTargets <- buildChromosomeTargets(chromosomes)

    # get wgbs which has coverage greater than 5
    wgbs <- wgbs %>% 
      select(-.data$X4, -.data$X7, -.data$X8, -.data$X9, -.data$X10) %>% 
      relocate(.data$X6, .after = .data$X3) %>%
      filter(.data$X5 > cvrgCut) %>%
      mutate(methylated = round(.data$X11/100*.data$X5))
 
    # get the exact position of methylated C in 1-base format 
    wgbsGrng <- wgbs %>% 
      select(-.data$X11) %>%
      filter(.data$X1 %in% chromTargets | length(chromTargets) == 0) %>%
      rename(coverage = .data$X5) %>%
      plyranges::as_granges(seqnames = .data$X1, start = (.data$X2 + 1), end = .data$X3, strand = .data$X6) %>%
      select(-.data$X2)

    GenomeInfoDb::seqlevelsStyle(wgbsGrng) <- 'NCBI'

    return(wgbsGrng)
}

# Preprocess splice junction features (ex. from ENCODE polyA-RNAseq PE data (whole cell) where read counts are aligned with star mapping converted to bedfile with junction position and strand information) by reading in the bed file, standardizing the chromosome style, and filtering by chromosomes, and returns GRanges object including junction position and strand information separated into 5' and 3' ends. GRanges object is 1-base resolution binary presence/absence of junction.
preprocessSpliceJunctionFeatures <- function(
  spliceJunctionsBedfile,
  chromosomes
) {
  sj <- readr::read_tsv(
    spliceJunctionsBedfile,
    col_names = FALSE,
    show_col_types = FALSE
  )
  if (ncol(sj) < 4) {
    stop("spliceJunctionsBedfile must contain at least 4 BED columns.")
  }

  strandCol <- if ("X6" %in% colnames(sj)) "X6" else "X4"
  strandValues <- as.character(sj[[strandCol]])
  strandValues <- dplyr::case_when(
    strandValues %in% c("1", "+", "plus", "PLUS") ~ "+",
    strandValues %in% c("-1", "-", "minus", "MINUS") ~ "-",
    TRUE ~ "*"
  )

  intronCut <- 180
  sjGrng <- sj %>%
    dplyr::mutate(strand = strandValues) %>%
    dplyr::transmute(
      seqnames = .data$X1,
      start = .data$X2 + 1L,
      end = .data$X3,
      strand = strand
    ) %>%
    plyranges::as_granges() %>%
    dplyr::filter(strand %in% c("+", "-")) %>%
    dplyr::filter(GenomicRanges::width(.) > intronCut)

  sjGrng <- filterByChromosomes(sjGrng, chromosomes)

  sj5 <- sjGrng %>%
    plyranges::anchor_5p() %>%
    dplyr::mutate(width = 1L) %>%
    unique()
  sj5$sj5 <- 1

  sj3 <- sjGrng %>%
    plyranges::anchor_3p() %>%
    dplyr::mutate(width = 1L) %>%
    unique()
  sj3$sj3 <- 1

  GenomeInfoDb::seqlevelsStyle(sj5) <- 'NCBI'
  GenomeInfoDb::seqlevelsStyle(sj3) <- 'NCBI'

  return(list(sj5 = sj5, sj3 = sj3))
}

getRPM <- function(bigwigPlus, bigwigMinus) {
  ## read in 
  bwPlus = rtracklayer::import.bw(bigwigPlus)
  bwMinus = rtracklayer::import.bw(bigwigMinus)
  
  plusSumCov = sum(bwPlus$score)
  minusSumCov = sum(bwMinus$score) 
  
  rpmBwPlus = bwPlus %>% 
    dplyr::mutate(score = score * 10^6 / plusSumCov) 
  
  rpmBwMinus = bwMinus %>% 
    dplyr::mutate(score = -1 * score * 10^6 / minusSumCov) 
  
  rpmBw <- c(rpmBwPlus, rpmBwMinus)
  rpmBw = rpmBw[rpmBw$score != 0]

  return(list(
    rpmBwPlus = rpmBwPlus,
    rpmBwMinus = rpmBwMinus
  ))
}

processBigWigs <- function(bigwigPlus, bigwigMinus) {
    commonSeqlevels <- intersect(
      GenomeInfoDb::seqlevels(bigwigPlus),
      GenomeInfoDb::seqlevels(bigwigMinus)
    )

    bigwigPlus <- GenomeInfoDb::keepSeqlevels(
      bigwigPlus,
      commonSeqlevels,
      pruning.mode = "coarse"
    )
    bigwigMinus <- GenomeInfoDb::keepSeqlevels(
      bigwigMinus,
      commonSeqlevels,
      pruning.mode = "coarse"
    )

    # merge plus and minus strand reads
    GenomicRanges::strand(bigwigPlus) <- "+"
    # clean up reads from minus strand, make negative score positive
    GenomicRanges::strand(bigwigMinus) <- "-"
    bigwigMinus$score <- abs(bigwigMinus$score)
    bigWigP3 <- c(bigwigPlus, bigwigMinus)
    GenomeInfoDb::seqlevels(bigWigP3) <- GenomeInfoDb::seqlevelsInUse(
      bigWigP3
    )

    return(bigWigP3)
}

makeGrangesBasepairResolution <- function(gr) {
  gr <- plyranges::as_granges(gr)
  if (length(gr) == 0L || all(GenomicRanges::width(gr) == 1L)) {
    return(gr)
  }
  w <- GenomicRanges::width(gr)
  idx <- rep.int(seq_along(gr), w)
  pos <- sequence(w, from = GenomicRanges::start(gr))
  out <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr)[idx],
    ranges = IRanges::IRanges(start = pos, width = 1L),
    strand = GenomicRanges::strand(gr)[idx]
  )
  S4Vectors::mcols(out) <- S4Vectors::mcols(gr)[idx, , drop = FALSE]
  GenomeInfoDb::seqinfo(out) <- GenomeInfoDb::seqinfo(gr)
  return(out)
}

# summarize raw reads count
summariseWdRc <- function(bw, grng) {
  rc <- grng %>%
    plyranges::find_overlaps_directed(bw) %>%
    dplyr::group_by(seqnames, start, end, strand, ensembl_gene_id) %>%
    dplyr::summarise(score = sum(score)) %>%
    tibble::as_tibble()
  
  # for those not overlapped with bw, means no reads count
  rc <- grng %>% 
    tibble::as_tibble() %>% 
    dplyr::left_join(rc, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
    tidyr::replace_na(list(score = 0)) 
  
  return(rc)
}

# Process PRO-seq features by reading in the plus and minus strand bigwig files, rpm normalizing the signal, and filtering by chromosomes. Create 10-bp tiled gene-body windows for smoothing. Create 1-bp backbone for feature merge. Summarize raw reads count in the gene-body windows
processPROseqFeatures <- function(bigwigPlus, bigwigMinus, geneBodyRegions, chromosomes) {

  rpmBw <- getRPM(bigwigPlus, bigwigMinus)
  processedBigWigs <- processBigWigs(rpmBw$rpmBwPlus, rpmBw$rpmBwMinus)

  # skip loess correction of bigwig signal

  chromTargets <- buildChromosomeTargets(chromosomes)

  filteredGeneBodies <- filterByChromosomes(geneBodyRegions, chromTargets)
  filteredBigWig <- filterByChromosomes(processedBigWigs, chromTargets)
  # Use 10-bp tiled gene-body windows for smoothing.
  binSize <- 10L
  gbWindows <- filteredGeneBodies %>%
    plyranges::as_granges() %>%
    plyranges::tile_ranges(width = binSize)
  if (!"ensembl_gene_id" %in% colnames(S4Vectors::mcols(filteredGeneBodies))) {
    stop("geneBodyRegions must include 'ensembl_gene_id' metadata.")
  }
  if ("partition" %in% colnames(S4Vectors::mcols(gbWindows))) {
    gbWindows$ensembl_gene_id <-
      filteredGeneBodies[gbWindows$partition]$ensembl_gene_id
  } else {
    stop("Tiled gene-body windows are missing 'partition' metadata.")
  }

  # Build 1-bp backbone for feature merge.
  gbOneBp <- makeGrangesBasepairResolution(filteredGeneBodies)
  gbForMerge <- summariseWdRc(filteredBigWig, gbOneBp)

  return(list(
    filteredGeneBodies = filteredGeneBodies,
    filteredBigWig = filteredBigWig,
    gbWindows = gbWindows,
    gbForMerge = gbForMerge
  ))
}

gaussianKernel <- function(bandwidth, y){
  y <- as.numeric(y)
  n <- length(y)
  if (n == 0L) {
    return(y)
  }

  maxR <- as.integer(floor((n - 1L) / 2L))
  if (maxR < 1L) {
    return(y)
  }

  effectiveBandwidth <- min(bandwidth, maxR / 4)
  if (effectiveBandwidth <= 0) {
    return(y)
  }

  r <- max(1L, as.integer(ceiling(4 * effectiveBandwidth)))
  r <- min(r, maxR)

  x <- seq.int(-r, r)
  gaussianTemp <- (1 / effectiveBandwidth) * exp(
    -0.5 * (x / effectiveBandwidth)^2
  )

  gaussianWeighted <- function(i, r, y, gaussianTemp) {
    sRange <- seq.int(i - r, i + r)
    z <- sum(gaussianTemp)
    sum(gaussianTemp * y[sRange]) / z
  }

  middleIdx <- seq.int(r + 1L, n - r)
  smoothedValues <- vapply(
    middleIdx,
    gaussianWeighted,
    numeric(1),
    r = r,
    y = y,
    gaussianTemp = gaussianTemp
  )

  leftIdx <- seq_len(r)
  rightIdx <- seq.int(n - r + 1L, n)
  smoothedY <- c(y[leftIdx], smoothedValues, y[rightIdx])
  return(smoothedY)
}

applyGaussianSmoothing <- function(gbwd, ctcf, histone) {
  binSize <- 10
  gbwdTbl <- tibble::as_tibble(gbwd)

    # get the ovp between gb and feature
    assignChiplikeScores <- function(ft, gbwd, rawBw, binSize, ftName){


        gbFt <- gbwd %>%
            plyranges::find_overlaps_directed(ft) %>%
            dplyr::group_by(seqnames, start, end, strand) %>%
            dplyr::summarise(ft = sum(score)) %>% #summarize bw signal in a bin
            tibble::as_tibble()
        
        gbFt <- gbwd %>%
            tibble::as_tibble() %>%
            dplyr::left_join(gbFt, by = c('seqnames', 'start', 'end', 'strand')) %>%
            tidyr::replace_na(list(ft = 0))  # replace NA as 0, handle missing values

        # pull out the vector of feature's raw data
        ftV <- gbFt %>% dplyr::pull(ft)
        
        # bandwidth selection based on raw bandwidth from meta plot and bin size of gb
        bandW <- rawBw / binSize # 100bp is the raw bandwidth (with unit of bp) selected from the meta-plot
        
        # smooth raw feature data with gaussian kernel
        ftSmoothed <- gaussianKernel(bandwidth = bandW, y = ftV)
        
        # compare smoothed data and raw data and check a narrowed region
        ftCompare <- tibble::tibble(
          index = seq_len(length(ftV)),
          original = ftV,
          smoothed = ftSmoothed
        )
        
        # assign missing values 0 and z-normalize smoothed profile
        gbwdTbl <- gbwd %>% tibble::as_tibble()
        joinBy <- intersect(
          c(
            "seqnames", "start", "end", "strand",
            "partition", "width", "ensembl_gene_id"
          ),
          colnames(gbwdTbl)
        )
        gbFt <- gbwdTbl %>%
            dplyr::left_join(gbFt, by = joinBy) %>%
            tidyr::replace_na(list(ft = 0)) %>%
            dplyr::mutate(ft = as.numeric(scale(ftSmoothed))) %>%
            dplyr::mutate(ft = tidyr::replace_na(ft, 0)) %>%
            dplyr::rename(!!ftName := ft)

        return(list(gbFt = gbFt, ftCompare = ftCompare))
    }

    # call and smooth 
    ctcfSmoothed <- assignChiplikeScores(
      ctcf,
      gbwd,
      rawBw = 100,
      binSize = 10,
      ftName = "ctcf"
    )
    smtf <- ctcfSmoothed$gbFt %>%
      dplyr::select(seqnames, start, end, strand, ensembl_gene_id, ctcf)

    histoneSmoothed <- gbwdTbl %>%
      dplyr::select(seqnames, start, end, strand, ensembl_gene_id)
    for (hsName in names(histone)) {
      hsSmoothed <- assignChiplikeScores(
        histone[[hsName]],
        gbwd,
        rawBw = 100,
        binSize = 10,
        ftName = hsName
      )$gbFt %>%
        dplyr::select(
          seqnames, start, end, strand, ensembl_gene_id,
          dplyr::all_of(hsName)
        )
      histoneSmoothed <- histoneSmoothed %>%
        dplyr::left_join(
          hsSmoothed,
          by = c(
            "seqnames", "start", "end", "strand",
            "ensembl_gene_id"
          )
        )
    }

    return(list(
      ctcfSmoothed = smtf,
      histoneSmoothed = histoneSmoothed
    ))
}

# Read in bed file of stem-loop coordinates already filtered by gini score, standardize the chromosome style, and filter by chromosomes. Return GRanges object of stem-loop coordinates with presence/absence of stem-loop indicator dms column.
preprocessStemLoopFeatures <- function(stemLoopsBedfile, chromosomes) {
    stemTbl <- readr::read_tsv(
      stemLoopsBedfile,
      col_names = FALSE,
      show_col_types = FALSE
    )
    if (ncol(stemTbl) < 4) {
      stop("stemLoopsBedfile must contain at least 4 BED columns.")
    }
    strandCol <- if ("X6" %in% colnames(stemTbl)) "X6" else "X4"
    stemLoops <- stemTbl %>%
      dplyr::mutate(strand = as.character(.data[[strandCol]])) %>%
      dplyr::mutate(
        strand = dplyr::case_when(
          strand %in% c("+", "1") ~ "+",
          strand %in% c("-", "-1") ~ "-",
          TRUE ~ "*"
        )
      ) %>%
      dplyr::transmute(
        seqnames = .data$X1,
        start = .data$X2 + 1L,
        end = .data$X3,
        strand = strand
      ) %>%
      plyranges::as_granges()

  stemLoops <- filterByChromosomes(stemLoops, chromosomes)
  stemLoops <- stemLoops %>%
    dplyr::mutate(dms = 1)
  GenomeInfoDb::seqlevelsStyle(stemLoops) <- 'NCBI'
  return(stemLoops)
}

# Read in bed file of low complexity sequence repeat coordinates from UCSC Table Browser, standardize the chromosome style, and filter by chromosomes. Return GRanges object of low complexity sequence repeat coordinates with presence/absence of repeat indicator rpts column.
preprocessRepeatsFeatures <- function(repeatsBedfile, chromosomes) {
    rpts = read_tsv(repeatsBedfile, col_names = FALSE)

    #select chromosome 1-22 and xy
    chrom <- paste0('chr', c(as.character(seq(1,22,1)), 'X', 'Y'))

    rpts <- rpts %>% 
    select(.data$X1, .data$X2, .data$X3, .data$X6) %>%
    plyranges::as_granges(seqnames = .data$X1, start = (.data$X2 + 1), end = .data$X3, strand = .data$X6) %>% # start + 1: change 0-base to 1-base
    filter(seqnames %in% chrom) %>% 
    select(-.data$X2) %>% 
    as_granges()

    GenomeInfoDb::seqlevelsStyle(rpts) <- 'NCBI'
    rpts <- filterByChromosomes(rpts, chromosomes)
    return(rpts)
}

# Produce feature matrix by expanding smoothed features to 1bp backbone and merging
produceFeatureMatrix <- function(
  ctcfSmoothed,
  histoneSmoothed,
  spliceJunctions,
  wgbs,
  stemLoops,
  repeats,
  gbsg,
  mergeGeneBodies
) {
    keyCols <- c("seqnames", "start", "end", "strand", "ensembl_gene_id")
    keyScoreCols <- c(keyCols, "score")

    safeZNorm <- function(x) {
      x <- as.numeric(x)
      x[is.na(x)] <- 0
      xSd <- stats::sd(x)
      if (is.na(xSd) || xSd == 0) {
        return(rep(0, length(x)))
      }
      (x - mean(x)) / xSd
    }

    baseTbl <- tibble::as_tibble(gbsg)
    if (!all(keyCols %in% colnames(baseTbl))) {
      stop("gbsg is missing merge key columns.")
    }
    if (!"score" %in% colnames(baseTbl)) {
      baseTbl$score <- 0
    }
    baseTbl <- baseTbl %>% dplyr::select(dplyr::all_of(keyScoreCols))
    baseGr <- plyranges::as_granges(baseTbl)
    mergeTbl <- tibble::as_tibble(mergeGeneBodies)
    if (!all(keyCols %in% colnames(mergeTbl))) {
      stop("mergeGeneBodies is missing merge key columns.")
    }
    mergeTbl <- mergeTbl %>% dplyr::select(dplyr::all_of(keyCols))
    mergeGr <- plyranges::as_granges(mergeTbl)

    sumByQuery <- function(queryIdx, values, nQuery) {
      if (length(queryIdx) == 0) {
        return(rep(0, nQuery))
      }
      summed <- rowsum(values, group = queryIdx, reorder = FALSE)
      out <- rep(0, nQuery)
      out[as.integer(rownames(summed))] <- as.numeric(summed[, 1])
      out
    }

    meanByQuery <- function(queryIdx, values, nQuery) {
      if (length(queryIdx) == 0) {
        return(rep(NA_real_, nQuery))
      }
      keep <- !is.na(values)
      if (!any(keep)) {
        return(rep(NA_real_, nQuery))
      }
      queryIdx <- queryIdx[keep]
      values <- values[keep]
      summed <- rowsum(values, group = queryIdx, reorder = FALSE)
      counted <- rowsum(
        rep(1, length(values)),
        group = queryIdx,
        reorder = FALSE
      )
      out <- rep(NA_real_, nQuery)
      idx <- as.integer(rownames(summed))
      out[idx] <- as.numeric(summed[, 1]) / as.numeric(counted[, 1])
      out
    }

    expandSmoothedToBackbone <- function(smoothedTbl, featureCols) {
      if (length(featureCols) == 0) {
        return(baseTbl[0])
      }

      smoothedSel <- smoothedTbl %>%
        dplyr::select(dplyr::all_of(c(keyCols, featureCols)))
      smoothedGr <- plyranges::as_granges(smoothedSel)
      hits <- GenomicRanges::findOverlaps(baseGr, smoothedGr)
      if (length(hits) == 0) {
        out <- baseTbl %>% dplyr::select(dplyr::all_of(keyCols))
        for (ftCol in featureCols) {
          out[[ftCol]] <- 0
        }
        return(out %>% dplyr::select(dplyr::all_of(featureCols)))
      }

      qIdx <- S4Vectors::queryHits(hits)
      sIdx <- S4Vectors::subjectHits(hits)
      qGene <- baseTbl$ensembl_gene_id[qIdx]
      sGene <- S4Vectors::mcols(smoothedGr)$ensembl_gene_id[sIdx]
      keep <- qGene == sGene
      qIdx <- qIdx[keep]
      sIdx <- sIdx[keep]

      out <- baseTbl %>%
        dplyr::select(dplyr::all_of(keyCols)) %>%
        dplyr::mutate(dummy = 0) %>%
        dplyr::select(-.data$dummy)
      nQuery <- nrow(out)
      for (ftCol in featureCols) {
        values <- S4Vectors::mcols(smoothedGr)[[ftCol]][sIdx]
        out[[ftCol]] <- meanByQuery(qIdx, values, nQuery)
      }
      fillList <- as.list(rep(0, length(featureCols)))
      names(fillList) <- featureCols
      out <- tidyr::replace_na(out, replace = fillList)
      out %>% dplyr::select(dplyr::all_of(featureCols))
    }

    summarizeOverlapWidthSignal <- function(featureGr, outCol) {
      hits <- GenomicRanges::findOverlaps(mergeGr, featureGr)
      if (length(hits) == 0) {
        return(baseTbl %>%
          dplyr::select(dplyr::all_of(keyCols)) %>%
          dplyr::mutate(!!outCol := 0) %>%
          dplyr::select(dplyr::all_of(outCol)))
      }

      qIdx <- S4Vectors::queryHits(hits)
      sIdx <- S4Vectors::subjectHits(hits)
      ovWidth <- GenomicRanges::width(
        GenomicRanges::pintersect(mergeGr[qIdx], featureGr[sIdx])
      )
      hitTbl <- mergeTbl[qIdx, keyCols, drop = FALSE] %>%
        dplyr::mutate(.ov = as.numeric(ovWidth)) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(keyCols))) %>%
        dplyr::summarise(sig = sum(.data$.ov), .groups = "drop")

      outTbl <- mergeTbl %>%
        dplyr::select(dplyr::all_of(keyCols)) %>%
        dplyr::left_join(hitTbl, by = keyCols) %>%
        tidyr::replace_na(list(sig = -1)) %>%
        dplyr::mutate(
          !!outCol := as.numeric(base::scale(sig, center = TRUE))
        ) %>%
        dplyr::select(dplyr::all_of(c(keyCols, outCol)))
      out1bpGr <- makeGrangesBasepairResolution(plyranges::as_granges(outTbl))
      out1bpTbl <- tibble::as_tibble(out1bpGr) %>%
        dplyr::select(dplyr::all_of(c(keyCols, outCol)))
      baseTbl %>%
        dplyr::select(dplyr::all_of(keyCols)) %>%
        dplyr::left_join(out1bpTbl, by = keyCols) %>%
        tidyr::replace_na(stats::setNames(list(0), outCol)) %>%
        dplyr::select(dplyr::all_of(outCol))
    }

    gbRpts <- summarizeOverlapWidthSignal(repeats, "rpts")

    summarizeMeanFeature <- function(featureGr, outCol, valueCol = "score") {
      hits <- GenomicRanges::findOverlaps(baseGr, featureGr)
      nQuery <- nrow(baseTbl)
      if (length(hits) == 0) {
        return(baseTbl %>%
          dplyr::select(dplyr::all_of(keyCols)) %>%
          dplyr::mutate(!!outCol := NA_real_))
      }
      qIdx <- S4Vectors::queryHits(hits)
      sIdx <- S4Vectors::subjectHits(hits)
      subjectTbl <- S4Vectors::mcols(featureGr)
      if (!valueCol %in% colnames(subjectTbl)) {
        valueCol <- "score"
      }
      if (!valueCol %in% colnames(subjectTbl)) {
        return(baseTbl %>%
          dplyr::select(dplyr::all_of(keyCols)) %>%
          dplyr::mutate(!!outCol := NA_real_))
      }
      values <- as.numeric(subjectTbl[[valueCol]][sIdx])
      meanVec <- meanByQuery(qIdx, values, nQuery)
      baseTbl %>%
        dplyr::select(dplyr::all_of(keyCols)) %>%
        dplyr::mutate(!!outCol := meanVec)
    }

    summarizeOverlapCount <- function(featureGr, outCol) {
      countVec <- as.numeric(GenomicRanges::countOverlaps(baseGr, featureGr))
      tibble::tibble(!!outCol := countVec)
    }

    summarizeWGBSFeature <- function(featureGr, outCol = "wgbs") {
      gbWgbs <- baseTbl %>%
        dplyr::select(dplyr::all_of(keyCols)) %>%
        plyranges::as_granges() %>%
        plyranges::find_overlaps_directed(featureGr) %>%
        dplyr::group_by(
          seqnames,
          start,
          end,
          strand,
          ensembl_gene_id
        ) %>%
        dplyr::summarise(
          coverage = sum(coverage),
          methylated = sum(methylated)
        ) %>%
        tibble::as_tibble()

      if (
        nrow(gbWgbs) == 0 ||
        !"coverage" %in% colnames(gbWgbs) ||
        !"methylated" %in% colnames(gbWgbs)
      ) {
        return(baseTbl %>%
          dplyr::select(dplyr::all_of(keyCols)) %>%
          dplyr::mutate(!!outCol := NA_real_))
      }

      gbWgbs %>%
        dplyr::mutate(
          !!outCol := as.numeric(
            scale(methylated / coverage, center = TRUE)
          )
        ) %>%
        dplyr::select(dplyr::all_of(keyCols), dplyr::all_of(outCol))
    }

    if (
      !is.list(spliceJunctions) ||
      !all(c("sj5", "sj3") %in% names(spliceJunctions))
    ) {
      stop("spliceJunctions must be a list with 'sj5' and 'sj3' GRanges.")
    }

    sj5Tbl <- summarizeOverlapCount(spliceJunctions$sj5, "sj5")
    sj3Tbl <- summarizeOverlapCount(spliceJunctions$sj3, "sj3")

    dmsKey <- c("seqnames", "start", "end", "strand")
    stemLoop1bp <- makeGrangesBasepairResolution(stemLoops)
    stemLoopTbl <- tibble::as_tibble(stemLoop1bp) %>%
      dplyr::select(dplyr::all_of(dmsKey)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(dms = 1)
    dmsTbl <- baseTbl %>%
      dplyr::select(dplyr::all_of(dmsKey)) %>%
      dplyr::left_join(stemLoopTbl, by = dmsKey) %>%
      tidyr::replace_na(list(dms = 0)) %>%
      dplyr::select(dms)

    wgbsTbl <- summarizeWGBSFeature(wgbs, "wgbs")
    wgbsTbl <- baseTbl %>%
      dplyr::select(dplyr::all_of(keyCols)) %>%
      dplyr::left_join(wgbsTbl, by = keyCols) %>%
      tidyr::replace_na(list(wgbs = 0)) %>%
      dplyr::select(wgbs)

    ctcfTbl <- expandSmoothedToBackbone(
      tibble::as_tibble(ctcfSmoothed),
      featureCols = c("ctcf")
    )

    histoneTbl <- tibble::as_tibble(histoneSmoothed)
    histoneCols <- setdiff(colnames(histoneTbl), keyScoreCols)
    wholeHs <- expandSmoothedToBackbone(histoneTbl, histoneCols)

    gbFt <- dplyr::bind_cols(
      baseTbl %>% dplyr::select(dplyr::all_of(keyScoreCols)),
      gbRpts,
      sj5Tbl,
      sj3Tbl,
      dmsTbl,
      wgbsTbl,
      ctcfTbl,
      wholeHs
    )
    return(gbFt)
}

#' @title Preprocess Epigenomic Features
#'
#' @description
#' Preprocessing
#'
#' @param bigwigPlus The path to a bigwig file for plus strand QC recording 
#' PRO-seq read counts. This can be generated with the proseq2.0 pipeline.
#' @param bigwigMinus The path to a bigwig file for minus strand QC recording #' PRO-seq read counts. This can be generated with the proseq2.0 pipeline.
#' @param geneBodyRegions A \code{\link[GenomicRanges]{GRanges-class}} object
#' that holds the gene body regions coordinates used for gene coordinate 
#' binning for feature smoothing
#' @param histoneInputFolder The path to a folder containing the histone bigwig #' input files
#' @param ctcfBigwig The path to a bigwig file containing the CTCF signal
#' @param wgbsBedfile The path to a bed file containing the WGBS signal
#' @param spliceJunctionsBedfile The path to a bed file containing the splice 
#' junction coordinates and strand information
#' @param stemLoopsBedfile The path to a bed file containing the stem-loop 
#' coordinates and strand information
#' @param repeatsBedfile The path to a bed file containing the repeat 
#' coordinates for low complexity regions
#' @param chromosomes a character vector of the chromosomes to process 
#' (default: # chr22). NULL means all chromosomes.
#' @return an \code{\link{EpigenomicFeatures-class}} object
#'
#' @examples
#'
#' @rdname EpigenomicFeatures-class
#' @export
methods::setGeneric(
    "preprocessEpigenomicFeatures",
    function(
        bigwigPlus,
        bigwigMinus,
        geneBodyRegions,
        histoneInputFolder,
        ctcfBigwig,
        wgbsBedfile,
        spliceJunctionsBedfile,
        stemLoopsBedfile,
        repeatsBedfile,
        chromosomes = c("22")
    ) {
        standardGeneric("preprocessEpigenomicFeatures")
    }
)

setMethod(
    "preprocessEpigenomicFeatures", "character",
    function(
      bigwigPlus,
      bigwigMinus,
      geneBodyRegions,
      histoneInputFolder,
      ctcfBigwig,
      wgbsBedfile,
      spliceJunctionsBedfile,
      stemLoopsBedfile,
      repeatsBedfile,
      chromosomes = c("22")
    ) {
        
    validateEpigenomicInputFiles(
        bigwigPlus = bigwigPlus,
        bigwigMinus = bigwigMinus,
        histoneInputFolder = histoneInputFolder,
        ctcfBigwig = ctcfBigwig,
        wgbsBedfile = wgbsBedfile,
        spliceJunctionsBedfile = spliceJunctionsBedfile,
        repeatsBedfile = repeatsBedfile
    )
    print("preprocessing epigenomic features")

    effectiveHistones <- preprocessHistoneFeatures(histoneInputFolder, chromosomes)
    print("preprocessed histone features")
    print(head(effectiveHistones, 5))
    ctcff <- preprocessCTCFFeature(ctcfBigwig, chromosomes)
    print("preprocessed ctcff features")
    wgbs <- preprocessWGBSFeature(wgbsBedfile, chromosomes)
    print("preprocessed wgbs features")
    spliceJunctions <- preprocessSpliceJunctionFeatures(spliceJunctionsBedfile, chromosomes) 
    print("preprocessed splice junctions features")
    stemLoops <- preprocessStemLoopFeatures(stemLoopsBedfile, chromosomes)
    print("preprocessed stem loops features")
    repeats <- preprocessRepeatsFeatures(repeatsBedfile, chromosomes)
    print("preprocessed repeats features")
    processedPROseq <- processPROseqFeatures(
      bigwigPlus,
      bigwigMinus,
      geneBodyRegions,
      chromosomes
    )
    filteredGeneBodies <- processedPROseq$filteredGeneBodies
    filteredBigWig <- processedPROseq$filteredBigWig
    gbWindows <- processedPROseq$gbWindows
    gbForMerge <- processedPROseq$gbForMerge
    # wgbs_handling.R
    print("preprocessed proseq features")
    smoothed <- applyGaussianSmoothing(
      gbWindows,
      ctcff,
      effectiveHistones
    )
    print("applied gaussian smoothing")
    ctcfSmoothed <- smoothed$ctcfSmoothed
    histoneSmoothed <- smoothed$histoneSmoothed

    featureMatrix <- produceFeatureMatrix(
      ctcfSmoothed,
      histoneSmoothed,
      spliceJunctions,
      wgbs,
      stemLoops,
      repeats,
      gbForMerge,
      filteredGeneBodies
    )
    print("produced feature matrix")
    return(featureMatrix)
})