#' @title Constructor for LocalElongationRates object
#'
#' @description
#' Class \code{LocalElongationRates}...
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
#' @slot spliceJunctionGranges a \code{\link[GenomicRanges]{GRanges-class}} #' object containing the splice junction coordinates and strand information
#' @slot stemLoopGranges a \code{\link[GenomicRanges]{GRanges-class}} object #' containing the stem-loop coordinates and strand information
#'
#' @name LocalElongationRates-class
#' @rdname LocalElongationRates-class
#' @importClassesFrom tibble tbl_df
#' @importFrom dplyr mutate select left_join rename %>%
#' @importFrom methods slot is slot<- validObject
#' @importFrom S4Vectors DataFrame splitAsList queryHits subjectHits mcols
#' @importFrom tibble as_tibble
#' @importFrom rtracklayer import.bw
#' @importFrom GenomicRanges width
#' @exportClass LocalElongationRates
methods::setClass("LocalElongationRates",
    slots = c(
        bigwigPlus = "character",
        bigwigMinus = "character",
        geneBodyRegions = "GRanges",
        histoneInputFolder = "character",
        ctcfBigwig = "character",
        wgbsBedfile = "character",
        spliceJunctionGranges = "GRanges",
        stemLoopGranges = "GRanges",
        repeatsBedfile = "character",
        chromosomes = "character"
    ),
    contains = "VIRTUAL"
)

# similar to estimateTranscriptionRates function, add inputValidationChecks like checkForOverlappingRegions function for geneBodyRegions

preprocessHistoneFeatures <- function(histoneInputFolder, chromosomes) {
    # List all .bw files in the directory
    bwFiles <- list.files(path = histoneInputFolder, pattern = "\\.bw$", full.names = TRUE)

    # Build one normalized GRanges per histone and merge into one object.
    histoneData <- lapply(bwFiles, function(file) {
        histoneName <- sub("\\.bw$", "", basename(file))
        data <- rtracklayer::import.bw(file)
        data$score <- data$score * 1e6 / sum(data$score)
        data$histoneType <- histoneName
        data
    })

    do.call(c, histoneData)

    histones = plyranges::as_granges(histoneData)
    seqlevelsStyle(histones) <- 'NCBI'

    chromosomes <- as.character(chromosomes)

    if (length(chromosomes) > 0) {
        histones <- histones[as.character(seqnames(histones)) %in% chromosomes]
    }

    # essentially at chr22_effective_histones_bigwig_rpm.Rdata
    # z standardization and smoothing

}

preprocessCTCFFeature <- function(ctcfBigwig, chromosomes) {

    data <- rtracklayer::import.bw(ctcfBigwig)
    data$score <- data$score * 1e6 / sum(data$score)
        
    ctcfGrng = plyranges::as_granges(data)
    seqlevelsStyle(ctcfGrng) <- 'NCBI'

    chromosomes <- as.character(chromosomes)

    if (length(chromosomes) > 0) {
        ctcfGrng <- ctcfGrng[as.character(seqnames(ctcfGrng)) %in% chromosomes]
    }

    return(ctcfGrng)
}
# Measures the DNA methylation (DNAm) feature 
# TODO: filter by chromosome?
preprocessWGBSFeature <- function(wgbsBedfile, chromosomes) {
    wgbs = readr::read_tsv(wgbsBedfile, col_names = F)
    ## set a cut for coverage
    cvrgCut <- 5
    selChr <- c(paste0('chr', seq(1,22)), 'chrX', 'chrY')

    # get wgbs which has coverage greater than 5
    wgbs <- wgbs %>%
    dplyr::select(-X4, -X7, -X8, -X9, -X10, -X12, -X13, -X14) %>%
    dplyr::relocate(X6, .after = X3) %>%
    dplyr::filter(X5 > cvrgCut) %>%
    dplyr::mutate(score = X11)

    # get the exact position of methylated C in 1-base format
    wgbsGrng <- wgbs %>%
    dplyr::select(-X11) %>%
    dplyr::filter(X1 %in% selChr) %>%
    dplyr::rename(coverage = X5) %>%
    plyranges::as_granges(seqnames = X1, start = (X2 + 1), end = X3, strand = X6) %>%
    dplyr::select(-X2, -coverage)

    # sanity check that all sites are 'CG'
    hsapiens <- BSgenome.Hsapiens.UCSC.hg38
    extendWgbs = wgbsGrng[1:1e3] %>%
    GenomicRanges::resize(width = 2, fix = "start")

    GenomeInfoDb::seqlevelsStyle(extendWgbs) = 'UCSC'
    allseq <- BSgenome::getSeq(hsapiens, extendWgbs) %>% as.character()

    allseq %>% unique() ## Yes, all is  'CG'

    bad <- allseq != "CG"
    if (any(bad)) {
        stop(
            paste0(
            "WGBS sanity check failed: ", sum(bad), " of ", length(allseq),
            " sampled sites are not 'CG'. Non-CG motifs: ",
            paste(unique(allseq[bad]), collapse = ", "),
            ". Check 0/1-based coordinates, strand assignment, and BED columns."
            )
        )
    }
}

preprocessSpliceJunctionFeatures <- function(spliceJunctionGranges, chromosomes) {
    spliceJunctionGranges <- spliceJunctionGranges[as.character(seqnames(spliceJunctionGranges)) %in% chromosomes]

    junctionCut <- 50

    # how to deal with duplicated splicing junctions
    up5 = promoters(spliceJunctionGrng, upstream = 0, downstream = 50) %>% unique
    down5 = spliceJunctionGrng %>% anchor_5p() %>% mutate(width = 1) %>% 
    shift_downstream(1) %>% resize(width = 50, fix = 'start') %>% unique
    width(down5) %>% summary

    up3 = spliceJunctionGrng %>% anchor_3p() %>% mutate(width = 1) %>% 
    shift_upstream(200)%>% resize(width = 200, fix = 'start') %>% unique
    down3 = spliceJunctionGrng %>% anchor_3p() %>% mutate(width = 1) %>%  
    shift_downstream(0) %>% # aviod overlap 1 base with up3 region 
    resize(width = 200, fix = 'start') %>% unique

    return(list(up5 = up5, down5 = down5, up3 = up3, down3 = down3))
}

get_rpm <- function(bigwigPlus, bigwigMinus) {
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
    # merge plus and minus strand reads
    strand(bigWigPlus) <- "+"
    # clean up reads from minus strand, make negative score positive
    strand(bigWigMinus) <- "-"
    bigWigMinus$score <- abs(bigWigMinus$score)
    bigWigP3 <- c(bigWigPlus, bigWigMinus)
    seqlevels(bigWigP3) <- seqlevelsInUse(bigWigP3)

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
    seqnames = as.character(GenomicRanges::seqnames(gr))[idx],
    ranges = IRanges::IRanges(start = pos, width = 1L),
    strand = as.character(BiocGenerics::strand(gr))[idx]
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

processPROseqFeatures <- function(bigwigPlus, bigwigMinus, geneBodyRegions, chromosomes) {

  rpm_bw <- get_rpm(bigwigPlus, bigwigMinus)
  processedBigWigs <- processBigWigs(rpm_bw$rpmBwPlus, rpm_bw$rpmBwMinus)

  # skip loess correction of bigwig signal for now

  filteredGeneBodies <- geneBodyRegions[geneBodyRegions$chromosome %in% chromosomes, ]
  filteredBigWig <- processedBigWigs[processedBigWigs$seqnames %in% chromosomes, ]

  return(list(
    filteredGeneBodies = filteredGeneBodies,
    filteredBigWig = filteredBigWig
  ))
}

gaussian_kernel <- function(bandwidth, y){
  r = 4 * bandwidth
  x = seq(-4*bandwidth, 4*bandwidth, 1)
  gaussian_temp <-  (1/bandwidth) * exp(-1/2 * (x / bandwidth)^2 )
  
  # first radius
  fst_r_range <- c(1:r)
  # # second radius
  snd_r_range <- c((length(y) - r + 1): length(y))
    
  # middle
  gaussian_weighted <- function(i, r, y, guassian_temp){
    s_range <- c((i - r) : (i + r))
    z <- sum(gaussian_temp)
    s_value <- sum(gaussian_temp * y[s_range]) / z
    
    return(s_value)
  }

  smoothed_values <-
    sapply(c((r + 1): (length(y)- r)), gaussian_weighted, 
           r = r, y = y, gaussian_temp)
  
  smoothed_y <- c(y[fst_r_range], smoothed_values, y[snd_r_range])
  
  return(smoothed_y)
}

applyGaussianSmoothing <- function(filteredGeneBodies, filteredBigWig, ctcf, histone) {
  bin_size <- 10
  gbwd <- filteredGeneBodies %>% 
    plyranges::as_granges() %>% 
    plyranges::tile_ranges(., width = bin_size)
  gbwd$ensembl_gene_id <- filteredGeneBodies[gbwd$partition]$ensembl_gene_id

    # get the ovp between gb and feature
    assign_chiplike_scores <- function(ft, gbwd, raw_bw, bin_size, ft_name){


        gb_ft <- gbwd %>%
            plyranges::find_overlaps_directed(ft) %>%
            dplyr::group_by(seqnames, start, end, strand) %>%
            dplyr::summarise(ft = sum(score)) %>% #summarize bw signal in a bin
            tibble::as_tibble()
        
        gb_ft <- gbwd %>%
            tibble::as_tibble() %>%
            dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand')) %>%
            tidyr::replace_na(list(ft = 0))  # replace NA as 0, handle missing values
        
        # pull out the vector of feature's raw data
        ft_v <- gb_ft %>% dplyr::pull(ft)
        
        # bandwidth selection based on raw bandwidth from meta plot and bin size of gb
        band_w <- raw_bw / bin_size # 100bp is the raw bandwidth (with unit of bp) selected from the meta-plot
        
        # smooth raw feature data with gaussian kernel
        ft_smoothed <- gaussian_kernel(bandwidth = band_w, y = ft_v)
        
        # compare smoothed data and raw data and check a narrowed region
        ft_compare <- tibble(index = seq(1, length(ft_v), 1), 
                            original = ft_v, smoothed = ft_smoothed)
        
        #get the gb-feature overlap with smoothed feature and normalize scores from 0 to 1
        # gb_ft <- gb_ft %>%
        #   dplyr::mutate(ft = ft_smoothed) %>%
        #   dplyr::mutate(ft = (ft - min(ft)) / (max(ft) - min(ft))) %>%
        #   dplyr::rename(!!ft_name := ft)

        # get the gb-features overlap with smoothed features and assigning missing values 0  before z normalization
        gb_ft <- gbwd %>%
            tibble::as_tibble() %>%
            dplyr::left_join(gb_ft, by = c('seqnames', 'start', 'end', 'strand', 'partition', 'width', 'ensembl_gene_id')) %>%
            tidyr::replace_na(list(ft = 0)) %>%
            dplyr::mutate(ft = scale(ft_smoothed)) %>% 
            dplyr::rename(!!ft_name := ft)

        return(list(gb_ft = gb_ft, ft_compare = ft_compare))
    }

    # call and smooth 
    ctcf_smoothed <- assign_chiplike_scores(ctcf, gbwd, raw_bw = 100, bin_size = 10, ft_name = 'ctcf')
    smtf <- tf_smoothed$gb_ft

    histone_smoothed <- assign_chiplike_scores(histone, gbwd, raw_bw = 100, bin_size = 10, ft_name = 'histone')
    smhistone <- histone_smoothed$gb_ft

    return(list(
      ctcf_smoothed = ctcf_smoothed,
      histone_smoothed = histone_smoothed
    ))
}

preprocessStemLoopFeatures <- function(stemLoopGranges, chromosomes) {
    stemLoopGranges <- stemLoopGranges[as.character(seqnames(stemLoopGranges)) %in% chromosomes]
    return(stemLoopGranges)
}

preprocessRepeatsFeatures <- function(repeatsBedfile, chromosomes) {
    rpts = read_tsv(repeatsBedfile, col_names = F)

    #select chromosome 1-22 and xy
    chrom = c(as.character(seq(1,22,1)), 'X', 'Y') %>% paste0('chr', .)

    rpts <- rpts %>% 
    select(X1, X2, X3, X6) %>%
    plyranges::as_granges(seqnames = X1, start = (X2 + 1), end = X3, strand = X6) %>% # start + 1: change 0-base to 1-base
    filter(seqnames %in% chrom) %>% 
    select(-X2) %>% 
    as_granges()

    seqlevelsStyle(rpts) <- 'NCBI'
    rpts <- rpts[as.character(seqnames(rpts)) %in% chromosomes]
    return(rpts)
}

produceFeatureMatrix <- function(ctcf_smoothed, histone_smoothed, spliceJunctions, wgbs, stemLoops, rpts, gbsg) {

    ctcf <- ctcf_smoothed %>% dplyr::select(-width)
    histone <- histone_smoothed %>% dplyr::select(-width)


    gb_rpts <- gbsg %>% 
    plyranges::find_overlaps_directed(repeats) %>% 
    unique() %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-width) %>% 
    tibble::add_column(rpts = 1)

    gb_rpts <- gbsg %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-width) %>% 
    dplyr::left_join(gb_rpts, by = c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id', 'score')) %>%
    tidyr::replace_na(list(rpts = 0)) %>% 
    dplyr::mutate(rpts = scale(rpts)) # z norm

    gb_ft <- gb_rpts %>%
    dplyr::inner_join(ctcf, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
    dplyr::inner_join(histone, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
    dplyr::inner_join(spliceJunctions, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>%
    dplyr::inner_join(stemLoops, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id')) %>% 
    dplyr::inner_join(wgbs, by =c('seqnames', 'start', 'end', 'strand', 'ensembl_gene_id'))

    key_cols <- c(
      'seqnames',
      'start',
      'end',
      'strand',
      'ensembl_gene_id',
      'score'
    )
    feature_cols <- setdiff(colnames(gb_ft), key_cols)
    feature_cols_to_process <- setdiff(feature_cols, 'rpts')

    safe_z_norm <- function(x) {
      x <- as.numeric(x)
      x[is.na(x)] <- 0
      x_sd <- stats::sd(x)
      if (is.na(x_sd) || x_sd == 0) {
        return(rep(0, length(x)))
      }
      (x - mean(x)) / x_sd
    }

    gb_ft <- gb_ft %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(feature_cols_to_process),
          ~ tidyr::replace_na(.x, 0)
        )
      ) %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(feature_cols_to_process),
          safe_z_norm
        )
      )

    return(gb_ft)
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
#' @param spliceJunctionGranges A \code{\link[GenomicRanges]{GRanges-class}} #' #' object containing the splice junction coordinates and strand information
#' @param stemLoopGranges A \code{\link[GenomicRanges]{GRanges-class}} object #' containing the stem-loop coordinates and strand information
#' @param repeatsBedfile The path to a bed file containing the repeat 
#' coordinates for low complexity regions
#' @param chromosomes a character vector of the chromosomes to process 
#' (default: # chr22). NULL means all chromosomes.
#' @return an \code{\link{LocalElongationRates-class}} object
#'
#' @examples
#'
#' @rdname LocalElongationRates-class
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
        spliceJunctionGranges,
        stemLoopGranges,
        repeatsBedfile,
        chromosomes = c("22")
    ) {
        standardGeneric("preprocessEpigenomicFeatures")
    }
)

setMethod(
    "preprocessEpigenomicFeatures", "character",
    function(bigwigPlus, bigwigMinus, geneBodyRegions, histoneInputFolder, ctcfBigwig, wgbsBedfile, spliceJunctionGranges, stemLoopGranges, repeatsBedfile, chromosomes = c("22")) { 
        
    browser()
    print("preprocessing epigenomic features")

    effectiveHistones <- preprocessHistoneFeatures(histoneInputFolder, chromosomes)
    print("preprocessed histone features")
    print(head(effectiveHistones, 5))
    ctcff <- preprocessCTCFFeature(ctcfBigwig, chromosomes)
    print("preprocessed ctcff features")
    wgbs <- preprocessWGBSFeature(wgbsBedfile, chromosomes)
    print("preprocessed wgbs features")
    # returns list of up5, down5, up3, down3
    spliceJunctions <- preprocessSpliceJunctionFeatures(spliceJunctionGranges, chromosomes) 
    print("preprocessed splice junctions features")
    stemLoops <- preprocessStemLoopFeatures(stemLoopGranges, chromosomes)
    print("preprocessed stem loops features")
    repeats <- preprocessRepeatsFeatures(repeatsBedfile, chromosomes)
    print("preprocessed repeats features")
    proseq_processed <- processPROseqFeatures(
      bigwigPlus,
      bigwigMinus,
      geneBodyRegions,
      chromosomes
    )
    filteredGeneBodies <- proseq_processed$filteredGeneBodies
    filteredBigWig <- proseq_processed$filteredBigWig
    print("preprocessed proseq features")
    smoothed <- applyGaussianSmoothing(
      filteredGeneBodies,
      filteredBigWig,
      ctcff,
      effectiveHistones
    )
    print("applied gaussian smoothing")
    ctcf_smoothed <- smoothed$ctcf_smoothed
    histone_smoothed <- smoothed$histone_smoothed

    featureMatrix <- produceFeatureMatrix(ctcf_smoothed, histone_smoothed, spliceJunctions, wgbs, stemLoops, repeats, filteredGeneBodies)
    print("produced feature matrix")
    return(featureMatrix)
})