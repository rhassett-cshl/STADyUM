#' @title Constructor for LocalElongationRates object
#'
#' @description
#' Class \code{LocalElongationRates}...
#'
#' @slot histoneInputFolder a path to folder containing the histone bigwig
#' input files
#' @slot chromosomes a character vector of the chromosomes to process (default: # chr22). NULL means all chromosomes.
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
        histoneInputFolder = "character"
        chromosomes = "character"
    ),
    contains = "VIRTUAL"
)

preprocessHistoneFeatures <- function(x, chromosomes) {
    # List all .bw files in the directory
    bw_files <- list.files(path = x, pattern = "\\.bw$", full.names = TRUE)

    # Build one normalized GRanges per histone and merge into one object.
    histone_data <- lapply(bw_files, function(file) {
        histone_name <- sub("\\.bw$", "", basename(file))
        data <- rtracklayer::import.bw(file)
        data$score <- data$score * 1e6 / sum(data$score)
        data$histone_type <- histone_name
        data
    })

    do.call(c, histone_data)

    histones = plyranges::as_granges(histone_data)
    seqlevelsStyle(histones) <- 'NCBI'

    chromosomes <- as.character(chromosomes)

    if (length(chromosomes) > 0) {
        histones <- histones[as.character(seqnames(histones)) %in% chromosomes]
    }

    # essentially at chr22_effective_histones_bigwig_rpm.Rdata
    # z standardization and smoothing

}

preprocessCTCFFeature <- function(ctcf_bwfile, chromosomes) {

    data <- rtracklayer::import.bw(ctcf_bwfile)
    data$score <- data$score * 1e6 / sum(data$score)
        
    ctcf_grng = plyranges::as_granges(data)
    seqlevelsStyle(ctcf_grng) <- 'NCBI'

    chromosomes <- as.character(chromosomes)

    if (length(chromosomes) > 0) {
        ctcf_grng <- ctcf_grng[as.character(seqnames(ctcf_grng)) %in% chromosomes]
    }

    return(ctcf_grng)
}
# Measures the DNA methylation (DNAm) feature 
# TODO: filter by chromosome?
preprocessWGBSFeature <- function(wgbs_bedfile, chromosomes) {
    wgbs = readr::read_tsv(wgbs_bedfile, col_names = F)
    ## set a cut for coverage
    cvrg_cut <- 5
    sel_chr <- c(paste0('chr', seq(1,22)), 'chrX', 'chrY')

    # get wgbs which has coverage greater than 5
    wgbs <- wgbs %>%
    dplyr::select(-X4, -X7, -X8, -X9, -X10, -X12, -X13, -X14) %>%
    dplyr::relocate(X6, .after = X3) %>%
    dplyr::filter(X5 > cvrg_cut) %>%
    dplyr::mutate(score = X11)

    # get the exact position of methylated C in 1-base format
    wgbs_grng <- wgbs %>%
    dplyr::select(-X11) %>%
    dplyr::filter(X1 %in% sel_chr) %>%
    dplyr::rename(coverage = X5) %>%
    plyranges::as_granges(seqnames = X1, start = (X2 + 1), end = X3, strand = X6) %>%
    dplyr::select(-X2, -coverage)

    # sanity check that all sites are 'CG'
    hsapiens <- BSgenome.Hsapiens.UCSC.hg38
    extend_wgbs = wgbs_grng[1:1e3] %>%
    GenomicRanges::resize(width = 2, fix = "start")

    GenomeInfoDb::seqlevelsStyle(extend_wgbs) = 'UCSC'
    allseq <- BSgenome::getSeq(hsapiens, extend_wgbs) %>% as.character()

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

# preprocess stem-loop, 3+5 splice sites, low complx, proseq data

#' @title Preprocess Epigenomic Features
#'
#' @description
#' Preprocessing
#'
#' @param x The path to a folder containing the histone bigwig input files
#' @param ctcf_bwfile The path to a bigwig file containing the CTCF signal
#' @param wgbs_bedfile The path to a bed file containing the WGBS signal
#' @param chromosomes a character vector of the chromosomes to process 
#' (default: # chr22). NULL means all chromosomes.
#' @return an \code{\link{LocalElongationRates-class}} object
#'
#' @examples
#'
#' @rdname LocalElongationRates-class
#' @export
setMethod(
    "preprocessEpigenomicFeatures", "character",
    function(x, chromosomes = c("22")) { 

    effectiveHistones <- preprocessHistoneFeatures(x, chromosomes)

    ctcff <- preprocessCTCFFeature(ctcf_bwfile, chromosomes)
    wgbs <- preprocessWGBSFeature(wgbs_bedfile, chromosomes)

})