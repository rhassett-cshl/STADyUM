####################################### for 1_EstimateRates part ####################################

# Select TSS with maximum score per gene
# If multiple TSSs with the same score, select the one with the most upstream TSS to break ties
# On "+" strand the most upstream TSS has the smallest coordinate;
# on "-" strand it has the largest coordinate.
keep_max_tsn <- function(tsn, tsn_cutoff = 5) {
  as.data.frame(tsn) %>%
    dplyr::filter(score >= tsn_cutoff) %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::arrange(
      dplyr::desc(score),
      ifelse(strand == "+", start, -start),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    as_granges()
}

# Select the single most upstream TSS per gene (strand-aware).
# On "+" strand the most upstream TSS has the smallest coordinate;
# on "-" strand it has the largest coordinate.
keep_upstream_tss <- function(tsn) {
  as.data.frame(tsn) %>%
    dplyr::group_by(ensembl_gene_id, strand) %>%
    dplyr::slice_min(order_by = ifelse(strand == "+", start, -start),
              n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    plyranges::as_granges()
}

# Build promoter-proximal pause regions and trimmed gene body regions
# for read counting. Filters TSNs by score, constructs pause windows of
# width `kmax`, and retains gene bodies longer than `gb_min_length`
# after trimming `trim_len` from each end.
build_readcount_regions <- function(bw_tsn,
                                    transcripts,
                                    tsn_cutoff    = 5,
                                    gb_min_length = 6000,
                                    gb_trim_len   = 1250,
                                    gb_max_length = 90000,
                                    kmax          = 200) {
  if (!inherits(bw_tsn, "GRanges"))
    stop("bw_tsn must be a GRanges.")
  if (!("ensembl_gene_id" %in% colnames(mcols(bw_tsn))))
    stop("bw_tsn must contain 'ensembl_gene_id'.")
  if (!("score" %in% colnames(mcols(bw_tsn))))
    stop("bw_tsn must contain 'score'.")
  if (!inherits(transcripts, "GRanges"))
    stop("transcripts must be a GRanges.")
  if (!("ensembl_gene_id" %in% colnames(mcols(transcripts))))
    stop("transcripts must contain 'ensembl_gene_id'.")

  # Filter low-signal TSNs and build pause windows
  bw_tsn   <- bw_tsn[bw_tsn$score >= tsn_cutoff]
  bw_pause <- promoters(bw_tsn, upstream = 0, downstream = kmax)
  rm_cols  <- intersect(c("score", "type"), colnames(mcols(bw_pause)))
  if (length(rm_cols) > 0) mcols(bw_pause)[rm_cols] <- NULL

  # Derive one TTS position (3' end, width = 1) per gene from the transcript annotation
  gngrng <- transcripts %>%
    dplyr::group_by(ensembl_gene_id) %>%
    plyranges::reduce_ranges_directed() %>%
    sort()
  bw_tts <- gngrng %>%
    plyranges::anchor_3p() %>%
    mutate(width = 1)

  # Build gene bodies spanning from TSN to TTS, apply length filter and trim
  generate_gene_body <- function(bw_tsn,
                               bw_tts,
                               bw_pause,
                               gb_min_length,
                               gb_trim_len,
                               gb_max_length) {
    idx_tts <- match(bw_pause$ensembl_gene_id, bw_tts$ensembl_gene_id)
    idx_tsn <- match(bw_pause$ensembl_gene_id, bw_tsn$ensembl_gene_id)

    keep <- !is.na(idx_tts) & !is.na(idx_tsn)
    if (!any(keep)) {
      stop("No overlapping genes between TSN/TTS and pause set.")
    }

    bw_tts_filtered <- bw_tts[idx_tts[keep]]
    bw_tsn_filtered <- bw_tsn[idx_tsn[keep]]

    gb_full <- punion(bw_tsn_filtered, bw_tts_filtered, fill.gap = TRUE)
    gb_full$gene_id <- bw_tsn_filtered$ensembl_gene_id

    min_full_width <- gb_min_length + 2 * gb_trim_len
    gb_full <- gb_full[width(gb_full) >= min_full_width]
    # Exclude 1,250 bp near TSS and 1,250 bp near transcript end
    gb_trimmed <- gb_full - gb_trim_len

    # Cap at 90 kb from the TSS-proximal side
    gb_capped <- resize(
      gb_trimmed,
      width = pmin(width(gb_trimmed), gb_max_length),
      fix = "start"
    )

    gb_capped
  }

  bw_gb_filtered <- generate_gene_body(bw_tsn, bw_tts, bw_pause, gb_min_length, gb_trim_len, gb_max_length)

  # Match pause regions to the filtered gene bodies
  match_idx         <- match(bw_gb_filtered$gene_id, bw_pause$ensembl_gene_id)
  keep2             <- !is.na(match_idx)
  bw_gb_filtered    <- bw_gb_filtered[keep2]
  bw_pause_filtered <- bw_pause[match_idx[keep2]]
  names(mcols(bw_pause_filtered))[names(mcols(bw_pause_filtered)) == "ensembl_gene_id"] <- "gene_id"

  stats <- list(
    tsn_after_cutoff = length(bw_tsn),
    gene_body_kept   = length(bw_gb_filtered)
  )
  message("Genes kept after length filter: ", stats$gene_body_kept)

  list(gene_body = bw_gb_filtered, pause = bw_pause_filtered, stats = stats)
}


merge_tss_from_two_samples <- function(tsn1, tsn2, tss_distance_threshold = 30) {
  # Align transcription start sites (TSSs) between two samples.
  # For genes present in both samples, nearby TSSs are merged to a common upstream
  # position (strand-aware), and labeled as identical or distinct.

  common <- intersect(tsn1$ensembl_gene_id, tsn2$ensembl_gene_id)
  tsn1i <- tsn1[match(common, tsn1$ensembl_gene_id)]
  tsn2i <- tsn2[match(common, tsn2$ensembl_gene_id)]

  # Calculate distances between TSSs of the two samples
  d <- abs(start(tsn1i) - start(tsn2i))
  names(d) <- tsn1i$ensembl_gene_id
  close <- d < tss_distance_threshold

  # Select the most upstream position based on strand
  # On the "+" strand: choose the smaller coordinate
  # On the "-" strand: choose the larger coordinate
  is_plus <- as.character(strand(tsn1i)) == "+"
  upstream_pos <- ifelse(is_plus,
                         pmin(start(tsn1i), start(tsn2i)),
                         pmax(start(tsn1i), start(tsn2i)))

  # For close TSSs, align both samples to the upstream position
  if (any(close)) {
    ranges(tsn1i)[close] <- IRanges(start = upstream_pos[close], end = upstream_pos[close])
    ranges(tsn2i)[close] <- IRanges(start = upstream_pos[close], end = upstream_pos[close])
  }

  # Add TSS labels to indicate if the TSS is distinct or identical
  tsn1$tss <- "Distinct TSS"
  tsn2$tss <- "Distinct TSS"
  i1 <- match(common, tsn1$ensembl_gene_id)
  i2 <- match(common, tsn2$ensembl_gene_id)

  tsn1$tss[i1][close] <- "Identical TSS"
  tsn2$tss[i2][close] <- "Identical TSS"

  # Update the ranges in the original objects with aligned positions
  ranges(tsn1)[i1] <- ranges(tsn1i)
  ranges(tsn2)[i2] <- ranges(tsn2i)

  # Return updated GRanges objects
  list(tsn1 = tsn1, tsn2 = tsn2, dist = d)
}


get_promoter_df_from_stadyum_obj <- function(stadyum_obj, promoter_len=500){
  tsn <- stadyum_obj@pauseRegions %>% anchor_5p() %>% mutate(width = 1) 
  tsn <- tsn[tsn$gene_id %in% rates(stadyum_obj)$geneId]
  promoter_grng <- resize(tsn, promoter_len, fix = "center")

  promoter_grng_df <- as.data.frame(promoter_grng)
  promoter_grng_df <- promoter_grng_df %>%
  left_join(rates(stadyum_obj), 
            by = c("gene_id" = "geneId"))

  promoter_grng_df$pauseCount <- sapply(promoter_grng_df$actualPauseSiteCounts, sum)
  promoter_grng_df <- promoter_grng_df %>% dplyr::select(-actualPauseSiteCounts)
  return(list(df = promoter_grng_df, grang=promoter_grng))
}



####################################### for 2_Nucleosome ####################################

get_promoter_motif_from_stadyum_obj <- function(stadyum_obj, promoter_df){
  # This function links pause regions from a STADyUM object to nearby gene promoters.
  # For each pause site, it identifies the closest TSS within 200 bp on the same strand.
  # Promoters are then classified into motif classes based on core promoter elements.

  tsn <- stadyum_obj@pauseRegions %>% anchor_5p() %>% mutate(width = 1)
  tsn_df <- as.data.frame(tsn) %>%
      as_tibble() %>%
      dplyr::rename(seqname = seqnames) %>%
      mutate(
          seqname = paste0("chr", seqname),
          tsn_pos = start,
          tsn_id  = dplyr::row_number())

  promoter_df <- promoter_df %>%
    dplyr::rename(gene_id = gene_ensembl_merged)

  pairs <- tsn_df %>%
    inner_join(
      promoter_df,
      by = c("seqname", "gene_id", "strand")
      ) %>%
    mutate(dist_bp = abs(tsn_pos - TSS))

  nearest <- pairs %>%
    dplyr::group_by(tsn_id) %>%
    slice_min(order_by = dist_bp, with_ties = FALSE) %>%
    ungroup() %>%
    filter(dist_bp <= 200)

  nearest <- nearest %>%
    mutate(
      motif_class = case_when(
        TATA.box == 1 ~ "TATA-box",
        CCAAT.box == 1 & GC.box == 1 & Inr == 0 & TATA.box == 0 ~ "GC-CCAAT",
        CCAAT.box == 1 & GC.box == 0 & Inr == 0 & TATA.box == 0 ~ "CCAAT only",
        GC.box == 1 & Inr == 0 & TATA.box == 0 & CCAAT.box == 0 ~ "GC-box only",
        GC.box == 1 & Inr == 1 & TATA.box == 0 & CCAAT.box == 0 ~ "GC-Inr",
        Inr == 1  & GC.box == 0 & TATA.box == 0 & CCAAT.box == 0 ~ "Inr only",
        TRUE ~ "Other"
      )
    )
  return(nearest)
}

get_1d_signal_from_contact_df <- function(contacts_df_path, seq_info) {
  # Reads a 1D signal BED file and converts it into a GRanges object 
  # with proper seqinfo.
    contacts_df <- read.table(contacts_df_path, sep = "\t", header = FALSE,
                 col.names = c("chrom","start","end","name","score","strand","tss"),
                 stringsAsFactors = FALSE)

    gr <- makeGRangesFromDataFrame(
        contacts_df,
        seqnames.field = "chrom",
        start.field    = "start",
        end.field      = "end",
        strand.field   = "strand",       
        keep.extra.columns = TRUE,
        starts.in.df.are.0based = TRUE
        )
    seqinfo(gr) <- seq_info[seqnames(seqinfo(gr))]
    return (gr)
}

loess_smooth <- function(y) {
  # Applies LOESS smoothing to a 2-kb signal profile centered at the TSS to reduce noise.
  x_pos <- seq(from = -1000, to = 1000, length.out = 2000)
  fit <- loess(y ~ x_pos, span = 0.05)
  predict(fit, x_pos)
}


get_tid_from_stadyum_obj <- function(stadyum_obj){
  # Extracts a 2-kb window centered on each pause region in a STADyUM object 
  # and aligns it to the corresponding genes.
    tsn <- stadyum_obj@pauseRegions %>% anchor_5p() %>% mutate(width = 1)
    tid <- resize(tsn, width = 2000, fix = "center")

    rate_df <- rates(stadyum_obj)
    matched_idx <- match(rate_df$geneId, tid$gene_id)
    tid <- tid[matched_idx[!is.na(matched_idx)], ]
    return (tid)
}

get_matrix_from_nucleosome_signal <- function(rate_object, 
                                              nucleosome_signal_in, 
                                              seq_info, 
                                              smooth_func=loess_smooth){
  # This function extracts nucleosome signal around transcription initiation regions.
  # It builds a position-by-position matrix aligned to TID windows from a STADyUM object.
  # Both the raw signal matrix and a smoothed version are returned.
  
  gr <- get_1d_signal_from_contact_df(nucleosome_signal_in, seq_info)
  tid <- get_tid_from_stadyum_obj(rate_object)

  score_matrix <- ScoreMatrix(gr, windows = tid, 
                    weight.col = "score", 
                    strand.aware = TRUE)
  score_matrix <- as.matrix(score_matrix@.Data)
  score_matrix_smooth <- t(apply(score_matrix, 1, smooth_func))

  return(list(raw = score_matrix, smooth = score_matrix_smooth, tid ))

}

get_histone_modification_sum <- function(promoter_df,promoter_grng, meta_sel){

  for (i in 1:nrow(meta_sel)){
    histone_bw <- import.bw(meta_sel$bw[i])
    hits <- findOverlaps(promoter_grng, histone_bw)

    query_hits <- queryHits(hits)
    subject_hits <- subjectHits(hits)
    df <- data.frame(
      promoter_idx = query_hits,
      histone_score = mcols(histone_bw)$score[subject_hits]
      )
    sum_scores <- df %>%
      dplyr::group_by(promoter_idx) %>%
      dplyr::summarize(histone_sum_scores = mean(histone_score))

    promoter_df[,meta_sel$Target[i]] <- sum_scores$histone_sum_scores
  }

  return(promoter_df)

}


get_match_matrix <- function(rate_df, matrix){
  rate_ids <- rate_df$geneId
  matrix_ids <- matrix[[1]]

  idx <- match(rate_ids, matrix_ids)
  matrix_ordered <- matrix[idx, ]
  matrix_ordered <- matrix_ordered[, -1]

  return(matrix_ordered)
}