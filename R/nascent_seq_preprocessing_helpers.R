####################################### for 1_EstimateRates part ####################################
build_readcount_regions <- function(bw_tsn,
                                    transcripts,
                                    tsn_cutoff = 5,
                                    gb_min_length = 6000,
                                    trim_len = 2000,
                                    kmax = 200) {
  # Construct genomic regions for read counting.
  # This function builds promoter-proximal pause regions and trimmed gene bodies
  # from TSN signals and transcript annotations, with basic filtering and QC stats.

  # ---- Step 1: Filter TSN ----
  bw_tsn <- bw_tsn[bw_tsn$score >= tsn_cutoff]

  # ---- Step 2: Build pause regions ----
  bw_pause <- promoters(bw_tsn, upstream = 0, downstream = kmax)
  rm_cols <- intersect(c("score", "type"), colnames(mcols(bw_pause)))
  if (length(rm_cols) > 0) mcols(bw_pause)[rm_cols] <- NULL

  # ---- Step 3: Build transcription ternimation site (3' ends, width = 1 bp)
  gngrng <- transcripts %>%
    group_by(ensembl_gene_id) %>%
    reduce_ranges_directed() %>%
    sort()
  bw_tts <- gngrng %>%
    plyranges::anchor_3p() %>%
    mutate(width = 1)
  # ---- Step 4: Build gene body regions ----
  generate_gene_body <- function(bw_tsn, bw_tts, bw_pause,
                                 gb_min_length, trim_len) {
    idx_tts <- match(bw_pause$ensembl_gene_id, bw_tts$ensembl_gene_id)
    idx_tsn <- match(bw_pause$ensembl_gene_id, bw_tsn$ensembl_gene_id)

    keep <- !is.na(idx_tts) & !is.na(idx_tsn)
    if (!any(keep)) stop("No overlapping genes between TSN/TTS and pause set.")

    bw_tts_filtered <- bw_tts[idx_tts[keep]]
    bw_tsn_filtered <- bw_tsn[idx_tsn[keep]]

    gb <- punion(bw_tsn_filtered, bw_tts_filtered, fill.gap = TRUE)
    gb$gene_id <- bw_tsn_filtered$ensembl_gene_id
    
    gb_filt <- gb[width(gb) > gb_min_length]
    gb_filt <- gb_filt - trim_len
    gb_filt
  }
  bw_gb_filtered <- generate_gene_body(
    bw_tsn, bw_tts, bw_pause, gb_min_length, trim_len
  ) 

  # ---- Step 5: Match pause to gene bodies ----
  match_idx <- match(bw_gb_filtered$gene_id, bw_pause$ensembl_gene_id)
  keep2 <- !is.na(match_idx)
  bw_gb_filtered <- bw_gb_filtered[keep2]
  bw_pause_filtered <- bw_pause[match_idx[keep2]]
  colnames(mcols(bw_pause_filtered)) <- c("gene_id", 'tss_type')
  # ---- Stats ----
  stats <- list(
    tsn_after_cutoff = length(bw_tsn),
    gene_body_kept = length(bw_gb_filtered)
  )
  message("Genes kept after length filter: ", stats$gene_body_kept)

  list(
    gene_body = bw_gb_filtered,
    pause = bw_pause_filtered,
    stats = stats
  )
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

process_sample_rate <- function(df, from = 20, to = 50, bins = 100) {
  # Add grouping variables based on rate-related metrics.
  # betaAdp, fkMean, and chi are split into quantile-based groups,
  # and fk variance is used to classify low and high SD genes.

  # Split betaAdp into five quantile-based groups (Q1–Q5)
  df$betaGroup <- cut(
    df$betaAdp,
    breaks = quantile(df$betaAdp, probs = seq(0, 1, 0.2), na.rm = TRUE),
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    include.lowest = TRUE
  )
  # Split fkMean into five quantile-based groups (Q1–Q5)
  df$fkMeanGroup <- cut(
    df$fkMean,
    breaks = quantile(df$fkMean, probs = seq(0, 1, 0.2), na.rm = TRUE),
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    include.lowest = TRUE
  )

  # Split chi into three quantile-based groups: Low, Medium, High
  df$chiGroup <- cut(
    df$chi,
    breaks = quantile(df$chi, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE)

  # Assign SD group based on the threshold
  df$fkSD <- sqrt(df$fkVar)
  threshold <- find_valley_threshold(df$fkSD, from, to, bins)
  print(threshold)
  df$sdGroup <- ifelse(df$fkSD <= threshold, "Sharp", "Broad")
  attr(df, "threshold") <- threshold
  return(df)
}
find_valley_threshold <- function(x, from, to, bins = 100) {
  hist_data <- hist(x, breaks = bins, plot = FALSE)
  mids <- hist_data$mids
  counts <- hist_data$counts
  candidate_idx <- which(mids > from & mids < to)
  # Find local minima in the histogram counts (valleys)
  minima_idx <- findpeaks(-counts[candidate_idx])[,2]
  mids[candidate_idx[minima_idx[which.min(counts[candidate_idx[minima_idx]])]]]
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
    group_by(tsn_id) %>%
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
      group_by(promoter_idx) %>%
      summarize(histone_sum_scores = mean(histone_score))

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