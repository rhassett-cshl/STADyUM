test_that("experiment_transcription_rates validation works", {
    # Create test data
    test_counts <- data.frame(
        gene_id = c("gene1", "gene2"),
        summarized_pause_counts = c(100, 200),
        pause_length = c(200, 200),
        summarized_gb_counts = c(1000, 2000),
        gb_length = c(1000, 2000)
    )

    test_rates <- tibble::tibble(
        gene_id = c("gene1", "gene2"),
        chi = c(1, 2),
        beta_org = c(0.1, 0.2),
        beta_adp = c(0.15, 0.25),
        fk_mean = c(50, 60),
        fk_var = c(10, 20)
    )

    test_pause_regions <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(1000, 2000), width = 200),
        gene_id = c("gene1", "gene2")
    )

    test_gene_body_regions <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(1200, 2200), width = c(1000, 2000)),
        gene_id = c("gene1", "gene2")
    )

    # Test valid object
    valid_obj <- methods::new("experiment_transcription_rates",
        counts = test_counts,
        bigwig_plus = "test_plus.bw",
        bigwig_minus = "test_minus.bw",
        pause_regions = test_pause_regions,
        gene_body_regions = test_gene_body_regions,
        gene_name_column = "gene_id",
        steric_hindrance = FALSE,
        omega_scale = NULL,
        rates = test_rates
    )

    expect_true(valid(valid_obj))

    # Test invalid counts
    invalid_counts <- test_counts
    invalid_counts$gene_id <- NULL
    invalid_obj <- methods::new("experiment_transcription_rates",
        counts = invalid_counts,
        bigwig_plus = "test_plus.bw",
        bigwig_minus = "test_minus.bw",
        pause_regions = test_pause_regions,
        gene_body_regions = test_gene_body_regions,
        gene_name_column = "gene_id",
        steric_hindrance = FALSE,
        omega_scale = NULL,
        rates = test_rates
    )
    expect_error(valid(invalid_obj))
})

test_that("accessor methods work", {
    # Create test object
    test_obj <- methods::new("experiment_transcription_rates",
        counts = data.frame(
            gene_id = "gene1", summarized_pause_counts = 100,
            pause_length = 200, summarized_gb_counts = 1000,
            gb_length = 1000
        ),
        bigwig_plus = "test_plus.bw",
        bigwig_minus = "test_minus.bw",
        pause_regions = GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = 1000, width = 200),
            gene_id = "gene1"
        ),
        gene_body_regions = GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = 1200, width = 1000),
            gene_id = "gene1"
        ),
        gene_name_column = "gene_id",
        steric_hindrance = FALSE,
        omega_scale = NULL,
        rates = tibble::tibble(
            gene_id = "gene1",
            chi = 1,
            beta_org = 0.1,
            beta_adp = 0.15,
            fk_mean = 50,
            fk_var = 10
        )
    )

    # Test getRates
    rates <- getRates(test_obj)
    expect_s3_class(rates, "tbl_df")
    expect_equal(nrow(rates), 1)

    # Test getCounts
    counts <- getCounts(test_obj)
    expect_s3_class(counts, "data.frame")
    expect_equal(nrow(counts), 1)

    # Test getRegions
    regions <- getRegions(test_obj)
    expect_type(regions, "list")
    expect_length(regions, 2)
    expect_s4_class(regions$pause_regions, "GRanges")
    expect_s4_class(regions$gene_body_regions, "GRanges")
})

test_that("filterByCounts works", {
    # Create test object with multiple genes
    test_obj <- methods::new("experiment_transcription_rates",
        counts = data.frame(
            gene_id = c("gene1", "gene2"),
            summarized_pause_counts = c(100, 5), # gene2 below threshold
            pause_length = c(200, 200),
            summarized_gb_counts = c(1000, 5), # gene2 below threshold
            gb_length = c(1000, 2000)
        ),
        bigwig_plus = "test_plus.bw",
        bigwig_minus = "test_minus.bw",
        pause_regions = GenomicRanges::GRanges(
            seqnames = c("chr1", "chr1"),
            ranges = IRanges::IRanges(start = c(1000, 2000), width = 200),
            gene_id = c("gene1", "gene2")
        ),
        gene_body_regions = GenomicRanges::GRanges(
            seqnames = c("chr1", "chr1"),
            ranges = IRanges::IRanges(
                start = c(1200, 2200),
                width = c(1000, 2000)
            ),
            gene_id = c("gene1", "gene2")
        ),
        gene_name_column = "gene_id",
        steric_hindrance = FALSE,
        omega_scale = NULL,
        rates = tibble::tibble(
            gene_id = c("gene1", "gene2"),
            chi = c(1, 2),
            beta_org = c(0.1, 0.2),
            beta_adp = c(0.15, 0.25),
            fk_mean = c(50, 60),
            fk_var = c(10, 20)
        )
    )

    # Filter with threshold of 10
    filtered_obj <- filterByCounts(test_obj, min_count = 10)

    # Check that gene2 was filtered out
    expect_equal(nrow(getCounts(filtered_obj)), 1)
    expect_equal(nrow(getRates(filtered_obj)), 1)
    expect_equal(getCounts(filtered_obj)$gene_id, "gene1")
})

test_that("exportToCSV works", {
    # Create test object
    test_obj <- methods::new("experiment_transcription_rates",
        counts = data.frame(
            gene_id = "gene1", summarized_pause_counts = 100,
            pause_length = 200, summarized_gb_counts = 1000,
            gb_length = 1000
        ),
        bigwig_plus = "test_plus.bw",
        bigwig_minus = "test_minus.bw",
        pause_regions = GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = 1000, width = 200),
            gene_id = "gene1"
        ),
        gene_body_regions = GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = 1200, width = 1000),
            gene_id = "gene1"
        ),
        gene_name_column = "gene_id",
        steric_hindrance = FALSE,
        omega_scale = NULL,
        rates = tibble::tibble(
            gene_id = "gene1",
            chi = 1,
            beta_org = 0.1,
            beta_adp = 0.15,
            fk_mean = 50,
            fk_var = 10
        )
    )

    # Create temporary file
    temp_file <- tempfile(fileext = ".csv")

    # Export to CSV
    exportToCSV(test_obj, temp_file)

    # Check that file was created and contains correct data
    expect_true(file.exists(temp_file))
    imported_data <- read.csv(temp_file)
    expect_equal(nrow(imported_data), 1)
    expect_equal(imported_data$gene_id, "gene1")

    # Clean up
    unlink(temp_file)
})
