test_that("runPORCUPINE completes full pipeline with minimal data", {
    # Create temporary test directory
    test_dir <- "tests/testdata"
    
    if (!dir.exists(test_dir)) {
        dir.create(test_dir, recursive = TRUE)
    }
    on.exit(unlink(test_dir, recursive = TRUE))
    # Create test network data - regulators x samples

    reg_net <- matrix(rnorm(1500), nrow = 30, ncol = 50)
    reg_net_file <- file.path(test_dir, "porcupine_network.RData")
    save(reg_net, file = reg_net_file)
    # Create edges file with proper structure (reg, tar, prior)
    edges_file <- file.path(test_dir, "porcupine_edges.txt")
    tfs <- rep(c("TF1", "TF2", "TF3"), each = 10)
    genes <- rep(paste0("GENE", 1:10), times = 3)
    priors <- rep(c(0, 1), length.out = 30)  # Alternate 0 and 1
    edges_data <- data.frame(
        reg = tfs,
        tar = genes,
        prior = priors
    )
    write.table(edges_data, edges_file, sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    # Create GMT file with pathways using target genes
    gmt_file <- file.path(test_dir, "porcupine_pathways.gmt")
    writeLines(c("PATHWAY1\tdescription\tGENE1\tGENE2\tGENE3\tGENE4\tGENE5",
                 "PATHWAY2\tdescription\tGENE6\tGENE7\tGENE8\tGENE9\tGENE10"),
               gmt_file)
    # Run pipeline with minimal parameters
    results <- runPORCUPINE(
        reg_net_file = reg_net_file,
        edge_file = edges_file,
        pathway_gmt_file = gmt_file,
        res_dir = test_dir,
        nperm = 10,  # Small for speed
        minSize = 3,
        maxSize = 10
    )
    # Check return structure
    expect_type(results, "list")
    expect_named(results, 
                 c("pathways_results", "random_results", "porcupine_results"))
    expect_s3_class(results$pathways_results, "data.frame")
    expect_s3_class(results$random_results, "data.frame")
    expect_s3_class(results$porcupine_results, "data.frame")
    # Check output files exist
    expect_true(file.exists(file.path(test_dir, "pathways_results.txt")))
    expect_true(file.exists(file.path(test_dir, 
                                      "pathways_results_random.txt")))
    expect_true(file.exists(file.path(test_dir, "porcupine_results.txt")))
})
