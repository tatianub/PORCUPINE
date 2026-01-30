# Helper function to create test data
create_test_universe <- function(size = 30) {
    paste0("GENE", 1:size)
}

# Helper function to create mock data for pca_random tests
create_mock_data <- function() {
    # Create deterministic edges ensuring all pathway genes are covered
    edges <- data.frame(
        reg = paste0("REG", 1:24),
        tar = c(
            # pathway1: 10 edges
            rep("GENE1", 3), rep("GENE2", 3), rep("GENE3", 2), 
            rep("GENE4", 2),  
            # pathway2: 6 edges
            rep("GENE5", 2), rep("GENE6", 2), rep("GENE7", 2),
            # pathway3: 6 edges
            rep("GENE8", 2), rep("GENE9", 2), rep("GENE10", 2),
            # extra gene: 2 edges
            rep("GENE11", 2)
        ),
        stringsAsFactors = FALSE
    )
    # Mock regulatory network matrix - match number of edges
    set.seed(123)  # For reproducible tests
    reg_net <- matrix(rnorm(24 * 20), nrow = 24, ncol = 20)
    rownames(reg_net) <- paste0("edge_", 1:24)
    colnames(reg_net) <- paste0("sample_", 1:20)
    # Mock target-to-rows mapping
    tar_to_rows <- split(seq_len(nrow(edges)), edges$tar)
    # Mock pathways - ensure they match the edges
    pathways_list <- list(
        # 4 genes, 10 edges
        "pathway1" = c("GENE1", "GENE2", "GENE3", "GENE4"),
        # 3 genes, 6 edges
        "pathway2" = c("GENE5", "GENE6", "GENE7"),
        # 3 genes, 6 edges
        "pathway3" = c("GENE8", "GENE9", "GENE10")
    )
    # Mock results from pca_pathway - match actual pathway sizes
    results_pca_pathways <- data.frame(
        pathway = c("pathway1", "pathway2", "pathway3"),
        pc1 = c(45.2, 38.7, 52.1),
        n_edges = c(10, 6, 6),        # Matches edge counts above
        pathway_size = c(4, 3, 3),    # Matches pathway lengths
        stringsAsFactors = FALSE
    )
    list(
        edges = edges,
        reg_net = reg_net,
        tar_to_rows = tar_to_rows,
        pathways_list = pathways_list,
        results_pca_pathways = results_pca_pathways
    )
}

# Tests for create_gene_set function
test_that("create_gene_set returns correct structure", {
    universe <- create_test_universe(10)
    result <- create_gene_set(universe, psize = 3, n_perm = 5)
    # Check output is a list
    expect_type(result, "list")
    # Check correct number of gene sets
    expect_equal(length(result), 5)
    # Check each gene set has correct size
    expect_true(all(sapply(result, length) == 3))
    # Check names are correctly formatted
    expected_names <- paste0("random_", 1:5)
    expect_equal(names(result), expected_names)
    # Check all genes are from universe
    all_genes <- unique(unlist(result))
    expect_true(all(all_genes %in% universe))
})



test_that("create_gene_set validates inputs correctly", {
    universe <- create_test_universe(10)
    # Test psize > universe size (should cause error in sample())
    expect_error(
        create_gene_set(universe, psize = 15, n_perm = 1),
        # R's sample() will give an error about insufficient elements
    )
    # Test empty universe
    expect_error(
        create_gene_set(character(0), psize = 1, n_perm = 1)
    )
    # Test psize = 0
    result_zero <- create_gene_set(universe, psize = 0, n_perm = 2)
    expect_true(all(sapply(result_zero, length) == 0))
})

test_that("create_gene_set handles default parameters", {
    universe <- create_test_universe(50)
    result <- 
        create_gene_set(universe, psize = 10)  # Using default n_perm = 1000
    expect_equal(length(result), 1000)
    expect_equal(names(result), paste0("random_", 1:1000))
})

test_that("pca_random validates function-specific inputs", {
    mock_data <- create_mock_data()
    # Test non-data.frame results_pca_pathways
    expect_error(
        pca_random(
            reg_net = mock_data$reg_net,
            edges = mock_data$edges,
            results_pca_pathways = "not_dataframe",
            pathways_list = mock_data$pathways_list,
            tar_to_rows = mock_data$tar_to_rows,
            n_perm = 10
        ),
        "'results_pca_pathways' must be a data frame"
    )
    # Test list instead of data.frame for results_pca_pathways
    expect_error(
        pca_random(
            reg_net = mock_data$reg_net,
            edges = mock_data$edges,
            results_pca_pathways = list(pathway_size = c(3, 4)),
            pathways_list = mock_data$pathways_list,
            tar_to_rows = mock_data$tar_to_rows,
            n_perm = 10
        ),
        "'results_pca_pathways' must be a data frame"
    )
})

test_that("pca_random returns correct output structure", {
    mock_data <- create_mock_data()
    # Use small n_perm for faster testing
    result <- pca_random(
        reg_net = mock_data$reg_net,
        edges = mock_data$edges,
        results_pca_pathways = mock_data$results_pca_pathways,
        pathways_list = mock_data$pathways_list,
        tar_to_rows = mock_data$tar_to_rows,
        n_perm = 5
    )
    # Check output is data.frame
    expect_s3_class(result, "data.frame")
    # Check required columns
    expected_cols <- c("pathway", "pc1", "n_edges", "pathway_size")
    expect_equal(colnames(result), expected_cols)
    # Check data types
    expect_type(result$pathway, "character")
    expect_type(result$pc1, "double")
    expect_type(result$n_edges, "integer")
    expect_type(result$pathway_size, "integer")
    # Check pathway names follow pattern
    expect_true(all(grepl("^random_", result$pathway)))
})

test_that("pca_random processes different pathway sizes correctly", {
    mock_data <- create_mock_data()
    result <- pca_random(
        reg_net = mock_data$reg_net,
        edges = mock_data$edges,
        results_pca_pathways = mock_data$results_pca_pathways,
        pathways_list = mock_data$pathways_list,
        tar_to_rows = mock_data$tar_to_rows,
        n_perm = 3
    )
    # Should have results for each unique pathway size
    unique_sizes <- unique(mock_data$results_pca_pathways$pathway_size)
    result_sizes <- unique(result$pathway_size)
    expect_setequal(result_sizes, unique_sizes)
    # Should have n_perm results for each pathway size
    for (size in unique_sizes) {
        count <- sum(result$pathway_size == size)
        expect_equal(count, 3)  
    }
})
