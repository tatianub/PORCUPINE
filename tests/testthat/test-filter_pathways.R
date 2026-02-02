# Test file for filter_pathways.R functions

test_that("filter_pathways works correctly", {
    # Create test data
    test_edges <- data.frame(
        reg = c("TF1", "TF2", "TF3", "TF1", "TF2"),
        tar = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
        prior = c(0, 1, 0, 1, 0)
    )
    test_pathways <- list(
        "PATHWAY1" = c("GENE1", "GENE2", "GENE6"),  # GENE6 not in edges
        "PATHWAY2" = c("GENE3", "GENE4", "GENE5"),
        "PATHWAY3" = c("GENE7", "GENE8")  # None in edges
    )
    # Test normal filtering
    result <- filter_pathways(test_pathways, test_edges)
    # Check that result is a list
    expect_type(result, "list")
    # Check that only genes present in edges are kept
    expect_equal(sort(result$PATHWAY1), c("GENE1", "GENE2"))
    expect_equal(sort(result$PATHWAY2), c("GENE3", "GENE4", "GENE5"))
    expect_false("PATHWAY3" %in% names(result))
    # Check that pathway names are preserved
    expect_true("PATHWAY1" %in% names(result))
    expect_true("PATHWAY2" %in% names(result))
    expect_false("PATHWAY3" %in% names(result))
})

test_that("filter_pathways handles situations when genes are not in edges", {
    # Create test data
    test_edges <- data.frame(
        reg = c("TF1", "TF2", "TF3", "TF1", "TF2"),
        tar = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
        prior = c(0, 1, 0, 1, 0)
    )
    test_pathways <- list(
        "PATHWAY1" = c("GENE7", "GENE8", "GENE9"), 
        "PATHWAY2" = c("GENE10", "GENE11", "GENE12")
    )
    result <- filter_pathways(test_pathways, test_edges)
    expect_type(result, "list")
    expect_length(result, 0)
})



test_that("filter_pathways_size works with default parameters", {
    test_pathways <- list(
        "SMALL" = c("GENE1", "GENE2", "GENE3", "GENE4"),
        "MEDIUM" = paste0("GENE", 1:50),
        "LARGE" = paste0("GENE", 1:200)
    )
    result <- filter_pathways_size(test_pathways)
    # should keep only MEDIUM pathway (size between 5 and 150)
    expect_type(result, "list")
    expect_length(result, 1)
    expect_true("MEDIUM" %in% names(result))
    expect_false("SMALL" %in% names(result))
    expect_false("LARGE" %in% names(result))
    expect_equal(length(result$MEDIUM), 50)
})

test_that("filter_pathways_size works with custom parameters", {
    test_pathways <- list(
        "PATHWAY1" = c("GENE1", "GENE2"),  # 2 genes
        "PATHWAY2" = c("GENE1", "GENE2", "GENE3"),  # 3 genes
        "PATHWAY3" = paste0("GENE", 1:8),  # 8 genes
        "PATHWAY4" = paste0("GENE", 1:12)  # 12 genes
    )
    # Filter with minSize = 3, maxSize = 10
    result <- filter_pathways_size(test_pathways, minSize = 3, maxSize = 10)
    expect_type(result, "list")
    expect_length(result, 2)
    expect_true("PATHWAY2" %in% names(result))
    expect_true("PATHWAY3" %in% names(result))
    expect_false("PATHWAY1" %in% names(result))  # Too small
    expect_false("PATHWAY4" %in% names(result))  # Too large
})
