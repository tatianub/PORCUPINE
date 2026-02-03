
# Helper functions to create mock data
create_mock_pca_pathway_result <- function(pathway = "pathway1",
    pathway_size = 10, pc1 = 25) {
    data.frame(
        pathway = pathway,
        pathway_size = pathway_size,
        pc1 = pc1,
        stringsAsFactors = FALSE
    )
}

create_mock_pca_random_result <- function(pathway_size = 10,
                                           n_samples = 100) {
  data.frame(
        pathway_size = rep(pathway_size, n_samples),
        pc1 = rnorm(n_samples, mean = 25, sd = 0.5),
        stringsAsFactors = FALSE
    )
}

create_mock_multiple_pathways <- function() {
    data.frame(
        pathway = c("pathway1", "pathway2", "pathway3", 
                    "pathway4"),
        pathway_size = c(10, 10, 15, 15),
        pc1 = c(25, 27, 23, 28),
        stringsAsFactors = FALSE
    )
}

create_mock_multiple_random <- function() {
    # Create random results for pathway sizes 10 and 15
    size_10 <- data.frame(
        pathway_size = rep(10, 200),
        pc1 = rnorm(200, mean = 25, sd = 0.5)
    )
    size_15 <- data.frame(
        pathway_size = rep(15, 200),
        pc1 = rnorm(200, mean = 25, sd = 0.5)
    )
    rbind(size_10, size_15)
}

# Tests for calculate_statistics function
test_that("calculate_statistics works with valid inputs", {
    # Create mock data
    pathway_result <- create_mock_pca_pathway_result()
    random_result <- create_mock_pca_random_result()
    # Test the function
    result <- calculate_statistics(pathway_result, random_result)
    # Check structure
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 1)
    expect_equal(ncol(result), 5)
    expect_named(result, 
                 c("pathway", "pathway_size", "pc1", "pval", "es"))
    # Check data types
    expect_type(result$pathway, "character")
    expect_type(result$pathway_size, "double")
    expect_type(result$pc1, "double")
    expect_type(result$pval, "double")
    expect_type(result$es, "double")
    # Check values
    expect_equal(result$pathway, "pathway1")
    expect_equal(result$pathway_size, 10)
    expect_equal(result$pc1, 25)
    expect_true(result$pval >= 0 && result$pval <= 1)
    expect_true(is.finite(result$es))
})




test_that("calculate_statistics fails with invalid inputs", {
    # Test with missing columns in pathway result
    invalid_pathway <- data.frame(pathway = "test", pc1 = 25)
    random_result <- create_mock_pca_random_result()
    expect_error(calculate_statistics(invalid_pathway, random_result))
    # Test with missing pc1 column in random result
    pathway_result <- create_mock_pca_pathway_result()
    invalid_random <- data.frame(pathway_size = rep(10, 100))
    expect_error(calculate_statistics(pathway_result, invalid_random),
                 "Missing required columns.")
})

# Tests for porcupine function
test_that("porcupine works with single pathway size", {
    # Create mock data with single pathway size
    pathways_result <- create_mock_pca_pathway_result()
    random_result <- create_mock_pca_random_result()
    # Suppress print output during testing
    result <- suppressMessages(porcupine(pathways_result, random_result))
    # Check structure
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 1)
    expect_equal(ncol(result), 5)
    expect_named(result, 
                 c("pathway", "pathway_size", "pc1", "pval", "es"))
})

test_that("porcupine works with multiple pathway sizes", {
    # Create mock data with multiple pathway sizes
    pathways_result <- create_mock_multiple_pathways()
    random_result <- create_mock_multiple_random()
    # Suppress print output during testing
    result <- suppressMessages(porcupine(pathways_result, random_result))
    # Check structure
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 4)  # Should have 4 pathways
    expect_equal(ncol(result), 5)
    # Check that all pathways are included
    expect_setequal(result$pathway, 
                    c("pathway1", "pathway2", "pathway3", "pathway4"))
    # Check that pathway sizes are preserved
    expect_setequal(result$pathway_size, c(10, 10, 15, 15))
    # Check that all p-values are valid
    expect_true(all(result$pval >= 0 & result$pval <= 1))
    # Check that all effect sizes are finite
    expect_true(all(is.finite(result$es)))
})

test_that("calculate_statistics handles empty inputs", {
    # Test with empty pathways
    empty_pathways <- data.frame(
        pathway = character(0),
        pathway_size = numeric(0),
        pc1 = numeric(0)
    )
    random_result <- create_mock_pca_random_result()
    # Expect error due to empty pathways
    expect_error(calculate_statistics(empty_pathways, random_result),
                 "Pathways results is empty.")
})

test_that("porcupine handles mismatched pathway sizes", {
    # Create pathways with sizes not present in random results
    pathways_result <- data.frame(
        pathway = "pathway1",
        pathway_size = 999,  # Size not in random results
        pc1 = 0.5,
        stringsAsFactors = FALSE
    )    
    random_result <- create_mock_pca_random_result(pathway_size = 10)
    # Expect error due to empty random results for the given size
    expect_error(porcupine(pathways_result, random_result),
                 "Permutation random results is empty.")
})

