test_that("run_pca validates matrix input correctly", {
    # Test non-matrix data
    expect_error(run_pca(c(1, 2, 3)), "Input 'data' must be a numeric matrix")
    expect_error(
        run_pca(data.frame(a = 1:100, b = 1:100)),
        "Input 'data' must be a numeric matrix"
    )
    # Test NA values
    na_matrix <- matrix(rnorm(20), nrow = 4, ncol = 5)
    na_matrix[1, 1] <- NA
    expect_error(run_pca(na_matrix), "Input data contains NA values.")
})

test_that("run_pca scaling and centering parameters affect results correctly", {
    # Create test data with different scales
    test_data <- 
        matrix(rnorm(20 * 20), ncol = 20, nrow = 20)
    result_scaled <- run_pca(
        test_data, 
        scale_data = TRUE,
        center_data = TRUE
    )
    result_unscaled <- run_pca(
        test_data, 
        scale_data = FALSE,
        center_data = FALSE
    )
    # Results should differ significantly due to scaling
    expect_false(identical(result_scaled$pc1, result_unscaled$pc1))
})

test_that("run_pca npcs parameter works correctly", {
    test_data <- matrix(rnorm(400), nrow = 20, ncol = 20)
    # Test that different npcs values work (only PC1 is returned)
    result_1pc <- run_pca(test_data, npcs = 1)
    result_5pc <- run_pca(test_data, npcs = 5)
    # Should return same structure regardless of npcs
    expect_equal(colnames(result_1pc), colnames(result_5pc))
})

test_that("run_pca returns correct output structure", {
    test_data <- matrix(rnorm(400), nrow = 20, ncol = 20)
    result <- run_pca(test_data, npcs = 1)
    # Check output is data.frame
    expect_s3_class(result, "data.frame")
    # Check required columns
    expected_cols <- c("pc1", "n_edges")
    expect_equal(colnames(result), expected_cols)
    # Check data types
    expect_type(result$pc1, "double")
    expect_type(result$n_edges, "integer")
})
