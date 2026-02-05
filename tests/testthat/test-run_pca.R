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
    set.seed(42)
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


test_that("run_pca errors on zero-variance features when scaling", {
  set.seed(1)
  x <- matrix(rnorm(5 * 10), nrow = 5, ncol = 10)
  x[1, ] <- 1  # feature 1 is constant -> sd = 0

  expect_error(
    run_pca(x, scale_data = TRUE, center_data = TRUE),
    "zero-variance"
  )
})

test_that("run_pca returns correct output structure", {
    test_data <- matrix(rnorm(400), nrow = 20, ncol = 20)
    result <- run_pca(test_data)
    # Check output is data.frame
    expect_s3_class(result, "data.frame")
    # Check required columns
    expected_cols <- c("pc1", "n_edges")
    expect_equal(colnames(result), expected_cols)
    # Check data types
    expect_type(result$pc1, "double")
    expect_type(result$n_edges, "integer")
})
