# Helper function to create mock data
create_mock_data <- function() {
    edges <- data.frame(
        reg = paste0("REG", 1:20),
        tar = c(rep("GENE1", 8), rep("GENE2", 7), 
                rep("GENE3", 5))
    )
    # Mock regulatory network matrix - match number of edges
    reg_net <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
    rownames(reg_net) <- paste0("edge_", 1:20)
    colnames(reg_net) <- paste0("sample_", 1:10)    
    # Mock target-to-rows mapping
    tar_to_rows <- split(seq_len(nrow(edges)), edges$tar)
    # Mock pathways
    pathways_list <- list(
        "pathway1" = c("GENE1", "GENE2"),
        "pathway2" = c("GENE2", "GENE3")
    )
    list(
        edges = edges,
        reg_net = reg_net,
        tar_to_rows = tar_to_rows,
        pathways_list = pathways_list
    )
}


test_that("pca_pathway handles invalid inputs correctly", {
    mock_data <- create_mock_data()
    # Test empty pathways_list
    expect_error(
        pca_pathway(list(), mock_data$reg_net, mock_data$edges, 
                    mock_data$tar_to_rows),
        "'pathways_list' must be a non-empty list"
    )
    # Test non-list pathways_list
    expect_error(
        pca_pathway("not_a_list", mock_data$reg_net, mock_data$edges, 
                    mock_data$tar_to_rows),
        "'pathways_list' must be a non-empty list"
    )
    # Test non-matrix reg_net
    expect_error(
        pca_pathway(mock_data$pathways_list, "not_matrix", mock_data$edges, 
                    mock_data$tar_to_rows),
        "'reg_net' must be a numeric matrix"
    )
    # Test non-numeric reg_net
    char_matrix <- matrix(letters[1:20], nrow = 5, ncol = 4)
    expect_error(
        pca_pathway(mock_data$pathways_list, char_matrix, mock_data$edges, 
                    mock_data$tar_to_rows),
        "'reg_net' must be a numeric matrix"
    )
    # Test NAs in reg_net
    na_matrix <- mock_data$reg_net
    na_matrix[1,2] <- NA
    expect_error(
        pca_pathway(mock_data$pathways_list, na_matrix, mock_data$edges, 
                    mock_data$tar_to_rows),
        "Input data contains NA values."
    )
})


test_that("pca_pathway returns correct output structure", {
    mock_data <- create_mock_data()
    result <- pca_pathway(
        mock_data$pathways_list,
        mock_data$reg_net,
        mock_data$edges,
        mock_data$tar_to_rows
    )
    # Check output is data.frame
    expect_s3_class(result, "data.frame")
    # Check required columns
    expected_cols <- c("pathway", "pc1", "n_edges", 
                       "pathway_size")
    expect_equal(colnames(result), expected_cols)
    # Check number of rows matches number of pathways
    expect_equal(nrow(result), length(mock_data$pathways_list))
    # Check data types
    expect_type(result$pathway, "character")
    expect_type(result$pc1, "double")
    expect_type(result$n_edges, "integer")
    expect_type(result$pathway_size, "integer")
})

test_that("pca_pathway considers scaling and centering parameters", {
    mock_data <- create_mock_data()
    # Test with different parameter combinations
    result_scaled <- pca_pathway(
        mock_data$pathways_list,
        mock_data$reg_net,
        mock_data$edges,
        mock_data$tar_to_rows,
        scale_data = TRUE,
        center_data = TRUE
    )
    result_unscaled <- pca_pathway(
        mock_data$pathways_list,
        mock_data$reg_net,
        mock_data$edges,
        mock_data$tar_to_rows,
        scale_data = FALSE,
        center_data = FALSE
    )
    # (assuming run_pca actually uses these parameters)
    expect_s3_class(result_scaled, "data.frame")
    expect_s3_class(result_unscaled, "data.frame")
    # Test that scaling/centering produces different PC1 values
    expect_false(identical(result_scaled$pc1, result_unscaled$pc1))
})

test_that("pca_pathway handles different npcs values", {
    mock_data <- create_mock_data()
    result <- pca_pathway(
        mock_data$pathways_list,
        mock_data$reg_net,
        mock_data$edges,
        mock_data$tar_to_rows,
        npcs = 2
    )
    expect_s3_class(result, "data.frame")
})
