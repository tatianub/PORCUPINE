#' Run PCA analysis on the data
#'
#' This function performs PCA analysis on the data. 
#'
#' @param data Numeric matrix with samples in columns, and features in rows
#' @param scale_data Logical, whether to scale the data (TRUE) or not (FALSE).
#'   Default is TRUE.
#' @param center_data Logical, whether to center the data (TRUE) or not (FALSE).
#'   Default is TRUE.
#' @param npcs Integer, number of principal components to compute. Default is 1.
#'
#' @return Data frame with PCA results containing
#' 
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' result <- run_pca(data, npcs = 1)
#' }
#'
#' @export
run_pca <- function(data, 
                    scale_data = TRUE, 
                    center_data = TRUE, 
                    npcs = 1) {
    # Input validation
    if (!is.numeric(data) || !is.matrix(data)) {
        stop("Input 'data' must be a numeric matrix")
    }
    if (any(is.na(data))) {
        stop("Input data contains NA values.")
    }
    # Transpose data for PCA (samples as rows, features as columns)
    data_t <- t(data)
    # Perform PCA using irlba
    res_pca <- irlba::prcomp_irlba(
        data_t, 
        n = npcs,
        scale. = scale_data,
        center = center_data
    )
    # Extract variance explained by the first PC
    pc1 <- summary(res_pca)$importance[2, 1] * 100
    # Create result data frame
    pca_result <- data.frame(
        pc1 = pc1,
        n_edges = ncol(data_t)
    )
    return(pca_result)
}