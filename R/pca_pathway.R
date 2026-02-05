#' Create Target-to-Rows Mapping
#' 
#' Creates a named list mapping each unique target node to the row indices 
#' in the edges data frame where that target appears. This is useful for 
#' quickly subsetting edges by target node.
#' 
#' @param edges A data frame containing network edges with a 'tar' column
#'              representing target nodes
#' @return A named list where names are unique target nodes and values are 
#'         integer vectors of row indices where each target appears
#' @examples 
#' # Example edges data frame
#' edges <- data.frame(
#'   src = c("A", "B", "C", "A"), 
#'   tar = c("X", "X", "Y", "Y")
#' )
#' target_mapping <- create_target_to_rows_mapping(edges)
#' # Result: list(X = c(1, 2), Y = c(3, 4))
create_target_to_rows_mapping <- function(edges) {
  if (!is.data.frame(edges)) {
    stop("edges must be a data frame")
  }
  if (!"tar" %in% colnames(edges)) {
    stop("edges must contain a 'tar' column")
  }
  tar_to_rows <- split(seq_len(nrow(edges)), edges$tar)
  return(tar_to_rows)
}



#' Run PCA analysis for a list of pathways on network edges
#'
#' This function performs Principal Component Analysis (PCA) for multiple 
#' biological pathways by extracting relevant network edges and computing
#' the first principal component to capture pathway-level variation.
#'
#' @param pathways_list Named list of character vectors, where each vector
#'   contains gene identifiers for a specific pathway
#' @param reg_net Numeric matrix with samples in columns and network edges
#'   in rows. Values represent regulatory relationships.
#' @param edges Data frame containing network edge information with columns
#'   "reg" (regulator) and "tar" (target)
#' @param tar_to_rows Named list providing fast lookup from target genes to
#'   row indices. Created as: split(seq_len(nrow(edges)), edges$tar)
#' @param ncores Integer, number of CPU cores to use for parallel processing.
#'   Default is 1.
#' @param scale_data Logical, whether to scale the data (TRUE) or not (FALSE).
#'   Default is TRUE. Recommended for different measurement scales.
#' @param center_data Logical, whether to center the data (TRUE) or not 
#'   (FALSE). Default is TRUE. Recommended for PCA analysis.
#'
#' @return Data frame with the following columns:
#'   \item{pathway}{Name of the pathway}
#'   \item{pc1}{Variance explained by the first principal component (%)}
#'   \item{n_edges}{Number of network edges used for the pathway}
#'   \item{pathway_size}{Total number of genes in the pathway}
#' @examples
#' \dontrun{
#' # Example pathway list
#' pathways <- list(
#'   "pathway1" = c("GENE1", "GENE2", "GENE3"),
#'   "pathway2" = c("GENE4", "GENE5", "GENE6")
#' )
#' 
#' # Run PCA analysis
#' results <- pca_pathway(
#'   pathways_list = pathways,
#'   reg_net = my_network,
#'   edges = edge_table,
#'   tar_to_rows = target_lookup,
#'   ncores = 4
#' )
#' }
#' @export
pca_pathway <- function(pathways_list,
                        reg_net,
                        edges,
                        tar_to_rows,
                        ncores = 1,
                        scale_data = TRUE,
                        center_data = TRUE) {
    # Input validation
    if (!is.list(pathways_list) || length(pathways_list) == 0) {
        stop("'pathways_list' must be a non-empty list")
    }
    if (!is.numeric(reg_net) || !is.matrix(reg_net)) {
        stop("'reg_net' must be a numeric matrix")
    }
    # Run PCA analysis for each pathway using parallel processing
    res <- parallel::mclapply(pathways_list, function(pathway) {
        # Fast lookup for target genes
        idx <- unlist(tar_to_rows[pathway], use.names = FALSE)
        # Remove NA indices that may occur if genes are not found
        idx <- idx[!is.na(idx)]
        # Skip pathways with insufficient edges
        if (length(idx) < 2) {
            return(data.frame(
                pc1 = NA_real_, 
                n_edges = 0
            ))
        }
        # Extract subnet for current pathway
        subnet <- reg_net[idx, , drop = FALSE]
        # Perform PCA analysis
        run_pca(
            subnet,
            scale_data = scale_data,
            center_data = center_data
        )
    }, mc.cores = ncores)
    # Combine results and add metadata
    res <- as.data.frame(do.call("rbind", res))
    # Handle cases where pathways_list might not have names
    pathway_names <- names(pathways_list)
    if (is.null(pathway_names)) {
        pathway_names <- paste0("pathway_", 1:length(pathways_list))
    }
    res$pathway <- pathway_names
    res$pathway_size <- lengths(pathways_list)
    res <- res[, c("pathway", "pc1", "n_edges", "pathway_size")]
    return(res)
}