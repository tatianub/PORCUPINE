#' Creates gene set from a universe of genes
#'
#' This function creates random gene sets
#' 
#' @param universe Gene universe to sample random gene set from
#' @param psize Number of genes in a gene set
#' @param n_perm Number of permutations to create a random gene set
#' (default: 1000)
#' 
#' @return Random gene set
#' @export

create_gene_set <- function(universe,
                    psize,
                    n_perm = 1000) {
    gene_set <- lapply(1:n_perm,
                    function(x) sample(universe, size = psize, replace = FALSE))
    names(gene_set) <- paste0("random_", 1:n_perm)
    return(gene_set)
}

#' Run PCA analysis on random gene sets for statistical comparison
#'
#' This function generates random gene sets of various sizes and performs 
#' Principal Component Analysis to create a null distribution for comparing
#' against real pathway PCA results. This enables statistical significance
#' testing of pathway-specific patterns.
#'
#' @param reg_net Numeric matrix with samples in columns and network edges
#'   (regulators) in rows. Values represent regulatory relationships.
#' @param edges Data frame containing network edge information with columns
#'   "reg" (regulator) and "tar" (target).
#' @param results_pca_pathways Data frame output from pca_pathway function,
#'   used to determine pathway sizes for random sampling.
#' @param pathways_list Named list of character vectors containing gene
#'   identifiers for each pathway.
#' @param tar_to_rows Named list providing fast lookup from target genes to
#'   row indices. Created as: split(seq_len(nrow(edges)), edges$tar).
#' @param n_perm Integer, number of random permutations to generate for each
#'   pathway size. Default is 1000.
#' @param ncores Integer, number of CPU cores to use for parallel processing.
#'   Default is 1.
#' @param scale_data Logical, whether to scale the data (TRUE) or not (FALSE).
#'   Default is TRUE. Recommended for different measurement scales.
#' @param center_data Logical, whether to center the data (TRUE) or not 
#'   (FALSE). Default is TRUE. Recommended for PCA analysis.
#'
#' @return Data frame with the following columns:
#'   \item{pathway}{Random gene set identifier}
#'   \item{pc1}{Variance explained by the first principal component (%)}
#'   \item{n_edges}{Number of network edges used for the gene set}
#'   \item{pathway_size}{Number of genes in the random gene set}
#'
#' @examples
#' \dontrun{
#' # Generate random gene sets and run PCA
#' random_results <- pca_random(
#'   reg_net = network_matrix,
#'   edges = edge_table,
#'   results_pca_pathways = pathway_pca_results,
#'   pathways_list = my_pathways,
#'   tar_to_rows = target_lookup,
#'   n_perm = 500,
#'   ncores = 4
#' )
#' }
#' @export
pca_random <- function(reg_net,
                       edges,
                       results_pca_pathways,
                       pathways_list,
                       tar_to_rows,
                       n_perm = 1000,
                       ncores = 1,
                       scale_data = TRUE,
                       center_data = TRUE) {
    # Input validation
    if (!is.numeric(reg_net) || !is.matrix(reg_net)) {
        stop("'reg_net' must be a numeric matrix")
    }
    if (!is.data.frame(results_pca_pathways)) {
        stop("'results_pca_pathways' must be a data frame")
    }
    # Extract unique pathway sizes and create gene universe
    pathways_size <- unique(results_pca_pathways$pathway_size)
    message("The total number of unique pathway sizes to process: ", 
            length(pathways_size))
    universe <- unique(unlist(pathways_list))

    # Pre-allocate result list for efficiency
    res_pca_random <- vector("list", length(pathways_size))
    total_sizes <- length(pathways_size)
    # Generate and analyze random gene sets for each pathway size
    for (m in seq_along(pathways_size)) {
        psize <- pathways_size[m]
        # Progress indication with percentage
        progress_pct <- round((m / total_sizes) * 100, 1)
        message("Processing pathways with size: ", psize, 
                " (", m, "/", total_sizes, " - ", progress_pct, "%)")
        # Create random gene sets of current size
        random_genes <- create_gene_set(
            universe = universe, 
            psize = psize, 
            n_perm = n_perm
        )
        # Run PCA analysis on random gene sets
        res_pca <- pca_pathway(
            pathways_list = random_genes,
            reg_net = reg_net,
            edges = edges,
            tar_to_rows = tar_to_rows,
            ncores = ncores,
            scale_data = scale_data,
            center_data = center_data
        )
        # Store results
        res_pca_random[[m]] <- res_pca
        # Show completion message for major progress milestones
        if (m %% max(1, round(total_sizes / 10)) == 0 || m == total_sizes) {
            completed_pct <- round((m / total_sizes) * 100, 1)
            message("Completed ", completed_pct,
                "% of pathway sizes (", m, "/", total_sizes, ")")
        }
    }
    # Combine all results into single data frame
    res_pca_random_all <- as.data.frame(do.call("rbind", res_pca_random))
    rownames(res_pca_random_all) <- NULL
    return(res_pca_random_all)
}