#' Run Complete PORCUPINE Analysis Pipeline
#' 
#' Performs the complete PORCUPINE analysis including network loading, pathway
#' filtering, PCA analysis, permutation testing, and result calculation.
#' Results are saved to specified directory.
#'
#' @param reg_net_file Character string, path to regulatory network file
#'   (.RData, .rda, .txt, .tsv, or .csv)
#' @param edge_file Character string, path to network edges file
#' @param pathway_gmt_file Character string, path to GMT file containing 
#'   pathway gene sets
#' @param res_dir Character string, directory path where results will be saved
#' @param ncores Integer, number of CPU cores to use for parallel 
#'   processing. Default is 1.
#' @param minSize Integer, minimum pathway size for filtering. Default is 5.
#' @param maxSize Integer, maximum pathway size for filtering. 
#'   Default is 150.
#' @param scale_data Logical, whether to scale data in PCA. Default is TRUE.
#' @param center_data Logical, whether to center data in PCA. Default is TRUE.
#' @param npcs Integer, number of principal components to compute. Default is 1.
#' @param nperm Integer, number of permutations for null model. 
#'   Default is 1000.
#' @param object_name Character string, specific object name to extract from
#'   RData files. If NULL, auto-detects. Default is NULL.
#' @param p_adjust_method Character string, multiple testing correction 
#'   method. Default is "fdr".
#' @param log_file Character string, path for log file. If NULL, creates
#'   timestamped log file in results directory. Default is NULL.
#'
#' @return List containing three data frames:
#'   \item{pathways_results}{PCA results for real pathways}
#'   \item{random_results}{PCA results for permuted data}
#'   \item{porcupine_results}{Final PORCUPINE results with p-values}
#'
#' @examples
#' \dontrun{
#' results <- runPORCUPINE(
#'   reg_net_file = "network.RData",
#'   edge_file = "edges.txt", 
#'   pathway_gmt_file = "pathways.gmt",
#'   res_dir = "results/",
#'   ncores = 4,
#'   nperm = 1000
#' )
#' }
#' @export
runPORCUPINE <- function(reg_net_file,
                        object_name = NULL,
                        edge_file,
                        pathway_gmt_file,
                        res_dir,
                        ncores = 1,
                        minSize = 5,
                        maxSize = 150,
                        scale_data = TRUE,
                        center_data = TRUE,
                        npcs = 1,
                        nperm = 1000,
                        p_adjust_method = "fdr",
                        p_adjust_threshold = 0.05,
                        log_file = NULL) {
    # Input validation
    if (!file.exists(reg_net_file)) {
        stop("Network file does not exist: ", reg_net_file)
    }
    if (!file.exists(edge_file)) {
        stop("Edge file does not exist: ", edge_file)
    }
    if (!file.exists(pathway_gmt_file)) {
        stop("Pathway GMT file does not exist: ", pathway_gmt_file)
    }
    # Setup logging
    if (is.null(log_file)) {
        log_file <- 
            file.path(res_dir, 
                paste0("porcupine_log_", 
                    format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    }
    if (!dir.exists(res_dir)) {
        dir.create(res_dir, recursive = TRUE)
    }
    # Create custom logging function for messages
    log_message <- function(...) {
        msg <- paste0(..., collapse = "")
        cat(msg, "\n")  # Display in console
        cat(msg, "\n", file = log_file, append = TRUE)  # Write to log
    }
    # Initialize log file
    cat("", file = log_file, append = FALSE)  # Create/clear log file
    on.exit({
        if (exists("start_time")) {
            cat("\nLog file closed at:", format(Sys.time()), "\n",
                file = log_file, append = TRUE)
        }
    })
    # Log session information
    log_message("=== PORCUPINE Analysis Started ===")
    log_message("Timestamp: ", format(Sys.time()))
    log_message("Log file: ", log_file)
    log_message("Parameters:")
    log_message("  Network file: ", reg_net_file)
    log_message("  Edge file: ", edge_file)
    log_message("  Pathway file: ", pathway_gmt_file)
    log_message("  Output directory: ", res_dir)
    log_message("  Cores: ", ncores)
    log_message("  Pathway size range: ", minSize, "-", maxSize)
    log_message("  Permutations: ", nperm)
    log_message("  P-value adjustment: ", p_adjust_method)
    log_message("  Scale data: ", scale_data)
    log_message("  Center data: ", center_data)
    if (!is.null(object_name)) {
        log_message("  Object name: ", object_name)
    }
    log_message("================================\n")
    start_time <- Sys.time()
    # Step 1: Load network data
    log_message("Step 1/6: Reading network file...")
    reg_net <- read_networks(reg_net_file, 
                object_name = object_name)
    log_message("Network dimensions: ", dim(reg_net)[1], " network edges ", 
            dim(reg_net)[2], " samples")
    # Step 2: Load edges
    log_message("Step 2/6: Reading edges file...")
    edges <- read_edges(edge_file)
    # Step 3: Load and filter pathways
    log_message("Step 3/6: Reading and filtering pathway file...")
    pathways <- load_gmt(pathway_gmt_file)
    pathways <- filter_pathways(pathways, edges)
    log_message("Number of pathways after edge filtering: ", 
            length(pathways))
    pathways_to_use <- filter_pathways_size(pathways, 
                                           minSize = minSize, 
                                           maxSize = maxSize)
    log_message("Number of pathways after size filtering (", minSize, "-", 
            maxSize, "): ", length(pathways_to_use))
    if (length(pathways_to_use) == 0) {
        stop("No pathways remain after filtering")
    }
    # Step 4: Run PCA analysis on pathways
    log_message("Step 4/6: Running PCA analysis on pathways...")
    tar_to_rows <- create_target_to_rows_mapping(edges)
    pca_res_pathways <- pca_pathway(pathways_to_use,
                                  reg_net,
                                  edges,
                                  tar_to_rows,
                                  npcs = npcs,
                                  ncores = ncores,
                                  scale_data = scale_data,
                                  center_data = center_data)
    # Save pathway results
    pathway_file <- file.path(res_dir, "pathways_results.txt")
    write.table(pca_res_pathways, pathway_file,
                col.names = TRUE, row.names = FALSE, 
                sep = "\t", quote = FALSE)
    log_message("Saved pathway results to: ", pathway_file)

    # Step 5: Run permutation analysis
    log_message("Step 5/6: Running permutation analysis (", nperm, 
            " permutations)...")
    log_message("This may take a long time...")
    pca_res_random <- pca_random(reg_net,
                               edges,
                               pca_res_pathways,
                               pathways_to_use,
                               tar_to_rows,
                               n_perm = nperm,
                               ncores = ncores,
                               scale_data = scale_data,
                               center_data = center_data)
    # Save random results
    random_file <- file.path(res_dir, "pathways_results_random.txt")
    write.table(pca_res_random, random_file,
                col.names = TRUE, row.names = FALSE, 
                sep = "\t", quote = FALSE)
    log_message("Saved permutation results to: ", random_file)
    # Step 6: Calculate final PORCUPINE results
    log_message("Step 6/6: Calculating final PORCUPINE results...")
    res_porcupine <- porcupine(pca_res_pathways, pca_res_random)
    res_porcupine$p.adjust <- p.adjust(res_porcupine$pval, 
                                     method = p_adjust_method)
    # Save final results
    porcupine_file <- file.path(res_dir, "porcupine_results.txt")
    write.table(res_porcupine, porcupine_file,
                col.names = TRUE, row.names = FALSE, 
                sep = "\t", quote = FALSE)
    log_message("Saved PORCUPINE results to: ", porcupine_file)
    # Summary
    end_time <- Sys.time()
    runtime <- end_time - start_time
    log_message("\nPORCUPINE analysis completed successfully!")
    log_message("Total runtime: ", round(runtime, 2), " ", 
            units(runtime))
    log_message("Results saved in directory: ", res_dir)
    
    significant_pathways <- sum(res_porcupine$p.adjust < p_adjust_threshold, 
                               na.rm = TRUE)
    log_message("Significant pathways (FDR < ", p_adjust_threshold, "): ", 
            significant_pathways, " out of ", nrow(res_porcupine))
    
    # Return results
    return(list(
        pathways_results = pca_res_pathways,
        random_results = pca_res_random,
        porcupine_results = res_porcupine
    ))
}
