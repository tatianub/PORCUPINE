#' Determine statistical significance of PCA score of a pathway
#'
#' Compares the observed PCA score for a pathway to a set of PCA scores of 
#' random gene sets of the same size as pathway. Uses a one-tailed t-test to
#' determine if the pathway PC1 score is significantly lower than random sets.
#' Also calculates Cohen's D effect size.
#' 
#' @param res_pca_pathway Data frame with pathway PCA results. Must contain
#'   columns: 'pathway' (character), 'pathway_size' (numeric), and 'pc1' 
#'   (numeric).
#' @param res_pca_rndm_set Data frame with random gene set PCA results. Must
#'   contain columns: 'pathway_size' (numeric) and 'pc1' (numeric).
#' @return Data frame with statistical results containing columns:
#'   \itemize{
#'     \item pathway: Pathway name (character)
#'     \item pathway_size: Number of genes in pathway (numeric)
#'     \item pc1: Observed PC1 score (numeric)
#'     \item pval: P-value from one-tailed t-test (numeric)
#'     \item es: Cohen's D effect size (numeric)
#'   }
#' @examples
#' \dontrun{
#' # Example pathway result
#' pathway_data <- data.frame(
#'   pathway = "KEGG_GLYCOLYSIS",
#'   pathway_size = 50,
#'   pc1 = 25
#' )
#' 
#' # Example random results
#' random_data <- data.frame(
#'   pathway_size = rep(50, 1000),
#'   pc1 = rnorm(1000, mean = 25, sd = 0.5)
#' )
#' 
#' stats <- calculate_statistics(pathway_data, random_data)
#' }
#' @export

calculate_statistics <- function(
        res_pca_pathway,
        res_pca_rndm_set) {
    if (nrow(res_pca_pathway) == 0) {
        stop("Pathways results is empty.")
    }
    if (nrow(res_pca_rndm_set) == 0) {
        stop("Permutation random results is empty.")
    }
    expected_columns_pca_pathway <- c("pathway", "pathway_size", "pc1")
    if (!all(expected_columns_pca_pathway %in% 
             colnames(res_pca_pathway))) {
        stop("Missing required columns.")
    }
    expected_columns_pca_rndm <- c("pathway_size", "pc1")
    if (!all(expected_columns_pca_rndm %in% 
             colnames(res_pca_rndm_set))) {
        stop("Missing required columns.")
    }
    pc1_rndm <- res_pca_rndm_set$pc1
    pathway <- res_pca_pathway[["pathway"]]
    path_size <- res_pca_pathway[["pathway_size"]]
    pc1_pathway <- as.numeric(res_pca_pathway[["pc1"]])
    pvalue <- t.test(pc1_rndm, mu = pc1_pathway, 
                     alternative = "less")$p.value
    effect_size <- lsr::cohensD(pc1_rndm, mu = pc1_pathway)
    stat_res <- data.frame(
        "pathway" = pathway,
        "pathway_size" = path_size,
        "pc1" = pc1_pathway,
        "pval" = pvalue,
        "es" = effect_size)
    return(stat_res)
}

#' Determine statistical significance of PCA scores of pathways
#'
#' This function compares PCA results for multiple pathways versus sets of
#' randomly sampled gene sets of matching sizes. For each pathway size group,
#' it calculates p-values and effect sizes by comparing observed pathway PC1
#' scores to the distribution of random PC1 scores.
#' 
#' @param res_pca_pathways Data frame with pathway PCA results. Must contain
#'   columns: 'pathway' (character), 'pathway_size' (numeric), and 'pc1'
#'   (numeric).
#' @param res_pca_rndm Data frame with random gene set PCA results. Must
#'   contain columns: 'pathway_size' (numeric) and 'pc1' (numeric).
#' @return Data frame with statistical results for all pathways, containing
#'   columns:
#'   \itemize{
#'     \item pathway: Pathway name (character)
#'     \item pathway_size: Number of genes in pathway (numeric)
#'     \item pc1: Observed PC1 score (numeric)
#'     \item pval: P-value from one-tailed t-test (numeric)
#'     \item es: Cohen's D effect size (numeric)
#'   }
#' @examples
#' \dontrun{
#' # Example pathway results
#' pathways <- data.frame(
#'   pathway = c("KEGG_GLYCOLYSIS", "KEGG_TCA_CYCLE"),
#'   pathway_size = c(50, 30),
#'   pc1 = c(25, 18)
#' )
#' 
#' # Example random results
#' random_results <- data.frame(
#'   pathway_size = c(rep(50, 500), rep(30, 500)),
#'   pc1 = c(rnorm(500, 25, 0.5), rnorm(500, 18, 0.4))
#' )
#' 
#' results <- porcupine(pathways, random_results)
#' }
#' @export

porcupine <- function(res_pca_pathways, res_pca_rndm) {
  pathways_size <- unique(res_pca_pathways$pathway_size)
  res_all <- list()
  for (k in 1:length(pathways_size)) {
      path_size <- pathways_size[k]
      print(path_size)
      # select pca results for the pathways of the size path_size
      res_pathway_set <- res_pca_pathways %>%
                dplyr::filter(pathway_size == path_size)
      # select pca results for random sets of the size of a pathway
      res_rndm_set <- res_pca_rndm %>%
                dplyr::filter(pathway_size == path_size)
      # calculate pvalue and effect size for each pathway versus random sets
      res <- list()
      for (i in 1:nrow(res_pathway_set)) {
          pathway_row <- res_pathway_set[i, , drop = FALSE]
          res[[i]] <- calculate_statistics(pathway_row, res_rndm_set)
      }
      names(res) <- NULL
      res_all <- c(res, res_all)
  }
  res_all <- do.call("rbind", res_all)
  return(res_all)
}
