#' Construct a Gene-by-Gene Correlation Matrix
#'
#' Generates a correlation matrix representing gene-gene relationships based on their expression across different conditions. 
#' This function computes a 'gene X gene' correlation matrix from an input 'gene X condition' matrix.
#'
#' @param expression A matrix of gene expression data. Rows represent genes and columns represent different conditions. The matrix should contain gene expression measurements across multiple conditions.
#' @param dropNArows Logical flag indicating whether to exclude rows where more than 50% of the observations are missing. Defaults to `TRUE`. When set to `TRUE`, rows with more than 50% missing values ('NA') are removed before computing the correlation matrix.
#' @param threads The number of threads to utilize for the correlation computation. This is particularly useful for handling large datasets. The correlation computation uses the `cor` function from the 'WGCNA' package.
#'
#' @return A square matrix representing the gene-by-gene correlations. The matrix is computed based on pairwise complete observations using Pearson's correlation method. 
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#' CorM1 <- get_corM(expression_matrix_1)
#'
#' @export
#'
#' @import WGCNA

get_corM <- function(expression, dropNArows = TRUE, threads = 8) {

  # library

  library("WGCNA")
  
  # Computes the correlation matrix

  rows_to_drop <- which(rowMeans(is.na(expression)) > 0.5)
  
  if (dropNArows == TRUE && length(rows_to_drop) > 0) {

  # 1. Remove rows for which more than 50% of the cells is "NaN"

  expr <- expression[-rows_to_drop, ]

  # 2. Compute 'gene X gene' correlation matrix of the transposed 'gene X condition' matrix

  corM <- WGCNA::cor(t(expr), y = NULL,
              use = "pairwise.complete.obs", method = c("pearson"), quick = 0, nThreads = threads)

  }

  else {

  # 3. Option to compute correlation matrix without removing rows with more than 50% "NaN"

    corM <- WGCNA::cor(t(expression), y = NULL,
                use = "pairwise.complete.obs", method = c("pearson"), quick = 0, nThreads = threads)

  }

  # 4. return

  corM

}
