#' Construct Correlation Matrix (corM)
#'
#'This function makes a 'gene X gene' correlation matrix from a 
#''gene X condition' matrix
#'
#'
#' @param expression A matrix consisting of gene expression 
#' measurements across multiple conditions (columns), for different genes (rows)
#' @param dropNArows option to indicate whether to remove rows of which more 
#' than 50% of the observations are missing. Defaults to "TRUE"
#' @param threads Number of threads to be used for the correlation calculation 
#' using the 'cor' function in the 'WGCNA' package
#'
#' @return Square gene X gene correlation matrix.
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
#' CorM1 <- correlationMatrix(expression_matrix_1)
#'
#' @export
#' 
#' @import WGCNA

get_corM_fast <- function(expression, dropNArows = TRUE, threads = 8) {

  # library

  library("WGCNA")
  
  # BiocManager::install("GO.db")
  # BiocManager::install("impute")
  # BiocManager::install("preprocessCore")
  
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
