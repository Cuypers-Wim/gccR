#' Remove Rows and Corresponding Columns with Zero Variance
#'
#' This internal function identifies and removes rows (and their corresponding columns)
#' with zero variance from two matrices, ensuring the matrices remain symmetrical. It is
#' particularly useful for preprocessing correlation matrices before further analysis.
#'
#' @param mat1 A numeric matrix for which rows and columns with zero variance are to be removed.
#' @param mat2 A numeric matrix for which rows and columns with zero variance are to be removed.
#' @return A named list with two elements, `matrixA` and `matrixB`, each being the input matrices
#'         with rows and columns of zero variance removed.
#' @examples
#' mat1 <- matrix(rnorm(100), 10)
#' mat2 <- matrix(rnorm(100), 10)
#' result <- internal_remove_zero_variance_rows(mat1, mat2)
#' @export
#' @keywords internal


internal_remove_zero_variance_rows <- function(mat1, mat2) {
  
  # Use apply to calculate variance for each row
  variances1 <- apply(mat1, 1, var)
  variances2 <- apply(mat2, 1, var)
  
  # Identify rows with zero variance
  rows_to_remove1 <- which(variances1 == 0)
  rows_to_remove2 <- which(variances2 == 0)
  
  rows_to_remove <- unique(c(rows_to_remove1, rows_to_remove2))
  
  # Print warnings for rows with zero variance
  if (length(rows_to_remove) > 0) {
    warnings <- paste("Warning: Row", rows_to_remove, "has zero variance and will be removed.")
    warning(paste(warnings, collapse = "\n"))
    
    # Remove the rows with zero variance
    # since these are symmetrical correlation matrices, also remove the corresponding columns
    
    mat_clean1 <- mat1[-rows_to_remove,-rows_to_remove, drop = FALSE]
    mat_clean2 <- mat2[-rows_to_remove,-rows_to_remove, drop = FALSE]
   
  } else {
    message("No zero variance rows were found.")
    mat_clean1 <- mat1
    mat_clean2 <- mat2
  }
  
  return(list(matrixA = mat_clean1, matrixB = mat_clean2))
  
}

#' Calculate EC for Subsets of an Expression Matrix Grouped per Sample
#'
#' This internal function calculates the EC for subsets of an
#' expression matrix. It groups the expression matrix per sample based on the specified
#' experiment combinations, computes correlation matrices, and then calculates the EC.
#' It is designed to work with subsets, particularly useful for analyses that require
#' partitioning of data into meaningful groups before EC calculation.
#'
#' @param exp_combo Vector of experiment identifiers to define subsets.
#' @param exprM Numeric matrix of expression data where columns represent experiments.
#' @param labelsM Matrix of labels corresponding to `exprM`.
#' @param experiments_all Vector of all possible experiment identifiers. If not NULL, used
#'        to filter `exprM` columns for subset EC calculation.
#' @param threads Integer, the number of threads to use for parallel computation (default is 1).
#' @param maxIter Maximum number of iterations for the EC computation (default is 200).
#' @param conv Convergence threshold for the EC computation (default is 0.001).
#' @return The final EC value
#' @details
#' The function first identifies columns in the expression matrix that match the specified
#' experiment combinations. It then splits the expression matrix into two halves based on these
#' combinations, computes correlation matrices for each half, and removes rows and columns
#' with zero variance. The EC is calculated from these processed correlation matrices.
#' @examples
#' exprM <- matrix(rnorm(1000), 100, 10)
#' labelsM <- matrix(sample(0:1, 1000, replace = TRUE), 100, 10)
#' experiments_all <- 1:10
#' exp_combo <- c(1, 5, 9)
#' ECresult <- internal_ec_subset(exp_combo, exprM, labelsM, experiments_all)
#' @export
#' @keywords internal

internal_ec_subset <- function(exp_combo, exprM, labelsM, 
                               experiments_all = NULL, threads = 1, 
                               maxIter = 200, conv = 0.001) {
  
  # get EC for subsets of an expression matrix grouped per sample
  
  
  columns_kr <- which(experiments_all %in% exp_combo)
  
  half_exprM1 <- exprM[ , columns_kr, drop = FALSE]
  half_exprM2 <- exprM[ , -columns_kr, drop = FALSE]
  
  # if the variance of a vector is zero, the correlation cannot be calculated
  # rows with zero variance can however be removed
  
  corM1 <- get_corM(half_exprM1, dropNArows = TRUE)
  corM2 <- get_corM(half_exprM2, dropNArows = TRUE)
  
  coma <- internal_remove_zero_variance_rows(corM1, corM2)
  
  subCorM1 = coma$matrixA[which(rownames(coma$matrixA) %in% rownames(labelsM)),
                   which(rownames(coma$matrixA) %in% rownames(labelsM))]
  subCorM2 = coma$matrixB[which(rownames(coma$matrixB) %in% rownames(labelsM)),
                   which(rownames(coma$matrixB) %in% rownames(labelsM))]
  
  subCorM <- extract_core_submatrix(subCorM1, subCorM2)
  csM2_ordered <- sort_matrix(subCorM$csM1, subCorM$csM2)
 
  is.matrix(subCorM$csM1)
  is.matrix(csM2_ordered)
  
  EC <- getEC(subCorM$csM1, csM2_ordered, 
              maxIter = maxIter, threads = threads, conv = conv)
  ECresult <- EC$ECfinal
  ECresult
  
} # end of function_ec_subset
