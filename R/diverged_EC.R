#' Calculate Diverged Background Expression Conservation (EC) Distribution
#'
#' Computes the 'diverged' background EC distribution by calculating EC values between permutated gene expression data and a given correlation matrix. This function helps in assessing the conservation of gene expression patterns in a context where gene relationships are expected to have diverged.
#'
#' @param expression The 'gene X condition' expression matrix for the first dataset.
#' @param csM1 The 'gene X gene' correlation matrix for core genes from the first dataset.
#' @param csM2 The 'gene X gene' correlation (sub)matrix from the second dataset. Gene IDs must correspond to those in the `expression` matrix.
#' @param weights Optional: Weights for calculating weighted correlations. If not provided, the function iteratively refines weights in each iteration to calculate EC.
#' @param conv Convergence criterion indicating the permissible difference in EC values between iterations for convergence.
#' @param maxIter The maximum number of iterations for EC calculation.
#' @param mc_cores Number of CPU cores to be used for parallel processing.
#'
#' @return A named vector of EC values representing a diverged or random distribution of gene expression conservation.
#'
#' @details The function permutes rows of the expression matrix and calculates the EC between each permutated row and the rows of the csM1 matrix. If weights are not provided, it employs an iterative approach to refine the weights in each iteration. This process aims to model the divergence in gene expression patterns.
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#' randomEC1 <- divergedEC(exprList$exprValues1, core_submatrix2_ordered,
#'                          ortho, EC$ECweights, conv = 0.001, maxIter = 200)
#'
#' @export

divergedEC <- function(expression = NULL, csM1 = NULL, 
                       csM2 = NULL, conv = 0.0001, maxIter = 100, 
                       threads = 1, weights_fixed = NULL) {


  # Calculates background distribution for diverged expression

  # library

  library("weights")
  library("WGCNA")
  library("parallel")


  # Checks & preprocessing -----------------------------------------------------
  
  if (is.matrix(csM1) == FALSE | is.matrix(csM2) == FALSE) stop("Argument 2 and/or 3 is not a matrix.")
  if (nrow(csM1) != ncol(csM2) | nrow(csM2) != ncol(csM2)) stop("Matrix 2 and/or matrix 3 is not a square matrix.")
  if(identical(dim(csM1), dim(csM2)) == FALSE) stop("The provided correlation matrices do not have the same dimensions.")
  

  # Check if rownames are identical between expression, csM1, and csM2
  if (!identical(rownames(expression), rownames(csM1)) || !identical(rownames(expression), rownames(csM2))) {
    stop("The row names of expression, csM1, and csM2 must be identical.")
  }
  
  # Check if column names are identical between csM1 and csM2
  if (!identical(colnames(csM1), colnames(csM2))) {
    stop("The column names of csM1 and csM2 must be identical.")
  }
  
  # retain rows of which the ID is in csM1 + reorder

  expr <- expression[which(rownames(expression) %in% rownames(csM1)), ]

  expr <- expr[rownames(csM1), ]
  
  if (nrow(expr) != nrow(csM1)) stop("the number of rows in expression matrix 1
                                     was not the same as the number of rows in 
                                     correlation matrix 1")
  
  # Preallocate -------------------------------------------------------------

  nRowCol <- nrow(csM2)
  corVec <- vector(length = nRowCol)
  ECVec <- vector(length = nRowCol)
  resetExpr <- matrix(nrow = nrow(expr), ncol = ncol(expr))

  # make list with permutated rows of expression values
  
  expression_list <- asplit(expr, 1)
  
  function_permutate <- function(x) {
    
    x[sample(1:length(x))]
    
  }
  
  resampled_rows <- lapply(expression_list, function_permutate)
  
  print("rows resampled, moving on to correlation calculations")
  
  # check if sums are still equal after resample
  # sum(expression_list[[3569]], na.rm = TRUE) == sum(resampled_rows[[3569]], na.rm = TRUE)
  

  # Cor permutated row and every other exprM row ----------------------------

  ## calculate correlation between each permutated row and every row in the original expression matrix
  
  function_row_cor <- function(x, resampled_row, threads) {
    
    # x = every row stored in expression list
    
    WGCNA::cor(x, resampled_row, use = "pairwise.complete.obs", method = c("pearson"), quick = 0, nThreads = threads)
    
  }
  
  row_cor_list <- vector("list", length = nrow(expr))
  
  for (index in 1:nrow(expr)) {
    
    print(paste("calculating correlation of row ", index, "with the rest of the compendium"))
    
    row_cor_list[[index]] <- unlist(mclapply(expression_list, function_row_cor, resampled_rows[[index]], threads = threads))
    
  }
  
  print("correlations with permutated rows calculated, moving on to EC calculations")
  

  # EC calculations ---------------------------------------------------------

  # now we have a row_cor_list (list of vectors) containing the correlation of 
  # every permutated row with all other rows in the other compendium
  
  # next, we loop over the row correlations in row_cor_list
  # each iteration, another row AND column in the provided correlation matrix will be replaced by a vector from rowcorrelations
  # this new corM will be used to calculate the random EC with the corM of compendium 2
  
  # run the entire EC functions for the most accurate weighting
  # Print initial state of variables
  
  print(paste("Initial state of csM1:", dim((csM1))))
  print(paste("Initial state of csM2:", dim((csM2))))
	  
  for (i in 1:nrow(expr)) {
    
    corVec <- row_cor_list[[i]]
	  
    print(paste("State of corVec at iteration", i, ":", length(corVec)))
	  
    corM_switched <- csM1
    corM_switched[i, ] <- corVec
    corM_switched[ ,i] <- corVec
    
    if(is.null(weights_fixed)) {
      
    # Code to execute when weights_fixed is NULL or not provided
    print("running in iterative mode calculating new weights")
      
    EC <- getEC(corM_switched, csM2, conv, maxIter, threads = threads)

    print(paste("Output of getEC at iteration", i, ":", length(EC$ECfinal)))

    ECVec[i] <- EC$ECfinal[i]
    
    } else {
      
      print("running in 'reuse weights' mode")
      
      # This option works when weights from 'get_EC' are provided
      # This will drastically speed up the calculations, as the EC does not need to be calculated until convergence anymore
      # The effect on the actual random distribution is expected to be minimal
      # The EC function should work with two vectors, and one vector as weights
      
      EC <- getEC(corM_switched, csM2, conv, maxIter, threads = threads, weights_fixed = weights_fixed)
      
      print(paste("Output of getEC at iteration", i, ":", length(EC$ECfinal)))
      
      ECVec[i] <- EC$ECfinal[i]
      
    }
    
    print(paste("number of iterations done:  ", i))
  }

  names(ECVec) <- rownames(csM2)
  ECVec
  


}
