#' diverged EC distribution
#'
#'   Calculates the 'diverged' background expression conservation (EC) distribution.
#'
#'
#' @param expression First 'gene X condition' expression matrix
#' @param csM1 Correlation matrix of core genes - sample 1 
#' @param csM2 Correlation (sub)matrix. Gene IDs must be the same as used in the exprM
#' @param ortho Two-column dataframe of gene IDs in species1/sample1,
#' and the corresponding orthologous gene IDs in species2/sample2
#' @param weights Optional argument: weights that will be used for calculating the weighted correlation.
#' If no weights are provided, the EC value will be calculated iteratively, refining the weights each iteration
#' @param conv Convergence criterion indicating how much the final EC value per gene can
#' differ from the result of the previous iteration
#' @param maxIter The maximum number of iterations
#' @param mc_cores Number of cpus to be used
#'
#' @return A named vector of EC values representing a random distribution
#
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
#' randomEC1 <- divergedEC(exprList$exprValues1, core_submatrix2_ordered,
#'                          ortho, EC$ECweights, conv = 0.001, maxIter = 200)
#'
#' @export

divergedEC_fast <- function(expression = exprList$exprValues1, csM1 = corM_ortho$csM1, 
                       csM2 = csM2_ordered, ortho, 
                       conv = 0.0001, maxIter = 100, threads = 1, weights = NA) {

	# toDO: change argument name 'expression' to something else
  # todo: add 'weights' option again to reuse weights

  # Calculates background distribution for diverged expression

  # library

  library("weights")
  library("WGCNA")
  library("parallel")

  # Remove rows for which more than 50% of the cells is "NaN" -> weglaten, ineens selectie en ordenen

  # expr <- expression[-which(rowMeans(is.na(expression)) > 0.5), ]

  # retain rows of which the ID is in csM1 + reorder

  expr <- expression[which(rownames(expression) %in% rownames(csM1)), ]

  expr <- expr[rownames(csM1), ]
  
  if (nrow(expr) != nrow(csM1)) stop("the number of rows in expression matrix 1
                                     was not the same as the number of rows in 
                                     correlation matrix 1")
  
  # preallocate

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
  
  # calculate correlation between each permutated row and every row in the original expression matrix
  
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
  
  # now we have a row_cor_list (list of vectors) containing the correlation of 
  # every permutated row with all other rows in the other compendium
  
  # next, we loop over the rowcorrelations in row_cor_list
  # each iteration, another row AND column will be replaced by a vector in rowcorrelations
  # this new corM will be used to calculate the random EC with the corM of compendium 2
  
  if(missing(weights)) {
    
  # run the entire EC functions for the most accurate weighting
  
  for (i in 1:nrow(expr)) {
    
    corVec <- row_cor_list[[i]]
    
    corM_switched <- csM1
    corM_switched[i, ] <- corVec
    corM_switched[ ,i] <- corVec
    
    EC <- getEC_fast(corM_switched, csM2, conv, maxIter, threads = threads)
    
    ECVec[i] <- EC$ECfinal[i]
    
    print(paste("number of iterations done:  ", i))
  }

  names(ECVec) <- rownames(csM2)
  ECVec
  
  } else {
    
    # UNDER CONSTRUCTION
    # to do: EC functie nodig die werkt met twee vectors en 1 waarde als gewicht
    
    csM2_list <- asplit(csM2, 1)
    
    EC_list <- vector("list", length = nrow(expr))
    
    for (i in 1:length(row_cor_list)) {
      
      EC_list[[i]] <- unlist(mclapply(csM2_list, function_row_cor, resampled_rows[[i]], threads = threads))
     
      print(paste("number of iterations done:  ", i))
    }
    
    
    names(ECVec) <- rownames(csM1)
    ECVec
    
    
  }


}
