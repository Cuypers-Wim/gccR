#' Get EC
#'
#'   Calculate the expression conservation (EC) between (orthologuous) genes from two correlation matrices.
#' This function iteratively calculates the weighted correlation of a given row of m1
#' with a row (with the same rowindex) in m2 until convergence.
#'
#' @param m1 First 'gene X gene' correlation matrix
#' @param m2 Second 'gene X gene' correlation matrix
#' @param conv Convergence criterion indicating how much the final EC value per gene can
#' differ from the result of the previous iteration. Default = 0.001
#' @param maxIter The maximum number of iterations. Default = 200
#' @param threads Number of threads to be used for WGCNA::cor and parallel:: mclapply. 
#' This will heavily speed up the analysis on linux and MAC computers, but 
#' defaults to 1 to support usage on windows computers.
#' 
#' @return List (EClist) with the following elements:
#'\enumerate{
#'   \item EClist$ECfinal = a vector of final converged EC scores per gene
#'   \item ECfinal$ECweights = a vector of the last set of weights used for the iterative, weighted correlation
#' }
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
#' EC <- getEC(corM_ortho$csM1, csM2_ordered, conv = 0.001, maxIter = 200, threads = 1)
#' View(ECfinal$ECfinal)
#'
#' @export
#' 


getEC_fast <- function(m1, m2, conv = 0.001, maxIter = 200, threads = 1) {

  # check arguments

  if (missing(m1) | missing(m2)) stop("Please specify m1 and m2.")
  if (is.matrix(m1) == FALSE | is.matrix(m2) == FALSE) stop("Argument 1 and/or 2 is not a matrix.")
  if (nrow(m1) != ncol(m2) | nrow(m2) != ncol(m2)) stop("Matrix 1 and/or matrix 2 is not a square matrix.")
  if(identical(dim(m1), dim(m2)) == FALSE) stop("The provided matrices do not have the same dimensions.")


  # library 
  
  library(parallel )
  library(WGCNA)
  library(weights)
  
  # helper functions
  
  function_geneCor <- function(x, matrix1 = m1, matrix2 = m2) {
    
    WGCNA::cor(matrix1[x[1], ], matrix2[x[2], ], use = "pairwise.complete.obs", method = c("pearson"), quick = 0, nThreads = threads)
    
  }
  
  function_w_geneCor <- function(x, matrix1 = m1, matrix2 = m2, weights = corVec) {
    
    wtd.cors(matrix1[x[1], ], matrix2[x[2], ], weights)
    
  }
  
  
  function_w_geneCor_iterative <- function(x, matrix1 = m1, matrix2 = m2, weights = weightsVec) {
    
    wtd.cors(matrix1[x[1], ], matrix2[x[2], ], weights)
    
  }

  # Iterative EC calculation
  #
  # part 1 --------------
  # EC first iteration, no weighing

  # echo status

  print("Calculating Iterative EC")

  # preallocate 

  nRowCol <- nrow(m1)
  
  i_vec1 <- vector(length = nRowCol)
  i_vec2 <- vector(length = nRowCol)
  
  corVec <- vector(length = nRowCol)
  weightsVec <- vector(length = nRowCol)
  
  # change name vector elements

  names(corVec) <- rownames(m1)

  # the piece of code below is to paralellise the cor calculation
  # it was adapted from this blogpost https://davetang.org/muse/2012/01/31/creating-a-correlation-matrix-with-r/
  
  
  # first, we fill preallocated vectors of length == nrow(m1) with rowindices, 
  # merge them into a matrix, and make a list of each matrix row.
  
  i_vec1 <- seq(from = 1, to = nrow(m1))
  i_vec2 <- seq(from = 1, to = nrow(m1))
  
  i_matrix <- cbind(i_vec1, i_vec2)
  
  combos <- asplit(i_matrix, 1)
  
  
  # now we can easily use mclapply() to parallelise the cor calculation per row
  
  corList <- mclapply(combos, function_geneCor, mc.cores = threads)
  
  corVec <- unlist(corList)

  #---------- EC first iteration done

  # part 2 --------------

  # calculated weighted EC iteratively
  
  # preallocate

  EClist <- vector("list", length = 2)
  ECscores <- vector("list", length = maxIter)

  ECscores[[1]] <- corVec

  # do second, weighted  iteration

  # set values lower than zero to zero

  corVec[corVec < 0] <- 0

  # weighted correlation
  
  ECscores[[2]] <- unlist(mclapply(combos, function_w_geneCor, mc.cores = threads))


  # now we keep looping until convergence
  # (i.e for example difference between EC calculations per gene < 0.001 compared to the previous iteration)
  # don't exceed number of maximum iterations

  i <- 2
  
  ## todo: minimum 10 iteraties -> dan check en dan pas while loop

  while ((all(ECscores[[i]] - ECscores[[i-1]] < conv) == "FALSE") && i < maxIter) {

    i <- i + 1

    # replace weighted correlation scores lower than zero by zero

    weightsVec <- ECscores[[i-1]]
    weightsVec[weightsVec < 0] <- 0

    # use the changed ECscores[[i-1]] vector as weight for the weighted correlation
    
    ECscores[[i]] <- unlist(mclapply(combos, function_w_geneCor_iterative, mc.cores = threads))

  }

  print(paste("Total number of iterations untill convergence:  ", i))
  
  names(ECscores[[i]]) <- rownames(m1)

  # the 'final weigth' is the final EC score with negative cor values set to 0

  finalWeightsVec <- ECscores[[i-1]]
  finalWeightsVec[finalWeightsVec < 0] <- 0
  
  # store vector with converged EC and final weight values in named list

  EClist <- list("ECfinal" = ECscores[[i]], "ECweights" = finalWeightsVec)

  # return

  EClist
}
