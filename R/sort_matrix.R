#' Align Rows and Columns of Correlation Matrices
#'
#' This function aligns the rows and columns of a second gene-gene correlation matrix (csM2) to match the order of another matrix (csM1). It ensures that the nth row (and column) in both matrices correspond to the same gene or its ortholog, maintaining the order of gene correlations.
#'
#' @param csM1 The reference 'gene X gene' correlation matrix.
#' @param csM2 The 'gene X gene' correlation matrix to be aligned with csM1.
#' @param ortho Optional: A two-column dataframe containing gene IDs in the first column and corresponding orthologous gene IDs in the second column. If provided, alignment is based on orthologous relationships. If omitted, the function assumes gene IDs in csM1 are identical to those in csM2.
#'
#' @return A matrix (csM2_ordered) that is reordered to match the row and column order of csM1. This reordering is done based on either identical gene IDs or ortholog gene IDs (if the 'ortho' argument is provided).
#'
#' @details The function reorders csM2 so that its genes (rows and columns) are in the same order as in csM1. If orthologous relationships are specified via the 'ortho' parameter, it uses this information to align genes correctly. Otherwise, it assumes identical gene IDs in both matrices and reorders accordingly.
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#' # Example when matrix 2 (csM2) has to be ordered according to matrix 1 (csM1) with orthology information:
#' csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, corM_ortho$orthologs)
#'
#' # Example when matrix 2 (csM2) has to be ordered according to matrix 1 (csM1) with identical gene IDs:
#' csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2)
#'
#' @export

sort_matrix <- function(csM1, csM2, ortho, rename = FALSE) {

  # checks

  if (missing(csM1) | missing(csM2)) stop("Please specify matrix1 and matrix2.")
  
  if (is.matrix(csM1) == FALSE) stop("Argument 1 is not a matrix.")
  if (is.matrix(csM2) == FALSE) stop("Argument 2 is not a matrix.")
  
  if (nrow(csM1) != ncol(csM2)) stop("Matrix 1 is not a square matrix.")
  if (nrow(csM2) != ncol(csM2)) stop("Matrix 2 is not a square matrix.")
  
  if(identical(dim(csM1), dim(csM2)) == FALSE) {
    stop("The provided matrices do not have the same dimensions.")
  }

  # preallocate vecs

  nRowCol <- nrow(csM1)
  rowVec <- vector(length = nRowCol)

  # preallocate matrix
  #  -- equal dimensions matrix 2

  csM2_ordered <- matrix(NA, nrow=nRowCol, ncol=nRowCol)

  #  -- make row and colnames equal to core_submatrix1 (csM1)

  rownames(csM2_ordered) <- rownames(csM1)
  colnames(csM2_ordered) <- colnames(csM1)

  if(missing(ortho)) {
    
    # row and col names are the same, so just reorder 

    rowVec <- rownames(csM1)
    csM2_ordered <- csM2[rowVec, rowVec]

  } else {
    if (!is.matrix(ortho) || ncol(ortho) != 2) {
      stop("Ortho is not a two-column matrix.")
    }

    for(i in 1:nRowCol) {
      
      # lookup orthologID and corresponding location in core submatrix 2
      
      rowGeneID <- rownames(csM1)[i]
      rowOrthologRowNumber <- which(ortho[ ,1] == rowGeneID)
      
      if (length(rowOrthologRowNumber) == 0) {
        warning(paste("No ortholog found in the provided ortholog list for", rowGeneID))
        next
      }
      
      
      rowOrthologID <- ortho[rowOrthologRowNumber, 2]
      orthoIndex <- which(rownames(csM2) == as.character(rowOrthologID))
      
      if (length(orthoIndex) == 0) {
        warning(paste("Ortholog ID not found in csM2 for", rowGeneID))
        next
      }
      
      rowVec[i] <- orthoIndex 
      
    }

  csM2_ordered <- csM2[rowVec, rowVec]

  }

  # return ordered matrix with original names

  csM2_ordered
  
  # Rename if the option is set to TRUE
  
  if (rename) {
    rownames(csM2_ordered) <- rownames(csM1)
    colnames(csM2_ordered) <- colnames(csM1)
  }
  
  # return renamed matrix 
  
  csM2_ordered

}
