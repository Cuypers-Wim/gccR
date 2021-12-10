#' Sort matrices
#'
#' This function sorts rows and columns of a symmetrical correlation matrix 
#'   (second argument; csM2) according to the other (first argument; csM1).
#' As a result, row n in Matrix 1 and row n in Matrix 2 correspond to the
#' same gene or ortholog, with the row values describing the correlation
#' of this gene with the other genes in the compendium in the exact same order.
#'
#' @param matrix1 First 'gene X gene' correlation matrix
#' @param matrix2 Second 'gene X gene' correlation matrix
#' @param ortho OPTIONAL argument. Two-column dataframe of gene IDs in species1/sample1,
#' and the corresponding orthologous gene IDs in species2/sample2. 
#'
#' @return A matrix (csM2_ordered) that is row- and column ordered according to matrix1
#' (either by identical gene IDs or ortholog gene IDs if the optional argument 'ortho' is provided)
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
#' # example when matrix 2 (csM2) has to be ordered according to matrix 1 (csM1)
#' # the orthology info (corM_ortho$orthologs) is taken into account
#'
#' csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, corM_ortho$orthologs)
#'
#' # example when matrix 2 (csM2) has to be ordered according to matrix 1 (csM1)
#' # bot matrices contain identical gene IDs
#'
#' csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2)
#'
#' @export

sort_matrix <- function(csM1, csM2, ortho) {

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

  for(i in 1:nRowCol) {

    # lookup orthologID and corresponding location in core submatrix 2

    rowGeneID <- (rownames(csM1))[i]
    rowOrthologRowNumber <- which(ortho[ ,1] == rowGeneID)
    rowOrthologID <- ortho[rowOrthologRowNumber, 2]
    rowVec[i] <- which(rownames(csM2) == toString(rowOrthologID))

  }

  csM2_ordered <- csM2[rowVec, rowVec]

  }

  # return ordered matrix

  csM2_ordered

}
