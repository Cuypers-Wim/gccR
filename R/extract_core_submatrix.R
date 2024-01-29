#' Extract Core Submatrices of Orthologous Genes
#'
#' Extracts submatrices from two 'gene X gene' correlation matrices. The function aligns these matrices by gene IDs or orthologous gene pairs, ensuring each submatrix contains genes present in the other. This is useful for comparative studies across species or conditions.
#' @param matrix1 The first 'gene X gene' correlation matrix.
#' @param matrix2 The second 'gene X gene' correlation matrix.
#' @param ortho Optional: A two-column data frame containing gene IDs in the first species/sample (first column) and their orthologous gene IDs in the second species/sample (second column). If not provided, the function assumes gene IDs in matrix1 are identical to those in matrix2.
#'
#' @return A list containing two core submatrices:
#'   \itemize{
#'     \item `csM1`: Core submatrix 1 derived from `matrix1`.
#'     \item `csM2`: Core submatrix 2 derived from `matrix2`.
#'   }
#' These submatrices contain only the genes (or orthologs) present in both original matrices (or in the ortho list, if provided).
#'
#' @details If the `ortho` parameter is not provided, the function matches genes based on identical IDs across the two input matrices. When `ortho` is used, it aligns the matrices based on the orthologous relationships defined in the `ortho` data frame. The function iteratively refines the matrices to ensure they contain only the corresponding genes or orthologs.
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#' corM_ortho <- extract_core_submatrix(corM1, corM2, singleCopyOrthologs)
#'
#' @export

extract_core_submatrix <- function(matrix1, matrix2, ortho) {

  # checks

  if (missing(matrix1) | missing(matrix2)) {
  stop("Please specify matrix1 and matrix2.")
  }
  
  if (is.matrix(matrix1) == FALSE) stop("Argument 1 is not a matrix.")
  if (is.matrix(matrix2) == FALSE) stop("Argument 2 is not a matrix.")
  
  if (nrow(matrix1) != ncol(matrix1)) stop("Matrix 1 is not a square matrix.")
  if (nrow(matrix2) != ncol(matrix2)) stop("Matrix 2 is not a square matrix.")
  
  # preallocate

  list <- vector(mode = "list", length = 2)

  # option in case gene IDs are identical in both matrices

  if(missing(ortho)) {

    # print lenghts before extraction

    print(paste("Number of rows matrix 1 before selection:  ", nrow(matrix1)))
    print(paste("Number of rows matrix 2 before selection:  ", nrow(matrix2)))

    # extract

    while (nrow(matrix1) != nrow(matrix2)) {

      print("Matrix1 and matrix2 still differ in length!")

      # retain only matrix rows and columns if their ID is in the other matrix

      matrix1 <- matrix1[which(rownames(matrix1) %in% rownames(matrix2)),
                         which(rownames(matrix1) %in% rownames(matrix2))]
      matrix2 <- matrix2[which(rownames(matrix2) %in% rownames(matrix1)),
                         which(rownames(matrix2) %in% rownames(matrix1))]

    }

  } else {

    # print lenghts before extraction

    print(paste("Number of rows matrix 1 before selection:  ", nrow(matrix1)))
    print(paste("Number of rows matrix 2 before selection:  ", nrow(matrix2)))
    print(paste("Number of rows singleCopyOrthologs before selection:  ", nrow(ortho)))

    # extract

    while (nrow(matrix1) != nrow(matrix2) && nrow(matrix2) != nrow(ortho)) {

      print("Matrix1, matrix2 and the list of orthologs still differ in length!")

      # extract

      # retain only ortholog IDs that exist as rowname in Matrix 1 and matrix 2

      ortho <- ortho[which(ortho[ ,1] %in% rownames(matrix1)), ]
      ortho <- ortho[which(ortho[ ,2] %in% rownames(matrix2)), ]

      # retain only matrix rows and columns if their ID is in the ortholog list

      matrix1 <- matrix1[which(rownames(matrix1) %in% ortho[ ,1]),
                         which(rownames(matrix1) %in% ortho[ ,1])]
      matrix2 <- matrix2[which(rownames(matrix2) %in% ortho[ ,2]),
                         which(rownames(matrix2) %in% ortho[ ,2])]

    }

  }

  # store the newly extracted matrices in a list if the number of rows is identical

  corM_ortho <- list("csM1" = matrix1, "csM2" = matrix2)

  # print the final number of core orthologs

  final = nrow(matrix1)

  print(paste("The number rows after selection:  ", final))

  # return list

  corM_ortho

}


