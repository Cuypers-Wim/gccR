#' Extract a submatrix of core genes
#'
#' This function extracts submatrices from two 'gene X gene' correlation matrices.
#' Each resulting submatrix will consist of genes or orthologs that are present in the other submatrix.
#' A list of orthologs can be provided (via the optional argument 'ortho').
#'  When using this option, bot matrices will consist exclusively of gene IDs that
#'  have an orthologuous counterpart in the other matrix.
#' If no list of orthologs is provided, the function will expect gene IDs in
#'  matrix 1 to be identical to gene IDs in matrix 2
#'
#'
#' @param matrix1 First 'gene X gene' correlation matrix
#' @param matrix2 Second 'gene X gene' correlation matrix
#' @param ortho OPTIONAL argument. Two-column data frame of gene IDs in species1/sample1,
#' and the corresponding orthologous gene IDs in species2/sample2. 
#'
#' @return A list (corM_ortho) consisting of 2 elements:
#' \enumerate{
#'   \item corM_ortho$csM1: core submatrix 1
#'   \item corM_ortho$csM2: core submatrix 2
#' }
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
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


