% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_core_submatrix.R
\name{extract_core_submatrix}
\alias{extract_core_submatrix}
\title{Extract Core Submatrices of Orthologous Genes}
\usage{
extract_core_submatrix(matrix1, matrix2, ortho)
}
\arguments{
\item{matrix1}{The first 'gene X gene' correlation matrix.}

\item{matrix2}{The second 'gene X gene' correlation matrix.}

\item{ortho}{Optional: A two-column data frame containing gene IDs in the first species/sample (first column) and their orthologous gene IDs in the second species/sample (second column). If not provided, the function assumes gene IDs in matrix1 are identical to those in matrix2.}
}
\value{
A list containing two core submatrices:
\itemize{
\item \code{csM1}: Core submatrix 1 derived from \code{matrix1}.
\item \code{csM2}: Core submatrix 2 derived from \code{matrix2}.
}
These submatrices contain only the genes (or orthologs) present in both original matrices (or in the ortho list, if provided).
}
\description{
Extracts submatrices from two 'gene X gene' correlation matrices. The function aligns these matrices by gene IDs or orthologous gene pairs, ensuring each submatrix contains genes present in the other. This is useful for comparative studies across species or conditions.
}
\details{
If the \code{ortho} parameter is not provided, the function matches genes based on identical IDs across the two input matrices. When \code{ortho} is used, it aligns the matrices based on the orthologous relationships defined in the \code{ortho} data frame. The function iteratively refines the matrices to ensure they contain only the corresponding genes or orthologs.
}
\examples{
corM_ortho <- extract_core_submatrix(corM1, corM2, singleCopyOrthologs)

}
\author{
Wim Cuypers, \email{wim.cuypers@uantwerpen.be}
}
