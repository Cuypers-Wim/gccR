% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_core_submatrix.R
\name{extract_core_submatrix}
\alias{extract_core_submatrix}
\title{Extract a submatrix of core genes}
\usage{
extract_core_submatrix(matrix1, matrix2, ortho)
}
\arguments{
\item{matrix1}{First 'gene X gene' correlation matrix}

\item{matrix2}{Second 'gene X gene' correlation matrix}

\item{ortho}{OPTIONAL argument. Two-column data frame of gene IDs in species1/sample1,
and the corresponding orthologous gene IDs in species2/sample2.}
}
\value{
A list (corM_ortho) consisting of 2 elements:
\enumerate{
\item corM_ortho$csM1: core submatrix 1
\item corM_ortho$csM2: core submatrix 2
}
}
\description{
This function extracts submatrices from two 'gene X gene' correlation matrices.
Each resulting submatrix will consist of genes or orthologs that are present in the other submatrix.
A list of orthologs can be provided (via the optional argument 'ortho').
When using this option, bot matrices will consist exclusively of gene IDs that
have an orthologuous counterpart in the other matrix.
If no list of orthologs is provided, the function will expect gene IDs in
matrix 1 to be identical to gene IDs in matrix 2
}
\examples{

corM_ortho <- extract_core_submatrix(corM1, corM2, singleCopyOrthologs)

}
\author{
Wim Cuypers, \email{wim.cuypers@uantwerpen.be}
}