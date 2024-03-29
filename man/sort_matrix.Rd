% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sort_matrix.R
\name{sort_matrix}
\alias{sort_matrix}
\title{Align Rows and Columns of Correlation Matrices}
\usage{
sort_matrix(csM1, csM2, ortho, rename = FALSE)
}
\arguments{
\item{csM1}{The reference 'gene X gene' correlation matrix.}

\item{csM2}{The 'gene X gene' correlation matrix to be aligned with csM1.}

\item{ortho}{Optional: A two-column dataframe containing gene IDs in the first column and corresponding orthologous gene IDs in the second column. If provided, alignment is based on orthologous relationships. If omitted, the function assumes gene IDs in csM1 are identical to those in csM2.}
}
\value{
A matrix (csM2_ordered) that is reordered to match the row and column order of csM1. This reordering is done based on either identical gene IDs or ortholog gene IDs (if the 'ortho' argument is provided).
}
\description{
This function aligns the rows and columns of a second gene-gene correlation matrix (csM2) to match the order of another matrix (csM1). It ensures that the nth row (and column) in both matrices correspond to the same gene or its ortholog, maintaining the order of gene correlations.
}
\details{
The function reorders csM2 so that its genes (rows and columns) are in the same order as in csM1. If orthologous relationships are specified via the 'ortho' parameter, it uses this information to align genes correctly. Otherwise, it assumes identical gene IDs in both matrices and reorders accordingly.
}
\examples{
# Example when matrix 2 (csM2) has to be ordered according to matrix 1 (csM1) with orthology information:
csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, corM_ortho$orthologs)

# Example when matrix 2 (csM2) has to be ordered according to matrix 1 (csM1) with identical gene IDs:
csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, as.matrix(singleCopyOrthologs), rename = TRUE)

}
\author{
Wim Cuypers, \email{wim.cuypers@uantwerpen.be}
}
