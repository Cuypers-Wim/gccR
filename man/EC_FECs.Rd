% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EC_FECs.R
\name{EC_FECs}
\alias{EC_FECs}
\title{EC per FEC}
\usage{
EC_FECs(m, clusters)
}
\arguments{
\item{exprM}{First 'gene X condition' expression matrix}

\item{labelsM}{Correlation (sub)matrix. Gene IDs must be the same as used in the exprM}
}
\value{
A  named vector of EC values representing perfect conservation, given technical variation
}
\description{
EC per FEC
}
\examples{

perfectEC1  <- perfect_EC(exprList$exprValues2, core_submatrix2_ordered, conv = 0.001, maxIter = 200)

}
\author{
Wim Cuypers, \email{wim.cuypers@uantwerpen.be}
}
