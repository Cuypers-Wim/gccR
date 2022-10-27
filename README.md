# gccR
gene co-expression conservation calculations in R


# Functionality
- Expression conservation (EC) calculations between two gene expression compendia
- Identification and analysis of functional expression classes (FECs)

# Installation

``` r
# Installation is currently only possible via the development version from GitHub:
# install.packages("devtools", auth_token="GHkeyHere")
devtools::install_github("Cuypers-Wim/CCR")
```

# Example usage

```R
library(CCR)

# path to two gene expression datasets (e.g downloaded from COLOMBOS.net) that you wish to compare

dataset1 <- file.path("/I/have/seen/the/bad/moon/rising/expressionCompendium1.txt")
dataset2 <- file.path("/I/wanna/know/have/you/ever/seen/the/rain/expressionCompendium2.txt")

# path to a list of orthologs (column1 = gene IDs dataset1, and column2 = orthologuous dataset2 IDs)

orthologs <- file.path("/this/path/is/truly/amazing/and/points/to/ortologs.txt")

# load gene expression data

exprList <- readData(dataset1, dataset2, NAstring = "NaN", source = "COLOMBOS")

# compute correlation matrices

corM1 <- get_corM(exprList$exprValues1, dropNArows = TRUE)
corM2 <- get_corM(exprList$exprValues2, dropNArows = TRUE)

# extract submatrices

corM_ortho <- extract_core_submatrix(corM1, corM2, singleCopyOrthologs)

# order matrix 2 according to matrix 1 (given a list of orthologs if IDs differ)

csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, singleCopyOrthologs)

# compute the expression conservation scores by means of iterative comparison of co-expression

EC <- getEC(corM_ortho$csM1, csM2_ordered, 0.001)

# estimate background distributions

perfectEC1  <- perfect_EC(exprList$exprValues2, csM2_ordered, conv = 0.001, maxIter = 200)
randomEC1 <- divergedEC(exprList$exprValues2, csM2_ordered, singleCopyOrthologs, EC$ECweights)

# These datasets can now be visualised and furter interrogated

```

# Credits
  Adapted from:
 - Sonego P. et al. (2015) 
 - Meysman P. et al. (2013)
