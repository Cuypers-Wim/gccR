##############
#            #
# R package  #
#            #
##############


# install package

# library("devtools")
# library(roxygen2)

# setwd("~/PhD/Data/GitHub")
# install("gccR")

library(gccR)

##############
#            #
# analysis   #
#            #
##############


##### LT2 vs SL1344 ####

dataset1 <- file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_lt2_exprdata_20151029.txt")

dataset2 <- file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_sl1344_exprdata_20151029.txt")

# path to a list of orthologs (column1 = gene IDs dataset1, and column2 = orthologuous dataset2 IDs)

# orthologs <- file.path("~/PhD/Data/p_expression/COLOMBOS/reference_genomes/singleCopyOrthologs.tab")

singleCopyOrthologs <- read.delim("~/PhD/Data/p_expression/COLOMBOS/reference_genomes/lt2_sl1344_singleCopyOrthologs.tab", header=FALSE)

singleCopyOrthologs <- singleCopyOrthologs[ ,c(2,1)]

# load gene expression data

exprList <- readData(dataset1, dataset2, NAstring = "NaN", source = "COLOMBOS")

# compute correlation matrices

corM1 <- get_corM_fast(exprList$exprValues1, dropNArows = TRUE, threads = 8)

corM2 <- get_corM_fast(exprList$exprValues2, dropNArows = TRUE, threads = 8)

# extract submatrices

corM_ortho <- extract_core_submatrix(corM1, corM2, singleCopyOrthologs)

# order matrix 2 according to matrix 1 (given a list of orthologs if IDs differ)

csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, singleCopyOrthologs)

# compute the expression conservation scores by means of iterative comparison of co-expression

EC <- getEC_fast(corM_ortho$csM1, csM2_ordered, 0.001, threads = 1)

# compute perfect EC

perfectEC <- perfect_EC_fast(exprM = exprList$exprValues1, labelsM = corM_ortho$csM1,
                             conv = 0.001, maxIter = 200, threads = 1)

# diverged EC

divergedEC <- divergedEC_fast(expression = exprList$exprValues1, csM1 = corM_ortho$csM1,
                              csM2 = csM2_ordered, ortho,
                              conv = 0.001, maxIter = 200, threads = 1)


#### plotting with function #####

# check whether vectors are named


random_EC <- read.csv("~/PhD/Data/GitHub/CCR_fast/random_EC.csv")

random_EC$id <- rownames(corM_ortho$csM1)
randomEC_vec <- c()
randomEC_vec <- c(random_EC$r_EC)
names(randomEC_vec) <- c(random_EC$id)

distr_plot <- plot_distributions(EC$ECfinal, perfectEC, randomEC_vec)
distr_plot

