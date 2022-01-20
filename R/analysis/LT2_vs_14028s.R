##############
#            #
# R package  #
#            #
##############


# library

library("devtools")
library(roxygen2)

# setwd("~/PhD/Data/GitHub")
# install("gccR")

library(gccR)

# process documentation

setwd("~/PhD/Data/GitHub/gccR")
document()

# re-install package

setwd("~/PhD/Data/GitHub")
install("gccR")

# or reload package after changes

devtools::load_all()


##############
#            #
# analysis   #
#            #
##############



##### LT2 VS 14028s #####

dataset1 <- file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_lt2_exprdata_20151029.txt")

dataset2 <- file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_14028s_exprdata_20151029.txt")

# path to a list of orthologs (column1 = gene IDs dataset1, and column2 = orthologuous dataset2 IDs)

# orthologs <- file.path("~/PhD/Data/p_expression/COLOMBOS/reference_genomes/singleCopyOrthologs.tab")

singleCopyOrthologs <- read.delim("~/PhD/Data/p_expression/COLOMBOS/reference_genomes/singleCopyOrthologs_reordered.tab", header=FALSE)

singleCopyOrthologs <- singleCopyOrthologs[ ,c(2,3)]

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

EC <- getEC_fast(corM_ortho$csM1, csM2_ordered, 0.0000001, threads = 1)

# compute perfect EC

# perfectEC <- perfect_EC_fast(exprM = exprList$exprValues1, labelsM = corM_ortho$csM1,
#                              conv = 0.001, maxIter = 200, threads = 1)
#

perfectEC <- perfect_EC_fast(exprM = exprList$exprValues2, labelsM = corM_ortho$csM2,
                             conv = 0.001, maxIter = 200, threads = 1, ortho = singleCopyOrthologs)


splits <- perfect_EC_fast(exprM = exprList$exprValues2, labelsM = corM_ortho$csM2, conv = 0.001,
                          maxIter = 200, threads = 1, ortho = singleCopyOrthologs,
                          multiple_splits = TRUE, experiment_info = exprList$headerInfo2)

splits
max_length <- max(lengths(splits))
add_na_function <- function(x, max_length) {
  length(x) = max_length
  x
}
adjusted <- lapply(splits, add_na_function, max_length)
df <- as.data.frame(do.call(cbind, adjusted))
molten_df <- melt(df)
df
molten_df
p <- ggplot(molten_df, aes(x=value, colour=variable )) +
  geom_density() +
  xlim(-0.8, 1.2) +
  theme_classic()
p
min(unlist(splits))
View(experiment_ids)




# diverged EC

divergedEC <- divergedEC_fast(expression = exprList$exprValues1, csM1 = corM_ortho$csM1,
                              csM2 = csM2_ordered, ortho,
                              conv = 0.001, maxIter = 200, threads = 1)


# perfect EC but with 10 subsets grouped per experiment

perfect_ECs <- perfect_EC_fast(exprM = exprList$exprValues2, labelsM = corM_ortho$csM2, conv = 0.001,
                               maxIter = 200, threads = 1, ortho = singleCopyOrthologs,
                               multiple_splits = TRUE, experiment_info = exprList$headerInfo2)

max_length <- max(lengths(splits))

min(unlist(splits))

add_na_function <- function(x, max_length) {

  length(x) = max_length

  x

}


adjusted <- lapply(splits, add_na_function, max_length)


df <- as.data.frame(do.call(cbind, adjusted))
View(df)


molten_df <- melt(df)


p <- ggplot(molten_df, aes(x=value, colour=variable )) +
  geom_density() +
  xlim(-0.8, 1.2) +
  theme_classic()


#### plotting LT2 vs 14028s distributions with function #####

# check whether vectors are named

random_EC <- read.csv("~/PhD/Data/GitHub/CCR_fast/lt2_14028s_random_EC.csv")
random_EC <- read.csv("~/PhD/Data/GitHub/CCR_fast/lt2_14028s_random_EC_NEW.csv")


random_EC$id <- rownames(corM_ortho$csM1)
randomEC_vec <- c()
randomEC_vec <- c(random_EC$r_EC)
names(randomEC_vec) <- c(random_EC$id)

distr_plot <- plot_distributions(EC$ECfinal, perfectEC, randomEC_vec)
distr_plot


distributions_list <- list(EC = EC$ECfinal, conserved = perfectEC, diverged = randomEC_vec)
names(distributions_list)


colnames_df <- c()
min_val <- c()
max_val <- c()

for (i in 1:length(distributions_list)) {

  colnames_df[i] <- names(distributions_list)[i]

  min_val[i] <- min(distributions_list[[i]])
  max_val[i] <- max(distributions_list[[i]])


}

cbind(names = colnames_df, min = min_val, max = max_val)


min(EC$ECfinal, perfectEC, randomEC_vec)


