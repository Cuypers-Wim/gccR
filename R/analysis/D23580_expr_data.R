
# Script info -------------------------------------------------------------

# parse and analyse expression data of D23580


# requirements per sample:

# -- accession per folder
# -- tx2geneFinal (links transcript to gene names)
# -- coldata = experimental design

# Library -----------------------------------------------------------------

library("tximportData")
library("tximport")
library("tximeta")
library("DESeq2")

# Import and parse --------------------------------------------------------

## tximport (to import abundance info) ------------------------------------

# path to file with id mapping

tx2genePath <- "D:/Documenten/PhD/Data/p_expression/typhimurium_D23580/R/tx2gene.txt"

tx2gene <- read.delim("~/PhD/Data/p_expression/typhimurium_D23580/R/tx2gene.txt", header=FALSE)

# paths with kallisto data files

paddo <- "D:/Documenten/PhD/Data/p_expression/typhimurium_D23580/kallisto_tphm_all"

# from tximportData Kallisto tutorial (https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto_with_TSV_files)

files <- list.files(paddo, full.names = TRUE)
names(files) <- sub('_abundance.tsv', '', basename(files))

txi_kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = FALSE)


## read dataset in DESeq and normalise ------------------------------------

# we consider every sample as if measured under the same condition

samplesVec <- colnames(txi_kallisto$counts)

conditions <- rep(c("A"), times = length(samplesVec))

samples <- cbind(samplesVec, conditions = "A")

rownames(samples) <- samples[ ,1]

samples <- samples[ , -1, drop = FALSE]


ddsTxi <- DESeqDataSetFromTximport(txi_kallisto,
                                   colData = samples,
                                   design = ~ 1)

# make the same dds object, but add normalised counts to a slot

dds_sizeF <- estimateSizeFactors(ddsTxi)

# get counts that are now normalised for library size

normalised_counts <- counts(dds_sizeF, normalized=TRUE)


## difference counts vs normalised ----------------------------------------

a <- txi_kallisto$counts
b <- normalised_counts

# quick check to see whether the dimension are still the same and the values different

dim(a)
dim(b)

colMeans(a) == colMeans(b)

## Construct compendia ----------------------------------------------------

measured_conditions <- read.csv("~/PhD/Data/p_expression/typhimurium_D23580/measured_conditions.csv", sep=";")

tpm_d23580_ids <- measured_conditions[which(measured_conditions$Strain == 'D23580'), 'SRA']

tpm_L474_ids <- measured_conditions[which(measured_conditions$Strain == 'L474'), 'SRA']

### Typhimurium D23580 compendium -------------------------------------------


if (all(tpm_d23580_ids %in% colnames(txi_kallisto$counts))) {

  tpm_d23580_nd <- normalised_counts[ , tpm_d23580_ids]

  if (length(tpm_d23580_ids) != ncol(tpm_d23580_nd))
    stop("check the number of provided IDs and normalised expression data columns")


}


# add pseudocounts of 0.001

tpm_d23580_nd[tpm_d23580_nd == 0] <- 0.001

# verify

min(tpm_d23580_nd)

# calculate log fold change

# tpm_d23580_nd is a matrix of normalised expression data from Typhimurium D23580
# the 'normal condition' has SRA id SRR7814100.
# the FC should therefore be calculated relative to SRR7814100.


function_calc_logFC <- function(row_x) {

  log2(row_x[2] / row_x[1])

}

# preallocate and set names

logFC_matrix <- matrix(, nrow = nrow(tpm_d23580_nd), ncol = (ncol(tpm_d23580_nd) -1))

rownames(logFC_matrix) <- rownames(tpm_d23580_nd)
colnames(logFC_matrix) <- 2:ncol(tpm_d23580_nd)

# add log FC

for (i in 2:ncol(tpm_d23580_nd)) {

  cols <- c(1, i)

  j <- i-1

  logFC_matrix[ ,j] <- apply(tpm_d23580_nd[, cols], 1, function_calc_logFC)

  colnames(logFC_matrix)[j] <- paste (
    colnames(tpm_d23580_nd)[1],
    colnames(tpm_d23580_nd)[i],
    sep = "-", collapse = NULL)


}

dim(logFC_matrix)
# -> 34 contrasts for 4521 genes


write.table(logFC_matrix,
            file="D:/Documenten/PhD/Data/p_expression/typhimurium_D23580/compendium/D23580_expr.txt",
            sep = "\t")


## Typhimurium L474 compendium -------------------------------------------


if (all(tpm_L474_ids %in% colnames(txi_kallisto$counts))) {

  tpm_L474_nd <- normalised_counts[ , tpm_L474_ids]

  if (length(tpm_L474_ids) != ncol(tpm_L474_nd))
    stop("check the number of provided IDs and normalised expression data columns")


}


# add pseudocounts of 0.001

tpm_L474_nd[tpm_L474_nd == 0] <- 0.001


min(tpm_L474_nd)


# calculate log fold change

# calculated relative to SRR7814133


function_calc_logFC <- function(row_x) {

  log2(row_x[2] / row_x[1])

}

logFC_matrix_L474 <- matrix(, nrow = nrow(tpm_L474_nd), ncol = (ncol(tpm_L474_nd) -1))

rownames(logFC_matrix_L474) <- rownames(tpm_L474_nd)
colnames(logFC_matrix_L474) <- 2:ncol(tpm_L474_nd)

for (i in 2:ncol(tpm_L474_nd)) {

  cols <- c(1, i)

  j <- i-1

  logFC_matrix_L474[ ,j] <- apply(tpm_L474_nd[, cols], 1, function_calc_logFC)

  colnames(logFC_matrix_L474)[j] <- paste (
    colnames(tpm_L474_nd)[1],
    colnames(tpm_L474_nd)[i],
    sep = "-", collapse = NULL)


}

dim(logFC_matrix_L474)
# -> 30 contrasts for 4521 genes


write.table(logFC_matrix_L474,
            file="D:/Documenten/PhD/Data/p_expression/typhimurium_D23580/compendium/L474_expr.txt",
            sep = "\t")


# Explore functional expression classes -----------------------------------


## Library ----------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)


## Clusters in the D23580 compendium --------------------------------------

compendium_D23580 <- logFC_matrix

corM = cor(t(compendium_D23580), y = NULL, use = "pairwise.complete.obs", method = c("pearson"))

# View(corM[1:100,1:100])

# get distance using the dist function in base R:
  # dist = dist(corM, method = "euclidian")
  # dendro <- hclust(dist, method="ward")
  # plot(dendro, cex = 0.1, main = "D23580 clusters")

# get distance using a parallel implementation of the dist function

library(parallelDist)
distPar <- parDist(corM, method = "euclidean")
dendro_par <- hclust(distPar, method="ward")
plot(dendro_par, cex = 0.1, main = "D23580 clusters")


# inspect cutoffs and select

abline(h=4000, col="red")

clusters = cutree(dendro, h=4000)

col_fun = colorRamp2(c(-1, -0.08, 0, 0.08, 1), c("green", "darkgreen", "black" ,"darkred", "red"))

is.matrix(corM)
rownames(corM)

ht1 = Heatmap(corM,
              col = col_fun,
              name = "Pearson correlation",
              cluster_rows = dendro_par,
              cluster_columns = dendro_par,
              show_row_names = FALSE,
              show_column_names = FALSE,
              width = unit(10, "cm"),
              height = unit(10, "cm"),
              row_split = 5,
              raster_device = "png",
              show_row_dend = FALSE)



## Clusters in the L474 compendium --------------------------------------

compendium_L474 <- logFC_matrix_L474
any(is.na(compendium_L474))

corML474 = cor(t(compendium_L474), y = NULL, use = "pairwise.complete.obs", method = c("pearson"))

df <- corML474[, sapply(corML474, function(x) { sd(x) == 0} )]


# get distance using a parallel implementation of the dist function

library(parallelDist)
distPar <- parDist(corML474, method = "euclidean")

corML474[,apply(corML474, 2, function(x) any(is.na(x)))]

View(corML474[1:100, 1:100])


dist = dist(corML474, method = "euclidian")

is.na(dist)

dist[is.infinite(dist)] <- 0

dendro <- hclust(dist, method="ward")
plot(dendro, cex = 0.1, main = "D23580 clusters")


dendro_L474 <- hclust(distPar, method="ward")
plot(dendro_L474, cex = 0.1, main = "D23580 clusters")


# inspect cutoffs and select

abline(h=4000, col="red")

clusters = cutree(dendro, h=4000)

col_fun = colorRamp2(c(-1, -0.08, 0, 0.08, 1), c("green", "darkgreen", "black" ,"darkred", "red"))

is.matrix(corM)
rownames(corM)

ht1 = Heatmap(corML474,
              col = col_fun,
              name = "Pearson correlation",
              cluster_rows = dendro_L474,
              cluster_columns = dendro_L474,
              show_row_names = FALSE,
              show_column_names = FALSE,
              width = unit(10, "cm"),
              height = unit(10, "cm"),
              row_split = 5,
              raster_device = "png",
              show_row_dend = FALSE)




