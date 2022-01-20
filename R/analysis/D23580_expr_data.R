
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

dim(a)
dim(b)

# the dimensions remained the same after normalisation

colMeans(a) == colMeans(b)

## select samples that originate from Typhimurium D23580 ------------------

measured_conditions <- read.csv("~/PhD/Data/p_expression/typhimurium_D23580/measured_conditions.csv", sep=";")

tpm_d23580_ids <- measured_conditions[which(measured_conditions$Strain == 'D23580'), 'SRA']

  

if (all(tpm_d23580_ids %in% colnames(txi_kallisto$counts))) {
  
  tpm_d23580_nd <- normalised_counts[ , tpm_d23580_ids]
  
  if (length(tpm_d23580_ids) != ncol(tpm_d23580_nd)) 
    stop("check the number of provided IDs and normalised expression data columns")

  
}


## add pseudocounts of 0.001 ---------------------------------------------

tpm_d23580_nd[tpm_d23580_nd == 0] <- 0.001


min(tpm_d23580_nd)


## calculate log fold change ---------------------------------------------

# tpm_d23580_nd is a matrix of normalised expression data from Typhimurium D23580
# the 'normal condition' has SRA id SRR7814100. 
# the FC should therefore be calculated relative to SRR7814100.


function_calc_logFC <- function(row_x) {
  
  log2(row_x[2] / row_x[1])
  
}

logFC_matrix <- matrix(, nrow = nrow(tpm_d23580_nd), ncol = (ncol(tpm_d23580_nd) -1))

rownames(logFC_matrix) <- rownames(tpm_d23580_nd)
colnames(logFC_matrix) <- 2:ncol(tpm_d23580_nd)

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
  

# Explore functional expression classes -----------------------------------


## Clusters in the D23580 compendium --------------------------------------

compendium_D23580 <- logFC_matrix

corM = cor(t(compendium_D23580), y = NULL, use = "pairwise.complete.obs", method = c("pearson"))

# View(corM[1:100,1:100])

dist = dist(corM, method = "euclidian")

dendro <- hclust(dist, method="ward")


plot(dendro, cex = 0.1, main = "D23580 clusters")


# 3. Add cutoff line (i.e expression classes)

abline(h=4000, col="red")

clusters = cutree(dendro, h=8000)



# Old Script ---------------------------------------------------------------

loopLength <- 6

listFC <- vector(mode = "list", length = loopLength)

for (i in (1:loopLength)) {

dir <- file.path(paste0(paddo, i))

#dir <- file.path("D:/Documenten/PhD/Data/rna-seq/kallisto/exp6")

colDataPath <- file.path(dir, "coldata_design.txt")
colData <- read.delim(colDataPath)
tx2geneFinal <- read.delim(tx2genePath, header=FALSE)

accessionsPath <- file.path(dir, "accessions.txt")
accessions <- read.table(accessionsPath, quote="\"", comment.char="")
# accessions per directory


accVec <- accessions$V1

# make vector of paths 

files <- file.path(dir, paste0(accVec, "_kallisto"), "abundance.tsv")

# give name to each vector element

fileNumber <- length(files)
names(files) <- paste0("sample", 1:fileNumber)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2geneFinal, ignoreAfterBar = FALSE)

#file.exists(files)
#read_counts = txi.kallisto.tsv$counts

#read_counts = as.integer(read_counts)

### DeSeq2 ###


library(DESeq2)

ddsTxi <- DESeqDataSetFromTximport(txi.kallisto.tsv,
                                   colData = colData,
                                   design = ~condition)

keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))

dds <- DESeq(dds)
res <- results(dds)

logFC <- res@listData[["log2FoldChange"]]

names(logFC) <- res@rownames

# add to list

listFC[[i]] <- logFC



n1 <- names(listFC[[1]])
n2 <- names(listFC[[2]])
n3 <- names(listFC[[3]])
n4 <- names(listFC[[4]])
n5 <- names(listFC[[5]])
n6 <- names(listFC[[6]])


allNames <- c(n1, n2, n3, n4, n5, n6)
unique_allNames <- unique(allNames)

#length(allNames)
#length(unique_allNames)

m1 <- matrix(, nrow = length(unique_allNames), ncol = 6)
rownames(m1) <- unique_allNames

i = 1
# fill matrix with values if rowname present, and otherwise leave "NA"

for (i in 1:6) {
for (j in 1:length(unique_allNames)) {
  if (rownames(m1)[j] %in% names(listFC[[i]])) {
    m1[j,i] <- listFC[[i]] [which(names(listFC[[i]]) == rownames(m1)[j])]
 }   
}
}

View(m1)

write.table(m1, file = "Typhi_gene_expression.txt")

nrow(m1)


View(df)



### !!!! #### nodig: log fold change per pairwise comp 


View(res)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

# calculate all possible FC within each experiment

is.matrix(read_counts)

df <- data.frame()

df <- combs1 

combsf = cbind(combs1, combs2, combs3, combs4)

combs <- combn(colnames(read_counts), 2)

combs1 <- combn(colnames(read_counts)[1:6], 2)
combs2 <- combn(colnames(read_counts)[7:18], 2)
combs3 <- combn(colnames(read_counts)[19:25], 2)
combs4 <- combn(colnames(read_counts)[26:31], 2)

merge(combs1,combs2)

# add pseudocounts
read_counts[read_counts == 0] <- 1

write.table(read_counts, file = "read_pseudo_counts")

View(read_counts)
View(combs)
View(combs1)

ncol(combs)
ncol(combsf)

logFC <- function(a, b) (log2(b) - log2(a))

logFC(10,100)


logfoldchanges <- apply(combsf, 2, function(col_names) logFC(read_counts[, col_names[1]], read_counts[, col_names[2]]))

View(logfoldchanges[1:100,1:100])
write.table(logfoldchanges, file = "logFC")

nrow(logfoldchanges)

dimnames(logfoldchanges)[[2]] <- apply(combsf, 2, paste, collapse = '_')

# do hclust on corm

#expr <- logfoldchanges[-which(rowMeans(is.na(logfoldchanges)) > 0.5), ]


View(logfoldchanges)

m2 = m1[-which(rowMeans(is.na(m1)) > 0.5), ]

nrow(m1)
ncol(m1)
nrow(m2)
ncol(m2)
          
          
          
corM = cor(t(m2), y = NULL, use = "pairwise.complete.obs", method = c("pearson"))
View(corM[1:100,1:100])

write.table(corM, file="TyphiCorM.txt")

# replace na by zero
corM[is.nan(corM)] <- 0
View(corM_noNaN)

rownames(corM)
colnames(corM)

corM[2,1:10]

dist = dist(corM_ortho$csM1, method = "euclidian")

dist[4211]

which(is.na(dist))

dist[is.na(dist)] <- 0

head(dist)
# 2. Hierarchical clustering

distanceM = as.matrix(dist)
distanceM[is.na(distanceM)] <- 0

ncol(distanceM)

dendro <- hclust(dist, method="ward")

title = "Typhi clusters"
plot(dendro, cex = 0.1, main = title)


# 3. Add cutoff line (i.e expression classes)

abline(h=8000, col="red")


clusters = cutree(dendro, h=8000)

# determine number of clusters given cutoff height (needed for Heatmap split visualisation)

numberOfClusters = max(clusters)

genes = which(clusters == 1)

length(genes)

library(topGO)

# 1. read/prepare mappings and cluster ids
goMappingPath <- file.path("D:/Documenten/PhD/Data/expression/geneid2go.map") #belangrijk!!
geneID2GO = readMappings(file = goMappingPath)
geneNames = names(geneID2GO)

# 2. initialise test set (gene cluster)

# get clusters
clusters = cutree(dendro, h=1000)

# make vector with "1" for each gene in cluster

myInterestingGenes <- names(which(clusters == 1))
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

# do GO enrichment

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultSumm <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 20)
  
# get genes (of test test) for a given overrepresented GO category (here GO:0035556)

allGO <- genesInTerm(GOdata)
all_genes_in_GO  <- allGO["GO:0035556"]
final_genes <- all_genes_in_GO[[1]][which(all_genes_in_GO[[1]] %in% names(which(clusters == 1)))]

  # 2. Make topGO object
  
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  # 3. run enrichment test
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # 4. get summary table
  
  resultSumm = GenTable(GOdata, classicFisher = resultFisher, topNodes = 20)
  
 
  
  
  savePath <- file.path(outputPath, paste(filename, "_GOenrichment", i, ".tab", sep = ""))
  write.table(resultSumm, savePath, sep="\t",row.names=FALSE)



# make heatmap, plot clustering and cluster

library(ComplexHeatmap)

# 1. color settings

library(circlize)
col_fun = colorRamp2(c(-1, -0.08, 0, 0.08, 1), c("green", "darkgreen", "black" ,"darkred", "red"))
cols1 = c("white", "red")
cols = c("darkred", "darkgreen", "white")
col_fun2 = colorRamp2(c(-0.05, -0.01, 0, 0.01, 0.05), c("green", "darkgreen", "black" ,"darkred", "red"))


#2. plot and save
outPath <- file.path("D:/Documenten/PhD/Data/rna-seq")


savePath <- file.path(outPath,paste("heatmapTyphi", ".pdf", sep = ""))

ht1 = Heatmap(corM, 
              col = col_fun, 
              name = "Pearson correlation", 
              cluster_rows = dendro, 
              cluster_columns = dendro,
              show_row_names = FALSE, 
              show_column_names = FALSE, 
              width = unit(10, "cm"), 
              height = unit(10, "cm"), 
              #row_split = numberOfClusters, 
              raster_device = "png", 
              show_row_dend = FALSE)




pdf(savePath, width = 20, height = 20)

draw(ht1)

dev.off()


exprM_corM_list1 = correlationMatrix(expression1)

corM


rownames(corM)

corM_ortho = extract_core_submatrix(exprM_corM_list1[[2]], corM, singleCopyOrthologs)
core_submatrix2_ordered = sort_matrix(corM_ortho$csM1, corM_ortho$csM2, corM_ortho$orthologs)
ECfinaleTemp = getEC_finalTemp(corM_ortho$csM1, core_submatrix2_ordered, 0.001)

ECfinaleTemp$ECfinal


plotECdistribution(as.data.frame(ECfinaleTemp$ECfinal))

min(ECfinaleTemp$ECfinal)
ECfinaleTemp$ECfinal[ECfinaleTemp$ECfinal < -0.5]











