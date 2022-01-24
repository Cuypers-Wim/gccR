
# EC calculations for Typhimurium D23580 ----------------------------------

library(gccR)

# Read data ---------------------------------------------------------------

## Orthologs ---------------------------------------------------------------

  # read data, set whitespace to NA

  orthologs_D23580_14028s_LT2 <- read.csv("~/PhD/Data/p_expression/typhimurium_D23580/D23580_14028s_orthologs.csv", sep=";", na.strings=c(""))

  D23580_14028s_orthologs <- orthologs_D23580_14028s_LT2[ ,c(2,1)]

  # View(D23580_14028s_orthologs)
  # dim(D23580_14028s_orthologs)

  # retain genes present in D23580

  orthologs <- D23580_14028s_orthologs[(!is.na(D23580_14028s_orthologs[ ,2])), ]

## D23580 ------------------------------------------------------------------

  tpm_D23580_expr <- read.delim("~/PhD/Data/p_expression/typhimurium_D23580/compendium/D23580_expr.txt")
  rownames(tpm_D23580_expr) <- tpm_D23580_expr[ ,1]
  tpm_D23580_expr <- tpm_D23580_expr[ ,-1]

  # View(tpm_D23580_expr)

## 14028s ------------------------------------------------------------------

  dataset1 <- file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_lt2_exprdata_20151029.txt")
  dataset2 <- file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_14028s_exprdata_20151029.txt")
  exprList <- readData(dataset1, dataset2, NAstring = "NaN", source = "COLOMBOS")

  tpm_14028s_expr <- exprList$exprValues2

  # View(tpm_14028s_expr)


# EC calculation ----------------------------------------------------------

  corM1 <- get_corM_fast(tpm_14028s_expr, dropNArows = TRUE, threads = 8)
  corM2 <- get_corM_fast(tpm_D23580_expr, dropNArows = TRUE, threads = 8)

  corM_ortho <- extract_core_submatrix(corM1, corM2, orthologs)

  csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, orthologs)

  EC <- getEC_fast(corM_ortho$csM1, csM2_ordered, 0.001, threads = 1)

  EC$ECfinal[EC$ECfinal < -0.49]

  EC$ECfinal[EC$ECfinal > 0.58]


# Co-expression - EC heatmap ----------------------------------------------


## heatmap - hclust --------------------------------------------------------

  # two corM need to be clustered
  # d23580 will be 'row-ordered' according to the 14028s compendium, but column-rdered accoding to hclust

  library(ComplexHeatmap)
  library(circlize)

### CorM 1 clustering -------------------------------------------------------

  dist = dist(corM_ortho$csM1, method = "euclidian")
  dendro <- hclust(dist, method="ward")

  plot(dendro)
  abline(h=1050, col="red")
  # define clusters based on hclust cutoff

  clusters <- cutree(dendro, h = 1050)
  number_of_cl <- max(clusters)


### CorM 2 clustering -------------------------------------------------------

  csM2_ordered_renamed <- csM2_ordered

  rownames(csM2_ordered_renamed) <- rownames(corM_ortho$csM1)
  colnames(csM2_ordered_renamed) <- colnames(corM_ortho$csM1)

  dist_c2 = dist(csM2_ordered_renamed, method = "euclidian")
  dendro_c2 <- hclust(dist_c2, method="ward")

  plot(dendro_c2)
  abline(h=1050, col="red")
  # define clusters based on hclust cutoff

  clusters_c2 <- cutree(dendro_c2, h = 1050)
  number_of_cl_c2 <- max(clusters_c2)


### Plot Heatmap ------------------------------------------------------------

  col_fun = colorRamp2(c(-1, -0.08, 0, 0.08, 1), c("green", "darkgreen", "black" ,"darkred", "red"))

  # checks

  is.matrix(corM_ortho$csM1) && is.matrix(csM2_ordered_renamed)
  all(rownames(corM_ortho$csM1) == rownames(csM2_ordered_renamed))

  #plot

  ht1 = Heatmap(corM_ortho$csM1,
                col = col_fun,
                name = "Pearson correlation",
                cluster_rows = dendro,
                cluster_columns = dendro,
                show_row_names = FALSE,
                show_column_names = FALSE,
                width = unit(10, "cm"),
                height = unit(10, "cm"),
                row_split = 6,
                raster_device = "png",
                show_row_dend = FALSE,
                column_title = "14028s")

  ht2 = Heatmap(EC$ECfinal,
                name = "EC",
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                width = unit(0.5, "cm"),
                raster_device = "png",
                column_title = "EC")

  ht3 = Heatmap(csM2_ordered_renamed,
                name = "D23580",
                col = col_fun,
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = dendro_c2,
                width = unit(10, "cm"),
                raster_device = "png",
                column_title = "D23580")

  final = (ht1 + ht2 + ht3)

  # export heatmap

  title <- "D23580_analysis"
  outPath <- "D:/Documenten/PhD/Data/p_expression/typhimurium_D23580/figures"

  savePath <- file.path(outPath,paste(title, "_heatmapEC1", ".pdf", sep = ""))

  pdf(savePath, width = 20, height = 20)

  draw(final)

  dev.off()


# Functional expression class annotation ----------------------------------


## Orthologs only ----------------------------------------------------------


  # GO enrichment
  # T3SS genes
  # nuccio & Baumler genes anaerobic metabolism
  # genes identified by wheeler et al.
  #

  # might be STM annotation and not STM14

### Gene ontology -----------------------------------------------------------

  BiocManager::install("topGO")

  library("topGO")

  goMappingPath <- file.path("D:/Documenten/PhD/Data/p_expression/expression/geneid2go.map") #belangrijk!!

  outputPath = file.path("D:/Documenten/PhD/Data/expression/output_final")
  filename = base1

  GO_clusters(dendrogram1sub, goMappingPath, outputPath, filename)

  # change labels of dendro (so that the labels of the clusters are from LT2, i.e matching GO mapping)

  dendro_ortho <- dendro

  i = 1000

  length(dendro_ortho$labels)

  for (i in 1:length(dendro_ortho$labels)) {

    row <- which(orthologs_D23580_14028s_LT2[ ,2] == dendro_ortho$labels[i])

    dendro_ortho$labels[i] <- orthologs_D23580_14028s_LT2[i,3]


  }

  # orthologs_D23580_14028s_LT2

  GO_clusters <- function(dendro, goMappingPath, outputPath, filename) {

    # carries out GO enrichment analysis on a set of genes (clusters)
    #
    # Args:
    #  goMappingPath: path to the custom GO mapping file for topGO
    #  nameVector: vector containing the names of the cluster objects
    #         cluster objects: list of gene identifiers from given expression class
    #
    # Return:
    #   List with Tables of enriched GO classes per gene cluster

    library(topGO)

    # 1. read/prepare mappings and cluster ids

    geneID2GO = readMappings(file = goMappingPath)
    geneNames = names(geneID2GO)

    # 2. initialise test set (gene cluster)

    clusters <- cutree(dendro_ortho, h = 1050)
    numberOfClusters = max(clusters)

    resultSumm_list <- vector(mode = "list", length = numberOfClusters)

    for(i in 1:max(clusters)) {
      j = names(which(clusters == i))
      myInterestingGenes = j
      geneList = factor(as.integer(geneNames %in% myInterestingGenes))
      names(geneList) = geneNames

      # 2. Make topGO object

      GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

      # 3. run enrichment test

      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

      # 4. get summary table

      resultSumm_list[[i]] <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 50)

      # savePath <- file.path(outputPath, paste(filename, "_GOenrichment", i, ".tab", sep = ""))
      # write.table(resultSumm, savePath, sep="\t",row.names=FALSE)
    }

  }



## All genes (supplementary data) ------------------------------------------




# EC per FEC --------------------------------------------------------------



