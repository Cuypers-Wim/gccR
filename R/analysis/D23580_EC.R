
# EC calculations for Typhimurium D23580 ----------------------------------

# Calculate the ...



# Packages ----------------------------------------------------------------

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

  # Check first in editor whether "ID" is included in the file.
  # If not, add this to avoid wrong column naming

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

  # EC

  corM1 <- get_corM_fast(tpm_14028s_expr, dropNArows = TRUE, threads = 8)
  corM2 <- get_corM_fast(tpm_D23580_expr, dropNArows = TRUE, threads = 8)

  corM_ortho <- extract_core_submatrix(corM1, corM2, orthologs)

  csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, orthologs)

  EC <- getEC_fast(corM_ortho$csM1, csM2_ordered, 0.0001, threads = 1)

  summary(EC$ECfinal)

  # conserved and diverged distributions

  ## conserved

  splits <- perfect_EC_fast(exprM = exprList$exprValues2, labelsM = corM_ortho$csM1, conv = 0.001,
                            maxIter = 200, threads = 1, splits = 30,
                            experiment_info = exprList$headerInfo2, tolerance_thresh = 0.25)

  split_halve <- perfect_EC_fast(exprM = exprList$exprValues2, labelsM = corM_ortho$csM1, conv = 0.001,
                            maxIter = 200, threads = 1, experiment_info = exprList$headerInfo2, tolerance_thresh = 0.025)


  summary(splits$ECfinal)

  # EC critical value based on perfect distribution == min(splits$ECfinal, na.rm = TRUE) (-0.38 in this case)
  # truly diverged genes:

  diverged_genes_vec <- EC$ECfinal[EC$ECfinal < -0.38]
  length(diverged_genes_vec)

  # inspect plot

  splits$split_plot

  # write genes and EC to table

  df <- data.frame(EC$ECfinal[EC$ECfinal < -0.38])

  write.table(df, file = "~/PhD/Data/p_expression/typhimurium_D23580/diverged.csv", append = FALSE, sep = "\t", dec = ".")


  # perform GO enrichment on these genes

    # custom mappings generated using gaf file and custom python script
    # were put in the 'example' folder of the package

  goMappingPath <- file.path("D:/Documenten/PhD/Data/p_expression/expression/geneid2go.map")
  output_path = file.path("D:/Documenten/PhD/Data/p_expression/go_D23580")
  goFilename <- "diverged_genes_D23580_GO"

  # we need to change the ids to LT2 IDs

  myInterestingGenes_D23580 <- names(diverged_genes_vec)
  myInterestingGenes <- orthologs_D23580_14028s_LT2[orthologs_D23580_14028s_LT2[ ,2] %in% myInterestingGenes_D23580,3]

  go_enrichment <- topGO_salmonella(goMappingPath, output_path,
                                    filename = goFilename, myInterestingGenes)


  ## diverged

  random_EC <- read.csv("~/PhD/Data/GitHub/CCR_fast/lt2_14028s_random_EC_NEW.csv")

  nrow(random_EC)
  length(rownames(corM_ortho$csM1))

  random_EC$id <- rownames(corM_ortho$csM1)

  randomEC_vec <- c()
  randomEC_vec <- c(random_EC$r_EC)
  names(randomEC_vec) <- c(random_EC$id)

  # plot distributions

 # distrib_names <- c(rep("EC", length(EC$ECfinal)),
 #                    rep("diverged", length(randomEC_vec)),
 #                    rep("conserved", length(final_perfect_EC)))

distrib_names <- c(rep("EC", length(EC$ECfinal)),
                   rep("diverged", length(randomEC_vec)),
                   rep("conserved", length(split_halve[[1]])))

 # values <- c(EC$ECfinal, randomEC_vec,final_perfect_EC)

 values <- c(EC$ECfinal, randomEC_vec, split_halve[[1]])


 df <- cbind(distrib_names, as.data.frame(values))

 p <- ggplot(df, aes(x=values, colour=distrib_names)) +
   geom_density() +
   xlim(-0.8, 1.2) +
   theme_classic()

# Co-expression - EC heatmap ----------------------------------------------


## heatmap - hclust --------------------------------------------------------

  # two corM need to be clustered
  # d23580 will be 'row-ordered' according to the 14028s compendium, but column-reordered according to hclust

  library(ComplexHeatmap)
  library(circlize)

### CorM 1 clustering -------------------------------------------------------

  dist = dist(corM_ortho$csM1, method = "euclidian")
  dendro <- hclust(dist, method="ward")

  plot(dendro)
  abline(h=1050, col="red")
  # define clusters based on hclust cutoff

  clusters_small <- cutree(dendro, h = 500)
  number_of_cl <- max(clusters_small)

  # cluster sizes

  for (i in 1:max(clusters_small)) {

    length_cl <- length(clusters_small[which(clusters_small == i)])

    print(paste0("the length of cluster ", i, " is ", length_cl))

  }

  # make a boxplot of the EC values per cluster

  EC_cl_df <- data.frame(clusters_small, EC$ECfinal)

  ggplot(EC_cl_df, aes(x=as.factor(clusters_small), y=EC.ECfinal)) +
    geom_boxplot()

 # toDo:
    # voeg naam of functie FEC toe, en geef titel met faceting
    # voeg lijnen toe met de EC cutoffs
    # voeg dots toe per meetpunt

  # facet_wrap(~treatment)


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

  cluster_color <- c("#a50026",
    "#d73027",
    "#f46d43",
    "#fdae61",
    "#fee090",
    "#ffffbf",
    "#e0f3f8",
    "#abd9e9",
    "#74add1",
    "#4575b4",
    "#313695")

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
                split = 11,
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

  ht3 = Heatmap(as.factor(clusters),
                col = cluster_color,
                name = "FEC",
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                width = unit(0.5, "cm"),
                raster_device = "png",
                column_title = "FEC")

  ht4 = Heatmap(csM2_ordered_renamed,
                name = "D23580",
                col = col_fun,
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = dendro_c2,
                width = unit(10, "cm"),
                raster_device = "png",
                column_title = "D23580")

  final = (ht1 + ht2+ ht3 + ht4)




  # export heatmap

  title <- "D23580_analysis"
  outPath <- "D:/Documenten/PhD/Data/p_expression/typhimurium_D23580/figures"

  savePath <- file.path(outPath,paste(title, "_heatmapEC2", ".pdf", sep = ""))

  pdf(savePath, width = 20, height = 20)

  draw(final)

  dev.off()

  citation("parallel")


# Functional expression class annotation ----------------------------------

### Gene ontology -----------------------------------------------------------

  # BiocManager::install("topGO")
  # BiocManager::install("Rgraphviz")

  library("topGO")
  library("Rgraphviz")

  goMappingPath <- file.path("D:/Documenten/PhD/Data/p_expression/expression/geneid2go.map") #belangrijk!!

  outputPath = file.path("D:/Documenten/PhD/Data/expression/output_final")
  filename = base1

  # GO_clusters(dendrogram1sub, goMappingPath, outputPath, filename)

  # change labels of dendro (so that the labels of the clusters are from LT2, i.e matching GO mapping)

  dendro_ortho <- dendro

  dendro$labels

  length(dendro_ortho$labels)


  for (i in 1:length(dendro_ortho$labels)) {

    row <- which(orthologs_D23580_14028s_LT2[ ,2] == dendro_ortho$labels[i])

    dendro_ortho$labels[i] <- orthologs_D23580_14028s_LT2[row,3]


  }

  # check
  length(dendro$labels[is.na(dendro_ortho$labels)])
  # 6 ids are 'NA', these are genes that are found in D23580 and 14028s, but not in LT2


    # carry out GO enrichment analysis on a set of genes (clusters)
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

    # plot(dendro_ortho)
    # abline(h=500, col="red")

    clusters <- cutree(dendro_ortho, h = 500)
    numberOfClusters = max(clusters)

    resultSumm_list <- vector(mode = "list", length = numberOfClusters)

    for(i in 1:max(clusters)) {
      j = names(which(clusters == i))
      myInterestingGenes = j
      geneList = factor(as.integer(geneNames %in% myInterestingGenes))
      names(geneList) = geneNames

      # 2. Make topGO object

      GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                    nodeSize = 10,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO)

      # 3. run enrichment test

      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

      # 4. get summary table

      resultSumm_list[[i]] <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 50, numChar = 1000)

      # showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
      # savePath <- file.path(outputPath, paste(filename, "_GOenrichment", i, ".tab", sep = ""))
      # write.table(resultSumm, savePath, sep="\t",row.names=FALSE)
    }

  }


  savePath <- file.path("D:/Documenten/PhD/Data/GitHub/gccR/results/go_cluster14028s")

  for (i in 1:length(resultSumm_list)) {

    filename <- paste0(savePath, "_cl_", i)
    write.table(resultSumm_list[[i]], filename, sep="\t",row.names=FALSE)

    }


  cluster_df <- data.frame(keyName=names(clusters), value=clusters, row.names=NULL)


    filename <- paste0(savePath, "_genes_per_cluster")
    write.table(cluster_df, filename, sep="\t",row.names=FALSE)




## Specialised annotation --------------------------------------------------

    # T3SS genes
    # nuccio & Baumler genes anaerobic metabolism
    # genes identified by wheeler et al.
    #
    # might be STM annotation and not STM14






## All genes (supplementary data) ------------------------------------------




# EC per FEC --------------------------------------------------------------



