
# EC calculations for Typhimurium D23580 ----------------------------------

# Calculate the ...
devtools::document()
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

  # get the number of orthologs identified

  orthologs_true <- orthologs[(!is.na(orthologs[ ,1])), ]

  nrow(orthologs_true)

## D23580 ------------------------------------------------------------------

  # Check first in editor whether "ID" is included in the file.
  # If not, add this to avoid wrong column naming

  tpm_D23580_expr <- read.delim("~/PhD/Data/p_expression/typhimurium_D23580/compendium/D23580_expr.txt")

  rownames(tpm_D23580_expr) <- tpm_D23580_expr[ ,1]
  tpm_D23580_expr <- tpm_D23580_expr[ ,-1]

  # View(tpm_D23580_expr)

## strain 4/74 (from the same study as D23580) ------------------------------

  tpm_474_expr <- read.delim("D:/Documenten/PhD/Data/p_expression/typhimurium_D23580/compendium/L474_expr.txt")

  ncol(tpm_474_expr)

## 14028s ------------------------------------------------------------------

  dataset1 <- file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_lt2_exprdata_20151029.txt")
  dataset2 <- file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_14028s_exprdata_20151029.txt")
  exprList <- readData(dataset1, dataset2, NAstring = "NaN", source = "COLOMBOS")

  tpm_14028s_expr <- exprList$exprValues2

  # View(tpm_14028s_expr)

df <- exprList$headerInfo2["Experiment_id", , drop = FALSE]

length(unique(as.character(df)))

# EC calculation ----------------------------------------------------------

  ## EC ----------------------------------------------------------------------

  corM1 <- get_corM_fast(tpm_14028s_expr, dropNArows = TRUE, threads = 8)
  corM2 <- get_corM_fast(tpm_D23580_expr, dropNArows = TRUE, threads = 8)

  corM_ortho <- extract_core_submatrix(corM1, corM2, orthologs)

  csM2_ordered <- sort_matrix(corM_ortho$csM1, corM_ortho$csM2, orthologs)

  EC <- getEC_fast(corM_ortho$csM1, csM2_ordered, 0.0001, threads = 1)

  summary(EC$ECfinal)

  length(EC$ECfinal)

  # shapiro.test(EC$ECfinal) -> not normally distributed


  ## EC and divergence -------------------------------------------------------

  # calculate conserved and diverged distributions and
  # use the min and max values to identify truly diverged and conserved genes


    ### Diverged ----------------------------------------------------------------

    splits <- perfect_EC_fast(exprM = exprList$exprValues2, labelsM = corM_ortho$csM1, conv = 0.001,
                              maxIter = 200, threads = 1, splits = 30,
                              experiment_info = exprList$headerInfo2, tolerance_thresh = 0.25)

    split_halve <- perfect_EC_fast(exprM = exprList$exprValues2, labelsM = corM_ortho$csM1, conv = 0.001,
                              maxIter = 200, threads = 1, experiment_info = exprList$headerInfo2, tolerance_thresh = 0.025)


    summary(splits$ECfinal)

    # EC critical value based on perfect distribution == min(splits$ECfinal, na.rm = TRUE) (-0.38 in this case)
    # truly diverged genes:

    diverged_genes_vec <- EC$ECfinal[EC$ECfinal < -0.35]
    length(diverged_genes_vec)

    # write genes and EC to table

    df <- data.frame(EC$ECfinal[EC$ECfinal < -0.35])

    write.table(df, file = "~/PhD/Data/p_expression/typhimurium_D23580/diverged.csv", append = FALSE, sep = "\t", dec = ".")

    read.table(file = "~/PhD/Data/p_expression/typhimurium_D23580/diverged.csv")

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
                                      filename = goFilename, myInterestingGenes, top_nodes = 40)

    ### Plot perfect EC splits + variance  -----------------------------------------

    # plot of all perfect distributions per split

    splits$split_plot

    # plot of the variance

    ec_variance_df <- cbind(splits$ECfinal, splits$variance)
    colnames(ec_variance_df) <- c("mean_EC_value","variance")

    library(ggplot2)
    ggplot(ec_variance_df, aes(x = mean_EC_value, y = variance)) +
      geom_point()

    plot(x=ec_variance_df[ , "mean_EC_value"],
         y=ec_variance_df[ , "variance"],
         xlab="Mean EC value", ylab="variance",
         main="EC versus variance")


    ### Conserved ---------------------------------------------------------------

    random_EC <- read.csv("~/PhD/Data/GitHub/CCR_fast/lt2_14028s_random_EC_NEW.csv")

    summary(random_EC)

    conserved_genes_vec <- EC$ECfinal[EC$ECfinal > 0.58  ]
    length(conserved_genes_vec)

    genenames_D23580 <- names(conserved_genes_vec)

    myInterestingGenes_conserved <- orthologs_D23580_14028s_LT2[orthologs_D23580_14028s_LT2[ ,2] %in% genenames_D23580,3]

    go_enrichment_conserved <- topGO_salmonella(goMappingPath, output_path,
                                      filename = "conserved_genes_go", myInterestingGenes = myInterestingGenes_conserved, top_nodes = 40)


    randomEC_vec <- c()
    randomEC_vec <- random_EC$r_EC
    names(randomEC_vec) <- random_EC$id

# Plot EC, conserved and diverged distributions --------------------------------


distrib_names <- c(rep("EC", length(EC$ECfinal)),
                   rep("diverged", length(randomEC_vec)),
                   rep("conserved", length(splits$ECfinal)))

values <- c(EC$ECfinal, randomEC_vec, splits$ECfinal)

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

  dist <- dist(corM_ortho$csM1, method = "euclidian")

  dendro <- hclust(dist, method="ward")

  plot(dendro)
  abline(h=600, col="darkred")
  # define clusters based on hclust cutoff

  superClusters <- cutree(dendro, h = 2500)
  clusters_small <- cutree(dendro, h = 600)
  number_of_cl <- max(clusters_small)

  # cluster sizes

  for (i in 1:max(clusters_small)) {

    length_cl <- length(clusters_small[which(clusters_small == i)])

    print(paste0("the length of cluster ", i, " is ", length_cl))

  }

  # make a boxplot of the EC values per cluster

  EC_cl_df <- data.frame(clusters_small, EC$ECfinal)

  ggplot(EC_cl_df, aes(x=as.factor(clusters_small), y=EC.ECfinal)) +
    geom_boxplot() +
    theme_classic() +
    labs(y = "EC", x = "Functional expression class")


  ggplot(EC_cl_df, aes(x=as.factor(clusters_small), y=EC.ECfinal)) +
    geom_violin() +
    geom_boxplot(width=0.2, fill="white") +
    theme_classic() +
    labs(y = "EC", x = "Functional expression class") +
    geom_hline(yintercept=0.58, linetype="dashed", color = "darkgreen") +
    geom_hline(yintercept=-0.35, linetype="dashed", color = "darkred")



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

  # palette from https://medialab.github.io/iwanthue/
  cluster_color <- c("#36dee6",
                      "#a27f38",
                      "#9750a1",
                      "#65913a",
                      "#ba496b",
                      "#6bc66f",
                      "#6777cf",
                      "#bb5437",
                      "#45bc8d",
                      "#c6a83e")

  names(cluster_color) <- c("1","2","3","4","5","6","7","8","9", "10")

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
                row_split = max(clusters_small),
                column_split = max(superClusters),
                row_gap = unit(0.5, "mm"),
                column_gap = unit(0.5, "mm"),
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

  ht3 = Heatmap(as.factor(clusters_small),
                col = cluster_color,
                name = "FEC",
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                width = unit(0.5, "cm"),
                raster_device = "png",
                column_title = "FEC")

  ht4 = Heatmap(as.factor(superClusters),
                name = "Cluster",
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                width = unit(0.5, "cm"),
                raster_device = "png",
                column_title = "Cluster")

  final = (ht1 + ht2+ ht3)

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


## EC and divergence  -> which FECls? --------------------------------------

  ec_divergence_clusters <- clusters_small[names(diverged_genes_vec)]

View(data.frame(ec_divergence_clusters))

ec_conserved_clusters <- clusters_small[names(conserved_genes_vec)]

table(ec_conserved_clusters)

length(conserved_genes_vec)

### Gene ontology -----------------------------------------------------------

  # GO enrichment on functional expression classes (FECls)

  goMappingPath <- file.path("D:/Documenten/PhD/Data/p_expression/expression/geneid2go.map") #belangrijk!!

  outputPath = file.path("D:/Documenten/PhD/Data/p_expression/221108_latest_results")

  # preallocate vector for the new ids

  clusters_newnames <- vector(length = length(clusters_small))

  for (i in 1:length(clusters_small)) {

    id_position <- which(orthologs_D23580_14028s_LT2[ , 2] == names(clusters_small)[i])

    clusters_newnames[i] <- orthologs_D23580_14028s_LT2[id_position, 3]

  }

  clusters_small_renamed <- clusters_small
  names(clusters_small_renamed) <- clusters_newnames

  # perform GO enrichment

    # results will be stored in a list (and written to a file - see function details)

  go_cluster_list <- vector(mode = "list", length = max(clusters_small))

  for (i in 1:max(clusters_small_renamed)) {

  fecl_n <- clusters_small_renamed[clusters_small_renamed == i]

  go_fecls <- topGO_salmonella(goMappingPath = goMappingPath, output_path = outputPath,
                                              filename = paste("go_fecl_small", i, sep=""),
                               myInterestingGenes = names(fecl_n),
                               top_nodes = 30)

  go_cluster_list[[i]] <- go_fecls

  } # end of for loop

  View(go_cluster_list[[10]])

# Specialised annotation --------------------------------------------------

    clusters_small

## T3SS genes --------------------------------------------------------------

    ids_T3SS1_df <- read.table("~/PhD/Data/sp_T3SS/SPI1_IDs.txt", quote="\"", comment.char="")
    ids_T3SS2_df <- read.table("~/PhD/Data/sp_T3SS/SPI2_IDs.txt", quote="\"", comment.char="")

    ids_T3SS1_vec <- ids_T3SS1_df[ ,1]
    ids_T3SS2_vec <- ids_T3SS2_df[ ,1]

## Anaerobic metabolism genes (nuccio & baumler) ---------------------------

    ids_anaerobic_df <- read.table("~/PhD/Data/GitHub/concord/invasiveness/gene_ids_anaerobic_n&b.txt", quote="\"", comment.char="")
    ids_anaerobic_vec <- ids_anaerobic_df[ ,1]

## Host adaptation genes (Wheeler et al.) ----------------------------------

    # all but 6 genes have a STM identifier

    ids_wheeler2018_df <- read.table("~/PhD/Data/GitHub/concord/invasiveness/topgenes_wheeler2018.txt", quote="\"", comment.char="")
    ids_wheeler2018_vec <- ids_wheeler2018_df[ ,1]



## Count occurence of gene sets across expression classes ------------------

    subset <- clusters_small_renamed[which(names(clusters_small_renamed) %in% ids_T3SS1_vec)]
    ids_T3SS2_vec

    length(ids_wheeler2018_vec)

    length(ids_anaerobic_vec)

    list_of_geneVecs <- list("T3SS1" = ids_T3SS1_vec,
                             "T3SS2" = ids_T3SS2_vec,
                             "anaerobic" = ids_anaerobic_vec,
                             "hostAdapt" = ids_wheeler2018_vec)


    for (geneSet in names(list_of_geneVecs)) {

      subset <- clusters_small_renamed[which(names(clusters_small_renamed) %in% list_of_geneVecs[[geneSet]])]

         for (clusterName in 1:10) {

          print(paste(geneSet,
                        "genes in cluster", clusterName, ":", length(which(subset == clusterName)), sep = " "))

            } # end of inner loop

    } # end of outer loop



## Get EC per extra annotation set -----------------------------------------

    # STM14 identifier were used in stead of STM (Typhimurium LT2)

    EC$ECfinal


## make binary annotation matrix for heatmap visualisation

    # get names core submatrix and change 14028s IDs to LT2 IDs

    allNames_vec <- orthologs_D23580_14028s_LT2[names(EC$ECfinal) %in% orthologs_D23580_14028s_LT2[ ,2] , 3]

    # preallocae matrix for binary values

    nrows <- length(allNames_vec)

    binary_anno_matrix <- matrix(0, ncol = 4, nrow = nrows)
    colnames(binary_anno_matrix) <- c("T3SS1", "T3SS2", "Anaerobic", "Extraintestinal")
    rownames(binary_anno_matrix) <- allNames_vec


    for (i in 1:length(allNamesVec)) {
      if (allNamesVec[i] %in% ids_T3SS1_vec) {
        binary_anno_matrix[i, 1] <- 1
      }
      else if (allNamesVec[i] %in% ids_T3SS2_vec) {
        binary_anno_matrix[i, 2] <- 1
      }
      else if (allNamesVec[i] %in% ids_anaerobic_vec) {
        binary_anno_matrix[i, 3] <- 1
      }
      else if (allNamesVec[i] %in% ids_wheeler2018_vec) {
        binary_anno_matrix[i, 4] <- 1
      }

      } # end of for loop

    # replace by STM14 identifiers and rmeove NA

    stm14rownames <- orthologs_D23580_14028s_LT2[rownames(binary_anno_matrix) %in% orthologs_D23580_14028s_LT2[ ,3] , 2]

    rownames(binary_anno_matrix) <- stm14rownames

    ba_fmatrix <- binary_anno_matrix[!is.na(rownames(binary_anno_matrix)), ]

    ba_fmatrix <- ba_fmatrix[rownames(ba_fmatrix) %in% names(EC$ECfinal), ]

    # binary_anno_matrix can now be added to the complex heatmap

    # problem, the binary matrix has STM IDs, while the other components of the heatmap have STM14 IDs

    ht4 = Heatmap(as.factor(ba_fmatrix[ ,1]),
                  col = c("0" = "#f1f1f1", "1" = "#a6cee3"),
                  name = "T3SS1",
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  cluster_rows = FALSE,
                  width = unit(0.5, "cm"),
                  raster_device = "png")
    ht5 = Heatmap(as.factor(ba_fmatrix[ ,2]),
                  col = c("0" = "#f1f1f1", "1" = "#1f78b4"),
                  name = "T3SS2",
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  cluster_rows = FALSE,
                  width = unit(0.5, "cm"),
                  raster_device = "png")
    ht6 = Heatmap(as.factor(ba_fmatrix[ ,3]),
                  col = c("0" = "#f1f1f1", "1" = "#b2df8a"),
                  name = "Anaerobic",
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  cluster_rows = FALSE,
                  width = unit(0.5, "cm"),
                  raster_device = "png")
    ht7 = Heatmap(as.factor(ba_fmatrix[ , 4]),
                  col = c("0" = "#f1f1f1", "1" = "#33a02c"),
                  name = "Extraintestinal",
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  cluster_rows = FALSE,
                  width = unit(0.5, "cm"),
                  raster_device = "png")



    final = (ht1 + ht2+ ht3 + ht4 + ht5 + ht6 + ht7)


    # test whether T3SS is in different lcusters indeed

    ids_T3SS1_vec
    ids_T3SS2_vec



    orthologs_D23580_14028s_LT2[names(clusters_small) %in% orthologs_D23580_14028s_LT2[ ,2] , 3]



## All genes (supplementary data) ------------------------------------------




# EC per FEC --------------------------------------------------------------



