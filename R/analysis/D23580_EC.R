
# EC calculations for Typhimurium D23580 ----------------------------------

library(gccR)

# Read data ---------------------------------------------------------------

## Orthologs ---------------------------------------------------------------

  # read data, set whitespace to NA

  D23580_14028s_orthologs <- read.csv("~/PhD/Data/p_expression/typhimurium_D23580/D23580_14028s_orthologs.csv", sep=";", na.strings=c(""))

  # View(D23580_14028s_orthologs)
  # dim(D23580_14028s_orthologs)

  # retain orthologs

  orthologs <- D23580_14028s_orthologs[(!is.na(D23580_14028s_orthologs[ ,2])), ]
  orthologs <- orthologs[ ,c(2,1)]

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

  library(ComplexHeatmap)
  library(circlize)

  dist = dist(corM_ortho$csM1, method = "euclidian")
  dendro <- hclust(dist, method="ward")

  plot(dendro)
  abline(h=1050, col="red")
  # define clusters based on hclust cutoff

  clusters <- cutree(dendro, h = 1050)
  number_of_cl <- max(clusters)


  # plot heatmap to inspect the clustering

  col_fun = colorRamp2(c(-1, -0.08, 0, 0.08, 1), c("green", "darkgreen", "black" ,"darkred", "red"))

  is.matrix(corM_ortho$csM1)
  rownames(corM_ortho$csM1)


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
                show_row_dend = FALSE)

  ht2 = Heatmap(EC$ECfinal,
                name = "EC",
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                width = unit(4, "cm"),
                raster_device = "png",
                column_title = "EC")

  ht3 = Heatmap(csM2_ordered,
                name = "D23580",
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                width = unit(10, "cm"),
                raster_device = "png",
                column_title = "EC")

  final = (ht1 + ht2 + ht3)

  nrow()
