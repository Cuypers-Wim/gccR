#' get functional expression classes
#'
#'
#' @param corM Correlation matrix
#' @param coreGenes Vector of core gene IDs 
#' @param height_cutoff 
#'
#' @return 
#
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
#' 
#'
#' @export


# toDo compare with k means or PCA
# todo add optional for core gene versus no core genes


# corM <- corM_ortho$csM1
# coreGenes <- rownames(corM)

# dim(corM)

# height_cutoff <- 800

get_expression_classes <- function(corM, coreGenes, height_cutoff = 1050) {
  
  # library
  
  library(ComplexHeatmap)
  library(circlize)
  
  # select core genes if required
  
  if(missing(coreGenes)) {
    
    corM <- corM
    
  } else {
    
    corM <- corM[which(rownames(corM) %in% coreGenes), which(rownames(corM) %in% coreGenes)]
    
  }
  
  # compute distances to obtain the distance matrix
  
  dist = dist(corM, method = "euclidian")
  
  # hierarchical clustering
  
  dendro <- hclust(dist, method="ward")
  
  # define clusters based on hclust cutoff
  
  clusters <- cutree(dendro, h = height_cutoff)
  
  number_of_cl <- max(clusters)

  print(paste("the number of clusters is", number_of_cl))
  
  # plot heatmap to inspect the clustering
  
  col_fun = colorRamp2(c(-1, -0.08, 0, 0.08, 1), c("green", "darkgreen", "black" ,"darkred", "red"))
  
  is.matrix(corM)
  rownames(corM)
  
  ht1 = Heatmap(corM, 
                col = col_fun, 
                name = "Pearson correlation", 
                cluster_rows = dendro, 
                cluster_columns = dendro,
                show_row_names = FALSE, 
                show_column_names = FALSE)
                
                
  ht1 = Heatmap(corM, 
                col = col_fun, 
                name = "Pearson correlation", 
                cluster_rows = dendro, 
                cluster_columns = dendro,
                show_row_names = FALSE, 
                show_column_names = FALSE, 
                width = unit(10, "cm"), 
                height = unit(10, "cm"), 
                row_split = clusters, 
                raster_device = "png", 
                show_row_dend = FALSE,
                column_title = title)
  
  outPath <- "~/PhD/Data/GitHub/CCR_fast/plots"
  title <- "typhimuriumLT2"
  savePath <- file.path(outPath,paste(title, "_heatmapEC1", ".pdf", sep = ""))
  
  pdf(savePath, width = 20, height = 20)
  
  draw(ht1)
  
  dev.off()
  

}
