#' perfect EC distribution
#'
#'   Wrapper function to perform GO enrichment for Salmonella genes using the TopGO package
#'
#' @param goMappingPath Path to the GO mapping file
#' (i.e. a file that associates each gene id with a set of GO terms)
#' @param outputPath Path for writing output files
#' @param myInterestingGenes Vector of gene IDs which represent the subset of genes of interest
#' @param ontology_type String indicating which of the three ontologies should be used.
#' (Options: "BP", "MF" and "CC", which refers to Biological Process, Molecular Function and Cellular Component)
#' @param top_nodes Number of top nodes that should be shown (Shows the 10 top scoring terms by default)
#'
#' @return A  named vector of EC values representing perfect conservation
#
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
#' go_enrichment <- topGO_salmonella(goMappingPath, output_path,
#' filename = "EC_diverged_GO",myInterestingGenes)
#'
#' @export
#'

topGO_salmonella <- function(goMappingPath = "examples/typhimurium_geneid2go.map",
                             output_path, filename = "go_output", myInterestingGenes, ontology_type = "BP", top_nodes = 10) {

  # library

  library("topGO")

  # make output file path

  outputPath <- file.path(output_path, paste(filename, ".txt", sep=""))

  # parse input for topGO object

  geneID2GO <- readMappings(goMappingPath)
  geneNames = names(geneID2GO)

  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  str(geneList)

  # BP = Biological process namespace

  GOdata <- new("topGOdata", ontology = ontology_type, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)

  # get stats and write output

  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultSumm = GenTable(GOdata, classicFisher = resultFisher, topNodes = top_nodes)

  write.table(resultSumm, outputPath, sep="\t",row.names=FALSE)

  # return output as df

  return(resultSumm)

}



