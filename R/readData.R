#' Read gene expression data
#'
#' This function can be used to read a gene expression compendium (a large table of log Fold Changes)
#'and store it in an R dataframe. If the data is retrieved from COLOMBOS, the option source = "COLOMBOS"
#'can be used. By default (for other data sources), a file tab delimited file must be supplied consisting of
#'gene expression measurements for different genes (rows), across multiple conditions (columns). The first
#'row must consist of the condition identifiers, while the first column must contain gene IDs.
#'
#' @param dataset1 Tab-delimited file of log Fold Changes for different genes (rows) measured under different conditions (columns)
#' @param dataset2 Tab-delimited file of log Fold Changes for different genes (rows) measured under different conditions (columns)
#' @param NAstring The NA string used (e.g "NaN")
#' @param source The data source of which the expression compendium was retrieved (e.g COLOMBOS).
#' Default: expression table (tab delimited file of gene expression measurements) containing row IDs in the first column, and experiment
#' IDs in the first row
#'
#' @return A list (exprList) of 6 named items:
#' \enumerate{
#'   \item $baseName1, $baseName2: basename of the file
#'   \item $exprValues1, $exprValues2: numeric dataframe or matrix of expression values
#'   \item $headerInfo1, $headerInfo2: dataframe containing additional headerlines explaning the experimental conditions
#' }
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#' # read data
#'
#' exprList <- readData(dataset1, dataset2, rownamePos1 = 1,
#' rownamePos2 = 1, NAstring = "NaN", source = "table")
#'
#' # inspect the dataframe of expression values
#'
#' View(head(exprList$exprValues1))
#'
#' @export

# dataset2 = file.path("~/PhD/Data/p_expression/COLOMBOS/exprData/Jan2020/colombos_sente_14028s_exprdata_20151029.txt")
# NAstring = "NaN"


readData <- function(dataset1, dataset2, NAstring = "NaN", source = "table") {
  
  # preallocate

  list1 <- vector(mode = "list", length = 6)

  # read in data depending on the source

  if (source == "table") {

    list1[[1]] <- gsub(pattern = "(.*)\\..*$", replacement = "\\1", basename(dataset1))
    list1[[2]] <- read.delim(dataset1, row.names=1, na.strings=NAstring)
    list1[[3]] <- "NoHeaderLines"
    list1[[4]] <- gsub(pattern = "(.*)\\..*$", replacement = "\\1", basename(dataset2))
    list1[[5]] <- read.delim(dataset2, row.names=1, na.strings=NAstring)
    list1[[6]] <- "NoHeaderLines"

  }

  if (source == "COLOMBOS") {

    list1[[1]] <- gsub(pattern = "(.*)\\..*$", replacement = "\\1", basename(dataset1)) # store the basename of the file (strip path and extension)
    list1[[2]] <- read.delim(dataset1, header=FALSE, na.strings=NAstring) # read the expression table in a dataframe

    list1[[2]] <- list1[[2]][ ,colSums(is.na(list1[[2]])) < nrow(list1[[2]])] # retain columns that are not exclusively 'NA'

    headerInfo1 <- list1[[2]][c(1:6),-c(1:3)] # store the header info for later use
    rownames(headerInfo1) <- c("Test-vs-Reference","Test_description","Reference_description","Experiment_id", "Data_source", "Platform")
    list1[[3]] <- headerInfo1

    list1[[2]] <- list1[[2]][-c(1:7),-c(2:3)]
    rownames1 <- list1[[2]][ ,1]
    list1[[2]] <- list1[[2]][ ,-1]

    list1[[2]] <- apply(list1[[2]], 2, as.numeric)

    colnames1 <- 1:ncol(list1[[2]])
    rownames(list1[[2]]) <- rownames1
    colnames(list1[[2]]) <- colnames1


    # dataset 2
    
    # test <- read.delim(dataset2, , header=FALSE, na.strings=NAstring)
    # View(test)
    
    list1[[4]] <- gsub(pattern = "(.*)\\..*$", replacement = "\\1", basename(dataset2)) # store the basename of the file (strip path and extension)
    list1[[5]] <- read.delim(dataset2, header=FALSE, na.strings=NAstring) # read the expression table in a dataframe
    
    list1[[5]] <- list1[[5]][ ,colSums(is.na(list1[[5]])) < nrow(list1[[5]])] # retain columns that are not exclusively 'NA'
  
    headerInfo2 <- list1[[5]][c(1:6),-c(1:3)] # store the header info for later use
    rownames(headerInfo2) <- c("Test-vs-Reference","Test_description","Reference_description","Experiment_id", "Data_source", "Platform")
    list1[[6]] <- headerInfo2

    list1[[5]] <- list1[[5]][-c(1:7),-c(2:3)]
    rownames2 <- list1[[5]][ ,1]
    list1[[5]] <- list1[[5]][ ,-1]

    list1[[5]] <- apply(list1[[5]], 2, as.numeric)

    colnames2 <- 1:ncol(list1[[5]])
    rownames(list1[[5]]) <- rownames2
    colnames(list1[[5]]) <- colnames2

  }

  # Name list elements

  names(list1) <- c("baseName1", "exprValues1", "headerInfo1", "baseName2", "exprValues2", "headerInfo2")

  # return

  list1

}


