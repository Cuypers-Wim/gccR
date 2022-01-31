#' perfect EC distribution
#'
#'   Calculate 'perfect' expression conservation (EC) distribution. This function splits a compendium
#'   in two parts, and calculates the EC for the two resulting matrices (i.e the weighted correlation of
#'   each row of matrix 1 (m1) with the corresponding row in matrix 2 (m2) is calculated).
#'
#' @param exprM First 'gene X condition' expression matrix
#' @param labelsM Correlation (sub)matrix. Gene IDs must be the same as used in the exprM
#' @param conv Convergence treshold for EC calculations
#' @param maxIter Maximum number of iterations for EC calculations
#' @param threads Number of threads to use for cor calculations and row-wise EC calculations
#' @param ortho OPTIONAL argument. Matrix consisting of IDs of genes used in the 
#' first compendium (first column) and gene IDs of the corresponding orthologs 
#' found in the secondcompendium
#'
#' @return A  named vector of EC values representing perfect conservation
#
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
#' perfectEC1  <- perfect_EC(exprList$exprValues2, core_submatrix2_ordered, 
#' conv = 0.001, maxIter = 200, threads = 8, ortho = singleCopyOrthologs_matrix)
#'
#' @export
#' 

# todo: update info function

perfect_EC_fast <- function(exprM, labelsM, conv = 0.001,
                            maxIter = 200, threads = 8, ortho,
                            multiple_splits = TRUE, experiment_info, tolerance_thresh = 0.05) {
  
  ##### helper function #####
  
  function_ec_subset <- function(exp_combo, exprM, labelsM) {
    
    # get EC for subsets of an expression matrix grouped per sample
    
    half_exprM1 <- exprM[ , which(experiments_all %in% exp_combo)]
    half_exprM2 <- exprM[ , -(which(experiments_all %in% exp_combo))]
    
    corM1 <- get_corM_fast(half_exprM1, dropNArows = TRUE, threads = 1)
    corM2 <- get_corM_fast(half_exprM2, dropNArows = TRUE, threads = 1)
    
    subCorM1 = corM1[which(rownames(corM1) %in% rownames(labelsM)),
                     which(rownames(corM1) %in% rownames(labelsM))]
    subCorM2 = corM2[which(rownames(corM2) %in% rownames(labelsM)),
                     which(rownames(corM2) %in% rownames(labelsM))]
    
    subCorM <- extract_core_submatrix(subCorM1, subCorM2)
    csM2_ordered <- sort_matrix(subCorM$csM1, subCorM$csM2)
    EC <- getEC_fast(subCorM$csM1, csM2_ordered, conv, maxIter, threads = threads)
    ECresult <- EC$ECfinal
    ECresult
    
  } # end of function_ec_subset

  
  # select relevant rows
  
  exprM <- exprM[which(rownames(exprM) %in% rownames(labelsM)), ]
  
  # unlist experiment info
  
  experiments_all <- unlist(experiment_info["Experiment_id", ])
  
  # multiple splits or single
  
  if (multiple_splits == TRUE) {
    
    # preallocate #
    
    exp_sizes_vec <- vector(length = ncol(exprM))
    exp_combos_list <- vector(mode = "list", length = 10)
    EC_list <- vector(mode = "list", length = 10)
    
   
    unique_experiment_ids <- unique(unlist(experiment_info["Experiment_id", ]))
    expression_cols <- colnames(experiment_info)
    
    if (length(experiments_all) != length(expression_cols)) stop("every row in the expression matrix should have an experiment id")
    
    for (i in 1:length(unique_experiment_ids)) {
      
      exp_sizes_vec[i] <- length(which(experiments_all == unique_experiment_ids[i]))
    }
    
    num_expr <- sum(exp_sizes_vec)
    combo_sum <- 0
    counter <- 0
    tolerance <- length(experiments_all)*tolerance_thresh
    
    while (counter < 10) {
      exp_selected_vec <- sample(exp_sizes_vec, (length(unique_experiment_ids)/2))
      combo_sum <- sum(exp_selected_vec)
      if (combo_sum > ((num_expr/2) - tolerance) && combo_sum < ((num_expr/2) + tolerance)) {
        counter <- counter + 1
        exp_combos_list[[counter]] <- unique_experiment_ids[which(exp_sizes_vec %in% exp_selected_vec)]
      } else {
        counter <- counter
      }
    }
    
    EC_list <- mclapply(exp_combos_list, function_ec_subset, exprM, labelsM, mc.cores = threads)
    
  } # end if 'multiple splits = TRUE'
  
  if (multiple_splits == FALSE) {
    
    # just split the matrix in half and calculate EC
    
    # preallocate
    
    ECresult <- vector(length = nrow(labelsM))
    names(ECresult) <- rownames(labelsM)
    exp_combos_list <- vector(mode = "list", length = 1)
    
    middle = round((length(experiments_all)/2), digits = 0)
    exp_combos_list[[1]] <- experiments_all[1:middle]
    EC_list <- mclapply(exp_combos_list, function_ec_subset, exprM, labelsM, mc.cores = 1)
    
  } # end if 'multiple splits = FALSE'
  
  ##### renaming, depending on the 'ortho' argument #####
  
  if(missing(ortho)) {
    
    # ECresult <- EC_list
    return(EC_list)
    
  } else {
    
    function_rename_EC <- function(EC_list_item, ortho) {
      
      for(i in 1:length(EC_list_item)) {
        
        # lookup orthologID and corresponding location in core submatrix 2
        # and replace the ID by the ID of the ortholog
        
        ECid <- names(EC_list_item)[i]
        orthoID <- ortho[which(ortho[ ,2] == ECid), 1]
        names(EC_list_item)[i] <- orthoID
        
      } # end of for loop
      
      EC_list_item
      
    } # end of 'function_rename_EC'
    
    ECresult_renamed <- lapply(EC_list, function_rename_EC, ortho)
    
    return(ECresult_renamed)
    
  } # end if else
  
  # return
  # ECresult
  
} # end of main function