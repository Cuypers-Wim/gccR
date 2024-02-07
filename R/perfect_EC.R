#' Calculate 'Perfect' Expression Conservation (EC) Distribution
#'
#' This function computes the 'perfect' expression conservation (EC) distribution by dividing a gene expression compendium into two parts and calculating the EC between the resulting matrices. The EC for each gene pair is obtained as a weighted correlation between corresponding rows in the two matrices.
#'
#' @param exprM The first 'gene X condition' expression matrix.
#' @param labelsM A correlation (sub)matrix with gene IDs matching those in `exprM`.
#' @param conv Convergence threshold for the EC calculation process.
#' @param maxIter The maximum number of iterations for the EC calculations.
#' @param threads Number of threads to utilize for correlation calculations and row-wise EC calculations.
#' @param ortho Optional: A matrix containing IDs of genes in the first compendium (first column) and corresponding orthologs in the second compendium (second column).
#' @param splits The number of ways to split the dataset for EC calculation. Default is 1. When more than one, the function iteratively calculates EC over multiple splits of the dataset.
#' @param experiment_info Information about the experiments, used to guide the splitting of the dataset.
#' @param tolerance_thresh Tolerance threshold for balancing the size of each split.
#'
#' @return If `splits` is 1, a named vector of EC values representing perfect conservation and a plot showing the EC values per split. If `splits` is more than 1, a list containing EC values for each split, the mean EC value across splits, and variance of EC values per gene.
#'
#' @details The function initially splits the expression matrix `exprM` into two halves and computes the EC for these two parts. When `splits` is more than 1, the function performs multiple splits based on `experiment_info`, calculates EC for each split, and computes the mean and variance of EC values across splits.
#'
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#' perfectEC1 <- perfect_EC(exprList$exprValues2, core_submatrix2_ordered,
#' conv = 0.001, maxIter = 200, threads = 8, ortho = singleCopyOrthologs_matrix)
#'
#' @export

# ToDo:
# - Update info function
# - add functionality to calculate mean per gene when iterating over multiple splits of the dataset

perfect_EC <- function(exprM = NULL, labelsM = NULL, conv = 0.001,
                            maxIter = 200, threads = 1, ortho = NULL,
                            splits = 1, experiment_info = NULL, tolerance_thresh = 0.05) {

  # libraries (for plotting)

  library("reshape2")
  library("ggplot2")
  
  # Define experiments_all based on experiment_info
  
  if (!is.null(experiment_info) && "Experiment_id" %in% colnames(experiment_info)) {
    experiments_all <- unlist(experiment_info["Experiment_id", ])
  } else {
    experiments_all <- colnames(exprM)
  }

  # Check if both inputs are matrices
  if (!is.matrix(exprM) || !is.matrix(labelsM)) {
    stop("Both inputs must be matrices.")
  }
  
  # stop if no similar rownames can be found
  
  common_row_names <- intersect(rownames(exprM), rownames(labelsM))
  
  # Stop if there are no common row names
  if (length(common_row_names) == 0) {
    stop("No common row names found.")
  }

  # select relevant rows

  exprM <- exprM[which(rownames(exprM) %in% rownames(labelsM)), ]

  # multiple splits or single

  if (splits > 1) {
    
    # preallocate #

    exp_sizes_vec <- vector(length = ncol(exprM))
    exp_combos_list <- vector(mode = "list", length = splits)
    EC_list <- vector(mode = "list", length = splits)


    unique_experiment_ids <- unique(experiments_all)
    
    # ToDo - check functionality COLOMBOS experiment info headers
    #expression_cols <- colnames(experiment_info)
    # if (length(experiments_all) != length(expression_cols)) stop("every column in the expression matrix should have an experiment id")
    
    
    for (i in 1:length(unique_experiment_ids)) {

      exp_sizes_vec[i] <- length(which(experiments_all == unique_experiment_ids[i]))
    }

    num_expr <- sum(exp_sizes_vec)
    combo_sum <- 0
    counter <- 0
    tolerance <- length(experiments_all)*tolerance_thresh
    
    while (counter <= splits) {

      # print iteration number

      print(paste("Iterations:  ", counter))

      # resample a set of experiment combinations within the defined boundaries (i.e more or less equal halves +/- tolerance treshold)

      exp_selected_vec <- sample(exp_sizes_vec, (length(unique_experiment_ids)/2))
      
      # Check if all values in exp_sizes_vec are 1
      
        if(all(exp_sizes_vec == 1)) {
          
          # All sizes are 1, handle this case specifically
          message("All experiment sizes are 1. Applying special handling.")
          
          # we want to split the unique_experiment_ids evenly
          half_length <- ceiling(length(unique_experiment_ids) / 2)
          
          # Sample half of the unique experiment IDs directly
          sampled_ids <- sample(unique_experiment_ids, half_length)
          
          # Assuming exp_combos_list and counter are already initialized
          counter <- counter + 1
          exp_combos_list[[counter]] <- sampled_ids
          
          # Since all experiments are of size 1, any split is valid, so we might skip further checks
         
        } else {
        
          combo_sum <- sum(exp_selected_vec)
        if (combo_sum > ((num_expr/2) - tolerance) && combo_sum < ((num_expr/2) + tolerance)) {
          counter <- counter + 1
          exp_combos_list[[counter]] <- unique_experiment_ids[which(exp_sizes_vec %in% exp_selected_vec)]
        } else {
          counter <- counter
        }
      }
      
    } # end while loop counter
    
    # calculate the EC for the resampled sets of experiments
    
    EC_list <- mclapply(exp_combos_list, internal_ec_subset, exprM = exprM, 
                        labelsM = labelsM, experiments_all = experiments_all, mc.cores = threads)

    # add list item with the mean EC value* per gene
    # (*) the mean of the EC values calculated for each of the splits

    # make a data frame first with all values per gene per split

    perfect_ec_df <- as.data.frame(do.call(cbind, EC_list))

    # make a plot of all the splits and add it to a slot

    molten_df <- melt(perfect_ec_df)
    colnames(molten_df) <- c("data_split","EC")

    p <- ggplot(molten_df, aes(x=EC)) +
      geom_density(aes(group=data_split)) +
      theme_classic()

    EC_list$split_plot <- p

    # get the mean per gene and add it to the list

    EC_list$ECfinal <- rowMeans(perfect_ec_df, na.rm = FALSE)

    # get the variance per row and add it to the list

    EC_list$variance <- apply(perfect_ec_df, 1, var)

  } # end if 'multiple splits = TRUE'
  

  if (splits == 1) {

    # split the matrix in half and calculate EC

    # preallocate

    ECresult <- vector(length = nrow(labelsM))
    names(ECresult) <- rownames(labelsM)
    EC_list <- vector(mode = "list", length = 1)
    names(EC_list) <- "EC_final"

    exp_combos_list <- vector(mode = "list", length = 1)

    middle = round((length(experiments_all)/2), digits = 0)
    exp_combos_list[[1]] <- experiments_all[1:middle]
    
    nested_list <- mclapply(exp_combos_list, internal_ec_subset, 
                            exprM = exprM, labelsM = labelsM, 
                            experiments_all = experiments_all, mc.cores = 1)
   
    # Flatten the list
    
    EC_list$EC_final <- do.call(c, nested_list)
    
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
