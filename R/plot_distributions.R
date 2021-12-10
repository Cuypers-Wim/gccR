#' Plot distributions
#'
#'
#'
#' @param EC Named vector of EC values
#' @param conserved Named vector of EC values representing conserved gene co-expression
#' @param diverged Named vector of EC values representing diverged gene co-expression

#' @return A  plot of distributions
#
#' @author Wim Cuypers, \email{wim.cuypers@@uantwerpen.be}
#'
#' @examples
#'
#' distr_plot <- plot_distributions(EC$ECfinal, perfectEC, randomEC_vec)
#' distr_plot
#'
#' @export

plot_distributions <- function(EC, conserved, diverged) {
 
  # checks
  
  combinedVec <- c(EC, conserved, diverged)
  
  if (is.vector(combinedVec) & 
      is.numeric(combinedVec) & 
      !is.null(names(combinedVec)) & 
      !any(is.na(names(combinedVec))) == FALSE) stop ("please provide named numeric vectors only")
  
  # define one plotting function to use for both possibilities in the if statement
  
  plotter <- function(df) {
    
    # library 
    
    library("reshape2")
    library("ggplot2")
    
    # plot
    
    molten_df <- melt(df, varnames=c('id', 'distribution'))
    
    p <- ggplot(molten_df, aes(x=value, colour=distribution)) + 
      geom_density() +
      xlim(-0.8, 1.2) +
      theme_classic()
    
    return(p)
    
  }
  

  if (length(EC) == length(conserved) &&
      length(EC) == length(diverged) &&  
      names(EC) == names(conserved) && 
      names(EC) == names(diverged)) {
    
    length(EC) == length(conserved)
    
    print("the provided vectors were equal in length and names")
    
    df <- cbind(EC = EC, conserved = conserved, diverged = diverged)
    
    distributions_plot <- plotter(df)
    
    distributions_plot
    
    # molten_df <- melt(df, varnames=c('id', 'distribution'))
    
    # print(ggplot(molten_df, aes(x=value, colour=distribution)) + geom_density())
    
    
    } else {
    
    print("the provided vectors were not equal in length and names")
    
    EC <- EC[which(names(EC) %in% names(conserved))]
    EC <- EC[which(names(EC) %in% names(diverged))]
    
    conserved <- conserved[which(names(conserved) %in% names(EC))]
    conserved <- conserved[which(names(conserved) %in% names(diverged))]
    
    diverged <- diverged[which(names(diverged) %in% names(conserved))]
    diverged <- diverged[which(names(diverged) %in% names(EC))]
    
    if (length(EC) != length(conserved) &&
        length(EC) != length(diverged) &&
        names(EC) != names(conserved) &&
        names(EC) != names(diverged)) stop ("check input vectors")
    
    df <- cbind(EC = EC, conserved = conserved, diverged = diverged)
    
    distributions_plot <- plotter(df)
    
    distributions_plot
    
    # molten_df <- melt(df, varnames=c('id', 'distribution'))
    
    # print(ggplot(molten_df, aes(x=value, colour=distribution)) + geom_density())
    
    }
  
}


