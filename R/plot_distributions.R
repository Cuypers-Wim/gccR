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
  
  
  df <- data.frame(value = c(EC, conserved, diverged),
                   distribution = factor(rep(c("EC", "Conserved", "Diverged"), 
                                             times = c(length(EC), length(conserved), length(diverged)))))
  
  # Plot
  ggplot(df, aes(x = value, colour = distribution)) + 
    geom_density() +
    xlim(-0.8, 1.2) +
    theme_classic()
}


