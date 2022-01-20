##############
#            #
# R package  #
#            #
##############


# install package

# library("devtools")
# library(roxygen2)

# setwd("~/PhD/Data/GitHub")
# install("gccR")

library(gccR)

##############
#            #
# analysis   #
#            #
##############


##### functional expression classes #####

# corM_ortho$csM1 contain genes that have an ortholog

get_expression_classes(corM_ortho$csM1, height_cutoff = 1050)


