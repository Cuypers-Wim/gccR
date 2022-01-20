
#### timing ####

# compute correlation matrices

system.time(corM1 <- get_corM(exprList$exprValues1, dropNArows = TRUE))

# user  system elapsed
# 21.86    0.09   22.11


system.time(corM2 <- get_corM(exprList$exprValues2, dropNArows = TRUE))

# user  system elapsed
# 26.84    0.11   27.07


system.time(corM1f <- get_corM_fast(exprList$exprValues1, dropNArows = TRUE))

# user  system elapsed
# 17.33    0.17   17.52
#

system.time(corM2f <- get_corM_fast(exprList$exprValues2, dropNArows = TRUE))

#  user  system elapsed
# 21.31    0.11   21.50
