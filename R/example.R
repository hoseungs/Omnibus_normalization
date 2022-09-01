# This example refers to the example in an R Package MiRKAT.

source('main.R')

data(throat.tree)
data(throat.otu.tab)
data(throat.meta)

# X: covariates
covar = cbind(throat.meta$Age, as.numeric(throat.meta$Sex == "Male"))

n = nrow(throat.meta)

# Y: taxa
taxa = throat.otu.tab 

# We want to include the first and second covariates in covar: i = c(1,2)
# The name(s) of clinical variable(s) of interest in covar is the second variable: k = 2
# We are interested in the sixth taxon in Y: j = 6
omni(X=covar, Y=taxa, i=c(1,2), j=6, k=2)
