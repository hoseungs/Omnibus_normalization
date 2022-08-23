# Omnibus_normalization

This repository provides the R function for the paper 

Accommodating Multiple Potential Normalizations in Microbiome Associations Studies (under review). \
Song, H., Ling, W., Zhao, N., Plantinga, A.M., Broedlow, C.A., Klatt, N.R., Hensley-McBain, T., Wu, M.C.

The main function, _omni_, provides the omnibus p-value by combining popular normalization strategies, \
such as none, rarefaction, TSS, CSS, and CLR, based on the Cauchy combination test (Liu et al., 2020).

The main function utilizes association testing methods: linear regression, ZINQ proposed by Ling et al. (2021), and QRank proposed by Song et al. (2017),\
which are have been shown to consistently protect type I error and do not depend on the particular choice of normalization.
