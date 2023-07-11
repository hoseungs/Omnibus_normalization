# Omnibus_normalization

This repository provides the _R_ function for the paper: 

**Song, H.**, Ling, W., Zhao, N., Plantinga, A.M., Broedlow, C.A., Klatt, N.R., \
  Hensley-McBain, T., Wu, M.C. (2023). \
  [Accommodating multiple potential normalizations in microbiome associations studies.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05147-w)   \
  **_BMC Bioinformatics, 24(1):22._**

* The main function, **omni**, provides the omnibus p-value by combining popular normalization strategies, such as none, rarefaction, TSS, CSS, and CLR, based on the Cauchy combination test (Liu et al., 2020).

* The main function utilizes association testing methods: linear regression, ZINQ proposed by Ling et al. (2021), and QRank proposed by Song et al. (2017), which are have been shown to consistently protect type I error and do not depend on the particular choice of normalization.


### References
* Liu, Y., Xie, J.: Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures. Journal of the American Statistical Association 115(529), 393–402 (2020)
* Song, X., Li, G., Zhou, Z., Wang, X., Ionita-Laza, I., Wei, Y.: Qrank: a novel quantile regression tool for eqtl discovery. Bioinformatics 33(14), 2123–2130 (2017)
* Ling, W., Zhao, N., Plantinga, A.M., Launer, L.J., Fodor, A.A., Meyer, K.A., Wu, M.C.: Powerful and robust non-parametric association testing for microbiome data via a zero-inflated quantile approach (zinq). Microbiome 9(1), 1–19 (2021)
