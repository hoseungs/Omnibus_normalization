library(GUniFrac) # for Rarefy
library(compositions) # for clr
library(metagenomeSeq) # for css
library(QRank) 
library(ZINQ)
source('util.R')

# This function provides the omnibus result (p-value) 
# based on some popular normalization methods: none, rarefaction, TSS, CSS, CLR,
# by the linear regression, ZINQ, and QRank association testing methods.

# X: covariates
# select interesting covariates through i, e.g., i = c(1, 2, 4)
# k: the index of clinical variable(s) of interest in X, e.g., k = 1 or 2 or 4
# Y: taxa
# select interesting taxa through j, e.g., j = 35

# use the linear regression, ZINQ, and QRank methods
# with normalization strategies: none, rarefaction, TSS, CSS, CLR

omni = function(X, Y, i, j, k) {
  Y_none = Y
  Y_rare = Rarefy(Y)$otu.tab.rff
  Y_tss = as.matrix(tss(Y))
  Y_css = t(as.matrix(css(Y))) 
  Y_clr = as.matrix(clr(Y+0.5))
  
  len = which(i==k)
  xnam <- paste("X", 1:ncol(X[,i]), sep="")
  
  # none
  Dat = data.frame(expression=Y_none[,j], X[,i])
  a = lm(expression~., data = Dat)
  lr_none = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len+1]) 
  
  result = ZINQ_tests(formula.logistic = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      formula.quantile = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      C = names(Dat)[len+1], y_CorD = "D", data = Dat, taus=tau2)
  zinq_none = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_none[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_none = QRank(gene=temp, snp=X[,k], cov=X[,setdiff(i,k)], tau=tau2)$composite.pvalue
  
  # rarefaction
  Dat = data.frame(expression=Y_rare[,j], X[,i])
  a = lm(expression~., data = Dat)
  lr_rare = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len+1]) 
  
  result = ZINQ_tests(formula.logistic = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      formula.quantile = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      C = names(Dat)[len+1], y_CorD = "D", data = Dat, taus=tau2)
  zinq_rare = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_rare[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_rare = QRank(gene=temp, snp=X[,k], cov=X[,setdiff(i,k)], tau=tau2)$composite.pvalue
  
  # TSS
  Dat = data.frame(expression=Y_tss[,j], X[,i])
  a = lm(expression~., data = Dat)
  lr_tss = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len+1]) 
  
  result = ZINQ_tests(formula.logistic = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      formula.quantile = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      C = names(Dat)[len+1], y_CorD = "D", data = Dat, taus=tau2)
  zinq_tss = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_tss[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_tss = QRank(gene=temp, snp=X[,k], cov=X[,setdiff(i,k)], tau=tau2)$composite.pvalue
  
  # CSS
  Dat = data.frame(expression=Y_css[,j], X[,i])
  a = lm(expression~., data = Dat)
  lr_css = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len+1]) 
  
  result = ZINQ_tests(formula.logistic = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      formula.quantile = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      C = names(Dat)[len+1], y_CorD = "D", data = Dat, taus=tau2)
  zinq_css = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_css[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_css = QRank(gene=temp, snp=X[,k], cov=X[,setdiff(i,k)], tau=tau2)$composite.pvalue
  
  # CLR
  Dat = data.frame(expression=Y_clr[,j], X[,i])
  a = lm(expression~., data = Dat)
  lr_clr = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len+1]) 
  
  result = ZINQ_tests(formula.logistic = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      formula.quantile = as.formula(paste("expression ~ ", paste(xnam, collapse= "+"))),
                      C = names(Dat)[len+1], y_CorD = "D", data = Dat, taus=tau2)
  zinq_clr = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_clr[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_clr = QRank(gene=temp, snp=X[,k], cov=X[,setdiff(i,k)], tau=tau2)$composite.pvalue
  
  # omni_lr
  x = c(lr_none, lr_rare, lr_tss, lr_css, lr_clr)
  cauchy.t = sum(tan((0.5-pmin(x,0.99))*pi))/length(x)
  lr_omni = 1 - pcauchy(cauchy.t)
  # omni_zinq
  x = c(zinq_none, zinq_rare, zinq_tss, zinq_css, zinq_clr)
  cauchy.t = sum(tan((0.5-pmin(x,0.99))*pi))/length(x)
  zinq_omni = 1 - pcauchy(cauchy.t)
  # omni_ qr
  x = c(qr_none, qr_rare, qr_tss, qr_css, qr_clr)
  cauchy.t = sum(tan((0.5-pmin(x,0.99))*pi))/length(x)
  qr_omni = 1 - pcauchy(cauchy.t)
  
  res = c('LR' = lr_omni, 'ZINQ' = zinq_omni, 'QRank' = qr_omni)
  
  return(res)
}
