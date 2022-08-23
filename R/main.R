library(GUniFrac) # for Rarefy
library(compositions) # for clr
library(metagenomeSeq) # for css
library(QRank) 
library(ZINQ)

tss = function(X) {
  temp = t( apply(X, 1, function(z){ z / sum(z)}) ) ### TSS sum=1
  return(temp)
}

css = function(X) {
  X = t(X)
  
  m = ncol(X) 
  ntaxa = nrow(X)
  
  colnames(X) = 1:m
  rownames(X) = 1:ntaxa
  
  OTUdata = AnnotatedDataFrame(data.frame(taxa=rownames(X)))
  obj = newMRexperiment(X, featureData=OTUdata)
  
  p = cumNormStatFast(obj)
  obj0 = cumNorm(obj, p = p)
  res = data.frame(MRcounts(obj0, norm=TRUE, log=TRUE))
  
  return(res)
}

tau1 = c(0.1, 0.25, 0.5, 0.75, 0.9)
tau2 = c(0.25, 0.5, 0.75)

# X: covariates. 
# select interesting covariates through i, e.g., i = c(1, 2, 4)
# k: the name(s) of clinical variable(s) of interest in X
# Y: taxa
# select interesting taxa through j, e.g., j = 35
# apply the linear regression, ZINQ, and QRank methods
# with normalization strategies: none, rarefaction, TSS, CSS, CLR
omni = function(X, Y, i, j, k) {
  Y_none = Y
  Y_rare = Rarefy(Y)$otu.tab.rff
  Y_tss = as.matrix(tss(Y))
  Y_css = t(as.matrix(css(Y)))
  Y_clr = as.matrix(clr(Y+0.5))
  
  len = ncol(X[,c(i,k)])
  
  # none
  Dat = data.frame(expression=Y_none[,j], X[,c(i,k)])
  a = lm(expression~., data = Dat)
  lr_none = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len]) 
  
  result = ZINQ_tests(formula.logistic = expression ~.,
                      formula.quantile = expression ~.,
                      C = names(Dat)[len], y_CorD = "D", data = Dat, taus=tau2)
  zinq_none = ZINQ_combination(result, method="Cauchy", taus=tau2)

  temp = Y_none[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_none = QRank(gene=temp, snp=X[,k], cov=X[,i], tau=tau2)$composite.pvalue
  
  # rarefaction
  Dat = data.frame(expression=Y_rare[,j], X[,c(i,k)])
  a = lm(expression~., data = Dat)
  lr_rare = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len]) 
  
  result = ZINQ_tests(formula.logistic = expression ~.,
                      formula.quantile = expression ~.,
                      C = names(Dat)[len], y_CorD = "D", data = Dat, taus=tau2)
  zinq_rare = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_rare[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_rare = QRank(gene=temp, snp=X[,k], cov=X[,i], tau=tau2)$composite.pvalue
  
  # TSS
  Dat = data.frame(expression=Y_tss[,j], X[,c(i,k)])
  a = lm(expression~., data = Dat)
  lr_tss = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len]) 
  
  result = ZINQ_tests(formula.logistic = expression ~.,
                      formula.quantile = expression ~.,
                      C = names(Dat)[len], y_CorD = "D", data = Dat, taus=tau2)
  zinq_tss = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_tss[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_tss = QRank(gene=temp, snp=X[,k], cov=X[,i], tau=tau2)$composite.pvalue
  
  # CSS
  Dat = data.frame(expression=Y_css[,j], X[,c(i,k)])
  a = lm(expression~., data = Dat)
  lr_css = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len]) 
  
  result = ZINQ_tests(formula.logistic = expression ~.,
                      formula.quantile = expression ~.,
                      C = names(Dat)[len], y_CorD = "D", data = Dat, taus=tau2)
  zinq_css = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_css[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_css = QRank(gene=temp, snp=X[,k], cov=X[,i], tau=tau2)$composite.pvalue
  
  # CLR
  Dat = data.frame(expression=Y_clr[,j], X[,c(i,k)])
  a = lm(expression~., data = Dat)
  lr_clr = as.numeric(summary(a)$coefficients[,"Pr(>|t|)"][len]) 
  
  result = ZINQ_tests(formula.logistic = expression ~.,
                      formula.quantile = expression ~.,
                      C = names(Dat)[len], y_CorD = "D", data = Dat, taus=tau2)
  zinq_clr = ZINQ_combination(result, method="Cauchy", taus=tau2)
  
  temp = Y_clr[,j]
  id = length(which(temp == 0)) 
  temp[temp == 0] =  temp[temp == 0] + rnorm(id, mean=0, sd=0.000000001)
  qr_clr = QRank(gene=temp, snp=X[,k], cov=X[,i], tau=tau2)$composite.pvalue
  
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
  
  res = c(lr_omni, zinq_omni, qr_omni)
  
  return(res)
}






