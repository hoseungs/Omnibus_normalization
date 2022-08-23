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