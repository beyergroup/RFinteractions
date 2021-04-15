######################################################################
### compile information needed for EPI() from an rfdev randomForest
### object
### rf=randomForest object (must be created with rfdev package)
### testset=logical indicating whether the testset (OOB) predictions
### or the training predictions should be used
######################################################################
EPI = function(x){
  d = abs(x$l-x$r)
  m = (abs(x$l)+abs(x$r))/2
  out = d-m
  out = out
  out = pmin(out,t(out))
  out[out<0] = 0
  return(out)
}
predsymm = function(rf,testset=TRUE){
  l = matrix(0,nrow=nrow(rf$importance), ncol=nrow(rf$importance))
  colnames(l) = rownames(rf$importance)
  rownames(l) = rownames(rf$importance)
  r = l
  ct = l
  xx = rf$forest$bestvar
  rr = rf$forest$rightDaughter
  ll = rf$forest$leftDaughter
  if(testset){
    pred = rf$forest$oobpred
  }else{
    pred = rf$forest$nodepred
  }
  xx[xx==0]=NA
  rr[rr==0]=NA
  ll[ll==0]=NA
  for(i in 1:ncol(xx)){
    if(all(is.na(xx[,i]))) next
    p = pred[ll[,i],i] - pred[rr[,i],i]
    ind = 1:nrow(xx)
    parents = xx[ifelse(ind%%2==0,match(ind,ll[,i]),match(ind,rr[,i])),i]
    ind = cbind(1:nrow(xx),parents,xx[,i],p)
    ind = ind[apply(ind,1,function(x) all(!is.na(x))),]
    if(is.null(nrow(ind))) next
    p = ind[,4]
    side = ind[,1]
    ind = ind[,2:3]
    il = side%%2!=0
    ir = !il
    l[ind[il,]] = l[ind[il,]] + p[il]
    r[ind[ir,]] = r[ind[ir,]] + p[ir]
    ct[ind] = ct[ind] + 1
  }
  return(list(l=l/rf$ntree,r=r/rf$ntree,counts=ct/sum(ct)))
}
