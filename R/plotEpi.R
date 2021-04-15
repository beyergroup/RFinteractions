######################################################################
### plot eQTL map with epistatic connections
### rf=randomForest object (must be created with rfdev package)
### P=matrix of P values
### cutoff=FDR threshold to be used
### main=title of plot
######################################################################
plotEpi = function(rf,P,cutoff=0.01,main=""){
  P[!upper.tri(P)] = NA
  P = as.data.frame(as.table(P))
  P = as.matrix(P[!is.na(P[,3]) & p.adjust(P[,3],'fdr')<cutoff,1:2])
  mode(P) = 'numeric'
  yl = max(rf$importance[,1])
  yl = c(-0.25*yl,1.1*yl)
  plot(rf$importance[,1],type='h',ylim=yl,axes=F,ylab='importance',
       main=main,col=grey(0.2))
  if(nrow(P)>0){
    apply(P,1,function(x){ 
      arc(min(x[1:2]),max(x[1:2]),0.7*yl[1],col='coral1')
    })
  }
  axis(1)
  axis(2,at=pretty(c(0,yl[2]),4))
}
## function to draw the connections
arc = function(x1,x2,depth=-1,col){
  x = c(x1,0.5*(x2-x1)+x1,x2)
  y = c(0,depth,0)
  lines(spline(x,y,x2-x1),col=col)
}