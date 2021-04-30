#' @title create a png and a pdf of a plot with one function
#'
#' @description
#' A convenience function that creates a pdf and a png with the same base
#' name for the  plotting command in
#' expr. Works for base R plotting as well as ggplot.
#'
#' @param expr An expression created with \code{expr()} with the plotting
#' command to be executed.
#' @param path The base file path of the figures that should be created
#' without file ending (".png" and ".pdf" are appended internally)
#' @param width width of the figure in inches.
#' @param heigth height of the figure in inches.
#' @param pointsize the pointsize of plotted text, passed on
#' to \code{png()} and \code{pdf()}
#' @param res The nominal resolution in ppi, passed on to \code{png()}
#' @return output of dev.off(). This function is mainly called for its
#' side effect of creating a pdf as well as a png file.
#'
#' @examples
pdfandpng = function(expr, path, width = 5, height = 5, pointsize = 12, res = 200){
  png(paste0(path, ".png"), width = width, height = height, units = "in",
      pointsize = pointsize, res = res)
  eval(expr)
  dev.off()
  pdf(paste0(path, ".pdf"), width = width, height = height,
      pointsize = pointsize)
  eval(expr)
  dev.off()
}

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
