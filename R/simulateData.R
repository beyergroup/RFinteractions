#' @title simulate Genotypes
#' 
#' @description
#' A function to simulate (0,1) genotype markers with tunable linkage 
#' disequilibrium
#' 
#' @param x A data frame or a matrix of predictors, or a formula describing 
#' the model to be fitted (for the print method, an randomForest object).
#' If x is a matrix it can contain either double or raw values. 
#' The use of raw values needs much less memory.
#' @param y A response vector. If a factor, classification is assumed, 
#' otherwise regression is assumed. If omitted, randomForest will 
#' run in unsupervised mode. 
#' @param n number of individuals
#' @param nmarkers number of markers
#' @param recombprob probablities of recombination; at each successive locus, 
#' one of these values is sampled randomly to determine how many of the
#'  individuals will "flip" their genotype; the more '0' values in this vector,
#'   the more LD there will be in the markers. Higher values will result 
#'   in more recombination.
#' @return A matrix with (0/1) genotypes.
#' @export
#'
#' @examples
simgeno = function(n,nmarkers=1000,recombprob=c(0,0,0,0,0.01,0.3)){
  geno = matrix(0,n,nmarkers)
  xn = sample(c(0,1),n,replace=T)
  geno[,1] = xn
  draw = ceiling(recombprob*n)
  for(i in 2:nmarkers){
    these = sample(1:n,sample(draw,1))
    xn[these] = !xn[these]
    geno[,i] = xn
  }
  return(geno)
}


#' @title simulate phenotypes
#' 
#' @description
#' A function to simulate traits, given an input logical vector and an effect
#' size for those "with" the trait (i.e. TRUE)
#' 
#' @param logic a logical vector indicating which samples "have" the 
#' trait (TRUE) and which do not (FALSE)
#' @param effect effect size (separation between those that have the trait
#'  and those that don't)
#' @param spread value that controls the spread (SD) of the distributions
#' @return A vector or numerical phenotypes
#' @export
#'
#' @examples
simtrait = function(logic,effect=1,spread=0.1){
  spread = spread*abs(effect)
  logic = as.logical(logic)
  out = numeric(length(logic))
  out[!logic] = rnorm(sum(!logic),0,spread)
  out[logic] = rnorm(sum(logic),effect,spread)
  return(out)
}
