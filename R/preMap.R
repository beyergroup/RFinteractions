#' @title Prepare mapping data
#' 
#' @description
#' Genotype and phenotype data are pre-processed and reformatted
#'  for QTL-mapping with random forest. 
#'
#' @details
#' The genotype matrix is checked and markers with too many missing
#' values are removed. Positions of remaining missing values and 
#' minor allele frequencies (MAF) are computed for later imputation
#' of missing genotypes. Highly correlated markers are collapsed to 
#' marker groups. The groupings are also stored and returned.
#' Population structure covariates, representing \code{propVar}%
#' of genetic variance are generated using emma. Phenotype values
#' are optionally scaled and centered. Genotype and phenotype matrices
#' are reordered according to the same sample order.
#' 
#' @param genotype A binary matrix with strain (rows) X marker (columns)
#'  entries. Do not include markers with MAF=1.
#' @param maxNAs Maximum number of missing values that are allowed for
#'  any marker. If this limit is exceeded, the marker is excluded 
#'  from the genotype matrix.
#' @param propVar Proportion of genotypic variance that is included
#'  in the population-structure covariates.
#' @param phenotype A vector or a matrix with numerical values. In 
#'  case of a matrix, each row should correspond to a trait, and each
#'  column to a strain. The matrix may contain NAs. 
#' @param sampleInfo A vector with integers indicating to which strain
#'  (row in the genotype) the measurement(s) in phenotype correspond.
#' @param scale Logical; Should centering and scaling of the phenotype be performed?
#'
#' @return List object containing the genotype, the groupings of the 
#'  markers, the mapping-covariates and the phenotype.
#' @export
#'
#' @examples
preMap <- function(genotype, maxNAs=floor(nrow(genotype)*0.5), propVar=0.75, phenotype, sampleInfo, scale=T){
  
  ###check input
  if(any(!is.matrix(genotype),!is.integer(genotype),genotype!=0&genotype!=1,na.rm = T)){
    stop("Genotype has to be a matrix containing integers (0/1)")
  }
  if(length(unique(sampleInfo))<nrow(genotype)){
    stop("Don't include genotype-data for strains for which no trait data is supplied. Details can be found in the manual.")
  }
  if(!is.integer(sampleInfo)){
    stop("sampleInfo has to be an integer vector.")
  }
  if(any(apply(genotype,2,FUN=function(m){length(unique(m))==1}))){
    stop("Don't include markers with MAF=1.")
  }
  if(any(!is.numeric(maxNAs),maxNAs>length(sampleInfo),maxNAs<0)){
    stop("maxNAs has to be a single integer>=0 and <= to the number of samples")
  }
  if(!is.numeric(propVar)){
    stop("propVar has to be numeric.")
  }
  if(propVar>=1|propVar<0){
    stop("propVar has to be at least zero and has to be smaller than one.")
  }
  if(any(!is.numeric(phenotype),(!is.vector(phenotype))&(!is.matrix(phenotype)))){
    stop("phenotype has to be numeric and either a vector or a matrix")
  }
  if(length(sampleInfo)!=ifelse(is.matrix(phenotype),ncol(phenotype),length(phenotype))){
    stop("sampleInfo has to be as long as there are samples")
  }
  if(ifelse(is.vector(phenotype),any(is.na(phenotype)),any(apply(is.na(phenotype),2,all)))){
    stop("Don't include samples without any trait data.")
  }
  
  
  ###reduce genotype while retaining information
  
  #find markers with too many NAs
  nNA <- apply(genotype,2,FUN=function(x){sum(is.na(x))})
  
  #group redundant markers into genotype groups, leave out markers with too many NAs
  genotype2group <- rep(NA,ncol(genotype))
  reversedGenotype <- abs(genotype-1)
  
  redGeno <- NULL
  nGroups <- 0
  
  for(i in 1:ncol(genotype)){
    if(nNA[i]>maxNAs){
      next
    }
    if(nGroups>0){
      id2G <- apply(redGeno,2,FUN=function(x){
        identical(x,genotype[,i])|identical(x,reversedGenotype[,i])
      })
      if(any(id2G)){
        genotype2group[i] <- which(id2G)
      }else{
        redGeno <- cbind(redGeno,genotype[,i])
        nGroups <- nGroups+1
        genotype2group[i] <- nGroups
      }
    }else{
      redGeno <- cbind(redGeno,genotype[,i])
      genotype2group <- 1
      nGroups <- 1
    }
  }
  group2genotype <- lapply(1:nGroups,FUN=function(x){which(genotype2group==x)})
  
  ###estimate the population-structure using the reduced genotype
  K<-emma.kinship(t(redGeno),use="pairwise.complete.obs")
  KPCA <- eigen(K)
  useEigen <- 1:which(cumsum(KPCA[[1]])>=sum(KPCA[[1]])*propVar)[1]
  populationStructure <- KPCA[[2]][,useEigen,drop=F]
  colnames(populationStructure) <- paste0("popStr",1:ncol(populationStructure))
  
  ###parse the genotype and the population-structure covariates according to sampleInfo
  mappingGenotype <- redGeno[sampleInfo,,drop=F]
  populationStructure <- populationStructure[sampleInfo,,drop=F]
  
  ###prepare a list with missing genotypes for imputation
  NApos <- extractNAs(mappingGenotype)
  
  ###scale the phenotype if needed
  if(scale){
    if(is.vector(phenotype)){
      phenotype <- scale(phenotype,center=T,scale=T)[,1]
    }else{
      phenotype <- t(apply(phenotype,1,FUN=function(x){
        out <- scale(x,center=T,scale=T)[,1]
        return(out)
      }))
    }
  }
  
  ###return output
  return(list(genotype=mappingGenotype,
              genotype2group=genotype2group,
              group2genotype=group2genotype,
              mappingCovariates=populationStructure,
              NAlist=NApos,
              phenotype=phenotype))
}

#' @title Allele frequencies and missing genotypes.
#' 
#' @description
#' Function that calculates the allele frequencies (1 or 0) for 
#' each marker (column) and identifies missing values.
#'
#' @param genotype Genotype matrix with one row per individual, one
#'  column per marker. Entries should be binary marker information (1,0,NA).
#'
#' @return List with one vector for each marker. The first element 
#'  is the frequency of allele 1. The following numbers indicate the
#'  indices of any missing values for the respective samples.
#' @export
#'
#' @examples
extractNAs <- function(genotype){
  
  apply(genotype,2,FUN=function(predictor){
    out <- length(which(predictor==1))/length(which(!is.na(predictor)))
    if(any(is.na(predictor))){
      out <- c(out,which(is.na(predictor)))
    }
    return(out)
  })
}

#' @title Impute missing genotypes
#' 
#' @description
#' Function that replaces missing genotypes according to the
#' local allele frequencies.
#'
#' @param genotype Genotype as a matrix with 1,0,NA.
#' @param NAlist List object with one vector per predictor, containing
#'  the frequency of the 1-allele, followed by the indices of 
#'  individuals with missing values. As created by preMap.
#'
#' @return Matrix with the same dimensions as the input-genotype.
#'  NAs are replaced with alleles (0,1) according to the allele
#'  frequencies in the input (NAlist).
#' @export
#'
#' @examples
replaceGenoNAs <- function(genotype,NAlist){
  sapply(1:ncol(genotype),FUN=function(x){ #per predictor
    lNA <- length(NAlist[[x]])-1
    if(lNA==0){
      return(genotype[,x]) #if no values are missing, return the original vector
    }
    genotype[NAlist[[x]][2:(lNA+1)],x] <- sample(c(1,0),size=lNA,prob = c(NAlist[[x]][1],1-NAlist[[x]][1]),replace=T) #sample NAs according to NAlist
    return(genotype[,x]) #return the extended genotype
  })
}