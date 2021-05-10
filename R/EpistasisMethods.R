# this is a deprecated function.
# it is taken from previous versions of RFepistasis, here
# for history/documentation reasons only.
RFepistasis <- function(mappingData,
                        markerInds1 = NULL, markerInds2 = NULL, pairMat = NULL,
                        mtry = NULL, ntree = 30000, npermut = 100,
                        nodesize = 5, minTest = 5, nthreads = 1,fullTable = F,
                        Rcutoff = 0.9){
  RFE = require(RandomForestExtended)
  if(!RFE){
    stop("RFepistasis depends on the package RandomForestExtended.
           It can be accessed at http://cellnet-sb.cecad.uni-koeln.de/resources/RandomForestExtended.")
  }
  missingG = any(is.na(mappingData$genotype))
  # check input
  if(any(!c("genotype","genotype2group","group2genotype","mappingCovariates",
            "NAlist","phenotype")%in%names(mappingData))){
    stop("mappingData has to be generated with preMap.")
  }
  if(!is.null(mtry) & !is.numeric(mtry)){
    stop("mtry has to be either NULL or a positive integer.")
  }
  if(is.numeric(mtry)){
    if(mtry < 0){
      stop("mtry has to be either NULL or a positive integer.")
    }
  }
  if(ntree<1 | !is.numeric(ntree)){
    stop("ntree has to be a positive integer.")
  }
  if(ntree<100 & is.numeric(ntree)){
    warning("number of trees is likely too small (<100)")
  }
  if(minTest<1 | !is.numeric(minTest)){
    stop("minTest has to be a positive integer.")
  }
  if(nodesize<1 | !is.numeric(nodesize)){
    stop("nodesize has to be a positive integer.")
  }
  if(nthreads<1 | !is.numeric(nthreads)){
    stop("nthreads has to be a positive integer.")
  }
  if(!is.logical(fullTable)){
    stop("fullTable has to be either TRUE or FALSE")
  }
  if(missingG){
    if(npermut < 1 | !is.numeric(npermut))
      stop("there are missing genotypes, but npermut is not a positive integer.")
    if(ntree %% npermut != 0)
      warning("there are missing genotypes, but ntree is not a multiple of npermut:\n final RF size may differ slightly from ntree")
  }
  if(Rcutoff<0 | Rcutoff>1 | !is.numeric(Rcutoff)){
    stop("Rcutoff has to be a number between 0 and 1.")
  }
  if(!is.null(markerInds1) & !is.integer(markerInds1)){
    stop("markerInds1 has to be either NULL or an integer vector")
  }
  if(!is.null(markerInds2) & !is.integer(markerInds2)){
    stop("markerInds2 has to be either NULL or an integer vector")
  }
  if(!is.null(pairMat) & !is.matrix(pairMat)){
    stop("pairMat has to be either NULL or an integer matrix with two columns")
  }
  if(is.null(markerInds1) & is.null(markerInds2) & is.null(pairMat)){
    message("no markerInds or pairMat submitted --> using all markers")
  }
  if(is.matrix(pairMat)){
    if(!is.integer(pairMat) & !ncol(pairMat) == 2){
      stop("the pairMat matrix has to be an integer matrix with two columns")
    }
  }

  # initialize some variables
  genotype = mappingData$genotype
  phenotypes = mappingData$phenotype
  if(is.vector(phenotypes))
    phenotypes = matrix(phenotypes, nrow = 1)
  popStr = mappingData$mappingCovariates
  NAlist = mappingData$NAlist

  if(is.null(pairMat)){
    if(is.null(markerInds1)){
      markerInds1 = 1:ncol(genotype)
    }
    if(is.null(markerInds2)){
      markerInds2 = 1:ncol(genotype)
    }
    #make a matrix with unique marker combinations
    pairMat <- expand.grid(markerInds1, markerInds2)
    pairMat <- t(apply(pairMat,1,function(x){ c(min(x), max(x)) }))
    pairMat <- pairMat[(pairMat[,1] != pairMat[,2]),]
    pairMat <- unique(pairMat)
  }else{
    markerInds1 = unique(pairMat[,1])
    markerInds2 = unique(pairMat[,2])
  }
  interestingMarkers <- union(markerInds1,markerInds2)
  interestingMarkers <- sort(interestingMarkers)
  cors <- cor(genotype[,interestingMarkers], use = "pairwise.complete.obs")
  ngenotype <- ncol(genotype)
  if(is.null(mtry)){
    mtry = floor(ncol(genotype)/3)
  }
  mHash <- rep(NA,ngenotype) # create a Hash in order to be able to re-map the original indices of the alleles
  ind <- 1
  for(i in 1:ngenotype){
    if(i %in% interestingMarkers){
      mHash[i] <- ind
      ind <- ind+1
    }
  }
  # rowHash <- matrix(NA, ncol = ngenotype, nrow = ngenotype)
  # for(i in 1:nrow(pairMat)){
  #    rowHash[pairMat[i,1],pairMat[i,2]] <- i
  #    rowHash[pairMat[i,2],pairMat[i,1]] <- i
  # }

  # iterate through the traits and call interactions
  phenoRes = lapply(seq_len(nrow(phenotypes)),function(p){
    phenotype = phenotypes[p,]
    # remove samples with missing phenotypes
    missingP <- !is.na(phenotype)
    phenotype <- phenotype[missingP]

    # grow RF
    if(missingG){
      #impute missing genotypes and grow the forest
      rf = lapply(1:npermut, FUN = function(x){
        imputed = replaceGenoNAs(genotype, NAlist)
        imputed = imputed[missingP,]
        rf <- RandomForestExtended::randomForest(x=imputed,y=phenotype,
                                                 ntree=ceiling(ntree/npermut),mtry=mtry,
                                                 nodesize=nodesize,keep.forest=T,
                                                 importance=F, nthreads=nthreads)
      })
      rf = combineForestsExtended(rf)
    } else{
      mappingGenotype = genotype[missingP,]
      rf <- RandomForestExtended::randomForest(x=mappingGenotype,y=phenotype,
                                               ntree=ntree,mtry=mtry,
                                               nodesize=nodesize,keep.forest=T,
                                               importance=F, nthreads=nthreads)
    }
    ntree = rf$ntree

    #organize the forest-matrices into an array
    forestArray <- array(0,dim = c(dim(rf$forest$nodestatus),5))
    forestArray[,,1] <- rf$forest$leftDaughter # rows: node, cols: tree
    forestArray[,,2] <- rf$forest$rightDaughter # rows: node, cols: tree
    forestArray[,,3] <- rf$forest$nodestatus # whether it is a terminal node or not (-3 means not terminal, -1 terminal, and 0 not used)
    forestArray[,,4] <- rf$forest$bestvar # index of the alleles that were used for the split at each node (and in each tree)
    forestArray[,,5] <- rf$forest$oobpred
    rm(rf)

    #create an object to store the distributions of side specific differences in means
    #maybe use sparse matrix representation here? for example library('Matrix') or library('slam')
    mdim <- length(interestingMarkers)
    left <- array(dim=c(mdim,mdim,3)) # dimensions: sum, sumofsquares, count
    left[,,3] <- 0
    right <- array(dim=c(mdim,mdim,3))
    right[,,3] <- 0
    counts <- array(0,dim=c(mdim,mdim,2))

    # walk through the forest and collect the distribution of side-specific means
    for (i in 1:ntree){
      tree <- forestArray[,i,]
      dim(tree) = c(length(tree)/5, 5)
      tS <- treeStructure(tree)
      ## for the paired frequency test:
      markers <- tree[,4]
      markers <- unique(markers[markers!=0]) # get all the markers that were used in this tree
      markers <- markers[markers %in% interestingMarkers]
      markers <- mHash[markers]
      counts[markers,markers,1] = counts[markers,markers,1] + 1 # increase the counter for the marker pairs that occured together
      counts[-markers,markers,2] = counts[-markers,markers,2] + 1 # increase the counter for the cases where only one of the markers was used (not used in row)

      ## left side:
      #build a matrix: col 1&2: nodes, col3&4: respective markers
      pairs <- cbind(tS$nl, matrix(mHash[tree[tS$nl,4]], ncol = 2,byrow = F))
      pairs = pairs[apply(pairs,1,function(x) !any(is.na(x))),,drop = F] # exclude pairs where one of the markers is not in interestingMarkers
      # for asymmetry tests
      if(nrow(pairs) > 0){
        for(j in 1:nrow(pairs)){
          stored <- left[pairs[j,3],pairs[j,4],3]
          childMeans <- tree[tree[pairs[j,2],1:2],5]
          if(any(childMeans==0)){ # if there are no oob samples for this node, the oob prediction is zero. Here these cases are excluded.
            next
          }
          if(stored >= 1){
            left[pairs[j,3],pairs[j,4],1] <- left[pairs[j,3],pairs[j,4],1] +
              (childMeans[2] - childMeans[1])
            left[pairs[j,3],pairs[j,4],2] <- left[pairs[j,3],pairs[j,4],2] +
              (childMeans[2] - childMeans[1]) * (childMeans[2] - childMeans[1])
          } else {
            left[pairs[j,3],pairs[j,4],1] <- childMeans[2] - childMeans[1]
            left[pairs[j,3],pairs[j,4],2] <- (childMeans[2] - childMeans[1]) * (childMeans[2] - childMeans[1])
          }
          left[pairs[j,3],pairs[j,4],3] <- stored + 1
        }
      }

      ## right side: (same stuff)
      pairs <- cbind(tS$nr, matrix(mHash[tree[tS$nr,4]], ncol = 2, byrow = F))
      pairs = pairs[apply(pairs,1,function(x) !any(is.na(x))),,drop = F] # exclude pairs where one of the markers is not in interestingMarkers
      if(nrow(pairs) > 0){
        for(j in 1:nrow(pairs)){
          stored <- right[pairs[j,3],pairs[j,4],3]
          childMeans <- tree[tree[pairs[j,2],1:2],5]
          if(any(childMeans == 0)){
            next
          }
          if(stored >= 1){
            right[pairs[j,3],pairs[j,4],1] <- right[pairs[j,3],pairs[j,4],1] +
              (childMeans[2] - childMeans[1])
            right[pairs[j,3],pairs[j,4],2] <- right[pairs[j,3],pairs[j,4],2] +
              (childMeans[2] - childMeans[1]) * (childMeans[2] - childMeans[1])
          } else{
            right[pairs[j,3],pairs[j,4],1] <- (childMeans[2] - childMeans[1])
            right[pairs[j,3],pairs[j,4],2] <- (childMeans[2] - childMeans[1]) * (childMeans[2] - childMeans[1])
          }
          right[pairs[j,3],pairs[j,4],3] <- stored+1
        }
      }

    }
    # get indices of marker pairs which have occured more than minTest times in the tree
    toTest <- which((left[,,3] >= minTest) | (right[,,3] >= minTest), arr.ind = T)
    temp = matrix(interestingMarkers[toTest],ncol = 2)
    # toTest = toTest[!is.na(rowHash[temp]),]

    nTrialsLeft <- rowSums(left[,,3]) # this is needed for the binom test
    nTrialsRight <- rowSums(right[,,3]) # for each marker, how many splits were there on the left/right side in total

    if(nrow(toTest) < 1){
      warning("no tests possible for this phenotype")
      return(NULL)
    } else {  # perform ttest, frequency test, and binomial test
      toTest <- cbind(toTest, NA, NA, NA, NA, NA, NA)
      toTest[,3:8] <- t(apply(toTest[,1:2], 1, FUN = function(row){
        if (abs(cors[row[1], row[2]]) >= Rcutoff) { # exclude markers in LD
          return(c(1, 0, 1, 0, 0, 1))
        }

        # paired frequency test:
        freqmat <- c(counts[row[1],row[2],1],
                     counts[row[1],row[2],2],
                     counts[row[2],row[1],2],
                     0)
        freqmat[4] <- ntree - sum(freqmat)
        dim(freqmat) <- c(2, 2)
        freq_res <- fisher.test(freqmat, alternative = "greater")$p.value

        # prepare the counts for sel. asymmetry:
        nleft <- left[row[1],row[2],3]
        nright <- right[row[1],row[2],3]

        # binomial test:
        ntrialsL <- nTrialsLeft[row[1]]
        ntrialsR <- nTrialsRight[row[1]]

        if(ntrialsL == 0 | ntrialsR == 0){
          return(c(NA, NA, NA, NA, NA, freq_res))
        }

        # t-test:
        if((nright >= minTest) & (nleft >= minTest)){ # we only want to compute the ttest if we have at least 5 slopes on the right and left side.
          n1 = left[row[1],row[2],3]
          m1 = left[row[1],row[2],1] / n1
          v1 = (left[row[1],row[2],2] / n1) - (m1 * m1)
          n2 = right[row[1],row[2],3]
          m2 = right[row[1],row[2],1] / n2
          v2 = (right[row[1],row[2],2] / n2) - (m2 * m2)
          se = sqrt((((n1 - 1) * v1) + ((n2 - 1) * v2)) / (n1 + n2 - 2))
          t  = sqrt(n1 * n2 / (n1 + n2)) * ((m1 - m2) / se)
          df <- n1 + n2 - 2
          tt_res = c(2 * pt(-abs(t), df = df), m1 - m2)
        } else
          tt_res <- c(NA,NA)

        chi_res <- tryCatch(unlist(prop.test(x = c(nleft,nright),
                                             n = c(ntrialsL,ntrialsR),
                                             correct = F,
                                             alternative = "two.sided")[3:4],
                                   use.names = FALSE),
                            warning=function(w) c(NA, NA, NA))

        return(c(tt_res, chi_res, freq_res))
      }))
    }

    # put all results in one table:
    toTest <- toTest[as.logical(rowSums(is.finite(toTest[,3:8]))),] # remove rows where no tests were performed
    toTest[,1] = interestingMarkers[toTest[,1]]
    toTest[,2] = interestingMarkers[toTest[,2]]

    pairMat <- cbind(pairMat, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    colnames(pairMat) <- c("markerA", "markerB",
                           "p_splitA_AbeforeB", "slope_diff_AbeforeB",
                           "p_splitA_BbeforeA", "slope_diff_BbeforeA",
                           "p_selA_AbeforeB", "prop_BleftOfA", "prop_BrightOfA",
                           "p_selA_BbeforeA", "prop_AleftOfB", "prop_ArightOfB",
                           "p_pairedSF_AbeforeB","p_pairedSF_BbeforeA")
    for(i in seq_len(nrow(toTest))){
      mA = toTest[i,1]
      mB = toTest[i,2]
      pairMatInd = which(min(mA,mB) == pairMat[,1] & max(mA,mB) == pairMat[,2])
      if(mA < mB){
        pairMat[pairMatInd,c(3,4,7,8,9,13)] <- toTest[i,3:8]
        # pairMat[rowHash[toTest[i,1],toTest[i,2]],3] <- toTest[i,3]
        # pairMat[rowHash[toTest[i,1],toTest[i,2]],4] <- toTest[i,4]
        # pairMat[rowHash[toTest[i,1],toTest[i,2]],7] <- toTest[i,5]
        # pairMat[rowHash[toTest[i,1],toTest[i,2]],8] <- toTest[i,6]
        # pairMat[rowHash[toTest[i,1],toTest[i,2]],9] <- toTest[i,7]
        # pairMat[rowHash[toTest[i,1],toTest[i,2]],13] <- toTest[i,8]
      }else{
        pairMat[pairMatInd,c(5,6,10,11,12,14)] <- toTest[i,3:8]
        # pairMat[rowHash[toTest[i,2],toTest[i,1]],5] <- toTest[i,3]
        # pairMat[rowHash[toTest[i,2],toTest[i,1]],6] <- toTest[i,4]
        # pairMat[rowHash[toTest[i,2],toTest[i,1]],10] <- toTest[i,5]
        # pairMat[rowHash[toTest[i,2],toTest[i,1]],11] <- toTest[i,6]
        # pairMat[rowHash[toTest[i,2],toTest[i,1]],12] <- toTest[i,7]
        # pairMat[rowHash[toTest[i,2],toTest[i,1]],14] <- toTest[i,8]
      }
    }
    res <- cbind(pairMat,
                 "splitA" = apply(pairMat[,c(3,5)],1,combinePvalues) ,
                 "selA" = apply(pairMat[,c(7,10)],1,combinePvalues),
                 "pairedSF" = apply(X = pairMat[,13:14],MARGIN = 1,FUN = function(x){
                   bool <- is.na(x)
                   ifelse(any(bool), yes = x[!bool], no = x[1])
                 }))
    res <- cbind(res,"ensemble" = apply(res[,15:17], 1, combinePvalues))
    if(fullTable)
      return(res)
    else
      return(res[,c(1:2,15:18)])
  })
  return(phenoRes)
}

#' @title Reformat RF tree structure
#'
#' @description
#' Function that reformats the structure of a RF tree for use in RFepistasis.
#'
#' @param treeMat Matrix representing the structure of a tree of a forest.
#'  Each row represents one node. The first two columns are the indices of
#'  the rows of the left and right child nodes, respectively.
#'  The third column indicates whether it is a terminal node
#'  or not (-3 means not terminal, -1 terminal, and 0 not used).
#'  The fourth column is the index in the genotype of the marker
#'  that was used for the split at this node. The last column is
#'  the mean out-of-bag (oob) trait value for this node.
#'
#' @return List of length 2. $nl: matrix with two columns containing
#'  indices of predictors. Predictors in the second column were used
#'  somewhere in the tree on the left side of the respective predictor
#'  in clumn a. $nr: same as $nl, but for pairs where predictors in
#'  column two were used somewhere on the right side of the predictors in column 1.
#' @export
#'
#' @examples
treeStructure <- function(treeMat){
  rownames(treeMat)=1:nrow(treeMat)
  treeMat <- treeMat[treeMat[,3] == -3,, drop = F] # get only nodes which have children
  L <- vector("list",nrow(treeMat))
  names(L) <- rownames(treeMat)
  R <- L
  nodes <- rev(rownames(treeMat))
  for(node in nodes){ # go through nodes and store all non-terminal left/right children in L/R
    dirChildren <- as.character.default(treeMat[node,1:2])
    logic <- (dirChildren %in% nodes)
    if(logic[1]){
      L[[node]] <- c(as.integer(dirChildren[1]), L[[dirChildren[1]]], R[[dirChildren[1]]])
    }
    if(logic[2]){
      R[[node]] <- c(as.integer(dirChildren[2]), L[[dirChildren[2]]], R[[dirChildren[2]]])
    }
  }

  ## transform above lists into matrices with valid marker combinations, first column mother, second column dauther node
  nl <- c(unlist(sapply(names(L), FUN = function(x){
    rep(as.numeric(x), length(L[[x]]))
  },USE.NAMES = F)), unlist(L, use.names = FALSE))
  if(is.null(nl))
    nl = matrix(data = NA, nrow = 0, ncol = 2)
  else
    dim(nl) <- c(length(nl)/2, 2)
  nr <- c(unlist(sapply(names(R), FUN = function(x){
    rep(as.numeric(x), length(R[[x]]))
  },USE.NAMES = F)), unlist(R, use.names = FALSE))
  if(is.null(nr))
    nr = matrix(data = NA, nrow = 0, ncol = 2)
  else
    dim(nr) <- c(length(nr)/2, 2)

  return(list(nl = nl, nr = nr))
}
