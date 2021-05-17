#' @title Normalized minimal depths extracted from maximal subtrees  
#'
#' @description
#' The function calculates the normalized minimal depths within all possible maximal subtrees from a random forest object (\code{rf}).
#' Here, rows represent the variables from which the maximal subtrees were extracted and the columns are showing the variables we are looking at within the subtrees. 
#' The diagonal entries are the normalized minimal depths of the variables within the whole tree. 
#' The resulting matrix is the sum of all min. depth matrices calculated for each tree within the forest. 
#' The min. depth values are normalized with the size of according the tree/ subtree.
#' In addition, the marginal effects of a variable (diagonal entry of the matrix) were subtracted from all minimal depths values of the variable for all subtrees. 
#'  
#' 
#' 
#' @rf Random forest object with keep.forest set to TRUE 
#' @markerInds Marker indices of interest. All possible combinations will be explored.If NULL, all markers will be selected.   
#' @return Matrix with minimal depth values for each subtree-variable pair (rows = subtree and cols = variable)
#'
#' @examples
maxsubtree_minDepth <- function(rf, markerInds = NULL){
  RFE = require(randomForest, Matrix)
  if(!RFE){
    stop("maxsubtree_minDepth depends on the package RandomForest.")
  }
  if(class(rf) != "randomForest"){
      stop("Input must be a rf object.")
  }
  
  
  ntree = rf$ntree
  
  if(!is.null(markerInds)){
    marker = names(rf$forest$xlevels)[unique(markerInds)]
  } else {
    marker = names(rf$forest$xlevels)
  }
  
  m = matrix(0, nrow = length(marker), ncol = length(marker), 
               dimnames = list(marker, marker)) # pxp matrix where p is the numbers of predictors 

  # Iterate through all trees
  for(t in 1:ntree){
    
    tree = getTree(rf, k=t, labelVar = T)
    sizeTree = nrow(tree[tree$status == 1,])

    
    ## Marginal effects = matrix diagonal ##
    marker_minDepth_marg <- sapply(marker, function(currentVar){ # Suppresses warnings bc NAs are produced 

      VarSplit = tree[which(tree$`split var` == currentVar),] # Get node where split with var x happens
      VarSplitNode = suppressWarnings(as.numeric(rownames(VarSplit[1,]))) # Number of the node (take the node closest to the root for minimal depth)
      minDepth = 0 # Start minimal depth counter
      
      
      if(is.na(VarSplitNode)){ # no split with marker x in this tree
        minDepth = 0 # No min. depth if not in the tree 

      } else if (VarSplitNode == 1 && !is.na(VarSplitNode)){ # first split at the root
        minDepth = 0 # Min. depth is zero

      } else {

        while (VarSplitNode != 1){ # count till we reach the root node 
          
          if ((VarSplitNode %% 2) == 0) { # if node is even, look at left parent node
            parent = tree[which(tree$`left daughter` == VarSplitNode),]
            
          } else { # if odd, look at right parent node 
            parent = tree[which(tree$`right daughter` == VarSplitNode),]
          }
          VarSplitNode = as.numeric(rownames(parent)) # Node number where separated by var x 
          minDepth = minDepth + 1 # Increase count until root is reached
        }
      }
      
      return(minDepth)
    })
    
    norm_minDepth_marg = marker_minDepth_marg / sizeTree # Normalize all depth by tree size 
    
    diag(m) = diag(m) + norm_minDepth_marg
  
   
    
     
    ## Min. depth within a subtree for all variables ##
    subtree_minDepth <- lapply(marker, function(subtreeVar){
      
      subtree = tree[which(tree$`split var` == subtreeVar),][1,] # Takes the root of the maximal subtree
      subtreeRoot = suppressWarnings(as.numeric(rownames(subtree))) 
      
      if(!is.na(subtreeRoot)){
        sizeSubtree = subtree_size(subtreeRoot, tree) # See function below | determines subtree size

        # Calculate minimal depth for each variable within the subtree
        minDepth_perSubtree <- sapply(marker[!marker == subtreeVar], function(currentVar){
          VarSplit = tree[which(tree$`split var` == currentVar),] # Get node where split with var x happens
          VarSplitNode = suppressWarnings(as.numeric(rownames(VarSplit[1,]))) # Number of the node (take the node closest to the root for minimal depth)
          minDepth = 0 # Start minimal depth counter
          
          while (VarSplitNode > subtreeRoot & !is.na(VarSplitNode)){ # count till we reach the root node of the subtree or above 
            if ((VarSplitNode %% 2) == 0) { # if node is even, look at left nodes
              parent = tree[which(tree$`left daughter` == VarSplitNode),]
              
            } else { # if odd, look at right nodes 
              parent = tree[which(tree$`right daughter` == VarSplitNode),]
            }
            VarSplitNode = as.numeric(rownames(parent)) # Node number where separated by var x 
            
            minDepth = minDepth + 1
          }
          
          if(VarSplitNode < subtreeRoot | is.na(VarSplitNode)){ # If we couldnt reach the subtree root -> var not in subtree OR variable is not in the tree
            minDepth = 0
            #### Here we could add the penalty equal to the max subtreesize +1 ####
            #minDepth = sizeSubtree + 1
          }
          
          return(minDepth)
        })
        
        # Add values to matrix 
        norm_minDepth_subtree = minDepth_perSubtree / sizeSubtree
      
        return(norm_minDepth_subtree)
        
      } else if (is.na(subtreeRoot)){ # no split with marker x in this tree -> no subtree possible
        no_subtree <- rep(0, length(marker)-1)
        names(no_subtree) <- marker[!marker == subtreeVar]
        return(no_subtree)
      }
    })
    
    # Add subtree min. depth to matrix 
    for (i in 1:length(subtree_minDepth)) {
      entry = subtree_minDepth[[i]]
      m[i,names(entry)] = m[i,names(entry)] + entry
    }
  }

  marg_effect = diag(m) # sum of all marginal effects
  m = m - marg_effect # remove marginal effects 
  diag(m) = marg_effect
  
  return(m)
}

  
    

#' @title Size of a subtree 
#'
#' @description
#' Calculates the size of a subtree (number of nodes,including the terminal nodes). 
#' The subtree is generated based on the input node number \code{subtree_root} which is going to be the root node of the generated subtree.
#' The number \code{tree_num} indicates the exact tree we are interested in among all trees of the random forest object. 
#' 
#' 
#' @subtree_root Node number of the full tree which should be the root of the subtree of interest 
#' @tree_num Number of the full tree of interest within the given random forest object 
#' @return Numeric size of a subtree
#'
#' @examples
# Calculates subtree size 
subtree_size = function(subtree_root, maintree){

  size1=0
  size2=0
  
  ld = maintree[subtree_root,1]
  
  if(ld != 0){
    size1 = subtree_size(ld, maintree) # if there is a daughter node, calculate its size by checking if it has further daughter nodes 
    counter = 1 # add +1 to size if node is non-terminal 
  } else if (ld == 0){
    counter = 0 # do not increase size if node is terminal (= no split)
  }
  
  rd = maintree[subtree_root,2]
  
  if(rd != 0){
    size2 = subtree_size(rd, maintree)
    counter = 1
  } else if (rd == 0){
    counter = 0
  }
  
  return(size1+size2+counter) # add node sizes of daughters + 1 for the parent node 
}
