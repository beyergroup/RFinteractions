## ToDo: ## 
## Final scores are the sum of all scores calculated for each tree. 
## A tree-wise score of zero occurs, if the variable is responsible for the first split at the root (indicates high importance) OR if there is a split but not within a subtree
## OR if there is no split with this variably it the whole tree. No splitting with this variable can be because of its unimportance or because it is just not sampled for the tree (mtry).
## Now how to decide whether a score is low because of importance or because it is so unimportant and therefore hardly ever appears in a tree. 
## Maybe normalize by the number a variable actually appears in trees? Or include something like a penalty score (maxtreesize +1 ?) for variables not in the tree but in mtry.
## Same penalty for variables not in the subtree but in the whole tree??
## 
## Normalize with treesize including terminal nodes or not?
## 
## How to compare the values? In general for the min depth in one tree -> the lower the min. depth the more important 
## Same variables split often close to each other (possible interaction) -> low min. depth of the variable within the subtree of the other variable in many trees
## The sum of these for all trees: 
## Often split close together -> addition of many close to zero values 
## Rarely split together -> addition of many zeros (here, the penalty would be good) or high values if the split is far down in the subtree 
## Subtraction of the diagonal entries will show, in how far the effect size could be due to the marginal effect of the variable 
## 
## 
##
## Compare row wise -> which variable has the largest effect within a subtree; then compare effect sizes
## Compare col wise -> within which subtree do we see the largest effect of variable x
## Look how far a values differs from zero?
## 
## 





#' @title Normalized minimal depths extracted from maximal subtrees  
#'
#' @description
#' The function calculates the normalized minimal depths within all possible maximal subtrees from a random forest object (\code{rf}).
#' Here, rows represent the variables from which the maximal subtrees were extracted and the columns are showing the variables we are looking at within the subtrees. 
#' The diagonal entries are the normalized minimal depths of the variables within the whole tree. 
#' The resulting matrix is the sum of all min. depth matrices calculated for each tree within the forest. 
#' The min. depth values are normalized with the size of according the tree/ subtree.
#' In addition, the marginal effects of a variable (diagonal entry of the matrix) was subtracted from all minimal depths values of the variable for all subtrees. 
#' #' 
#' 
#' 
#' @rf Random forest object generated with the rf function from RandomForestExtended with keep.forest set to TRUE 
#' @return A dataframe with minimal depth values for each subtree-variable pair(rows = subtree and cols = variable)
#'
#' @examples
maxsubtree_minDepth <- function(rf, markerInds = NULL){
  RFE = require(RandomForestExtended)
  if(!RFE){
    stop("maxsubtree_minDepth depends on the package RandomForestExtended.
           It can be accessed at http://cellnet-sb.cecad.uni-koeln.de/resources/RandomForestExtended.")
  }
  if(any(!c("call","type","predicted","mse","rsq","oob.times", "importance", "ntree", "mtry", "importanceSD", "localImportance", "proximity", "coefs", 
            "forest","y", "test", "inbag")%in%names(rf))){
    stop("Random Forest object has to be generated with the rf function from RandomForestExtended with keep.forest set to TRUE.")
  }
  ntree = rf$ntree
  
  if(!is.null(markerInds)){
    marker = names(rf$forest$xlevels)[unique(markerInds)]
  } else {
    marker = names(rf$forest$xlevels)
  }

  # Iterate through all trees
  matrix_minDepth <- lapply(1:ntree, function(t){
    m = matrix(NA, nrow = length(marker), ncol = length(marker), 
               dimnames = list(marker, marker)) # pxp matrix where p is the numbers of predictors 
    
    sizeTree = treesize(rf, terminal =  F)[t] 
    tree = getTree(rf, k=t, labelVar = T)
    
    
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
    diag(m) = norm_minDepth_marg #Add normalized values to matrix diagonal
  
   
    
     
    ## Min. depth within a subtree for all variables ##
    subtree_minDepth <- lapply(marker, function(subtreeVar){
      
      subtree = tree[which(tree$`split var` == subtreeVar),][1,] # Takes the maximal subtree
      subtreeRoot = suppressWarnings(as.numeric(rownames(subtree))) 
      
      if(!is.na(subtreeRoot)){
        sizeSubtree = subtree_size(subtreeRoot, tree_num = t) # See function below

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
          }
          
          return(minDepth)
        # min depth im subtee  
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
      m[i,names(entry)] = entry
    }
    
    return(m)
  })

  matrix_minDepth = Reduce("+", matrix_minDepth) # Sum of all min. depth matrices of all rf trees 
  
  # Subtract marginal effect of the marker from the maxsubtree min. depth values
  for(col in 1:ncol(matrix_minDepth)){
    marg_effect = matrix_minDepth[col, col]
    matrix_minDepth[-col, col] = matrix_minDepth[-col, col] - marg_effect
  }
  
  return(as.data.frame(matrix_minDepth))
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
subtree_size = function(subtree_root, tree_num){
  tree = getTree(rf, k=tree_num, labelVar = T)
  
  size1=0
  size2=0
  ld = tree[subtree_root,1]
  
  if(ld == 0){
    size1=0  # if there is no further daughter node -> size does not increase
  } else {
    size1 = subtree_size(ld, tree_num) # if there is a daughter node, calculate its size by checking if it has further daughter nodes 
  }
  
  rd = tree[subtree_root,2]
  
  if(rd == 0){
    size2=0
  } else {
    size2 = subtree_size(rd, tree_num)
  }
  return(size1+size2+1) # add node sizes of daughters + 1 for the parent node 
}
