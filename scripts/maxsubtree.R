## ToDo scores werden später über alle Bäume aufsummieret, wenn einer aber nicht vorkommt ist der score desw niedriger nicht weil er wichtiger ist
## Kommt er nicht vor, weil unwichtig oder weil nicht gesampled für den tree ??

maxsubtree_minDepth <- function(rf){
  ntree = rf$ntree
  marker = names(rf$forest$xlevels)

  # Iterate through all trees
  matrix_minDepth <- lapply(1:ntree, function(t){
    m = matrix(NA, nrow = length(marker), ncol = length(marker), 
               dimnames = list(marker, marker)) # pxp matrix where p is the numbers of predictors 
    
    sizeTree = treesize(rf, terminal =  F)[t] #TODO: normalize with treesize including terminal nodes or not?
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
    
    # Subtract marginal effect of the marker from the maxsubtree min. depth values
    for(col in 1:ncol(m)){
      marg_effect = m[col, col]
      m[-col, col] = m[-col, col] - marg_effect
    }
    
    
    return(m)
  })

  matrix_minDepth = Reduce("+", matrix_minDepth) # Sum of all min. depth matrices of all rf trees 
  return(as.data.frame(matrix_minDepth))
}

  
    



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
