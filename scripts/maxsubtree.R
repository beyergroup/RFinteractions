# Work in progress ...

matrix_minDepth <- function(rf){
  ntree = rf$ntree
  marker = names(rf$forest$xlevels)
  #interactionMatrix = matrix(0, nrow = length(marker), ncol = length(marker), 
   #                          dimnames = list(marker, marker))
  
  mat <- lapply(1:ntree, function(t){
    m = matrix(NA, nrow = length(marker), ncol = length(marker), 
               dimnames = list(marker, marker)) # pxp Matrix where p is the numbers of predictors 
    
    sizeTree = treesize(rf, terminal =  F)[t]
    tree = getTree(rf, k=t, labelVar = T)
    
    
    ## Marginal effects ##
    marker_minDepth_marg <- sapply(marker, function(currentVar){ # Suppresses warnings bc NAs are produced 

      VarSplit = tree[which(tree$`split var` == currentVar),] # Get node where split with var x happens
      VarSplitNode = suppressWarnings(as.numeric(rownames(VarSplit[1,]))) # Number of the node (take the node closest to the root for minimal depth)
      minDepth = 0 # Start minimal depth counter
      
      
      if(is.na(VarSplitNode)){ # no split with marker x in this tree
        minDepth = NA # No min. depth if not in the tree 

                ## ToDo scores werden später über alle Bäume aufsummieret, wenn einer aber nicht vorkommt ist der score desw niedriger nicht weil er wichtiger ist
        ## Kommt er nicht vor, weil unwichtig oder weil nicht gesampled für den tree ??
        
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
          minDepth = minDepth + 1
        }
      }
      
      return(minDepth)
    })
    
    norm_minDepth_marg = marker_minDepth_marg / sizeTree
    
    diag(m) = norm_minDepth_marg
    
    
    ## Min. depth within a subtree ##
    sapply(marker, function(subtreeVar){
      
      subtree = tree[which(tree$`split var` == subtreeVar),][1,] # Take the maximal subtree
      subtreeRoot = suppressWarnings(as.numeric(rownames(subtree))) 
      
      if(!is.na(subtreeRoot)){
        sizeSubtree = subtree_size(subtreeRoot) # See function below
      
      
        minDepth_perSubtree <- sapply(marker[!marker == subtreeVar], function(currentVar){
          VarSplit = tree[which(tree$`split var` == currentVar),] # Get node where split with var x happens
          VarSplitNode = suppressWarnings(as.numeric(rownames(VarSplit[1,]))) # Number of the node (take the node closest to the root for minimal depth)
          minDepth = 0 # Start minimal depth counter
          
          while (VarSplitNode > subtreeRoot & !is.na(VarSplitNode)){ # count till we reach the root node of the subtree or above 
            if ((VarSplitNode %% 2) == 0) { # if node is even, look at left parent node
              parent = tree[which(tree$`left daughter` == VarSplitNode),]
              
            } else { # if odd, look at right parent node 
              parent = tree[which(tree$`right daughter` == VarSplitNode),]
            }
            VarSplitNode = as.numeric(rownames(parent)) # Node number where separated by var x 
            
            minDepth = minDepth + 1
          }
          
          if(VarSplitNode < subtreeRoot | is.na(VarSplitNode)){ # If we couldnt reach the subtree root -> var not in subtree OR variable is not in the tree
            minDepth = NA
          }
          
          return(minDepth)
          
        })

        norm_minDepth_subtree = minDepth_perSubtree / sizeSubtree
      
        m[subtreeVar,marker[!marker == subtreeVar]] =  norm_minDepth_subtree
    
    } else {
      
    }
      }
      
  
  return(mat)  
  })
  
    
}




# Calculates subtree size 
subtree_size = function(subtree_root){
  size1=0
  size2=0
  ld = tree[subtree_root,1]
  
  if(ld == 0){
    size1=1
  } else {
    size1 = subtree_size(ld)
  }
  
  rd = tree[subtree_root,2]
  
  if(rd == 0){
    size2=1
  } else {
    size2 = subtree_size(rd)
  }
  return(size1+size2+1)
}
