#' Remove elements from a conditioning set by using conditional independence
#'
#' @param DAG Directed Acyclic Graph
#' @param node node
#' @param cond_set vector of nodes in conditioning set
#'
#' @returns a vector containing the nodes that cannot be removed from the
#' conditioning set.
#'
#'
#' @examples
#'
#' DAG = create_empty_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' remove_CondInd(DAG = DAG, node = "U1", cond_set = c("U2"))
#' remove_CondInd(DAG = DAG, node = "U3", cond_set = c("U1"))
#'
#' @export
#'
remove_CondInd <- function(DAG, node, cond_set){
  new_cond_set = cond_set
  for (i in cond_set){
    # check if dsep(node, i | cond_set\{i} )
    if (bnlearn::dsep(DAG, node,i,setdiff(cond_set,i))){
      new_cond_set = new_cond_set[-which(new_cond_set == i)]
    }

  }
  return(new_cond_set)
}

#' Create empty DAG
#'
#' This function creates a directed graph with a total of `N_nodes` nodes and
#' no arcs. The nodes are named `U1`, `U2`, etc.
#'
#' @param N_nodes An integer equal to the number of nodes
#'
#' @returns A \code{bnlearn} graph object with `N_nodes` nodes and no arcs
#'
#' @seealso [bnlearn::empty.graph()] which this function wraps.
#'
#' @examples
#' create_empty_DAG(6)
#' create_empty_DAG(10)
#'
#' @export
create_empty_DAG <- function(N_nodes){
  node.names = c()
  for (i in 1:N_nodes){
    node.names[i] = paste("U",as.character(i),sep="")
  }
  DAG = bnlearn::empty.graph(node.names)
  adj = matrix(0L, ncol = N_nodes, nrow = N_nodes,
               dimnames = list(node.names, node.names))
  bnlearn::amat(DAG) = adj
  return(DAG)
}


#' D-separation of two nodes given a set in a DAG
#'
#' @param DAG Directed Acyclic Graph
#' @param X node
#' @param Y node
#' @param Z set
#'
#' @return \code{TRUE} if the sets are d-separated and \code{FALSE} if not
#'
#' @examples
#'
#' DAG = create_empty_DAG(5)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' dsep_set(DAG, 'U1', 'U5')
#'
#' @export
#'
dsep_set <- function(DAG, X, Y, Z = NULL){
  if (length(Z) == 0){
    Z = character(0)
  }
  for (x in X){
    for (y in Y){
      if (!bnlearn::dsep(DAG, x, y, Z)){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}
