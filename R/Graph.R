#' Remove elements from a conditioning set by using conditional independence
#'
#' @param DAG Directed Acyclic Graph
#' @inheritParams create_empty_DAG
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
#' DAG = set_arc(DAG, 'U1', 'U3')
#' DAG = set_arc(DAG, 'U2', 'U3')
#'
#' remove_CondInd(DAG = DAG, node = "U1", cond_set = c("U2"))
#' remove_CondInd(DAG = DAG, node = "U3", cond_set = c("U1"))
#'
#' @export
#'
remove_CondInd <- function(DAG, node, cond_set, graph_backend = NULL){
  DAG = .prepare_graph(DAG, graph_backend)
  return (.remove_CondInd(DAG, node, cond_set))
}

.remove_CondInd <- function(DAG, node, cond_set){
  new_cond_set = cond_set
  for (i in cond_set){
    # check if dsep(node, i | cond_set\{i} )
    if (.graph_dsep(DAG, node,i,setdiff(cond_set,i))){
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
#' @param N_nodes A positive integer equal to the number of nodes
#' @param graph_backend Graph backend, either \code{"bnlearn"} or
#' \code{"caugi"}. For functions accepting a graph, \code{NULL} (the
#' default) uses its backend. An explicit backend converts the input graph.
#' Creating a new graph defaults to \code{"bnlearn"}. The caugi backend
#' requires the optional package \pkg{caugi} (version 1.3.0 or later).
#'
#' @returns A native graph for the selected backend, with `N_nodes` nodes and no arcs.
#'
#' @seealso \code{bnlearn::empty.graph()} and \code{caugi::caugi()}.
#'
#' @examples
#' create_empty_DAG(6)
#' create_empty_DAG(10)
#' if (requireNamespace("caugi", quietly = TRUE)) {
#'   DAG = create_empty_DAG(3, graph_backend = "caugi")
#'   DAG = set_arc(DAG, "U1", "U3")
#'   dsep_set(DAG, "U1", "U2")
#' }
#'
#' @export
create_empty_DAG <- function(N_nodes, graph_backend = "bnlearn"){
  graph_backend = .graph_backend(graph_backend = graph_backend)
  if (!is.numeric(N_nodes) || length(N_nodes) != 1L || is.na(N_nodes) ||
      !is.finite(N_nodes) || N_nodes < 1 || N_nodes != floor(N_nodes)){
    stop("'N_nodes' must be a positive integer.", call. = FALSE)
  }
  node.names = paste0("U", seq_len(N_nodes))
  if (graph_backend == "caugi"){
    .require_caugi()
    return(caugi::caugi(nodes = node.names, class = "DAG"))
  }
  return(bnlearn::empty.graph(node.names))
}


#' D-separation of two nodes given a set in a DAG
#'
#' @param DAG Directed Acyclic Graph
#' @inheritParams create_empty_DAG
#' @param X node
#' @param Y node
#' @param Z set
#'
#' @return \code{TRUE} if the sets are d-separated and \code{FALSE} if not
#'
#' @examples
#'
#' DAG = create_empty_DAG(5)
#' DAG = set_arc(DAG, 'U1', 'U3')
#' DAG = set_arc(DAG, 'U2', 'U3')
#' DAG = set_arc(DAG, 'U1', 'U4')
#' DAG = set_arc(DAG, 'U2', 'U4')
#' DAG = set_arc(DAG, 'U3', 'U4')
#'
#' dsep_set(DAG, 'U1', 'U5')
#'
#' @export
#'
dsep_set <- function(DAG, X, Y, Z = NULL, graph_backend = NULL){
  DAG = .prepare_graph(DAG, graph_backend)
  return (.dsep_set(DAG, X, Y, Z))
}

.dsep_set <- function(DAG, X, Y, Z = NULL){
  if (length(Z) == 0){
    Z = character(0)
  }
  for (x in X){
    for (y in Y){
      if (!.graph_dsep(DAG, x, y, Z)){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}
