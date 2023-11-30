#' D-separation of two nodes given a set in a DAG
#'
#' @param DAG Directed Acyclic Graph
#' @param X node
#' @param Y node
#' @param Z set
#'
#' @return TRUE if the sets are d-separated and FALSE if not
#'
dsep_set <- function(DAG, X, Y, Z=c()){
  dseparated = TRUE
  for (x in X){
    for (y in Y){
      if (length(Z)==0){
        if (!bnlearn::dsep(DAG,x,y)){
          dseparated = FALSE
        }
      } else{
        if (!bnlearn::dsep(DAG,x,y,Z)){
          dseparated = FALSE
        }
      }
    }
  }
  return(dseparated)
}


#' Turns an order hash into a matrix
order_hash_to_mat <- function(DAG, order_hash) {
  adj.mat = bnlearn::amat(DAG)
  for (v in bnlearn::nodes(DAG)) {
    order = order_hash[[v]]
    k = 1
    for (w in order) {
      adj.mat[w, v] = k
      k = k + 1
    }
  }
  return(adj.mat)
}

#' Turns a matrix of orders into an order hash
mat_to_order_hash <- function(DAG, order_mat) {
  order_hash = r2r::hashmap()
  for (v in nodes(DAG)) {
    order_hash[[v]] = NULL
    if (sum(order_mat[, v] != 0) > 0) {
      order = c()
      for (k in 1:sum(order_mat[, v] != 0)) {
        w = nodes(DAG)[which(order_mat[, v] == k)]
        order = append(order, w)
      }
      order_hash[[v]] = order
    }
  }
  return(order_hash)
}

#' Makes a copy of a hashmap
#'
#' @param hash hashmap
#'
#' @returns hashmap
copy_hash <- function(hash) {
  new = r2r::hashmap()
  if (length(r2r::keys(hash)) > 0) {
    for (key in r2r::keys(hash)) {
      new[[key]] = hash[[key]]
    }
  }
  return(new)
}

#' Take symmetric difference of two vectors
#'
#' @param a Vector 1
#' @param b Vector 2
#'
#' @returns Symmetric difference of a and b
sym_diff <- function(a,b) setdiff(union(a,b), intersect(a,b))



