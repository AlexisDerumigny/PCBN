
# Turns an order hash into a matrix
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

# Turns a matrix of orders into an order hash
mat_to_order_hash <- function(DAG, order_mat) {
  order_hash = r2r::hashmap()
  for (v in bnlearn::nodes(DAG)) {
    order_hash[[v]] = NULL
    if (sum(order_mat[, v] != 0) > 0) {
      order = c()
      for (k in 1:sum(order_mat[, v] != 0)) {
        w = bnlearn::nodes(DAG)[which(order_mat[, v] == k)]
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






####################################################
##### For data with standard normal margins  #######
####################################################

### Fits normal margins and returns a list with the mean and sd for each node
fit_normal_margins <- function(data){
  nodes = colnames(data)
  fit_list = list()
  for (node in nodes){
    fit = MASS::fitdistr(data[[node]],"normal")
    estimates = list( mean = fit$estimate[[1]], sd = fit$estimate[[2]])
    fit_list[[node]] =  estimates
  }
  return(fit_list)
}


### Scale data-frame with uniform margins to standard normal margins
to_normal_scale <- function(data){
  data_new = data
  for (i in 1:length(data)){
    data_new[[i]] = stats::qnorm(data[[i]], 0, 1)
  }
  return(data_new)
}

### Scale data-frame with uniform margins to standard normal margins
to_uniform_scale <- function(data){
  data_new = data
  for (i in 1:length(data)){
    data_new[[i]] = VineCopula::pobs(data[[i]])
  }
  return(data_new)
}






###################################
###### not used? ##################
##################################


# Takes in two PCBNs and returns TRUE if they have the same ordering
check_correct_order <- function(PCBN1, PCBN2) {
  order1 = PCBN1$order_hash
  order2 = PCBN2$order_hash
  DAG1 = PCBN1$DAG
  DAG2 = PCBN2$DAG

  # TODO: Graphs should have distance equal to 0

  # For now assume they do
  node.names = bnlearn::nodes(DAG1)
  # Loop over all v-structures
  for (v in node.names){
    if (length(DAG1$nodes[[v]]$parents) > 1){
      # The orders should be the same everywhere
      if (!min(order1[[v]] == order2[[v]])) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

