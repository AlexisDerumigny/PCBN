

#' Compute the distance between two DAGs or two PCBN structures
#'
#' @param DAG1,DAG2 two DAGS with the same set of nodes
#' @noRd
distance_DAGs = function(DAG1, DAG2, typeDistance)
{
  node_set = bnlearn::nodes(DAG1)

  switch(
    typeDistance,

    "distance_arcs" = {

      distance = sum(abs(DAG1$arcs - DAG2$arcs))
    },

    "distance_vstructs" = {
      distance = 0
      for (v in node_set){
        parents1 = DAG1$nodes[[v]]$parents
        parents2 = DAG2$nodes[[v]]$parents

        distance = distance + length(sym_diff(parents1, parents2))
      }
    }
  )

  return(distance)
}

# symmetric difference function
sym_diff <- function(a,b) setdiff(union(a,b), intersect(a,b))


# Takes in two PCBNs and returns TRUE if they have the same ordering
have_2PCBN_theSameOrdering <- function(PCBN1, PCBN2) {
  order1 = PCBN1$order_hash
  order2 = PCBN2$order_hash
  DAG1 = PCBN1$DAG
  DAG2 = PCBN2$DAG

  # Assume they have the same nodes
  node.names = bnlearn::nodes(DAG1)
  # TODO: Graphs should have distance equal to 0

  # Loop over all v-structures
  for (v in node.names){
    if (length(DAG1$nodes[[v]]$parents) > 1){
      # The orders should be the same everywhere
      if (!min(order1[[v]] == order2[[v]])) {
        return(FALSE)
      }
    }

    if (length(DAG1$nodes[[v]]$parents) > 1){
      # The orders should be the same everywhere
      if (!min(order1[[v]] == order2[[v]])) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}


#
#
# ### Computes the Kullback-Leibler divergence of a estimated PCBN to the true PCBN
# KL_divergence_PCBN <- function(data, PCBN, margins, PCBN_fit, margins_fit) {
#   KL = 0
#   L = length(data[[1]])
#   KL_vec = rep(0, L)
#   data_uniform = to_uniform_scale(data)
#
#   log_lik_true = logLik_PCBN(data, PCBN, margins)
#   log_lik_est = logLik_PCBN(data, PCBN_fit, margins_fit)
#   KL = (log_lik_true - log_lik_est)/L
#   return(KL)
# }
#
# ### Computes the Kullback-Leibler divergence of a estimated GBN to the true PCBN
# KL_divergence_GBN <- function(data, PCBN, margins, GBN) {
#   KL = 0
#   L = length(data[[1]])
#   KL_vec = rep(0, L)
#   data_uniform = to_uniform_scale(data)
#
#   log_lik_true = logLik_PCBN(data, PCBN, margins)
#   log_lik_est = stats::logLik(GBN, data)
#   KL = (log_lik_true - log_lik_est)/L
#   return(KL)
# }
#
#

