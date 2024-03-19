
#' Turns a general graph into a restricted graph.
#'
#' @param DAG Directed Acyclic Graph.
#'
#' @returns Restricted DAG.
#'
#' @examples
#' # DAG with an active cycle at node 5
#' DAG = create_DAG(5)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U4', 'U5')
#'
#' # Fixed graph has extra arcs 1 -> 5, 2 -> 5
#' fixed_DAG = DAG_to_restricted(DAG)
#'
#'
#' # DAG with an interfering v-structures node 3
#' DAG = create_DAG(5)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
#'
#' # Fixed graph has extra arc 1 -> 5
#' fixed_DAG = DAG_to_restricted(DAG)
#'
#' @export
DAG_to_restricted <- function(DAG) {
  DAG_copy = DAG

  # TODO: This should be a while loop, since fixing active cycles and int_vs can
  # create more of them

  # Remove active cycles
  active_cycles = active_cycles(DAG)
  if (length(active_cycles) > 0) {
    # Point arcs from all nodes to the v-structure
    for (ac in active_cycles) {
      vstruc = ac[[1]]
      rest = ac[which(ac != vstruc)]
      for (node in rest) {
        DAG_copy = bnlearn::set.arc(DAG_copy, node, vstruc)
      }
    }
  }

  # Remove interfering v-structures
  res = find_B_sets(DAG)
  if (res$has_interfering_vstrucs) {
    for (node in res$nodes_with_inter_vs) {
      B_sets = res$B_sets[[node]]

      # B_sets consists of empty B-set, full B-set, at least two interfering
      # B-sets
      # Empty and full B-set never provide problems, it remains to check the others
      N_rows_B_sets = dim(B_sets)[1]
      for (i in 2:(N_rows_B_sets - 1)) {
        for (j in (i + 1):(N_rows_B_sets - 1)) {
          # B_sets[i] not in B_sets[j]
          Bset_i = B_sets[i, ]
          Bset_j = B_sets[j, ]
          increasing = all(Bset_i <= Bset_j)

          if (! increasing)
          {
            # Add arcs from all nodes in Bset_i not in Bset_j to the bq of Bset_j
            nodes_in_i_not_in_j = names(which(Bset_i > Bset_j))
            bq_j = rownames(B_sets)[j]
            for (node in nodes_in_i_not_in_j){
              DAG_copy = bnlearn::set.arc(DAG_copy, node, bq_j)
            }
          }
        }
      }
    }
  }

  return(DAG_copy)
}
