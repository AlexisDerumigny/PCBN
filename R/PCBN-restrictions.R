
#' Does a DAG satisfy the restrictions of no active cycle and no
#' interfering v-structures
#'
#' @param DAG the DAG object
#' @param verbose if \code{verbose} is \code{2}, details are printed.
#' If \code{verbose} is \code{1}, details are printed only if an active cyle
#' or an interfering v-structure is found.
#' If \code{verbose} is \code{0} the function does not print and only returns
#' the result.
#' @param check_both if \code{TRUE}, both v-structures and active cycles are
#' checked anyway. If \code{FALSE}, the function stops early if it already found
#' any v-structures.
#'
#' @return a Boolean indicating whether or not the PCBN satisfies the
#' restrictions.
#'
#' @examples
#'
#' DAG = create_empty_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' is_restrictedDAG(DAG)  # 1 active cycle
#'
#' @export
is_restrictedDAG <- function(DAG, verbose = 2, check_both = TRUE)
{
  has_vstructs = has_interfering_vstrucs(DAG)
  if (has_vstructs && verbose >= 1){
    cat("At least one v-structure was found.\n")
  } else if (!has_vstructs && verbose >= 2){
    cat("No v-structures were found.\n")
  }

  if (!check_both && has_vstructs){
    # Early stopping: we know that the conditions are not satisfied
    return (FALSE)
  }

  active_cycle = active_cycles(DAG = DAG, early.stopping = TRUE)
  has_active_cycles = (length(active_cycle) > 0)

  if (has_active_cycles && verbose >= 1){
    cat("At least one active cycle was found.\n")
  } else if (!has_active_cycles && verbose >= 2){
    cat("No active cycle were found.\n")
  }

  is_restricted = !has_vstructs && !has_active_cycles
  return (is_restricted)
}

# This returns TRUE only if the DAG is indeed acyclic
# and if it satisfies the restrictions.
# It has early stopping, so it returns FALSE if it is not a DAG
# and only then does it check the restrictions.
is_DAG_and_restricted <- function(DAG, verbose = 0)
{
  is_DAG = bnlearn::acyclic(DAG)
  if (!is_DAG){
    return (FALSE)
  }

  return (is_restrictedDAG(DAG, verbose = verbose, check_both = FALSE))
}


#' Turns a general graph into a restricted graph.
#'
#' @param DAG Directed Acyclic Graph.
#'
#' @returns Restricted DAG.
#'
#' @examples
#' # DAG with an active cycle at node 5
#' DAG = create_empty_DAG(5)
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
#' DAG = create_empty_DAG(5)
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

  repeat {

    # Remove active cycles
    active_cycles = active_cycles(DAG)
    hasActiveCycles = length(active_cycles) > 0

    if (hasActiveCycles) {
      DAG = fix_active_cycles(DAG, active_cycles)
    }

    # Remove interfering v-structures
    res = find_B_sets(DAG)
    hasInterferingVStructs = res$has_interfering_vstrucs

    if (hasInterferingVStructs) {
      DAG = fix_interfering_vstructs(DAG, res)
    }

    if (!hasActiveCycles && !hasInterferingVStructs){
      break
    }
  }

  return(DAG)
}


#' @rdname DAG_to_restricted
#' @export
fix_active_cycles <- function(DAG, active_cycles) {

  # Point arcs from all nodes to the v-structure
  for (ac in active_cycles) {
    vstruc = ac[[1]]
    rest = ac[which(ac != vstruc)]
    for (node in rest) {
      DAG = bnlearn::set.arc(DAG, node, vstruc)
    }
  }
  return (DAG)
}

#' @rdname DAG_to_restricted
#' @export
fix_interfering_vstructs <- function(DAG, all_B_sets){

  for (v in all_B_sets$nodes_with_inter_vs) {
    B_sets = all_B_sets$B_sets[[v]]
    N_rows_B_sets = dim(B_sets)[1]
    if (N_rows_B_sets <= 3) next

    # B_sets now consists of empty B-set, full B-set, and at least two
    # interfering B-sets. Empty and full B-sets never provide problems,
    # it remains to check the others
    for (i in 2:(N_rows_B_sets - 2)) {
      for (j in (i + 1):(N_rows_B_sets - 1)) {
        # B_sets[i] not in B_sets[j]
        Bset_i = B_sets[i, ]
        Bset_j = B_sets[j, ]
        increasing = all(Bset_i <= Bset_j)

        if (increasing) next

        # Add arcs from all nodes in Bset_i not in Bset_j to the bq of Bset_j
        nodes_in_i_not_in_j = names(which(Bset_i > Bset_j))
        bq_j = rownames(B_sets)[j]

        for (node in nodes_in_i_not_in_j){
          DAG = bnlearn::set.arc(DAG, node, bq_j)
        }
      }
    }
  }
  return (DAG)
}
