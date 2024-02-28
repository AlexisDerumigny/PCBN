#' Checks if graph has interfering v-structures
#'
#' @param DAG Directed Acyclic Graph
#'
#' @returns TRUE if graph contains interfering vs and FALSE if not
#'
interfering_vstrucs_check <- function(DAG){
  return(find_B_sets(DAG)$has_interfering_vstrucs)
}


#' Find all the B-sets of a given DAG
#'
#' @param DAG A bnlearn graph object
#' @param v node at which we want to find the B-sets
#'
#' @returns \code{find_B_sets} list with three elements \itemize{
#'    \item \code{B_sets} list of B-sets for each node
#'
#'    \item \code{has_interfering_vstrucs} a boolean specifying if the graph
#'    contains interfering v-structures or not
#'
#'    \item \code{nodes_with_inter_vs} a list containing nodes forming the
#'    interfering v-structures.
#' }
#' \code{find_B_sets_v} returns a boolean matrix with \code{(2 + length(children))}
#' columns and \code{length(parents)} rows.
#'
#' @examples
#' DAG = create_DAG(6)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U4', 'U5')
#'
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U6')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U6')
#' DAG = bnlearn::set.arc(DAG, 'U5', 'U6')
#'
#' find_B_sets_v(DAG, v = 'U5')
#' B_sets = find_B_sets(DAG)
#' B_sets$B_sets
#'
#' @export
#'
find_B_sets <- function(DAG)
{
  has_interfering_vstrucs = FALSE
  nodes_with_inter_vs = list()
  node.names = bnlearn::nodes(DAG)
  adj.mat = bnlearn::amat(DAG)
  B_set_list = list()

  for (v in node.names) {
    B_set = find_B_sets_v(DAG = DAG, v = v)
    B_set_list[[v]] = B_set

    # TODO: check code for interfering v-structures

    parents = DAG$nodes[[v]]$parent
    if (length(parents) > 0) {
      if (!B_sets_are_increasing(B_set)) {
        has_interfering_vstrucs = TRUE
        nodes_with_inter_vs = append(nodes_with_inter_vs, v)
      }
    }
  }
  return( list( B_sets = B_set_list,
                has_interfering_vstrucs = has_interfering_vstrucs,
                nodes_with_inter_vs = nodes_with_inter_vs
  ) )
}


#' Find the B sets for a given node v
#'
#' @rdname find_B_sets
#' @export
find_B_sets_v <- function(DAG, v)
{
  parents = DAG$nodes[[v]]$parent
  nparents = length(parents)
  children = DAG$nodes[[v]]$children

  if (nparents == 0) {
    all_B_sets = matrix(nrow = 2 + length(children), ncol = 0)
    return (all_B_sets)
  }

  if (length(children) == 0) {
    all_B_sets = rbind(
      rep(FALSE, nparents) ,
      rep(TRUE, nparents)
    )
  } else if (length(children) == 1) {
    all_B_sets = rbind(
      rep(FALSE, nparents) ,
      parents %in% DAG$nodes[[children]]$parents,
      rep(TRUE, nparents)
    )
  } else { # We know now that length(children) > 1

    # We put all the B-sets in a matrix
    all_B_sets = vapply(X = 1:length(children),
                        FUN.VALUE = rep(TRUE, nparents),
                        FUN = function(i){
      # This returns a vector of booleans of the same size as `parents`
      parents %in% DAG$nodes[[children[i]]]$parents
    })

    # We transpose this because `vapply` makes it a column for each child
    # instead of the desired 1 row for each child
    all_B_sets = t(all_B_sets)

    # We add the trivial B-sets
    all_B_sets = rbind(
      rep(FALSE, nparents) ,
      all_B_sets,
      rep(TRUE, nparents)
    )
  }

  # # We remove duplicates
  # all_B_sets = unique(all_B_sets)

  rownames(all_B_sets) <- c("Empty B-set", children, "Full B-set")
  colnames(all_B_sets) <- parents

  if (length(children) <= 1){
    return (all_B_sets)
  }

  # We sort the B-sets by size
  B_sets_sizes = rowSums(all_B_sets)
  order_B_sets_sizes = order(B_sets_sizes, decreasing = FALSE)
  all_B_sets = all_B_sets[order_B_sets_sizes , ]

  return (all_B_sets)
}


#' Checks if the B-sets for a particular node form an increasing sequence.
#'
#' @param B_sets a boolean matrix with \code{(2 + length(children))}
#' columns and \code{length(parents)} rows.
#' They are assumed to be sorted in increasing order of row sums,
#' i.e. by increasing order of set cardinality.
#' Typically, this will be the output of \code{find_B_sets_v}.
#'
#' @returns TRUE if the list forms an ordered sequence, FALSE if not.
#'
#' @examples
#' B_sets = matrix(c(FALSE, FALSE, FALSE, FALSE,
#'                   TRUE , FALSE, FALSE, FALSE,
#'                   TRUE , TRUE , FALSE, FALSE,
#'                   TRUE , TRUE , TRUE ,  TRUE),
#'                 nrow = 4, byrow = TRUE)
#'
#' B_sets_are_increasing(B_sets)
#'
#' B_sets = matrix(c(FALSE, FALSE, FALSE, FALSE,
#'                   TRUE , FALSE, TRUE , FALSE,
#'                   TRUE , TRUE , FALSE, FALSE,
#'                   TRUE , TRUE , TRUE ,  TRUE),
#'                 nrow = 4, byrow = TRUE)
#'
#' B_sets_are_increasing(B_sets)
#'
#' @export
B_sets_are_increasing <- function(B_sets){
  n_Bsets = nrow(B_sets)

  if (n_Bsets <= 3 || ncol(B_sets) == 0){
    return (TRUE)
  }

  for (i in 2:n_Bsets){
    Bset_i_1 = B_sets[i - 1,]
    Bset_i = B_sets[i,]
    increasing = all(Bset_i_1 <= Bset_i)
    if ( ! increasing ){
      return (FALSE)
    }
  }

  return (TRUE)
}

