
#' Checks if a given (conditional) copula has already been specified
#'
#' @param DAG Directed Acyclic Graph object corresponding to the model
#' @param order_hash hashmap of orders of the parental sets
#' @param w node in DAG
#' @param v node in DAG
#' @param cond vector of nodes in DAG.
#' It is assumed to have been already sorted.
#'
#' @returns \code{TRUE} if the conditional copula \eqn{C_{w, v | cond}}
#' has been specified in the model
#'
#' @examples
#'
#' DAG = create_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#'
#' is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
#'                     w = "U1", v = "U3", cond = c())
#' # returns TRUE because the copula c_{1,3} is known
#'
#' is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
#'                     w = "U2", v = "U3", cond = c())
#' # returns FALSE because the copula c_{2,3} is not known
#'
#' is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
#'                     w = "U2", v = "U3", cond = c("U1"))
#' # returns TRUE because the copula c_{2,3 | 1} is known
#'
#' @export
#'
is_cond_copula_specified <- function(DAG, order_hash, w, v, cond){
  if (dsep_set(DAG, w, v, cond)){
    return(TRUE)
  } else {
    if(length(cond) == 0){ # We transform NULL to character(0)
      cond = character(0)
    }
  }

  parents_v = order_hash[[v]]
  if (w %in% parents_v){
    index_w_in_parents = which(parents_v == w)
    parents_up_to_w = if(index_w_in_parents == 1) { character(0) } else {
      parents_v[1:(index_w_in_parents - 1)] }

    if (identical(sort(parents_up_to_w), cond) ){
      return(TRUE)
    }
  }

  parents_w = order_hash[[w]]
  if (v %in% parents_w){
    index_v_in_parents = which(parents_w == v)
    parents_up_to_v = if(index_v_in_parents == 1) { character(0) } else {
      parents_w[1:(index_v_in_parents - 1)] }

    if (identical(sort(parents_up_to_v), cond)){
      return(TRUE)
    }
  }
  return(FALSE)
}


#' Find among parents of a node, the one that has a conditional copula specified
#'
#' @param DAG Directed Acyclic Graph object corresponding to the model
#' @param order_hash hashmap of orders of the parental sets
#' @param v node in DAG
#' @param cond vector of nodes in DAG. This must not be empty.
#' It is assumed that conditionally independent nodes have already been
#' removed by the function \code{\link{remove_CondInd}}.
#' It is assumed to have been already sorted.
#'
#' @returns a list with \itemize{
#'    \item a node \code{w} such that the conditional copula
#'    \eqn{C_{w, v | cond[-v]}} has been specified in the model.
#'
#'    If no such node can be found, an error message is raised.
#'
#'    \item the set \code{cond[-v]}
#' }
#'
#' @examples
#'
#' DAG = create_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#'
#' find_cond_copula_specified(DAG = DAG, order_hash = order_hash,
#'                            v = "U3", cond = c("U1"))
#' # returns "U1" because the copula c_{1,3} is known
#'
#' find_cond_copula_specified(DAG = DAG, order_hash = order_hash,
#'                            v = "U3", cond = c("U1", "U2"))
#' # returns "U2" because the copula c_{2,3|1} is known
#'
#'
#' @export
#'
find_cond_copula_specified <- function(DAG, order_hash, v, cond)
{
  # Find specified c_{wv|cond_set_minus_w}
  w = NULL
  for (i_w in 1:length(cond)) {
    w_proposed = cond[i_w]
    cond_set_minus_w = cond[-i_w]
    if (is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                 w = w_proposed, v = v,
                                 cond = cond_set_minus_w)) {
      w = w_proposed
      break
    }
  }
  if (is.null(w)) {
    if (length(cond) == 0){
      stop("'find_cond_copula_specified' should not be called ",
           "with an empty conditioning set.")
    }
    stop("no specified conditional copula found.\n",
         "We are at node: ", v, " and the conditioning set is: ",
         dputCharacterVec(cond), "\n",
         "Check that the PCBN satisfies the restrictions ",
         "and that the orders of the parents are all compatible.")
  }

  return (list(w = w, cond_set_minus_w = cond_set_minus_w))
}

