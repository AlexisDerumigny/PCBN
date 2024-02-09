#' Finds all possible copula assignments given a DAG
#'
#' @param DAG Directed Acyclic Graph
#'
#' @returns a list of hashmaps containing the possible orders
find_all_orders <- function(DAG) {
  # Start with empty order
  order_hash = r2r::hashmap()
  all_orders = list(order_hash)

  well_ordering = bnlearn::node.ordering(DAG)
  # Add all possible orders for node v to each order_hash
  for (node in well_ordering) {
    all_orders = extend_orders(DAG, all_orders, node)
  }
  return(all_orders)
}

#' Fills in all possible orders for the next node for each possible order
#'
#' @param DAG Directed Acyclic Graph
#' @param all_orders list of orders
#' @param node node
#'
#' @returns list of order hashmaps
#'
extend_orders <- function(DAG, all_orders, node) {
  extended = list()
  # For each order hashmap we add all possible orders for v
  for (order in all_orders) {
    copy_order = copy_hash(order)
    # If node has no parents do add NULL
    if (length(DAG$nodes[[node]]$parents) == 0) {
      copy_order[[node]] = NULL
      extended[[length(extended) + 1]] = copy_order
    } else{
      # Else add all possible orders for v to the list
      orders_node = find_all_orders_v(DAG, node, order)
      for (node_order in orders_node) {
        copy_order = copy_hash(order)
        copy_order[[node]] = node_order
        extended[[length(extended) + 1]] = copy_order
      }
    }
  }
  return(extended)
}

#' Finds all possible orders of node v given previous copula assignments
#'
#' @param DAG Directed Acyclic Graph
#' @param v node
#' @param order_hash hashmap of orders
#'
#' @returns list of vectors containing all possible orders for v
#'
find_all_orders_v <- function(DAG, v, order_hash) {
  parents = DAG$nodes[[v]]$parents
  B_sets = find_B_sets(DAG)$B_sets
  B_sets_v = B_sets[[v]]

  # Order_list contains the partial orders (starting with empty)
  order_list = list(NULL)
  for (i in 1:length(parents)) {
    new_order_list = list()
    for (order in order_list) {
      B_minus_O = find_B_minus_O(B_sets_v, order)
      # Each possible candidate results in a different order
      for (w in possible_candidates(DAG, v, order, order_hash, B_minus_O)) {
        new_order = append(order, w)
        new_order_list[[length(new_order_list) + 1]] = new_order
      }
    }
    order_list = new_order_list
  }
  return(order_list)
}

#' Finds the smallest B-set of v larger than the current partial order
#'
#' @param B_sets list of B-sets for a particular node
#' @param partial_order order list of parents for particular node
#'
#' @returns Smallest B-set strictly larger than the partial order
#'
find_B_minus_O <- function(B_sets, partial_order){
  for (q in 1:length(B_sets)){
    if (sets::as.set(partial_order)<sets::as.set(B_sets[[q]])){
      return(setdiff(B_sets[[q]], partial_order))
    }
  }
}

#' Gives possible candidates to be added to a partial order
#'
#' @param DAG Directed Acyclic Graph.
#' @param v node in DAG
#' @param order_v partial order for node v.
#' @param order_hash hashmap of parental orders
#' @param B_minus_O B-set setminus partial order
#'
#' @returns vector of possible candidates
#'
possible_candidates <- function(DAG, v, order_v, order_hash, B_minus_O){
  Poss.Cand = c()

  for (w in B_minus_O){
    # Independence
    if (dsep_set(DAG, w, order_v)){
      Poss.Cand = append(Poss.Cand, w)
    } else if (incoming_arc(DAG, w, v, order_v, order_hash)){
      Poss.Cand = append(Poss.Cand, w)
    } else if (outgoing_arc(DAG, w, v, order_v, order_hash)){
      Poss.Cand = append(Poss.Cand, w)
    }
  }
  return(Poss.Cand)
}


#' Checks if w can be added by an incoming arc
#'
#' @param DAG Directed Acyclic Graph.
#' @param w node in DAG
#' @param v node in DAG
#' @param order_v partial order for node v
#' @param order_hash hashmap of parental orders
#'
#' @returns TRUE if w is a possible candidate, FALSE if not
#'
incoming_arc <- function(DAG, w, v, order_v, order_hash){
  adj.mat = bnlearn::amat(DAG)

  for (o in order_v){
    if (adj.mat[w,o]==1){
      order_o = order_hash[[o]]
      if ((which(order_o == w)-1)>0){
        pa_o_up_to_w = order_o[1:which(order_o == w)-1]
      } else{
        pa_o_up_to_w = c()
      }

      if (length(setdiff(order_v, union(o, pa_o_up_to_w)))==0){
        return(TRUE)
      } else if (sets::as.set(pa_o_up_to_w) < sets::as.set(order_o)){
        if (dsep_set(DAG, w, setdiff(order_v, union(o, pa_o_up_to_w)), union(o, pa_o_up_to_w))){
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}

#' Checks if w can be added by an outgoing arc
#'
#' @param DAG Directed Acyclic Graph.
#' @param w node in DAG
#' @param v node in DAG
#' @param order_v partial order for node v
#' @param order_hash hashmap of parental orders
#'
#' @returns TRUE if w is a possible candidate, FALSE if not
#'
outgoing_arc <- function(DAG, w, v, order_v, order_hash){
  adj.mat = bnlearn::amat(DAG)

  for (o in order_v){
    if (adj.mat[o,w]==1){
      order_w = order_hash[[w]]
      if ((which(order_w == o)-1)>0){
        pa_w_up_to_o = order_w[1:which(order_w == o)-1]
      } else{
        pa_w_up_to_o = c()
      }

      if (length(setdiff(order_v, union(o, pa_w_up_to_o)))==0){
        return(TRUE)
      } else if (sets::as.set(pa_w_up_to_o) < sets::as.set(order_w)){
        if (dsep_set(DAG, w,setdiff(order_v, union(o, pa_w_up_to_o)), union(o, pa_w_up_to_o))){
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}


#' Checks if a given (conditional) copula has already been specified
#'
#' @param DAG Directed Acyclic Graph object corresponding to the model
#' @param order_hash hashmap of orders of the parental sets
#' @param w node in DAG
#' @param v node in DAG
#' @param cond vector of nodes in DAG
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
  }
  parents_v = bnlearn::parents(x = DAG, node = v)
  if (w %in% parents_v){
    index_w_in_parents = which(parents_v == w)
    parents_up_to_w = if(index_w_in_parents == 1) { c() } else {
      parents_v[1:(index_w_in_parents - 1)] }

    if (sets::as.set(parents_up_to_w) == sets::as.set(cond)){
      return(TRUE)
    }
  }
  parents_w = bnlearn::parents(x = DAG, node = w)
  if (v %in% parents_w){
    index_v_in_parents = which(parents_w == v)
    parents_up_to_v = if(index_v_in_parents == 1) { c() } else {
      parents_w[1:(index_v_in_parents - 1)] }

    if (sets::as.set(parents_up_to_v) == sets::as.set(cond)){
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
         "We are at node: ", v, " and the conditioning set is: ", cond, "\n",
         "Check that the PCBN satisfies the restrictions ",
         "and that the orders of the parents are all compatible.")
  }

  return (list(w = w, cond_set_minus_w = cond_set_minus_w))
}

