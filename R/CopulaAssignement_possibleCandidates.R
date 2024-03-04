
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
    } else if (!is.null(incoming_arc(DAG, w, v, order_v, order_hash))){
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
#' @param w,v nodes in DAG. \code{w} is assumed to be a parent of \code{v}.
#'
#' @param order_v partial order for node v
#' @param order_hash hashmap of parental orders
#'
#' @returns a node \code{o} such \code{w} is a parent of \code{o},
#' and \code{w} can be used as an incoming arc to \code{v} by the node \code{o}.
#' If no such \code{o} can be found, \code{w} cannot be used as a potential
#' candidate for the order of \code{v} by incoming arc.
#' This function returns then \code{NULL}.
#'
#' @export
#'
incoming_arc <- function(DAG, w, v, order_v, order_hash)
{
  adj.mat = bnlearn::amat(DAG)

  order_v_sorted = sort(order_v)
  # We inspect each node in the order
  for (o in order_v){
    # We must have w -> o for w to be an incoming arc
    if (adj.mat[w, o] == 1){
      # We look at the order of parents of o (assumed to have been already ordered)
      order_o = order_hash[[o]]
      index_w_in_parents_o = which(order_o == w)
      if (index_w_in_parents_o == 1){
        # If w is the first parent of o,
        # then there is no other conditioning to be added
        pa_o_up_to_w = c()
      } else{
        # Else we take all the parents of o that are before w in the order
        pa_o_up_to_w = order_o[1:index_w_in_parents_o - 1]
      }

      # This vector is the vector of nodes containing
      # o, and its parents up to w (not including w)
      pa_o_up_to_w_and_o = sort(union(o, pa_o_up_to_w))

      if (identical(order_v_sorted, pa_o_up_to_w_and_o) ) {
        # If Ovk is the same set as o union pa_o_up_to_w
        # then we can compute the margin u_{w | Ovk}
        # using the known copula C_{w, o | pa_o_up_to_w} (already estimated at node o)

        return(o)
      } else if ( all(pa_o_up_to_w_and_o %in% order_v) )
      {
        # If all the pa_o_up_to_w are in order_v
        # we check whether the d-separation holds:
        if (dsep_set(DAG = DAG, X = w,
                     Y = setdiff(order_v, pa_o_up_to_w_and_o),
                     Z = pa_o_up_to_w_and_o) ) {
          return(o)
        }
      }
    }
  }

  return(NULL)
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
#' @export
outgoing_arc <- function(DAG, w, v, order_v, order_hash){
  adj.mat = bnlearn::amat(DAG)

  order_v_sorted = sort(order_v)

  for (o in order_v){
    if (adj.mat[o,w]==1){
      order_w = order_hash[[w]]
      index_o_in_parents_w = which(order_w == o)

      if (index_o_in_parents_w == 1){
        pa_w_up_to_o = c()
      } else{
        pa_w_up_to_o = order_w[1:which(order_w == o)-1]
      }

      pa_w_up_to_o_and_o = sort(union(o, pa_w_up_to_o))
      if (identical(order_v_sorted, pa_w_up_to_o_and_o)){
        return(o)
      } else if (setdiff(pa_w_up_to_o, order_w) ){
        if (dsep_set(DAG = DAG, X = w, Y = setdiff(order_v, pa_w_up_to_o_and_o),
                     Z = pa_w_up_to_o_and_o)) {
          return(o)
        }
      }
    }
  }
  return(NULL)
}

