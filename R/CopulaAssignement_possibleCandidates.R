
#' Possible candidates to be added to a partial order
#'
#' When given a partial order of a PCBN, one can complete it by adding one
#' of the parents' node to the partial order. Some nodes can be added; they are
#' then called "possible candidates".
#'
#' \code{possible_candidate_incoming_arc} returns a node \code{o} such \code{w}
#' is a parent of \code{o}, and \code{w} can be used as an incoming arc to
#' \code{v} by the node \code{o}. If no such \code{o} can be found, \code{w}
#' cannot be used as a potential candidate for the order of \code{v} by incoming
#' arc. Then, the function \code{possible_candidate_incoming_arc} returns
#' \code{NULL}.
#'
#' In the same way, \code{possible_candidate_outgoing_arc} returns a node
#' \code{o} such \code{o} is a parent of \code{w}, and \code{w} can be used as
#' an outgoing arc to \code{v} by the node \code{o}.
#'
#' @param DAG Directed Acyclic Graph.
#' @param w,v nodes in DAG. \code{w} is assumed to be a parent of \code{v}.
#' @param order_v partial order for node v.
#' @param order_hash hashmap of parental orders
#' @param B_minus_O this is the current B-set, without the elements of
#' \code{order_v}, i.e. this is the set of elements that could be considered
#' possible candidates.
#'
#' @returns \code{possible_candidates} returns a vector of possible candidates,
#' potentially empty.
#' Both \code{possible_candidate_incoming_arc} and
#' \code{possible_candidate_outgoing_arc} return either a node \code{o}, or
#' \code{NULL} if they could not find such a node.
#'
#' @seealso \code{\link{dsep_set}} for checking whether two sets of nodes are d-separated
#' by another set.
#' \code{\link{find_B_sets}} to find the B-sets.
#'
#'
#' @examples
#'
#' DAG = create_empty_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#'
#' # Node of interest
#' v = "U4"
#'
#' # If we start by 1, then the arc 1 -> 3 cannot be used as an incoming arc
#' # (it is actually an outgoing arc)
#' possible_candidate_incoming_arc(
#'   DAG = DAG, w = "U3", v = v, order_v = c("U1"), order_hash = order_hash)
#' possible_candidate_outgoing_arc(
#'   DAG = DAG, w = "U3", v = v, order_v = c("U1"), order_hash = order_hash)
#' possible_candidates(
#'   DAG = DAG, v = v, order_v = c("U1"), order_hash = order_hash, B_minus_O = "U2")
#'
#'
#' @export
possible_candidates <- function(DAG, v, order_v, order_hash, B_minus_O)
{
  Poss.Cand = c()
  # We check every element in `B_minus_O` to see whether it can be a
  # possible candidate.
  for (w in B_minus_O){
    # Independence
    if (dsep_set(DAG, w, order_v)){
      Poss.Cand = append(Poss.Cand, w)
    } else if (!is.null(possible_candidate_incoming_arc(DAG, w, v,
                                                        order_v, order_hash))){
      Poss.Cand = append(Poss.Cand, w)
    } else if (!is.null(possible_candidate_outgoing_arc(DAG, w, v,
                                                        order_v, order_hash))){
      Poss.Cand = append(Poss.Cand, w)
    }
  }
  return(Poss.Cand)
}


#' @rdname possible_candidates
#' @export
possible_candidate_incoming_arc <- function(DAG, w, v, order_v, order_hash)
{
  adj.mat = bnlearn::amat(DAG)

  order_v_sorted = sort(order_v)
  # We inspect each node in the order
  for (o in order_v){
    # We must have w -> o for w to be an incoming arc
    if (adj.mat[w, o] == 1){
      # We look at the order of parents of o
      # (assumed to have been already ordered)
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
        # using the known copula C_{w, o | pa_o_up_to_w}
        # (already estimated at node o)

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

#' @rdname possible_candidates
#' @export
#'
possible_candidate_outgoing_arc <- function(DAG, w, v, order_v, order_hash){
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
      } else if ( all(pa_w_up_to_o %in% order_v) ){
        if (dsep_set(DAG = DAG, X = w, Y = setdiff(order_v, pa_w_up_to_o_and_o),
                     Z = pa_w_up_to_o_and_o)) {
          return(o)
        }
      }
    }
  }
  return(NULL)
}

