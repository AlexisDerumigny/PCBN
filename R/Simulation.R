
#' Samples from a specified PCBN
#'
#'
#' @param object PCBN object to sample from.
#' \bold{This does not work if the PCBN does not abide by the B-sets.
#' And in general, it does not work if the PCBN is outside of
#' the class of restricted PCBNs.}
#'
#' @param N sample size
#'
#' @return a data frame of N samples
#'
#' @examples
#' DAG = create_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#'
#' fam = matrix(c(0, 1, 1,
#'                0, 0, 1,
#'                0, 0, 0), byrow = TRUE, ncol = 3)
#'
#' # FIXME
#' # automatically put the dimnames
#' # when creating PCBN objects
#'
#' rownames(fam) <- c("U1", "U2", "U3")
#' colnames(fam) <- c("U1", "U2", "U3")
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = sample_PCBN(my_PCBN, N = 5)
#'
#'
#' @export
#'
sample_PCBN <- function(object, N) {

  # Initialize data frame
  nodes = bnlearn::nodes(object$DAG)
  data = data.frame(matrix(ncol = length(nodes), nrow = N))
  colnames(data) <- nodes

  well_ordering = bnlearn::node.ordering(object$DAG)
  for (node in well_ordering) {
    parents = object$order_hash[[node]]

    # Simulating is analogous to regular vine
    if (length(parents) > 0) {
      for (i_parent in length(parents):1) {
        fam = object$copula_mat$fam[parents[i_parent], node]
        tau = object$copula_mat$tau[parents[i_parent], node]
        par = VineCopula::BiCopTau2Par(fam, tau)

        # We must compute the conditional margin parent|lower using a proper recursion
        lower = if (i_parent == 1) {c()
          } else {parents[1:(i_parent - 1)]}
        parent_given_lower = compute_sample_margin(object = object, data = data,
                                                   v = parents[i_parent],
                                                   cond_set = lower)
        data[, node] = VineCopula::BiCopHinv1(u1 = parent_given_lower,
                                              u2 = stats::runif(N, 0, 1),
                                              family = fam,
                                              par = par)
      }
    } else { # if there are no parents
      data[, node] = stats::runif(N, 0, 1)
    }
  }
  return(data)
}


#' Computes a conditional margin during sampling
#'
#' @param object PCBN object to sample from.
#' \bold{This does not work if the PCBN does not abide by the B-sets.
#' And in general, it does not work if the PCBN is outside of
#' the class of restricted PCBNs.}
#'
#' @param data data frame of observations of size \code{n}
#' @param v name of the node
#' @param cond_set conditioning set.
#' This is a vector containing the names of all the nodes in the conditioning set.
#'
#' @return a vector of size \eqn{n} of realizations \eqn{u_{i, v | cond_set}}
#' for \eqn{i = 1, \dots, n}.
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
#' fam = matrix(c(0, 1, 1,
#'                0, 0, 1,
#'                0, 0, 0), byrow = TRUE, ncol = 3)
#'
#' rownames(fam) <- c("U1", "U2", "U3")
#' colnames(fam) <- c("U1", "U2", "U3")
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' # Initialize data frame
#' N = 100
#' nodes = bnlearn::nodes(my_PCBN$DAG)
#' data = data.frame(matrix(ncol = length(nodes), nrow = N))
#' colnames(data) <- nodes
#'
#' data[, "U1"] = stats::runif(N)
#' data[, "U2"] = stats::runif(N)
#' u_1_given2 = compute_sample_margin(object = my_PCBN, data = data,
#'                                    v = "U1", cond_set = c("U2"))
#'
#' identical(data[, "U1"], u_1_given2)
#'
#' @export
#'
compute_sample_margin <- function(object, data, v, cond_set) {

  # Unpack PCBN object
  DAG = object$DAG
  order_hash = object$order_hash
  copula_mat = object$copula_mat

  # Remove nodes by conditional independence
  cond_set_dependent = remove_CondInd(DAG, v, cond_set)

  # If there is no conditioning set, returns the unchanged observations of v
  if (length(cond_set_dependent) == 0) {
    return (data[, v])
  }

  # Find specified c_{wv|cond_set_minus_w}
  w = NULL
  for (i_w in 1:length(cond_set_dependent)) {
    w_proposed = cond_set_dependent[i_w]
    cond_set_minus_w = cond_set_dependent[-i_w]
    if (is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                 w = w_proposed, v = v,
                                 cond = cond_set_minus_w)) {
      w = w_proposed
      break
    }
  }
  if (is.null(w)) {
    stop("no specified conditional copula found in `compute_sample_margin`.\n",
         "Check that the PCBN satisfies the restrictions ",
         "and that the orders of the parents are all compatible.")
  }

  # we must have w->v or v<-w
  if (copula_mat$fam[w, v] != 0) {
    fam = copula_mat$fam[w, v]
    tau = copula_mat$tau[w, v]
    par = VineCopula::BiCopTau2Par(fam, tau)

    # FIXME: the performance can potentially be improved
    # by caching these sample margins to reuse them whenever needed

    w_given_rest = compute_sample_margin(object, data, w, cond_set_minus_w)
    v_given_rest = compute_sample_margin(object, data, v, cond_set_minus_w)

    v_given_cond = VineCopula::BiCopHfunc1(u1 = w_given_rest,
                                           u2 = v_given_rest,
                                           family = fam,
                                           par = par)
  } else {
    fam = copula_mat$fam[v, w]
    tau = copula_mat$tau[v, w]
    par = VineCopula::BiCopTau2Par(fam, tau)

    w_given_rest = compute_sample_margin(object, data, w, cond_set_minus_w)
    v_given_rest = compute_sample_margin(object, data, v, cond_set_minus_w)

    v_given_cond = VineCopula::BiCopHfunc2(u1 = v_given_rest, u2 = w_given_rest,
                                           family = fam, par = par)
  }

  return(v_given_cond)
}
