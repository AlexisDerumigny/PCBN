#' Simulate data from a specified PCBN
#'
#'
#' @param object PCBN object to simulate from.
#' \bold{This does not work if the PCBN does not abide by the B-sets.
#' And in general, it does not work if the PCBN is outside of
#' the class of restricted PCBNs.}
#'
#' @param check_PCBN check whether the given PCBN satisfies the restrictions.
#' If this is set to \code{FALSE}, no checking is performed.
#' This means that the error due to the a non-restricted PCBN object (if this
#' is the case) will occur later in the computations (and may not be so clear -
#' typically it is because of failing to find a given conditional copula).
#' Nevertheless, even with `check_PCBN = TRUE` if could be that some error happen
#' later if the parental orderings are not compatible with each other.
#'
#' @param N sample size
#'
#' @param verbose if \code{0}, don't print anything.
#' If \code{verbose >= 1}, print information about the simulation procedure.
#'
#' @return a data frame of N samples
#'
#' @examples
#' DAG = create_empty_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#'
#' fam = matrix(c(0, 0, 1,
#'                0, 0, 1,
#'                0, 0, 0), byrow = TRUE, ncol = 3)
#'
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = PCBN_sim(my_PCBN, N = 5)
#'
#'
#' @export
#'
PCBN_sim <- function(object, N, check_PCBN = TRUE, verbose = 1)
{
  if (check_PCBN){
    .checkPCBNobject_for_simulation(object, verbose = verbose)
  }

  # Initialize data frame
  nodes = bnlearn::nodes(object$DAG)
  data = data.frame(matrix(ncol = length(nodes), nrow = N))
  colnames(data) <- nodes

  well_ordering = bnlearn::node.ordering(object$DAG)
  for (node in well_ordering) {
    parents = object$order_hash[[node]]

    # Simulating is analogous to regular vine
    marginal = stats::runif(N, 0, 1) # Start with uniform
    if (length(parents) > 0) {
      # If node has parents, we apply recursion of h-functions
      # The innermost part of the recursion starts with the biggest conditioning
      for (i_parent in length(parents):1) {
        fam = object$copula_mat$fam[parents[i_parent], node]
        tau = object$copula_mat$tau[parents[i_parent], node]
        par = VineCopula::BiCopTau2Par(fam, tau)

        # We must compute the conditional margin parent|lower using a proper recursion
        lower = if (i_parent == 1) {c()} else {parents[1:(i_parent - 1)]}

        parent_given_lower = compute_sample_margin(object = object, data = data,
                                                   v = parents[i_parent],
                                                   cond_set = lower,
                                                   check_PCBN = FALSE,
                                                   verbose = verbose)

        if (verbose > 1){
          cat("Using the invert h function of ", node, " given ", lower, "\n")
        }

        marginal = VineCopula::BiCopHinv1(u1 = parent_given_lower,
                                          u2 = marginal,
                                          family = fam,
                                          par = par)
      }
    }
    data[, node] = marginal
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
#' @param check_PCBN check whether the given PCBN satisfies the restrictions.
#' If this is set to \code{FALSE}, no checking is performed.
#' This means that the error due to the a non-restricted PCBN object (if this
#' is the case) will occur later in the computations (and may not be so clear -
#' typically it is because of failing to find a given conditional copula).
#'
#' @param verbose if \code{0}, don't print anything.
#' If \code{verbose >= 1}, print information about the fitting procedure.
#'
#' @return a vector of size \eqn{n} of realizations \eqn{u_{i, v | cond\_set}}
#' for \eqn{i = 1, \dots, n}.
#'
#' @examples
#'
#' DAG = create_empty_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#'
#' fam = matrix(c(0, 0, 1,
#'                0, 0, 1,
#'                0, 0, 0), byrow = TRUE, ncol = 3)
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
compute_sample_margin <- function(object, data, v, cond_set, check_PCBN = TRUE,
                                  verbose = 1)
{
  if (check_PCBN){
    .checkPCBNobject_for_simulation(object, verbose = verbose)
  }

  # Unpack PCBN object
  DAG = object$DAG
  order_hash = object$order_hash
  copula_mat = object$copula_mat

  if (verbose > 1){
    cat("Trying to estimated the cond cdf of ", v, " given ", cond_set, "...\n")
  }

  # Remove nodes by conditional independence
  cond_set_dependent = remove_CondInd(DAG, v, cond_set)

  if (verbose > 1 & !identical(cond_set, cond_set_dependent)){
    cat("Trying to estimated the cond cdf of ", v, " given ", cond_set_dependent, "...\n")
  }

  # If there is no conditioning set, returns the unchanged observations of v
  if (length(cond_set_dependent) == 0) {
    return (data[, v])
  }

  # Find specified c_{wv|cond_set_minus_w}
  cop_specified = find_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                             v = v, cond = cond_set_dependent)
  w = cop_specified$w
  cond_set_minus_w = cop_specified$cond_set_minus_w

  # we must have w->v or v<-w
  if (copula_mat$fam[w, v] != 0) {
    fam = copula_mat$fam[w, v]
    tau = copula_mat$tau[w, v]
    par = VineCopula::BiCopTau2Par(fam, tau)

    # FIXME: the performance can potentially be improved
    # by caching these sample margins to reuse them whenever needed

    w_given_rest = compute_sample_margin(object, data, w, cond_set_minus_w,
                                         check_PCBN = FALSE, verbose = verbose)

    v_given_rest = compute_sample_margin(object, data, v, cond_set_minus_w,
                                         check_PCBN = FALSE, verbose = verbose)

    v_given_cond = VineCopula::BiCopHfunc1(u1 = w_given_rest,
                                           u2 = v_given_rest,
                                           family = fam,
                                           par = par)
  } else {
    fam = copula_mat$fam[v, w]
    tau = copula_mat$tau[v, w]
    par = VineCopula::BiCopTau2Par(fam, tau)

    w_given_rest = compute_sample_margin(object, data, w, cond_set_minus_w,
                                         check_PCBN = FALSE, verbose = verbose)

    v_given_rest = compute_sample_margin(object, data, v, cond_set_minus_w,
                                         check_PCBN = FALSE, verbose = verbose)

    v_given_cond = VineCopula::BiCopHfunc2(u1 = v_given_rest, u2 = w_given_rest,
                                           family = fam, par = par)
  }

  return(v_given_cond)
}


# This is the internal function used by the functions above to check that the
# PCBN object indeed satisfies the restrictions of no active cycle nor interfering
# v-structure.
#
# We use `verbose = 1` so that it only prints a message if it finds anything.
# Else it prints nothing.
# We use `check_both = FALSE` to save time so that we do not check a second time
# if the first condition is already not satisfied.
#
# If the conditions are not satisfied, this function raises a classed error.
# The class of this error condition is "UnRestrictedPCBNError". This helps for
# the unit test (written in the corresponding test page) so that we can detect
# that the error in the simulation function indeed comes from this test
# (and is therefore user-friendly, as opposed to an error that would come later
# and would be hard to understand for the user - typically that a given
# conditional copula has not been specified).
#
# In the same way, if (at least) one of the ordering does not abide by the B sets,
# this also raises a classed error.
#
# Nevertheless, it does not raises an error if the orderings are not compatible
# with each other (this can happen even if all the orderings abide by the B-sets).
# Testing for this would amount to pre-compute the whole computation tree.
#
# FIXME: maybe we want to do this anyway?
#
.checkPCBNobject_for_simulation <- function(PCBN, verbose)
{
  is_restricted = is_restrictedDAG(PCBN$DAG, verbose = verbose, check_both = FALSE)
  if (!is_restricted){
    stop(errorCondition(
      message = paste0("The DAG does not satisfy the restrictions. ",
                       "Therefore, simulation is not possible."),
      class = "UnRestrictedPCBNError"))
  }

  is_order_abiding_Bsets = is_order_abiding_Bsets(DAG = PCBN$DAG,
                                                  order_hash = PCBN$order_hash)
  if (!is_order_abiding_Bsets){
    stop(errorCondition(
      message = paste0("The parental orderings do not abide by the B-sets. ",
                       "Therefore, simulation is not possible."),
      class = "ParentalOrderingsBsetsError"))
  }
}

