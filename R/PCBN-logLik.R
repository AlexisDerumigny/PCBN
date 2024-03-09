

#' Log-likelihood of a PCBN object
#'
#' This function computes the log-likelihood of the PCBN model given a dataset
#' of i.i.d. observations uniformly (or approximatively uniformly) distributed
#' on \eqn{[0,1]}. This is the same as the logarithm of the density of the PCBN
#' at the observations.
#'
#' @param PCBN the PCBN object
#' @param data_uniform the dataset for which the log-likelihood is computed.
#' It must have already been standardized to uniform margins.
#' @param ... other arguments, ignored for the moment
#'
#' @return the log-likelihood of the PCBN model for the given dataset
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
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = sample_PCBN(my_PCBN, N = 10)
#'
#' logLik(my_PCBN, mydata)
#'
#' @export
#'
logLik.PCBN <- function(PCBN, data_uniform, ...){

  # Unpack PCBN object
  DAG = PCBN$DAG
  order_hash = PCBN$order_hash
  copula_mat = PCBN$copula_mat

  log_lik = 0

  well_ordering = bnlearn::node.ordering(DAG)
  # For every node v, and for every parent w,
  # we need to compute the density of arc c_{ wv | pa(v; w)\{w} }
  for (v in well_ordering) {
    parents = order_hash[[v]]
    if (length(parents) > 0) {
      for (i_parent in 1:length(parents)) {
        w = parents[i_parent]
        fam = copula_mat$fam[w, v]
        tau = copula_mat$tau[w, v]
        par = VineCopula::BiCopTau2Par(fam, tau)

        # Compute the required margins
        lower = if (i_parent == 1) {c()
        } else {parents[1:(i_parent - 1)]}
        v_given_lower = compute_sample_margin(PCBN, data_uniform, v, lower)
        w_given_lower = compute_sample_margin(PCBN, data_uniform, w, lower)

        log_lik_arc_w_to_v = sum(log(
          VineCopula::BiCopPDF(w_given_lower, v_given_lower, family = fam, par = par)))

        log_lik = log_lik + log_lik_arc_w_to_v
      }
    }
  }
  return(log_lik)
}


#' PDF of a PCBN model
#'
#' This function computes the Probability Density Function of a PCBN model.
#'
#' @param PCBN PCBN object
#' @param newdata new data on which the PDF should be computed
#'
#' @return the probability density at newdata.
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
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = sample_PCBN(my_PCBN, N = 10)
#'
#' PCBN_PDF(my_PCBN, mydata)
#'
PCBN_PDF <- function(PCBN, newdata)
{
  return (exp(logLik(PCBN, data_uniform = newdata)))
}

