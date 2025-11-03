

#' Log-likelihood of a PCBN object
#'
#' This function computes the log-likelihood of the PCBN model given a dataset
#' of i.i.d. observations uniformly (or approximatively uniformly) distributed
#' on \eqn{[0,1]}. This is the same as the logarithm of the density of the PCBN
#' at the observations.
#'
#' @param object the PCBN object
#' @param data_uniform the dataset for which the log-likelihood is computed.
#' It must have already been standardized to uniform margins.
#' @param ... other arguments, ignored for the moment
#'
#' @return the log-likelihood of the PCBN model for the given dataset
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
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = PCBN_sim(my_PCBN, N = 10)
#'
#' logLik(my_PCBN, mydata)
#'
#' @export
#'
logLik.PCBN <- function(object, data_uniform, ...){

  # Unpack PCBN object
  DAG = object$DAG
  order_hash = object$order_hash
  copula_mat = object$copula_mat

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
        v_given_lower = compute_sample_margin(object, data_uniform, v, lower)
        w_given_lower = compute_sample_margin(object, data_uniform, w, lower)

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
#' This is a wrapper to \code{\link{logLik.PCBN}}.
#'
#' @param PCBN PCBN object
#' @param newdata new data on which the PDF should be computed
#'
#' @return the probability density at newdata.
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
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = PCBN_sim(my_PCBN, N = 10)
#'
#' PCBN_PDF(my_PCBN, mydata)
#'
#' @export
#'
PCBN_PDF <- function(PCBN, newdata)
{
  return (exp(stats::logLik(PCBN, data_uniform = newdata)))
}



# # Old functions relating to log-likelihood with normal margins ===============
#
#
#
# ### Fits normal margins and returns a list with the mean and sd for each node
# fit_normal_margins <- function(data){
#   nodes = colnames(data)
#   fit_list = list()
#   for (node in nodes){
#     fit = MASS::fitdistr(data[[node]],"normal")
#     estimates = list( mean = fit$estimate[[1]], sd = fit$estimate[[2]])
#     fit_list[[node]] =  estimates
#   }
#   return(fit_list)
# }
#
# ### Scale data-frame with uniform margins to standard normal margins
# to_normal_scale <- function(data){
#   data_new = data
#   for (i in 1:length(data)){
#     data_new[[i]] = stats::qnorm(data[[i]], 0, 1)
#   }
#   return(data_new)
# }
#
# ### Scale data-frame with uniform margins to standard normal margins
# to_uniform_scale <- function(data){
#   data_new = data
#   for (i in 1:length(data)){
#     data_new[[i]] = VineCopula::pobs(data[[i]])
#   }
#   return(data_new)
# }
#
# # Log-likelihood assuming normal margins
# #
# # Computes the marginal part of the log-likelihood
# logLik_margins_PCBN <- function(data_normal, margins){
#   log_lik = 0
#   nodes = colnames(data_normal)
#
#   for (v in nodes) {
#     mean_v = margins[[v]]$mean
#     sd_v = margins[[v]]$sd
#     log_lik = log_lik + sum(log(stats::dnorm(data_normal[[v]], mean = mean_v, sd = sd_v)))
#   }
#   return(log_lik)
# }
#
# # Compute full log-likelihood of PCBN with normal margins
# logLik_PCBN <- function(data, PCBN, margins){
#   data_uniform = to_uniform_scale(data)
#   log_lik = logLik_margins_PCBN(data, margins) + stats::logLik(PCBN, data_uniform)
#   return(log_lik)
# }
#
# # Computes several performance metrics of the GBN:
# # logLik, AIC, BIC, KL divergence and CvM distance
# metrics_GBN <- function(data, PCBN, margins, GBN_fit, N_monte_carlo){
#   metrics =list()
#
#   log_lik = stats::logLik(GBN_fit, data)
#   AIC = -2*stats::BIC(GBN_fit, data)
#   BIC = -2*stats::AIC(GBN_fit, data)
#
#
#   # Generate data for KL and CvM
#   data_help_uniform = PCBN_sim(PCBN, N_monte_carlo)
#   data_help = to_normal_scale(data_help_uniform)
#
#   KL = KL_divergence_GBN(data_help, PCBN, margins, GBN_fit)
#
#   return(list(logLik = log_lik, AIC = AIC, BIC = BIC, KL = KL))
# }
#
#
# # Computes several performance metrics of the PCBN:
# # logLik, AIC, BIC, KL divergence and CvM distance
# metrics_PCBN <- function(data, PCBN, margins, PCBN_fit, margins_fit, N_monte_carlo){
#   metrics =list()
#
#   log_lik = PCBN_fit$metrics$logLik
#   AIC = -PCBN_fit$metrics$AIC
#   BIC = -PCBN_fit$metrics$BIC
#
#   # To compute the logLik, AIC and BIC we need to add the information from the margins
#   log_lik_margins = 0
#   nodes = colnames(data)
#   for (v in nodes) {
#     mean_v = margins[[v]]$mean
#     sd_v = margins[[v]]$sd
#     log_lik_margins = log_lik_margins +
#       sum( log( stats::dnorm(data[[v]], mean = mean_v, sd = sd_v) ) )
#   }
#
#   log_lik = log_lik + log_lik_margins
#   AIC = AIC + (2*(2 * length(margins_fit)) - 2*log_lik_margins)
#   BIC = BIC + (log(length(data))*(2 * length(margins_fit)) - 2*log_lik_margins)
#
#
#   # Generate data for KL and CvM
#   data_help_uniform = PCBN_sim(PCBN, N_monte_carlo)
#   data_help = to_normal_scale(data_help_uniform)
#
#   KL = KL_divergence_PCBN(data_help, PCBN, margins, PCBN_fit, margins_fit)
#
#   return(list(logLik = log_lik, AIC = AIC, BIC = BIC, KL = KL))
# }
#
#
# ### Computing the marginal part of the PDF of a PCBN given normal data
# PDF_margins_PCBN <- function(data_normal, margins){
#   likelihood = 1
#   nodes = colnames(data_normal)
#
#   for (v in nodes) {
#     mean_v = margins[[v]]$mean
#     sd_v = margins[[v]]$sd
#     likelihood = likelihood * prod(stats::dnorm(data_normal[[v]], mean = mean_v, sd = sd_v))
#   }
#   return(likelihood)
# }
#
#
