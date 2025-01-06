
# Log-likelihood assuming normal margins
#
# Computes the marginal part of the log-likelihood
logLik_margins_PCBN <- function(data_normal, margins){
  log_lik = 0
  nodes = colnames(data_normal)

  for (v in nodes) {
    mean_v = margins[[v]]$mean
    sd_v = margins[[v]]$sd
    log_lik = log_lik + sum(log(stats::dnorm(data_normal[[v]], mean = mean_v, sd = sd_v)))
  }
  return(log_lik)
}

# Compute full log-likelihood of PCBN with normal margins
logLik_PCBN <- function(data, PCBN, margins){
  data_uniform = to_uniform_scale(data)
  log_lik = logLik_margins_PCBN(data, margins) + stats::logLik(PCBN, data_uniform)
  return(log_lik)
}

###################################################
########  Performance metrics #####################
##################################################

# Computes several performance metrics of the GBN:
# logLik, AIC, BIC, KL divergence and CvM distance
metrics_GBN <- function(data, PCBN, margins, GBN_fit, N_monte_carlo){
  metrics =list()

  log_lik = stats::logLik(GBN_fit, data)
  AIC = -2*stats::BIC(GBN_fit, data)
  BIC = -2*stats::AIC(GBN_fit, data)


  # Generate data for KL and CvM
  data_help_uniform = PCBN_sim(PCBN, N_monte_carlo)
  data_help = to_normal_scale(data_help_uniform)

  KL = KL_divergence_GBN(data_help, PCBN, margins, GBN_fit)

  return(list(logLik = log_lik, AIC = AIC, BIC = BIC, KL = KL))
}


# Computes several performance metrics of the PCBN:
# logLik, AIC, BIC, KL divergence and CvM distance
metrics_PCBN <- function(data, PCBN, margins, PCBN_fit, margins_fit, N_monte_carlo){
  metrics =list()

  log_lik = PCBN_fit$metrics$logLik
  AIC = -PCBN_fit$metrics$AIC
  BIC = -PCBN_fit$metrics$BIC

  # To compute the logLik, AIC and BIC we need to add the information from the margins
  log_lik_margins = 0
  nodes = colnames(data)
  for (v in nodes) {
    mean_v = margins[[v]]$mean
    sd_v = margins[[v]]$sd
    log_lik_margins = log_lik_margins +
      sum( log( stats::dnorm(data[[v]], mean = mean_v, sd = sd_v) ) )
  }

  log_lik = log_lik + log_lik_margins
  AIC = AIC + (2*(2 * length(margins_fit)) - 2*log_lik_margins)
  BIC = BIC + (log(length(data))*(2 * length(margins_fit)) - 2*log_lik_margins)


  # Generate data for KL and CvM
  data_help_uniform = PCBN_sim(PCBN, N_monte_carlo)
  data_help = to_normal_scale(data_help_uniform)

  KL = KL_divergence_PCBN(data_help, PCBN, margins, PCBN_fit, margins_fit)

  return(list(logLik = log_lik, AIC = AIC, BIC = BIC, KL = KL))
}


### Computing the marginal part of the PDF of a PCBN given normal data
PDF_margins_PCBN <- function(data_normal, margins){
  likelihood = 1
  nodes = colnames(data_normal)

  for (v in nodes) {
    mean_v = margins[[v]]$mean
    sd_v = margins[[v]]$sd
    likelihood = likelihood * prod(stats::dnorm(data_normal[[v]], mean = mean_v, sd = sd_v))
  }
  return(likelihood)
}


