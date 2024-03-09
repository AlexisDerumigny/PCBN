
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
  data_help_uniform = sample_PCBN(PCBN, N_monte_carlo)
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
  data_help_uniform = sample_PCBN(PCBN, N_monte_carlo)
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

### Computes the Kullback-Leibler divergence of a estimated PCBN to the true PCBN
KL_divergence_PCBN <-
  function(data, PCBN, margins, PCBN_fit, margins_fit) {
    KL = 0
    L = length(data[[1]])
    KL_vec = rep(0, L)
    data_uniform = to_uniform_scale(data)

    log_lik_true = logLik_PCBN(data, PCBN, margins)
    log_lik_est = logLik_PCBN(data, PCBN_fit, margins_fit)
    KL = (log_lik_true - log_lik_est)/L
    return(KL)
  }

### Computes the Kullback-Leibler divergence of a estimated GBN to the true PCBN
KL_divergence_GBN <- function(data, PCBN, margins, GBN) {
  KL = 0
  L = length(data[[1]])
  KL_vec = rep(0, L)
  data_uniform = to_uniform_scale(data)

  log_lik_true = logLik_PCBN(data, PCBN, margins)
  log_lik_est = stats::logLik(GBN, data)
  KL = (log_lik_true - log_lik_est)/L
  return(KL)
}

# Returns the distance between two graphs
distance_DAGs = function(DAG1, DAG2){
  distance = 0
  node_set =  bnlearn::nodes(DAG1)
  amat1 = bnlearn::amat(DAG1)
  amat2 = bnlearn::amat(DAG2)

  ### Compute difference in the skeleton differences in the skeleton

  # Turn both into undirected graphs
  g1 = igraph::as.undirected(bnlearn::as.igraph(DAG1))
  g2 = igraph::as.undirected(bnlearn::as.igraph(DAG2))

  # Sum all arcs and substract the intersection
  int <- igraph::graph.intersection(g1,g2)
  distance = distance +
    igraph::ecount(g1) + igraph::ecount(g2) - 2 * igraph::ecount(int)


  ### Compute correct v-structures
  for (v in node_set){
    parents1 = DAG1$nodes[[v]]$parents
    parents2 = DAG2$nodes[[v]]$parents

    # One must be part of a v-structure
    if ((length(parents1)>0) | (length(parents2)>0)){
      distance = distance + length(sym_diff(parents1, parents2))
    }
  }
  return(distance)
}

