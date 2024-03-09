

####################################################
########## Computes density/CDF/log-lik of PCBN ####
####################################################

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

# Computes the copula part of the log-likelihood
logLik_copulas_PCBN <- function(data_uniform, PCBN){
  # Unpack PCBN object
  DAG = PCBN$DAG
  order_hash = PCBN$order_hash
  copula_mat = PCBN$copula_mat

  log_lik = 0

  well_ordering = bnlearn::node.ordering(DAG)
  # For every node v for every parent w, we need to compute the density of arc c_{ wv|pa(v; w)\{w} }
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

        log_lik_arc_w_to_v = sum(log(VineCopula::BiCopPDF(w_given_lower, v_given_lower, family = fam, par = par)))

        log_lik = log_lik + log_lik_arc_w_to_v
      }
    }
  }
  return(log_lik)
}

# Compute full log-likelihood of PCBN with margins
logLik_PCBN <- function(data, PCBN, margins){
  data_uniform = to_uniform_scale(data)
  log_lik = logLik_margins_PCBN(data, margins) + logLik_copulas_PCBN(data_uniform, PCBN)
  return(log_lik)
}

