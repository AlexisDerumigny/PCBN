#' Samples from a specified PCBN
#' 
#' @param object PCBN object
#' @param N sample size
#' 
#' @return a data frame of N samples
#' 
sample_PCBN <- function(object, N) {
  # Unpack PCBN object
  DAG = object$DAG
  order_hash = object$order_hash
  copula_mat = object$copula_mat
  
  # Initialize data frame
  nodes = bnlearn::nodes(DAG)
  data = data.frame(matrix(ncol = length(nodes), nrow = N))
  colnames(data) <- nodes
  
  well_ordering = bnlearn::node.ordering(DAG)
  for (node in well_ordering) {
    V = runif(N, 0, 1)
    parents = order_hash[[node]]
    # Simulating is analogous to regular vine
    if (length(parents) > 0) {
      for (parent in rev(parents)) {
        fam = copula_mat$fam[parent, node]
        tau = copula_mat$tau[parent, node]
        par = VineCopula::BiCopTau2Par(fam, tau)
        
        # We must compute the conditional margin parent|lower using a proper recursion
        lower = parents[0:(which(parents == parent) - 1)]
        parent_given_lower = compute_sample_margin(object, data, parent, lower)
        V = VineCopula::BiCopHinv1(parent_given_lower,
                                   V,
                                   family = fam,
                                   par = par)
      }
    }
    data[node] = V
  }
  return(data)
}

#'Computes the conditional margins during sampling
#'
#' @param object PCBN object
#' @param data data frame
#' @param v node
#' @param cond_set conditioning set
#' 
#' @return a vector of realizations \eqn{u_{v|cond_set}}
#' 
compute_sample_margin <- function(object, data, v, cond_set) {
  # Unpack PCBN object
  DAG = object$DAG
  order_hash = object$order_hash
  copula_mat = object$copula_mat
  
  # Remove nodes by conditional independence
  cond_set = PCBN:::remove_CondInd(DAG, v, cond_set)
  
  if (length(cond_set) == 0) {
    return(data[[v]])
  } else{
    # Find specified c_{wv|cond_set_minus_w}
    for (w in cond_set) {
      cond_set_minus_w = cond_set[!cond_set == w]
      if (copula_is_specified(DAG, order_hash, w, v, cond_set_minus_w)) {
        break
      }
    }
    # we must have w->v or v<-w
    if (copula_mat$fam[w, v] != 0) {
      fam = copula_mat$fam[w, v]
      tau = copula_mat$tau[w, v]
      par = VineCopula::BiCopTau2Par(fam, tau)
      
      w_given_rest = compute_sample_margin(object, data, w, cond_set_minus_w)
      v_given_rest = compute_sample_margin(object, data, v, cond_set_minus_w)
      
      v_given_cond = VineCopula::BiCopHfunc1(w_given_rest,
                                 v_given_rest,
                                 family = fam,
                                 par = par)
    } else{
      fam = copula_mat$fam[v, w]
      tau = copula_mat$tau[v, w]
      par = VineCopula::BiCopTau2Par(fam, tau)
      
      w_given_rest = compute_sample_margin(object, data, w, cond_set_minus_w)
      v_given_rest = compute_sample_margin(object, data, v, cond_set_minus_w)
      
      v_given_cond = VineCopula::BiCopHfunc2(v_given_rest, w_given_rest, family = fam, par = par)
    }
  }
  return(v_given_cond)
}