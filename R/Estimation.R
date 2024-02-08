
#' Fits the copula joining w and v given cond_set abiding
#' by the conditional independencies of the graph
#'
#' @param data data frame
#' @param DAG Directed Ayclic Graph
#' @param v,w nodes of the graph
#' @param cond_set vector of nodes of DAG in the conditioning set
#' @param familyset vector of copula families
#' @param order_hash hashmap of parental orders
#'
#' @returns copula object
#'
BiCopCondFit <- function(data, DAG, w, v, cond_set, familyset, order_hash){
  if (is.null(copula_hash[[create_copula_tag(DAG, order_hash, w, v, cond_set)]])){
    w_given_cond = ComputeCondMargin(data, DAG, w, cond_set, familyset, order_hash)
    v_given_cond = ComputeCondMargin(data, DAG, v, cond_set, familyset, order_hash)
    C_wv = VineCopula::BiCopSelect(w_given_cond, v_given_cond, familyset = familyset)
    copula_hash[[create_copula_tag(DAG, order_hash, w, v, cond_set)]] = C_wv
  }
  return(copula_hash[[create_copula_tag(DAG, order_hash, w, v, cond_set)]])
}

#' Computes the margin v given cond_set abiding by the conditional independenciese of the graph
#'
#' @param data data frame
#' @param DAG Directed Ayclic Graph
#' @param v node of the graph
#' @param cond_set vector of nodes of DAG in the conditioning set
#' @param familyset vector of copula families
#' @param order_hash hashmap of parental orders
#'
ComputeCondMargin <- function(data, DAG, v, cond_set, familyset, order_hash){
  # remove elements by conditional independence
  cond_set = remove_CondInd(DAG, v, cond_set)
  if (length(cond_set)==0){
    margin_hash[[create_margin_tag(DAG, order_hash, v, cond_set)]] = data[[v]]
  } else{


    # Check if we already computed the margin
    if (is.null(margin_hash[[create_margin_tag(DAG, order_hash, v, cond_set)]])){
      # To compute we need C_{w,v|cond_set\{w}} and use the h-function
      # Pick a w such that this copula is specified
      for (w in cond_set){
        cond_set_minus_w = cond_set[!cond_set==w]
        if (is_cond_copula_specified(DAG, order_hash, w, v, cond_set_minus_w)){
          break
        }
      }
      # If our copula is not in the hash map -> fit it!
      if (is.null(copula_hash[[create_copula_tag(DAG, order_hash, w, v, cond_set_minus_w)]])){
        C_wv = BiCopCondFit(data, DAG, w, v, cond_set_minus_w, familyset, order_hash)
      }
      else{
        C_wv = copula_hash[[create_copula_tag(DAG, order_hash, w, v, cond_set_minus_w)]]
      }
      # Compute v|cond_set_minus_w and w|cond_set_minus_w with the h-functions
      w_given_rest = ComputeCondMargin(data, DAG, w, cond_set_minus_w, familyset, order_hash)
      v_given_rest = ComputeCondMargin(data, DAG, v, cond_set_minus_w, familyset, order_hash)

      # Compute v|cond_set
      v_given_cond = VineCopula::BiCopHfunc1(w_given_rest, v_given_rest, obj = C_wv)
      margin_hash[[create_margin_tag(DAG, order_hash, v, cond_set)]] = v_given_cond
    }
  }
  return(margin_hash[[create_margin_tag(DAG, order_hash, v, cond_set)]])
}


#' Fit all possible orders given a DAG
#'
#' @param data data frame
#' @param DAG Directed Acyclic Graph
#' @param familyset vector of copula families
#'
#' @returns list containing best fit and all fitted models
#'
fit_all_orders <- function(data, DAG, familyset = c(1, 3, 4, 5, 6), reuse_hash = FALSE) {
  if (!reuse_hash){
    assign("copula_hash", r2r::hashmap(), envir = .GlobalEnv)
    assign("margin_hash", r2r::hashmap(), envir = .GlobalEnv)
  }

  all_orders = find_all_orders(DAG)
  fitted_list = list()
  for (order in all_orders) {
    fitted_PCBN = fit_copulas(data, DAG, order, familyset)
    fitted_list[[length(fitted_list) + 1]] = fitted_PCBN
  }

  best_fit = fitted_list[[1]]
  for (i in 1:length(fitted_list)) {
    fit = fitted_list[[i]]
    if (fit$metrics$BIC > best_fit$metrics$BIC) {
      best_fit = fit
    }
  }

  return(list(best_fit = best_fit, fitted_list = fitted_list))
}

#' Fit the copulas of a PCBN given data
#'
#' @param DAG Directed Acyclic Graph
#' @param order_hash hashmap of parental orders
#' @param familyset vector of copula families
#'
#' @returns all fitted copulas
#'
#' @seealso [BiCopCondFit] which this function wraps.
#'
fit_copulas <-
  function(data,
           DAG,
           order_hash,
           familyset = c(1, 3, 4, 5, 6)) {
    tau = bnlearn::amat(DAG)
    fam = bnlearn::amat(DAG)

    logLik = 0
    BIC = 0
    AIC = 0

    node.names = bnlearn::node.ordering(DAG)
    for (v in node.names) {
      parents = order_hash[[v]]
      if (length(parents) > 0) {
        for (i_parent in 1:length(parents)) {
          w = parents[i_parent]
          parents_up_to_w = if (i_parent == 1) {c()
          } else {parents[1:(i_parent - 1)]}

          C = BiCopCondFit(data, DAG, w, v, parents_up_to_w, familyset, order_hash)
          from = i_parent

          # FIXME: we can get a small performance speedup by removing
          # this `which` command, and changing the loop in v above.
          to = which(bnlearn::nodes(DAG) == v)
          tau[from, to] = C$tau
          fam[from, to] = C$family

          logLik = logLik + C$logLik
          BIC = BIC - C$BIC
          AIC = AIC - C$AIC
        }
      }
    }
    copula_mat = list(tau = tau, fam = fam)
    metrics = as.data.frame(list(
      logLik = logLik,
      BIC = BIC,
      AIC = AIC
    ))
    PCBN = new_PCBN(
      DAG = DAG,
      order_hash = order_hash,
      copula_mat = copula_mat
    )
    PCBN$metrics = metrics
    return(PCBN)
  }
