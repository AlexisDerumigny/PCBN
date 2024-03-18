
# Default starting environment
#'
#'
#' @rdname BiCopCondFit
#' @export
#'
default_envir <- function(){
  e = new.env()
  e$copula_hash = r2r::hashmap()
  e$margin_hash = r2r::hashmap()
  e$keychain = r2r::hashmap()
  return (e)
}


#' Fits the copula joining w and v given cond_set abiding
#' by the conditional independencies of the graph
#'
#' @param data data frame
#' @param DAG Directed Ayclic Graph
#' @param v,w nodes of the graph
#' @param cond_set vector of nodes of DAG. They should all be parents of v.
#' They should be ordered from the smallest to the biggest.
#' @param familyset vector of copula families
#' @param order_hash hashmap of parental orders
#' @param e environment containing all the hashmaps
#' @param verbose if \code{0}, don't print anything.
#' If \code{verbose = 1}, print information about the fitting procedure.
#'
#' @return \code{default_envir} returns an environment to be passed
#' to \code{BiCopCondFit} or to \code{ComputeCondMargin}. \code{BiCopCondFit}
#' returns the fitted copula object of \code{v}, \code{w} given \code{cond_set}.
#' \code{ComputeCondMargin} returns the fitted conditional margins of \code{v}
#' given \code{cond_set}.
#'
#' Both functions store all intermediary results in \code{e} to save computation
#' time.
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
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = PCBN_sim(my_PCBN, N = 5)
#' e = default_envir()
#' ls(e)
#' C_13 = BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U3",
#'                     cond_set = c(), familyset = 1, order_hash = order_hash,
#'                     e = e)
#'
#' C_23_1 = BiCopCondFit(data = mydata, DAG = DAG, v = "U2", w = "U3",
#'                       cond_set = "U1", familyset = 1, order_hash = order_hash,
#'                       e = e)
#'
#' U_2_13 = ComputeCondMargin(data = mydata, DAG = DAG,
#'                            v = "U2", cond_set = c("U1", "U3"),
#'                            familyset = 1, order_hash = order_hash, e = e)
#'
#' @export
#'
BiCopCondFit <- function(data, DAG, v, w, cond_set, familyset, order_hash, e,
                         verbose = 1)
{
  if (v > w){
    # We switch them. From now on, v < w
    return (BiCopCondFit(data = data, DAG = DAG, v = w, w = v,
                         cond_set = cond_set, familyset = familyset,
                         order_hash = order_hash, e = e))
  }

  # If the conditioning set is null,
  # we assume it means that there is no conditioning, i.e. a vector of size 0
  if (is.null(cond_set)) {cond_set = character(0)}

  # We look for the copula in the hash.
  copula_key = e$keychain[[list(margins = c(v, w), cond = cond_set)]]
  if (!is.null(copula_key))
  {
    C_wv = e$copula_hash[[copula_key]]
  } else {
    # We look for the keys of the two conditional margins in the hash
    # We try to simplify the conditioning set first
    cond_set_v = remove_CondInd(DAG, v, cond_set)
    v_key = e$keychain[[list(margin = v, cond = cond_set_v)]]

    cond_set_w = remove_CondInd(DAG, w, cond_set)
    w_key = e$keychain[[list(margin = w, cond = cond_set_w)]]

    if (!is.null(v_key) && !is.null(w_key)){
      copula_key <- make_and_store_keyCopula(v = v, w = w, cond = cond_set,
                                             v_key = v_key, w_key = w_key, e = e)

      C_wv = e$copula_hash[[copula_key]]
    } else {
      # key not available yet
      C_wv = NULL
    }
  }

  if (!is.null(C_wv) ){
    # We have already the copula, so we can just return it
    return (C_wv)
  }

  # We now need to estimate the (conditional) copula
  # so we first get the two margins
  v_given_cond = ComputeCondMargin(data, DAG, v, cond_set, familyset, order_hash,
                                   e = e, verbose = verbose)

  w_given_cond = ComputeCondMargin(data, DAG, w, cond_set, familyset, order_hash,
                                   e = e, verbose = verbose)

  # We can now estimate the (simplified) conditional copula
  if (verbose > 0){
    cat("Estimating the copula of ", v, " and ", w,
        if (length(cond_set)) {c(" given ", cond_set)} else {c()}, "\n")
  }
  C_wv = VineCopula::BiCopSelect(w_given_cond, v_given_cond, familyset = familyset)

  if (is.null(copula_key)){
    v_key = e$keychain[[list(margin = v, cond = cond_set_v)]]
    w_key = e$keychain[[list(margin = w, cond = cond_set_w)]]

    # The key for this conditional copula is not present
    # so we rebuild it ourselves, from the two (conditional) marginal keys
    copula_key <- make_and_store_keyCopula(v = v, w = w, cond = cond_set,
                                           v_key = v_key, w_key = w_key, e = e)
  }

  # We finally store the copula in the hash
  e$copula_hash[[copula_key]] = C_wv

  # and we announce that two new conditional margins are available in the keychain
  make_and_store_keyMargin(v = v, cond = sort(c(w, cond_set)),
                           copula_key = copula_key, e = e)

  make_and_store_keyMargin(v = w, cond = sort(c(v, cond_set)),
                           copula_key = copula_key, e = e)

  return(C_wv)
}


# Computation of conditional margins
#' @rdname BiCopCondFit
#'
#' @export
#'
ComputeCondMargin <- function(data, DAG, v, cond_set, familyset, order_hash,
                              e, verbose = 1)
{
  # Remove as much elements as possible by conditional independence
  cond_set = remove_CondInd(DAG, v, cond_set)

  # 1- We see if the result is already available ===============================

  # If there are no more elements in cond_set
  # this means that we are in the case of an unconditional margin
  if (length(cond_set) == 0){
    v_key_result = make_and_store_keyMargin(v = v, cond = character(0),
                                            copula_key = NULL, e = e)
    # We can just save this in the hashmap
    e$margin_hash[[v_key_result]] = data[, v]
    # and the margin information in the keychain
    e$keychain[[list(margin = v, cond = character(0))]] = v_key_result

    return (data[, v] )
  }

  # Check if we already computed this margin
  v_key_result = e$keychain[[list(v, cond_set)]]
  if ( !is.null(v_key_result) )
  {
    return ( e$margin_hash[[v_key_result]] )
  }

  # 2- We find a good `w` ======================================================

  # To compute we need C_{w,v|cond_set\{w}} and use the h-function
  # Pick a `w` in `cond_set` such that this copula is specified
  cop_specified = find_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                             v = v, cond = cond_set)
  w = cop_specified$w
  cond_set_minus_w = cop_specified$cond_set_minus_w

  # We try to simplify this new conditioning set
  cond_set_minus_w_simpV = remove_CondInd(DAG, v, cond_set_minus_w)


  # 3- We find all necessary keys ==============================================

  # We look for the two conditional margins in the hash
  v_key = e$keychain[[list(margin = v, cond = cond_set_minus_w_simpV)]]
  if (is.null(v_key)){
    stop("The conditional margin ", v, " | ", cond_set_minus_w_simpV,
         " is not available.")
  }

  # We try to simplify the conditioning set for w
  cond_set_minus_w_simpW = remove_CondInd(DAG, w, cond_set_minus_w)
  w_key = e$keychain[[list(margin = w, cond = cond_set_minus_w_simpW)]]
  if (is.null(w_key)){
    stop("The conditional margin ", w, " | ", cond_set_minus_w_simpW,
         " is not available.")
  }

  # We look for the copula in the hash.
  copula_key = e$keychain[[list(margins = sort(c(v, w)), cond = cond_set_minus_w)]]
  if (is.null(copula_key)){
    # The key for this conditional copula is not present
    # so we rebuild it ourselves, from the two (conditional) marginal keys
    if (v < w){
      copula_key <- make_and_store_keyCopula(v = v, w = w, cond = cond_set_minus_w,
                                             v_key = v_key, w_key = w_key, e = e)
    } else {
      copula_key <- make_and_store_keyCopula(v = w, w = v, cond = cond_set_minus_w,
                                             v_key = w_key, w_key = v_key, e = e)
    }
  }

  # 4- We get the copula =======================================================

  # We try to get the copula from the hash
  C_wv = e$copula_hash[[copula_key]]

  # If our copula is not in the hash map -> fit it!
  if (is.null(C_wv)){
    C_wv = BiCopCondFit(data, DAG, w, v, cond_set_minus_w, familyset, order_hash,
                        e = e, verbose = verbose)
  }

  # 5- We get the two conditional margins ======================================

  # Get v|cond_set_minus_w and w|cond_set_minus_w on the hash or recompute them
  v_given_rest = e$margin_hash[[v_key]]
  w_given_rest = e$margin_hash[[w_key]]
  if (is.null(v_given_rest)){
    v_given_rest <- ComputeCondMargin(data = data, DAG = DAG, v = v,
                                      cond_set = cond_set_minus_w,
                                      familyset = familyset,
                                      order_hash = order_hash, e = e, verbose = verbose)
  }
  if (is.null(w_given_rest)){
    w_given_rest <- ComputeCondMargin(data = data, DAG = DAG, v = w,
                                      cond_set = cond_set_minus_w_simpW,
                                      familyset = familyset,
                                      order_hash = order_hash, e = e, verbose = verbose)
  }

  # 6- We compute the new conditional margin ===================================

  # Compute v|cond_set under the simplifying assumption
  if (verbose > 0){
    cat("Estimating the cond cdf of ", v, " given ", cond_set, "\n")
  }
  v_given_cond = VineCopula::BiCopHfunc1(w_given_rest, v_given_rest,
                                         obj = C_wv)

  # 7- We save and return the result ===========================================

  # We can just save this in the hashmap
  v_key_result = make_and_store_keyMargin(v = v, cond = cond_set,
                                          copula_key = copula_key, e = e)
  e$margin_hash[[v_key_result]] = v_given_cond

  return( v_given_cond )
}


#' Fit the copulas of a PCBN given data
#'
#' @param data data frame
#' @param DAG Directed Acyclic Graph
#' @param order_hash hashmap of parental orders
#' @param familyset vector of copula families
#' @param e environment containing all the hashmaps
#' @param score_metric name of the metric used to choose the best order.
#' Possible choices are \code{logLik}, \code{AIC} and \code{BIC}.
#'
#' @returns \code{fit_copulas} returns the fitted PCBN, with an additional
#' element called \code{metrics} which is a named vector of length 3 with names
#' \code{c("logLik", "BIC", "AIC")}, where
#' \eqn{AIC = - 2 * logLik + 2 * nparameters}
#' and \eqn{BIC = - 2 * logLik + log(n) * nparameters},
#' for a sample size \code{n} and \code{nparameters} is the number of parameters.
#'
#' \code{fit_all_orders} returns a list containing: \itemize{
#'   \item \code{best_fit} the PCBN which is the best according to the metric
#'   \code{score_metric}.
#'
#'   \item \code{fitted_list} the list of all fitted PCBNs.
#'
#'   \item \code{metrics} the matrix of metrics (logLik, BIC, AIC).
#'   Each row \code{i} of this matrix corresponds to a PCBN with a different set
#'   of parents' orderings, and corresponds to element \code{i} of
#'   \code{fitted_list}.
#' }
#'
#' @seealso [BiCopCondFit] which this function wraps.
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
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = PCBN_sim(my_PCBN, N = 5)
#' e = default_envir()
#'
#' result = fit_copulas(data = mydata, DAG = DAG,
#'                      order_hash = order_hash,
#'                      familyset = 1, e = e)
#'
#' result_all_orders = fit_all_orders(data = mydata, DAG = DAG,
#'                                    familyset = 1, e = e)
#'
#' # The two fitted PCBNs are:
#' print(result_all_orders$fitted_list[[1]])
#' print(result_all_orders$fitted_list[[2]])
#' # and the metrics matrix is:
#' print(result_all_orders$metrics)
#'
#' # The PCBN corresponding to the true order U1 < U2 is usually better
#' # than the second one. This Will be more clear with a bigger sample size.
#'
#'
#' @export
#'
fit_copulas <- function(data,
                        DAG,
                        order_hash,
                        familyset = c(1, 3, 4, 5, 6),
                        e) {
  tau = bnlearn::amat(DAG)
  fam = bnlearn::amat(DAG)

  logLik = 0
  BIC = 0
  AIC = 0

  node.names = bnlearn::node.ordering(DAG)
  for (i_v in 1:length(node.names)) {
    v = node.names[i_v]
    parents = order_hash[[v]]

    # We only have a copula to estimate if `v` has any parent(s).
    if (length(parents) > 0) {
      for (i_parent in 1:length(parents)) {
        w = parents[i_parent]
        parents_up_to_w = if (i_parent == 1) {c()
        } else {parents[1:(i_parent - 1)]}

        C = BiCopCondFit(data, DAG, w, v, parents_up_to_w, familyset, order_hash, e = e)

        tau[i_parent, i_v] = C$tau
        fam[i_parent, i_v] = C$family

        logLik = logLik + C$logLik
        BIC = BIC + C$BIC
        AIC = AIC + C$AIC
      }
    }
  }
  copula_mat = list(tau = tau, fam = fam)
  metrics = c(logLik = logLik, BIC = BIC, AIC = AIC)

  PCBN = new_PCBN(
    DAG = DAG,
    order_hash = order_hash,
    copula_mat = copula_mat
  )
  PCBN$metrics = metrics
  return(PCBN)
}



#' Fit all possible orders given a DAG
#'
#' @rdname fit_copulas
#' @export
#'
fit_all_orders <- function(data, DAG, familyset = c(1, 3, 4, 5, 6),
                           e, score_metric = "BIC")
{
  all_metrics = c("logLik", "BIC", "AIC")
  if (! (score_metric %in% all_metrics) ){
    stop("Invalid 'score_metric': ", score_metric,
         "\nPossible choices are: ",
         paste0(paste0("'", all_metrics, "'", sep = ""),
                collapse = ", "), ".")
  }
  all_orders = find_all_orders(DAG)
  fitted_list = list()

  metrics = matrix(nrow = length(all_orders), ncol = 3)
  colnames(metrics) <- all_metrics

  for (i_order in 1:length(all_orders)) {
    order = all_orders[[i_order]]
    fitted_PCBN = fit_copulas(data, DAG, order, familyset, e = e)

    fitted_list[[i_order]] = fitted_PCBN
    metrics[i_order, ] = fitted_PCBN$metrics
  }

  if (score_metric == "logLik"){
    i_best_fit = which.max(metrics[, score_metric])[1]
  } else {
    i_best_fit = which.min(metrics[, score_metric])[1]
  }

  best_fit = fitted_list[[i_best_fit]]

  return(list(best_fit = best_fit, fitted_list = fitted_list,
              metrics = metrics))
}

