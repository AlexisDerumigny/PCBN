
######################################################
###########     Hill climbing   ###################
###################################################

#' Hill climbing algorithm for restricted PCBNs
#'
#' @param data data frame
#'
#' @param start starting Directed Acyclic Graph. If \code{NULL} (the default),
#' we start with the empty graph.
#'
#' @param familyset vector of copula families
#'
#' @param score_metric name of the metric used to choose the best order.
#' Possible choices are \code{logLik}, \code{AIC} and \code{BIC}.
#'
#' @param verbose if \code{0}, don't print anything.
#' If \code{verbose >= 1}, print information about the procedure.
#'
#' @param e environment containing all the hashmaps
#'
#'
#' @returns DAG which locally maximizes BIC based score function
#'
#' @examples
#'
#' DAG = create_empty_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U2", "U1")
#' order_hash[['U4']] = c("U2", "U3")
#'
#' fam = matrix(c(0, 0, 1, 0,
#'                0, 0, 1, 1,
#'                0, 0, 0, 1,
#'                0, 0, 0, 0), byrow = TRUE, ncol = 4)
#' tau = 0.8 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = PCBN_sim(my_PCBN, N = 5)
#'
#' result = hill.climbing.PCBN(data = mydata, familyset = 1)
#'
#' @export
#'
hill.climbing.PCBN <- function(data, familyset = c(1, 3, 4, 5, 6), verbose = 2,
                               start = NULL, e = NULL, score_metric = "BIC")
{
  if (verbose > 0) {
    cat("----------------------------------------------------------------\n")
  }

  if (is.null(start)){
    nodes = colnames(data)
    DAG = bnlearn::empty.graph(nodes)
    if (verbose > 0){
      cat("* Starting from the empty graph.\n")
    }
  } else{
    DAG = start
  }
  if (is.null(e)){
    e = default_envir()
  }
  fitted = fit_all_orders(data = data, DAG = DAG, familyset = familyset, e = e,
                          score_metric = score_metric, verbose = verbose - 2)
  reference = fitted$best_fit$metrics[score_metric]

  if (verbose > 0) {
    cat("* current score:", reference, "\n")
  }

  iter = 1
  repeat{
    if (verbose > 0){
      cat(paste0("* Starting iteration ", iter, "\n"))
    }
    allowed.operations = allowed_operations_fromDAG(DAG)

    # Compute data frame with the score delta of each operation
    df = operation_score_deltas(data, DAG, familyset, allowed.operations, e = e,
                                score_metric = score_metric,
                                verbose = verbose - 1)

    # Select best operation based on the score
    if (score_metric == "logLik"){
      whichBest = which.max(df$score)
      bestop = df[whichBest, ]
      improvement = bestop$score - reference
    } else {
      whichBest = which.min(df$score)
      bestop = df[whichBest, ]
      improvement = reference - bestop$score
    }
    if (verbose == 1){
      cat("Number of possible operations:" , nrow(allowed.operations), "\n")
      cat("*best operation:\n")
      print(bestop)

    } else if (verbose >= 2){
      cat("----------------------------------------------------------------\n")
      cat("* possible operations:\n")
      df$best = ""
      df$best[whichBest] = "best"
      print(df)
    }

    if (improvement > 0){
      DAG = operation_do(DAG, op = bestop)
      reference = bestop$score

      # if (verbose > 0){
      #   print(DAG)
      # }
    }
    else{
      if (verbose > 0){
        cat("No improvement found. Stopping Hill Climbing algorithm.\n")
      }
      break;
    }

    if (verbose > 0) {
      cat("----------------------------------------------------------------\n")
    }
    iter = iter + 1
  }

  fitted = fit_all_orders(data = data, DAG = DAG, familyset = familyset, e = e,
                          verbose = verbose - 2)

  if (verbose > 0){
    cat("----------------------------------------------------------------\n")
    cat("Best PCBN:\n")
    print(fitted$best_fit)
  }

  return(fitted)
}

# Computes the score delta for all allowed operations
operation_score_deltas = function(data, DAG, familyset, allowed.operations,
                                  e, score_metric, verbose = 1)
{
  # Loop over all allowed operations
  for (i in 1:nrow(allowed.operations)){
    op = allowed.operations[i,]
    DAG_new = operation_do(DAG, op)

    # Fit all possible orders
    fitted = fit_all_orders(data, DAG_new, familyset, e = e,
                            verbose = verbose - 1)
    score = fitted$best_fit$metrics[score_metric]

    allowed.operations$score[i] = score
  }
  return(allowed.operations)
}


#' Finds all operations resulting in a restricted DAG
#'
#' @param DAG the current DAG.
#'
#' @param op the operation to be applied
#' It should be a list with at least the elements \code{from}, \code{to}, and
#' \code{operation}. \code{From} and \code{to} are the nodes (characters) and
#' \code{operation} is a character describing the operation to be applied.
#'
#' @return \code{allowed_operations_fromDAG} returns a \code{data.frame}
#' with 3 columns: \code{from}, \code{to}, \code{operation}.
#' Possible operations are \code{"set"}, \code{"drop"} and \code{"reverse"}.
#'
#' \code{operation_do} and \code{operation_undo} return the modified DAG after
#' applying or undoing the operation. Note that this does not modify the
#' original DAG.
#'
#' @examples
#'
#' # We create an empty DAG with 4 nodes
#' DAG1 = create_empty_DAG(4)
#'
#' # Which kind of operations are possible?
#' operations = allowed_operations_fromDAG(DAG1)
#'
#' # We apply the first possible operation
#' op = operations[1,]
#' DAG2 = operation_do(DAG1, op)
#'
#' # and then undo it
#' DAG3 = operation_undo(DAG2, op)
#'
#' # We come back to the original DAG
#' stopifnot(identical(DAG1, DAG3))
#'
#' @export
#'
allowed_operations_fromDAG <- function(DAG){
  nodes = bnlearn::nodes(DAG)
  n.nodes = length(nodes)
  adj.mat = bnlearn::amat(DAG)

  # Create data frame to store all operations
  list_operations = list()
  i_operations = 1

  # Loop over all edges
  for (i in 1:nrow(adj.mat)){
    for (j in 1:nrow(adj.mat)){
      # We could potentially change this to a `for` loop as:
      # for (i in 1:(nrow(adj.mat)-1)){
      #   for (j in (i+1):nrow(adj.mat)){
      # to avoid the `if` statement below. This could also avoid
      # checking twice for the non-zer adjancency list. But this would
      # mean that part of th code is duplicated for (i,j) and for (j,i).
      # For the moment, let's keep it like this.

      if (i != j){
        from = nodes[i]
        to = nodes[j]

        # Addition
        if (adj.mat[i,j] == 0 && adj.mat[j,i] == 0){
          # If you set check.cycles to TRUE it stops your code I think
          DAG_new = bnlearn::set.arc(DAG, from, to, check.cycles = FALSE)
          if (is_DAG_and_restricted(DAG_new)){
            list_operations[[i_operations]] = c("from" = from, "to" = to,
                                                operation = "set")
            i_operations = i_operations + 1
          }
        } else if (adj.mat[i,j] == 1){
          # Removal
          DAG_new = bnlearn::drop.arc(DAG, from, to)
          if (is_DAG_and_restricted(DAG_new)){
            list_operations[[i_operations]] = c("from" = from, "to" = to,
                                                operation = "drop")
            i_operations = i_operations + 1
          }

          # Reversal
          DAG_new = bnlearn::reverse.arc(DAG, from, to, check.cycles = FALSE)
          if (is_DAG_and_restricted(DAG_new)){
            list_operations[[i_operations]] = c("from" = from, "to" = to,
                                                operation = "reverse")
            i_operations = i_operations + 1
          }
        }
      }
    }
  }

  # We now make it a data.frame
  df = do.call(what = rbind, args = list_operations) |>
    as.data.frame.matrix()

  return(df)
}


#' @rdname allowed_operations_fromDAG
#' @export
operation_do <- function (DAG, op){
  DAG_new = switch (
    op$operation,
    'set' = bnlearn::set.arc(DAG, op$from, op$to),
    'drop' = bnlearn::drop.arc(DAG, op$from, op$to),
    'reverse' = bnlearn::reverse.arc(DAG, op$from, op$to),
    {stop(errorCondition(
      message = paste0("Operation '", op$operation, "' is not supported. ",
                       "Possible types are: 'set', 'drop' and 'reverse'."),
      class = "UnknownOperationError"))}
  )
  return (DAG_new)
}

#' @rdname allowed_operations_fromDAG
#' @export
operation_undo <- function (DAG, op){
  DAG_new = switch (
    op$operation,
    'set' = bnlearn::drop.arc(DAG, op$from, op$to),
    'drop' = bnlearn::set.arc(DAG, op$from, op$to),
    'reverse' = bnlearn::reverse.arc(DAG, op$from, op$to),
    {stop(errorCondition(
      message = paste0("Operation '", op$operation, "' is not supported. ",
                       "Possible types are: 'set', 'drop' and 'reverse'."),
      class = "UnknownOperationError"))}
  )
  return (DAG_new)
}

