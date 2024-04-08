
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
#' @param debug to print debug information
#'
#' @param e environment containing all the hashmaps
#'
#'
#' @returns DAG which locally maximizes BIC based score function
#'
#' @examples
#'
#' DAG = create_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U2", "U1")
#' order_hash[['U4']] = c("U2", "U3")
#'
#' fam = matrix(c(0, 1, 1, 1,
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
hill.climbing.PCBN <- function(data, familyset = c(1, 3, 4, 5, 6), verbose = 1,
                               start = NULL, e = NULL, score_metric = "BIC")
{
  if (is.null(start)){
    nodes = colnames(data)
    DAG = bnlearn::empty.graph(nodes)
  } else{
    DAG = start
  }
  if (is.null(e)){
    e = default_envir()
  }
  fitted = fit_all_orders(data = data, DAG = DAG, familyset = familyset, e = e,
                          score_metric = score_metric)
  reference = fitted$best_fit$metrics[score_metric]

  if (verbose > 0) {
    cat("----------------------------------------------------------------\n")
    cat("* starting from the following network:\n")
    print(DAG)
    cat("* current score:", reference, "\n")
  }

  iter = 1
  repeat{
    allowed.operations = allowed.operations.general(DAG)
    if (verbose > 0){
      cat("Number of possible operations:" , nrow(allowed.operations), "\n")
    }

    # Compute data frame with the score delta of each operation
    df = operation_score_deltas(data, DAG, familyset, allowed.operations, e = e,
                                score_metric = score_metric)

    # Select best operation based on the score
    if (score_metric == "logLik"){
      bestop = df[which.max(df$score), ]
      improvement = bestop$score - reference
    } else {
      bestop = df[which.min(df$score), ]
      improvement = reference - bestop$score
    }

    if (improvement > 0){
      if (verbose > 0){
        cat("----------------------------------------------------------------\n")
        cat("* possible operations:\n")
        print(df)
        cat("*best operation:\n")
        print(bestop)
      }

      DAG = switch (
        bestop$operation,
        'set' = bnlearn::set.arc(DAG, bestop$from, bestop$to),
        'drop' = bnlearn::drop.arc(DAG, bestop$from, bestop$to),
        'reverse' = bnlearn::reverse.arc(DAG, bestop$from, bestop$to)
      )
      reference = bestop$score

      if (verbose > 0){
        print(DAG)
      }
    }
    else{
      break;
    }
    iter = iter + 1
  }

  fitted = fit_all_orders(data = data, DAG = DAG, familyset = familyset, e = e)

  return(fitted)
}

# Computes the score delta for all allowed operations
operation_score_deltas = function(data, DAG, familyset, allowed.operations,
                                  e, score_metric)
{
  # Loop over all allowed operations
  for (i in 1:nrow(allowed.operations)){
    op = allowed.operations[i,]
    DAG_new = switch (
      op$operation,
      'set' = bnlearn::set.arc(DAG, op$from, op$to),
      'drop' = bnlearn::drop.arc(DAG, op$from, op$to),
      'reverse' = bnlearn::reverse.arc(DAG, op$from, op$to)
    )

    # Fit all possible orders
    fitted = fit_all_orders(data, DAG_new, familyset, e = e)
    score = fitted$best_fit$metrics[score_metric]

    allowed.operations$score[i] = score
  }
  return(allowed.operations)
}


#' Finds all operations resulting in a restricted DAG
#'
#' @param DAG the current DAG
#'
#' @returm a matrix with 3 columns `from`, `to`, `operation`.
#' Possible operations are "set", "drop" and "reverse".
#'
#' @noRd
allowed.operations.general <- function(DAG){
  nodes = bnlearn::nodes(DAG)
  n.nodes = length(nodes)
  adj.mat = bnlearn::amat(DAG)

  # Create data frame to store all operations
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("from", "to", "operation")

  # Loop over all edges
  for (i in 1:nrow(adj.mat)){
    for (j in 1:nrow(adj.mat)){
      if (i != j){
        from = nodes[i]
        to = nodes[j]

        # Addition
        if (adj.mat[i,j] == 0 & adj.mat[j,i] == 0){
          # If you set check.cycles to TRUE it stops your code I think
          DAG_new = bnlearn::set.arc(DAG, from, to, check.cycles = FALSE)
          if (bnlearn::acyclic(DAG_new)){
            if (!(has_interfering_vstrucs(DAG_new))){
              if (length(active_cycles(DAG_new, early.stopping = TRUE)) == 0){
                df = rbind(df, list(from = from, to = to, operation = "set"))
              }
            }
          }
        } else if (adj.mat[i,j] == 1){
          # Removal
          DAG_new = bnlearn::drop.arc(DAG, from, to)
          if (bnlearn::acyclic(DAG_new)){
            if (!(has_interfering_vstrucs(DAG_new))){
              if (length(active_cycles(DAG_new, early.stopping = TRUE)) == 0){
                df = rbind(df, list(from = from, to = to, operation = "drop"))
              }
            }
          }

          # Reversal
          DAG_new = bnlearn::reverse.arc(DAG, from, to, check.cycles = FALSE)
          if (bnlearn::acyclic(DAG_new)){
            if (!(has_interfering_vstrucs(DAG_new))){
              if (length(active_cycles(DAG_new, early.stopping = TRUE)) == 0){
                df = rbind(df, list(from = from, to = to, operation = "reverse"))
              }
            }
          }
        }
      }
    }
  }
  return(df)
}


