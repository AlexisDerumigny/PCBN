
######################################################
###########     Hill climbing   ###################
###################################################

#' Hill climbing algorithm for restricted PCBNs
#'
#' @param data data frame
#' @param start starting Directed Acyclic Graph
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
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#' mydata = sample_PCBN(my_PCBN, N = 5)
#'
#' # Does not work yet
#' # TODO: fix storing of the trees in the hash thing
#' #
#' # result = hill.climbing.PCBN(data = mydata, start = create_DAG(4),
#' #                             familyset = 1, debug=FALSE)
#'
#' @export
#'
hill.climbing.PCBN <- function(data, start, familyset, debug=FALSE,
                               copula_hash = r2r::hashmap(),
                               margin_hash = r2r::hashmap())
{

  nodes = names(data)
  n.nodes = length(nodes)
  adj.mat = bnlearn::amat(start)
  nparents = colSums(adj.mat)
  iter = 1
  DAG = start
  fitted = fit_all_orders(data, DAG)
  reference = fitted$best_fit$metrics$BIC

  if (debug) {
    cat("----------------------------------------------------------------\n")
    cat("* starting from the following network:\n")
    print(DAG)
    cat("* current score:", reference, "\n")
  }

  repeat{
    if (iter>1){
      cat("----------------------------------------------------------------\n")
      cat("* current network:\n")
      print(DAG)
      cat("* current score:", reference, "\n")
    }

    allowed.operations = allowed.operations.general(DAG)

    # Compute data frame with the score delta of each operation
    df = operation_score_deltas(data, DAG, familyset, reference, allowed.operations,
                                copula_hash = copula_hash,
                                margin_hash = margin_hash)

    ## Select best operation based on column order
    bestop = df[which.max(df$score.delta),]

    ## Select best operation at random (function slice_max from tidyverse library)
    # maxima = slice_max(df, order_by = score.delta)
    # bestop = sample_n(maxima, 1)

    if (bestop$improve){
      if (debug){
        cat("----------------------------------------------------------------\n")
        cat("* possible operations:\n")
        print(df)
        cat("*best operation:\n")
        print(bestop)
      }

      # There is a function for this below
      if (bestop$operation == 'set'){
        DAG = bnlearn::set.arc(DAG, bestop$from, bestop$to)
      }
      if (bestop$operation == 'drop'){
        DAG = bnlearn::drop.arc(DAG, bestop$from, bestop$to)
      }
      if (bestop$operation == 'reverse'){
        DAG = bnlearn::reverse.arc(DAG, bestop$from, bestop$to)
      }
      reference = reference + bestop$score.delta
    }
    else{
      break;
    }
    iter = iter + 1
  }
  return(DAG)
}

# Computes the score delta for all allowed operations
operation_score_deltas = function(data, DAG, familyset, reference, allowed.operations,
                                  copula_hash, margin_hash){
  nodes = names(data)
  n.nodes = length(nodes)
  adj.mat = bnlearn::amat(DAG)

  # Create dataframe to store all operations
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(df) <- c("from", "to", "operation", "score.delta", "improve")

  # Loop over all allowed operations
  for (i in 1:nrow(allowed.operations)){
    op = allowed.operations[i,]
    DAG_new = apply.operation(DAG, op)

    # Fit all possible orders
    fitted = fit_all_orders(data, DAG_new, familyset,
                            copula_hash = copula_hash,
                            margin_hash = margin_hash)
    score.delta = fitted$best_fit$metrics$BIC - reference


    if (score.delta>0){improve = TRUE}
    else{improve=FALSE}
    df = rbind(df, data.frame(list(from=op$from, to=op$to, operation=op$operation,
                                   score.delta=score.delta, improve=improve)))
  }
  return(df)
}

# Applies operation op to DAG
apply.operation <- function(DAG, op){
  if (op$operation == 'set'){
    DAG_new = bnlearn::set.arc(DAG, op$from, op$to, check.cycles = FALSE)
  }
  if (op$operation == 'drop'){
    DAG_new = bnlearn::drop.arc(DAG, op$from, op$to)
  }
  if (op$operation == 'reverse'){
    DAG_new = bnlearn::reverse.arc(DAG, op$from, op$to, check.cycles = FALSE)
  }
  return(DAG_new)
}

# Finds all operations resulting in a restricted DAG
allowed.operations.general <- function(DAG){
  nodes = bnlearn::nodes(DAG)
  n.nodes = length(nodes)
  adj.mat = bnlearn::amat(DAG)

  # Create data frame to store all operations
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(df) <- c("from", "to", "operation")

  # Loop over all edges
  for (i in 1:nrow(adj.mat)){
    for (j in 1:nrow(adj.mat)){
      if (i!=j){
        from = nodes[i]
        to = nodes[j]

        # Addition
        if (adj.mat[i,j]==0 & adj.mat[j,i]==0){
          # If you set check.cycles to TRUE it stops your code I think
          DAG_new = bnlearn::set.arc(DAG, from, to, check.cycles = FALSE)
          if (bnlearn::acyclic(DAG_new)){
            if (!(interfering_vstrucs_check(DAG_new))){
              if (!(active_cycle_check(DAG_new)$active_cycles)){
                df = rbind(df, list(from = from, to = to, operation = "set"))
              }
            }
          }
        }
        # Removal
        if (adj.mat[i,j]==1){
          DAG_new = bnlearn::drop.arc(DAG, from, to)
          if (bnlearn::acyclic(DAG_new)){
            if (!(interfering_vstrucs_check(DAG_new))){
              if (!(active_cycle_check(DAG_new)$active_cycles)){
                df = rbind(df, list(from = from, to = to, operation = "drop"))
              }
            }
          }
        }
        # Reversal
        if (adj.mat[i,j]==1){
          DAG_new = bnlearn::reverse.arc(DAG, from, to, check.cycles = FALSE)
          if (bnlearn::acyclic(DAG_new)){
            if (!(interfering_vstrucs_check(DAG_new))){
              if (!(active_cycle_check(DAG_new)$active_cycles)){
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


