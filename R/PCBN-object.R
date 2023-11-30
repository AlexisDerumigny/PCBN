
#' Initializes PCBN class
#'
#' @param DAG the corresponding DAG (a `bn` object)
#' @param order_hash a hashmap of character vectors
#' Each character vector corresponds
#' to the order of the parents for the given node.
#'
#' @param copula_mat list with the matrix of families
#' and the matrix of estimated parameters
#' (assuming a 1-dimensional family)
#'
#' @return the new PCBN
#'
#' @examples
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
#' tau = 0.2 * fam
#'
#' my_PCBN = new_PCBN(
#'   DAG, order_hash,
#'   copula_mat = list(tau = tau, fam = fam))
#'
#'
new_PCBN <- function(DAG, order_hash, copula_mat){
  PCBN = list(DAG = DAG,
              order_hash = order_hash,
              copula_mat = copula_mat)

  class(PCBN) <- "PCBN"
  return (PCBN)
}

# Plot function
plot.PCBN <- function(object){
  plot(object$DAG)
}

# Print function
print.PCBN <- function(object){
  object$copula_mat
}

