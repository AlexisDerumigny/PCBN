
#' Initializes PCBN class
#'
#' @param DAG the corresponding DAG (a `bn` object)
#' @param order_hash a hashmap of character vectors
#' Each character vector corresponds
#' to the order of the parents for the given node.
#'
#' @param copula_mat a list with at least two components: \itemize{
#'   \item \code{fam} the matrix of families
#'   \item \code{tau} the matrix of Kendall's tau
#' }
#' They both should be matrices of size \code{d * d},
#' where \code{d} is the number of nodes in the graph \code{DAG}.
#'
#' @param verbose If \code{0}, no message is printed.
#' If \code{1} (recommended),
#' information is printed during the checking of the PCBN.
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
#' @export
#'
new_PCBN <- function(DAG, order_hash, copula_mat, verbose = 0)
{
  # Checking the components of copula_mat
  if (!is.list(copula_mat) | length(copula_mat) < 2){
    stop("'copula_mat' must be a list with at least 2 elements. ",
         "See 'help(new_PCBN)'.")
  }
  if (! ("tau" %in% names(copula_mat)) |
      ! ("fam" %in% names(copula_mat)) ){
    stop("'copula_mat' must include at least elements named 'tau' and 'fam'.")
  }
  d = bnlearn::nnodes(DAG)

  if (!is.matrix(copula_mat$tau) |
      nrow(copula_mat$tau) != ncol(copula_mat$tau) |
      nrow(copula_mat$tau) != d){
    stop("'copula_mat$tau' should be a squared matrix of size d*d, ",
         "where d is the number of nodes of DAG.")
  }
  if (!is.matrix(copula_mat$fam) |
      nrow(copula_mat$fam) != ncol(copula_mat$fam) |
      nrow(copula_mat$fam) != d){
    stop("'copula_mat$tau' should be a squared matrix of size d*d, ",
         "where d is the number of nodes of DAG.")
  }

  # Check 'order_hash' and complete if necessary
  complete_and_check_orders(DAG, order_hash)

  # Checking the matching of colnames/rownames to the names of the nodes
  nodes_names = bnlearn::nodes(DAG)
  if (length(setdiff(colnames(copula_mat$fam), nodes_names))){
    stop("The colnames of copula_mat$fam do not match with the names in the DAG")
  }
  if (length(setdiff(rownames(copula_mat$fam), nodes_names))){
    stop("The rownames of copula_mat$fam do not match with the names in the DAG")
  }
  if (length(setdiff(colnames(copula_mat$tau), nodes_names))){
    stop("The colnames of copula_mat$tau do not match with the names in the DAG")
  }
  if (length(setdiff(rownames(copula_mat$tau), nodes_names))){
    stop("The rownames of copula_mat$tau do not match with the names in the DAG")
  }

  # Automatically add names if there is none.
  if (is.null(colnames(copula_mat$fam)) | is.null(rownames(copula_mat$fam)) ){
    if (verbose > 0){
      cat("Missing colnames and/or rownames in 'copula_mat$fam'.\n",
          "They are automatically replaced by 'bnlearn::nodes(DAG)', ",
          "assuming that the column and the rows of 'copula_mat$fam' ",
          "are in the right order.\n")
    }
    colnames(copula_mat$fam) <- rownames(copula_mat$fam) <- nodes_names
  }
  if (is.null(colnames(copula_mat$tau)) | is.null(rownames(copula_mat$tau)) ){
    if (verbose > 0){
      cat("Missing colnames and/or rownames in 'copula_mat$tau'.\n",
          "They are automatically replaced by 'bnlearn::nodes(DAG)', ",
          "assuming that the column and the rows of 'copula_mat$tau' ",
          "are in the right order.\n")
    }
    colnames(copula_mat$tau) <- rownames(copula_mat$tau) <- nodes_names
  }

  PCBN = list(DAG = DAG,
              order_hash = order_hash,
              copula_mat = copula_mat)

  class(PCBN) <- "PCBN"
  return (PCBN)
}

#' Plot function
#'
#' @param x PCBN object
#' @param ... other arguments, unused
#'
#' @export
plot.PCBN <- function(x, ...){
  plot(x$DAG)
}

#' Print function
#'
#' @param x PCBN object
#'
#' @param print.orders if \code{all}, print all orders.
#' If \code{non-empty}, this only prints the non-empty ones.
#'
#' @param ... other arguments, unused
#'
#' @export
print.PCBN <- function(x, print.orders = "non-empty", ...){
  cat("Copula matrix:\n")
  print(x$copula_mat)

  cat("Parental orderings:\n")
  well_ordering = bnlearn::node.ordering(DAG)
  # Show the ordering for each node
  for (node in well_ordering) {
    if (length(order_hash[[node]]) > 0 || print.orders == "all")
    cat(paste0(node, ":", sep = ""),
        paste0(order_hash[[node]], collapse = " < "),
        "\n")
  }
}

