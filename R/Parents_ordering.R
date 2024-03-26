#' Finds all possible copula assignments given a DAG
#'
#' @param DAG Directed Acyclic Graph
#'
#' @returns a list of hashmaps containing the possible orders
#'
#' @examples
#' DAG = create_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#' all_orders = find_all_orders(DAG)
#' length(all_orders)
#' # 8 orders
#' for (i in 1:length(all_orders)){
#'   cat("Order ", i, ": \n")
#'   cat("U3:", all_orders[[i]][['U3']])
#'   cat(" ; U4:", all_orders[[i]][['U4']], "\n")
#' }
#'
#' @export
find_all_orders <- function(DAG) {
  # Start with empty order
  order_hash = r2r::hashmap()
  all_orders = list(order_hash)

  well_ordering = bnlearn::node.ordering(DAG)
  # Add all possible orders for node v to each order_hash
  for (node in well_ordering) {
    all_orders = extend_orders(DAG, all_orders, node)
  }
  return(all_orders)
}


#' Fills in all possible orders for the next node for each possible order
#'
#' @param DAG Directed Acyclic Graph
#' @param all_orders list of orders
#' @param node node
#'
#' @returns list of order hashmaps
#'
#' @examples
#' DAG = create_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' # Start with empty order
#' order_hash = r2r::hashmap()
#'
#' all_orders_3 = find_all_orders_v(DAG, v = "U3", order_hash = order_hash)
#' print(all_orders_3)
#'
#' # Two possible choices for node 3, let's use the first
#' order_hash[['U3']] = all_orders_3[[1]]
#'
#' extended_orders = extend_orders(DAG, list(order_hash), node = 'U4')
#' length(extended_orders)
#' # We can extend this order in 4 ways:
#' for (i in 1:length(extended_orders)){
#'   print(extended_orders[[i]][['U4']])
#' }
#' # We never pick U2 and U3 first, because their copula is not specified
#'
#' @export
extend_orders <- function(DAG, all_orders, node)
{
  parents = DAG$nodes[[node]]$parents
  # If the node has no parents then there is no order to be added.
  if (length(parents) == 0) {
    return (all_orders)
  }
  if (length(parents) == 1) {
    # If there are no orders yet, create a default one
    if (length(all_orders) == 0){
      order = r2r::hashmap()
      order[[node]] = parents
      all_orders = list(order)
      return (all_orders)
    }
    # Else, update directly all existing orders and return them
    for (order in all_orders){
      order[[node]] = parents
    }
    return (all_orders)
  }
  # We now know that `node` has at least 2 parents.

  extended = list()
  # For each order hashmap we add all possible orders for v
  for (order in all_orders) {
    copy_order = copy_hash(order)
    # If node has no parents do add NULL
    if (length(DAG$nodes[[node]]$parents) == 0) {
      copy_order[[node]] = NULL
      extended[[length(extended) + 1]] = copy_order
    } else{
      # Else add all possible orders for v to the list
      orders_node = find_all_orders_v(DAG, node, order)
      for (node_order in orders_node) {
        copy_order = copy_hash(order)
        copy_order[[node]] = node_order
        extended[[length(extended) + 1]] = copy_order
      }
    }
  }
  return(extended)
}

#' Finds all possible orders of node v given previous copula assignments
#'
#' @param DAG Directed Acyclic Graph
#' @param v node
#' @param order_hash hashmap of orders
#'
#' @returns list of vectors containing all possible orders for v
#'
#' @examples
#' DAG = create_DAG(5)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U4', 'U5')
#'
#' # Start with empty order
#' order_hash = r2r::hashmap()
#'
#' all_orders_3 = find_all_orders_v(DAG, v = "U3", order_hash = order_hash)
#' print(all_orders_3)
#'
#' # Two possible choices for node 3, let's use the first
#' order_hash[['U3']] = all_orders_3[[1]]
#'
#' all_orders_4 = find_all_orders_v(DAG, v = "U4", order_hash = order_hash)
#' print(all_orders_4)
#'
#' # Four possible choices for node 4, let's use the third
#' order_hash[['U4']] = all_orders_4[[3]]
#'
#' all_orders_5 = find_all_orders_v(DAG, v = "U5", order_hash = order_hash)
#' print(all_orders_5)
#'
#' # Eight possible orders for node 5; let's use the fourth
#' order_hash[['U5']] = all_orders_5[[4]]
#'
#'
#' @export
find_all_orders_v <- function(DAG, v, order_hash)
{
  parents = DAG$nodes[[v]]$parents
  B_sets_v = find_B_sets_v(DAG = DAG, v = v)
  B_sets_v = unique(B_sets_v)
  delta_B_sets = B_sets_cut_increments(B_sets_v)

  # Order_list contains the partial orders (starting with empty)
  order_list = list(NULL)
  for (i_delta_B_sets in 1:length(delta_B_sets)){
    delta_B_set = delta_B_sets[[i_delta_B_sets]]
    for (i in 1:length(delta_B_set)) {
      new_order_list = list()
      for (order in order_list)
      {
        B_minus_O = setdiff(delta_B_set, order)

        # Each possible candidate results in a different order
        for (w in possible_candidates(DAG, v, order, order_hash, B_minus_O)) {
          new_order = append(order, w)
          new_order_list[[length(new_order_list) + 1]] = new_order
        }
      }
      order_list = new_order_list
    }
  }
  return(order_list)
}


#' Complete an order and check whether these are valid orders on parents sets
#'
#' @param DAG the DAG
#' @param order_hash the hashmaps of orders
#'
#' @returns \code{NULL}. This function has only side-effects,
#' and modifies \code{order_hash}. It stops if the orders are not valid orders
#' on the parents sets.
#'
#'
#' @examples
#'
#' DAG = create_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' order_hash = r2r::hashmap()
#' try({complete_and_check_orders(DAG, order_hash)})
#' # Error because the order of the parents on "U3" should be specified.
#'
#' order_hash[['U3']] = c("U1", "U2")
#' complete_and_check_orders(DAG, order_hash)
#' r2r::keys(order_hash)
#' # We obtain "U3" and "U4" because they both have parents
#'
#' @export
#'
complete_and_check_orders <- function(DAG, order_hash)
{
  node.names = bnlearn::nodes(DAG)
  for (i_node in 1:length(node.names)){
    node = node.names[[i_node]]

    parents = bnlearn::parents(x = DAG, node = node)

    if (is.null(order_hash[[node]]))
    {
      if (length(parents) > 0){
        if (length(parents) == 1){
          # This is easy, there is only one parent so we add it to the hash
          order_hash[[node]] <- parents
        } else {
          stop("Order missing for node '", node, "'.\n",
               "You need to provide an order for the set of its parents, for example ",
               dputCharacterVec(sort(parents)), ".\n",
               "Remember that this order has to be compatible with all ",
               "of the other orders in 'order_hash'.")
        }
      }
    } else {
      is_valid = identical( sort(parents), sort(order_hash[[node]]) )
      if (! is_valid){
        stop("Bad set of parents for node '", node, "'\n",
             "Its parents are: ", dputCharacterVec(sort(parents)), "\n",
             "But the order given by 'order_hash[[", node, "]] is", order_hash[[node]], ".")
      }
    }
  }
}


#' Check whether a certain order abides by the B-sets
#'
#' @param DAG the considered DAG
#' @param order_hash the hashmaps of parents ordering
#'
#' @param B_sets matrix of B-sets, assumed to be increasing.
#' This can be the output of \code{\link{find_B_sets_v}}
#' or of \code{\link{B_sets_make_unique}}.
#'
#' @param orderParents a vector of characters, interpreted as the ordered
#' parents
#'
#' @return It returns `TRUE` if the order abides by the B-sets, and `FALSE` else.
#'
#' @examples
#'
#' DAG = create_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#' order_hash[['U4']] = c("U1", "U3")
#' is_order_abiding_Bsets(DAG, order_hash)
#' order_hash[['U3']] = c("U2", "U1")
#' is_order_abiding_Bsets(DAG, order_hash)
#'
#' @export
is_order_abiding_Bsets <- function(DAG, order_hash)
{
  nodes = bnlearn::nodes(DAG)
  for (i_node in 1:length(nodes)){
    node = nodes[i_node]
    B_sets = find_B_sets_v(DAG, node)
    # We test only if there are at least 2 parents
    # otherwise this is always TRUE.
    if (ncol(B_sets) >= 2){
      is_abiding_v = is_order_abiding_Bsets_v(B_sets = B_sets,
                                              orderParents = order_hash[[node]])
      if (!is_abiding_v){
        # Early stopping
        return (FALSE)
      }
    }
  }
  return (TRUE)
}

#' @rdname is_order_abiding_Bsets
#' @export
is_order_abiding_Bsets_v <- function(B_sets, orderParents)
{
  B_sets_increments = B_sets_cut_increments(B_sets)

  # We start looking at 0
  position_in_ordering = 0
  for (i in 1:length(B_sets_increments)){
    nodes = B_sets_increments[[i]]
    # We skip empty increments
    if (length(nodes) > 0){
      indices_ordering = position_in_ordering + (1:length(nodes))
      if ( ! identical(sort(nodes) ,
                       sort(orderParents[indices_ordering]) ) ){
        return (FALSE)
      }
      # We update the reference in the ordering vector
      position_in_ordering = position_in_ordering + length(nodes)
    }
  }
  return (TRUE)
}


# Nice printing of a character vector
dputCharacterVec <- function (vec){
  return (c("c(",
            paste0(paste0("'", vec, "'"), collapse = ","),
            ")"))
}

