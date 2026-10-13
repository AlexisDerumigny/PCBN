.caugi_available <- function() {
  requireNamespace("caugi", quietly = TRUE) &&
    getNamespaceVersion("caugi") >= "1.3.0"
}

.require_caugi <- function() {
  if (!.caugi_available()) {
    stop(
      "The caugi backend requires 'caugi' >= 1.3.0. ",
      "Install it with install.packages('caugi').",
      call. = FALSE
    )
  }
}

.graph_backend <- function(DAG = NULL, graph_backend = NULL) {
  if (!is.null(graph_backend)) {
    if (
      !is.character(graph_backend) ||
        length(graph_backend) != 1L ||
        is.na(graph_backend) ||
        !graph_backend %in% c("bnlearn", "caugi")
    ) {
      stop("'graph_backend' must be 'bnlearn' or 'caugi'.", call. = FALSE)
    }
  }
  if (is.null(DAG) || inherits(DAG, "bn")) {
    input_backend = "bnlearn"
  } else if (
    requireNamespace("caugi", quietly = TRUE) && caugi::is_caugi(DAG)
  ) {
    input_backend = "caugi"
  } else {
    stop("'DAG' must be a bnlearn or caugi graph.", call. = FALSE)
  }
  if (is.null(graph_backend)) {
    return(input_backend)
  }
  return(graph_backend)
}

.prepare_graph <- function(DAG, graph_backend = NULL) {
  backend = .graph_backend(DAG, graph_backend)
  if (backend == "bnlearn" && inherits(DAG, "bn")) {
    return(DAG)
  }

  .require_caugi()
  if (inherits(DAG, "bn")) {
    return(caugi::as_caugi(DAG, class = "DAG"))
  }
  if (!caugi::is_dag(DAG)) {
    stop(
      "'DAG' must contain only directed edges and be acyclic.",
      call. = FALSE
    )
  }
  if (DAG@graph_class != "DAG") {
    DAG = caugi::mutate_caugi(DAG, class = "DAG")
  }
  if (backend == "bnlearn") {
    return(caugi::as_bnlearn(DAG))
  }
  return(DAG)
}

.graph_nodes <- function(DAG) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::nodes(DAG))
  }
  return(caugi::nodes(DAG)$name)
}

.graph_parents <- function(DAG, node) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::parents(DAG, node))
  }
  return(as.character(caugi::parents(DAG, nodes = node)))
}

.graph_children <- function(DAG, node) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::children(DAG, node))
  }
  return(as.character(caugi::children(DAG, nodes = node)))
}

.graph_adjacency <- function(DAG) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::amat(DAG))
  }
  return(caugi::as_adjacency(DAG))
}

.graph_node_ordering <- function(DAG) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::node.ordering(DAG))
  }
  return(caugi::topological_sort(DAG))
}

.graph_acyclic <- function(DAG) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::acyclic(DAG))
  }
  return(caugi::is_acyclic(DAG))
}

.graph_dsep <- function(DAG, X, Y, Z = character(0)) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::dsep(DAG, X, Y, Z))
  }
  return(caugi::d_separated(DAG, X, Y, Z))
}

#' Add a directed arc to a graph
#'
#' Adds an arc from one node to another, using the graph's backend by default.
#' An existing arc in the same direction is left unchanged. Both nodes must
#' already exist in the graph. Arcs that would create a directed cycle are
#' rejected.
#'
#' @param DAG A bnlearn or caugi directed acyclic graph.
#' @param from Name of the node where the arc starts.
#' @param to Name of the node where the arc ends.
#' @inheritParams create_empty_DAG
#' @return The updated graph in the resolved backend.
#' @examples
#' DAG = create_empty_DAG(3)
#' DAG = set_arc(DAG, "U1", "U3")
#' DAG = set_arc(DAG, "U2", "U3")
#'
#' if (requireNamespace("caugi", quietly = TRUE)) {
#'   DAG = create_empty_DAG(3, graph_backend = "caugi")
#'   DAG = set_arc(DAG, "U1", "U3")
#'   DAG = set_arc(DAG, "U2", "U3")
#' }
#' @export
set_arc <- function(DAG, from, to, graph_backend = NULL) {
  DAG = .prepare_graph(DAG, graph_backend)
  nodes = .graph_nodes(DAG)

  for (node in list(from, to)) {
    if (!is.character(node) || length(node) != 1L || is.na(node) ||
        !node %in% nodes) {
      stop("'from' and 'to' must each name an existing node.", call. = FALSE)
    }
  }

  return(.graph_set_arc(DAG, from, to))
}

.graph_set_arc <- function(DAG, from, to) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::set.arc(DAG, from, to))
  }

  if (from %in% .graph_parents(DAG, to)) {
    return(DAG)
  }

  return(caugi::add_edges(DAG, from = from, edge = "-->", to = to))
}

.graph_as_igraph <- function(DAG) {
  if (inherits(DAG, "bn")) {
    return(bnlearn::as.igraph(DAG))
  }

  return(caugi::as_igraph(DAG))
}

.graph_plot <- function(DAG, ...) {
  if (inherits(DAG, "bn")) {
    plot(DAG, ...)
  } else {
    # Caugi constructs a plot object; printing it actually draws the graph.
    print(plot(DAG, ...))
  }

  invisible(NULL)
}

.graph_plot_highlight <- function(DAG, arcs) {
  if (inherits(DAG, "bn")) {
    bnlearn::graphviz.plot(
      DAG,
      highlight = list(arcs = arcs, col = "red", lwd = 3)
    )
  } else {
    by_edge = list()
    for (i in seq_len(nrow(arcs))) {
      by_edge[[arcs$from[i]]][[arcs$to[i]]] = list(col = "red", lwd = 3)
    }
    .graph_plot(DAG, edge_style = list(by_edge = by_edge))
  }
  invisible(NULL)
}
