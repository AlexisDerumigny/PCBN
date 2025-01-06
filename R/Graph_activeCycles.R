#' Checks if a graph contains active cycles
#'
#' @param DAG Directed Acyclic
#'
#' @param early.stopping if \code{TRUE}, stop at the first active cycle that is
#' found.
#'
#' @param active_cycles_list a list of active cycles as given by
#' \code{active_cycles}. If this is \code{NULL}, the function
#' \code{active_cycles} is run on \code{DAG} to find the active cycles to be
#' displayed.
#'
#' @returns \code{active_cycles} returns a list containing the active cycles.
#' Each active cycle is a character vector of the name of the nodes involved in
#' the active cycle. The first element of this vector is the converging node of
#' the active cycle.
#'
#' \code{plot_active_cycles} is called for its side-effects only. It plots the
#' active cycles if any, and else prints a message.
#'
#' @examples
#'
#' DAG = create_empty_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' active_cycles(DAG)  # no active cycle
#'
#' DAG = create_empty_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' active_cycles(DAG)  # 1 active cycle
#'
#' DAG = create_empty_DAG(5)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
#'
#' active_cycles(DAG)  # 2 active cycles
#' active_cycles(DAG, early.stopping = TRUE)  # The first active cycle
#'
#' # Plotting the active cycles
#' plot_active_cycles(DAG)
#' # which is the same as
#' plot_active_cycles(DAG, active_cycles_list = active_cycles(DAG))
#'
#' fixedDAG = fix_active_cycles(DAG)
#' plot_active_cycles(fixedDAG)
#'
#' @export
active_cycles <- function(DAG, early.stopping = FALSE)
{
  node.names = bnlearn::nodes(DAG)
  adj.mat = bnlearn::amat(DAG)
  active_cycle_list = list()

  # Turn DAG into undirected graph
  DAG_igraph = igraph::as_undirected(bnlearn::as.igraph(DAG))

  for (v in node.names){
    parents = DAG$nodes[[v]]$parent
    children = DAG$nodes[[v]]$children
    # We look for active cycles of the form w -> v <- z
    # where there is a trail between w and z that does not involve v
    # such that this trail has no converging connection and no chords

    if (length(parents) < 2){
      next # go to next node `v`
    } # We now know v has at least two parents.

    for (i_w in 2:length(parents)) {
      w = parents[i_w]
      # We can only restrict ourselves to the parents up to w
      # so that we do not count the couples (w,z) and (z,w) twice.
      for (i_z in 1:(i_w - 1)) {
        z = parents[i_z]

        # If parents are adjacent then there is a chord directly
        if (adj.mat[w,z] == 1 || adj.mat[z,w] == 1){
          next # go to the next parent `i_z`
        } # We now know that the considered parents are not adjacent

        # Parents in an active cycle are joined by a trail with no chords and
        # with no converging connection consisting of nodes not adjacent to v
        # (this means that there are no chords in the undirected cycle).

        # Therefore, no node in pa(v) nor in ch(v) can be on the trail

        # 1: Remove all nodes in pa(v)\{w,z}, ch(v) and v
        toRemove = c(v, parents[-c(i_z, i_w)], children)
        DAG_igraph_rem = DAG_igraph - toRemove

        # Find all undirected paths between w and z
        paths = igraph::all_simple_paths(DAG_igraph_rem, w, z)
        if (length(paths) == 0){
          next # there are no path to consider
        }

        # 2: Check for converging connections and chords
        for (i in 1:length(paths)){
          names_path_nodes = names(paths[[i]])

          if (path_hasConvergingConnections(DAG, names_path_nodes) ||
              path_hasChords(DAG, names_path_nodes)){
            next # not an active cycle
          } # we now know that this is an active cycle

          # We store the new active cycle (node v, path) as a new element
          # of the list
          active_cycle_list[[length(active_cycle_list) + 1
                             ]] <- c(v, names_path_nodes)

          if (early.stopping){
            return (active_cycle_list = active_cycle_list)
          }
        }
      }
    }
  }
  return(active_cycle_list)
}


#' Checks a path for converging connections and chords.
#'
#' @param DAG Directed Acyclic Graph.
#' @param path character vector of nodes in the DAG forming a trail.
#'
#' @return \code{path_hasConvergingConnections} returns \code{TRUE}
#' if the path contains a converging connection.
#' \code{path_hasChords} returns \code{TRUE}
#' if the path contains a chord..
#'
#' @examples
#'
#' DAG = create_empty_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#' path_hasConvergingConnections(DAG, c('U1', 'U2', 'U3', 'U4') ) # FALSE
#' path_hasChords(DAG, c('U1', 'U2', 'U3', 'U4') ) # FALSE
#'
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' path_hasConvergingConnections(DAG, c('U1', 'U2', 'U3', 'U4') ) # FALSE
#' path_hasChords(DAG, c('U1', 'U2', 'U3', 'U4') ) # TRUE: has a chord
#'
#' DAG = create_empty_DAG(4)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U4', 'U3')
#' path_hasConvergingConnections(DAG, c('U1', 'U2', 'U3', 'U4') )
#' # TRUE: has a converging connection
#' path_hasChords(DAG, c('U1', 'U2', 'U3', 'U4') ) # FALSE
#'
#' @export
path_hasConvergingConnections <- function(DAG, path){
  node.names = bnlearn::nodes(DAG)
  adj.mat = bnlearn::amat(DAG)
  N = length(path)

  for (i in 2:(N - 1)){
    # Check v-structures at i
    if ((adj.mat[path[i-1], path[i]] == 1) &&
        (adj.mat[path[i+1], path[i]] == 1)){
      return(TRUE)
    }
  }
  return(FALSE)
}

#' @rdname path_hasConvergingConnections
#' @export
path_hasChords <- function(DAG, path){
  node.names = bnlearn::nodes(DAG)
  adj.mat = bnlearn::amat(DAG)
  N = length(path)

  for (i in 1:(N - 2)){
    for(j in (i+2):N){ # Ensures that we don't check arcs twice
      if (adj.mat[path[i], path[j]] ||
          adj.mat[path[j], path[i]] > 0){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}



#' @rdname active_cycles
#' @export
#'
plot_active_cycles = function(DAG, active_cycles_list = NULL){
  if (!requireNamespace("Rgraphviz", quietly = TRUE)) {
    warning("The package 'Rgraphviz' needs to be installed for this function ",
            "to work.\nYou can download it using the following command:\n",
            "install.packages('BiocManager')\nBiocManager::install()\n",
            "BiocManager::install('Rgraphviz')")
    return ()
  }
  if (is.null(active_cycles_list)){
    active_cycles_list = active_cycles(DAG)
  }

  if (length(active_cycles_list)==0){
    cat("The list of active cycles is empty. \n")
    return ()
  }

  adj.mat = bnlearn::amat(DAG)
  L = length(active_cycles_list)
  for (i in 1:L){
    if (i > 1){
      more.plots <- readline(prompt = "Plot next active cycle? (Y/N): ")
      if (substring(more.plots, 1, 1) %in% c("N", "n")){
        break
      }
    }

    cat("Plotting active cycle ", i, "of", L, "\n")
    active_cycle = active_cycles_list[[i]]

    # Graphviz requires a dataframe of the arcs to highlight them
    # So, vector active_cycle -> dataframe of arcs along this active cycle df
    df <- data.frame(matrix(ncol = 2, nrow = 0))
    for (j in 1:(length(active_cycle))){
      node1 = active_cycle[j]

      if (j<length(active_cycle)){
        node2 = active_cycle[j+1]
      } else{
        node2 = active_cycle[1]
      }

      if (adj.mat[node1,node2]==1){ # node1 -> node2
        df = rbind(df, data.frame(list(from=node1, to=node2)))
      }
      if (adj.mat[node2,node1]==1){ # node2 -> node1
        df = rbind(df, data.frame(list(from=node2, to=node1)))
      }
    }
    bnlearn::graphviz.plot(DAG, highlight = list(arcs = df, col = "red", lwd = 3))
  }
}


