#' Remove elements from a conditioning set by using conditional independence
#'
#' @param DAG Directed Acyclic Graph
#' @param node node
#' @param cond_set vector of nodes in conditioning set
#'
#' @returns a vector containing the nodes that cannot be removed from the
#' conditioning set.
#'
#'
#' @examples
#'
#' DAG = create_empty_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' remove_CondInd(DAG = DAG, node = "U1", cond_set = c("U2"))
#' remove_CondInd(DAG = DAG, node = "U3", cond_set = c("U1"))
#'
#' @export
#'
remove_CondInd <- function(DAG, node, cond_set){
  new_cond_set = cond_set
  for (i in cond_set){
    # check if dsep(node, i | cond_set\{i} )
    if (bnlearn::dsep(DAG, node,i,setdiff(cond_set,i))){
      new_cond_set = new_cond_set[-which(new_cond_set == i)]
    }

  }
  return(new_cond_set)
}

#' Create empty DAG
#'
#' This function creates a directed graph with a total of `N_nodes` nodes and
#' no arcs. The nodes are named `U1`, `U2`, etc.
#'
#' @param N_nodes An integer equal to the number of nodes
#'
#' @returns A \code{bnlearn} graph object with `N_nodes` nodes and no arcs
#'
#' @seealso [bnlearn::empty.graph()] which this function wraps.
#'
#' @examples
#' create_empty_DAG(6)
#' create_empty_DAG(10)
#'
#' @export
create_empty_DAG <- function(N_nodes){
  node.names = c()
  for (i in 1:N_nodes){
    node.names[i] = paste("U",as.character(i),sep="")
  }
  DAG = bnlearn::empty.graph(node.names)
  adj = matrix(0L, ncol = N_nodes, nrow = N_nodes,
               dimnames = list(node.names, node.names))
  bnlearn::amat(DAG) = adj
  return(DAG)
}


#' D-separation of two nodes given a set in a DAG
#'
#' @param DAG Directed Acyclic Graph
#' @param X node
#' @param Y node
#' @param Z set
#'
#' @return \code{TRUE} if the sets are d-separated and \code{FALSE} if not
#'
#' @examples
#'
#' DAG = create_empty_DAG(5)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#'
#' dsep_set(DAG, 'U1', 'U5')
#'
#' @export
#'
dsep_set <- function(DAG, X, Y, Z = NULL){
  if (length(Z) == 0){
    Z = character(0)
  }
  for (x in X){
    for (y in Y){
      if (!bnlearn::dsep(DAG, x, y, Z)){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}


#' Checks if a graph contains active cycles
#'
#' @param DAG Directed Acyclic
#' @param early.stopping if \code{TRUE}, stop at the first active cycle that is
#' found.
#'
#' @returns a list containing the active cycles.
#' Each active cycle is a character vector of the name of the nodes involved in
#' the active cycle. The first element of this vector is the converging node of
#' the active cycle.
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
