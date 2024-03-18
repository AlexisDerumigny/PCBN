#' Remove elements from a conditioning set by using conditional independence
#'
#' @param DAG Directed Acyclic Graph
#' @param node node
#' @param cond_set vector of nodes in conditioning set
#'
#' @returns a vector containing the nodes not removable from the conditioning set
#'
#'
#' @examples
#'
#' DAG = create_DAG(3)
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

#' Create empty graph
#'
#' `create_DAG()` creates a directed graph with a total of `N_nodes` nodes and no arcs.
#' The graph is a bnlearn object.
#' The nodes are named `U1`, `U2`, etc.
#'
#' @param N_nodes An integer equal to the number of nodes
#' @returns A bnlearn graph object with `N_nodes` nodes and no arcs
#' @seealso [bnlearn::empty.graph()] which this function wraps.
#' @examples
#' create_DAG(6)
#' create_DAG(10)
#' @export
create_DAG <- function(N_nodes){
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
#' DAG = create_DAG(5)
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
#'
active_cycle_check <- function(DAG, early.stopping = TRUE)
{
  node.names = bnlearn::nodes(DAG)
  adj.mat = bnlearn::amat(DAG)
  active_cycle_list = list()

  # Turn DAG into undirected graph
  DAG_igraph = igraph::as.undirected(bnlearn::as.igraph(DAG))

  for (v in node.names){
    parents = DAG$nodes[[v]]$parent
    children = DAG$nodes[[v]]$children
    # We look for active cycles of the form w -> v <- z
    # where there is a trail between w and z that does not involve v
    # such that this trail has no converging connection and no chords

    if (length(parents) > 1)
    {
      # We can only restrict ourselves to the parents up to w
      # so that we do not count the couples (w,z) and (z,w) twice.
      for (i_w in 2:length(parents))
      {
        w = parents[i_w]
        for (i_z in 1:(i_w - 1))
        {
          z = parents[i_z]

          # Parents must be non-adjacent, otherwise there is a chord directly
          if (adj.mat[w,z] == 0 & adj.mat[z,w] == 0){

            # Parents in an active cycle are joined by a trail
            # with no chords and with no converging connection
            # consisting of nodes not adjacent to v
            # (the undirected cycle has no chords)

            # This means that no node in pa(v) or ch(v) can be on the trail

            # 1: Remove all nodes in pa(v)\{w,z}, ch(v) and v
            toRemove = c(v, parents[-c(i_z, i_w)], children)
            DAG_igraph_rem = DAG_igraph - toRemove

            # Find all undirected paths between w and z
            paths = igraph::all_simple_paths(DAG_igraph_rem, w, z)

            # 2: Check for converging connections and chords
            if (length(paths) > 0){
              for (i in 1:length(paths)){
                if (path_check(DAG,paths[i])){ # Checks for v-strucs and chords
                  L = length(active_cycle_list)
                  path_vec = c(v, names(paths[[i]])) # v + path = active cycle
                  active_cycle_list[[L + 1]] = path_vec

                  if (early.stopping){
                    return (active_cycle_list = active_cycle_list)
                  }
                }
              }
            }
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
#' @param path list of nodes in DAG forming a trail.
#'
#' @returns TRUE if the path contains no converging connections and chords.
#'
path_check <- function(DAG, path){
  node.names = bnlearn::nodes(DAG)
  adj.mat = bnlearn::amat(DAG)
  no_chords_vstrucs = TRUE

  path.nodes = names(unlist(path))
  N = length(path.nodes)

  for (i in 1:N){
    # Check vstrucs
    if (i > 1 & i < N-1){
      if ((adj.mat[path.nodes[i-1], path.nodes[i]]==1) &
          (adj.mat[path.nodes[i+1], path.nodes[i]]==1)){
        no_chords_vstrucs = FALSE
      }
    }

    # Check chords
    if(i < length(path.nodes)-1){
      for(j in (i+2):N){ # Ensures that we don't check arcs twice
        if (adj.mat[path.nodes[i],path.nodes[j]] +
            adj.mat[path.nodes[j],path.nodes[i]] > 0){
          no_chords_vstrucs = FALSE
        }
      }
    }
  }
  return(no_chords_vstrucs)
}


#' Turns a general graph into a restricted graph.
#'
#' @param DAG Directed Acyclic Graph.
#'
#' @returns Restricted DAG.
#'
#' @examples
#' # DAG with an active cycle at node 5
#' DAG = create_DAG(5)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U4', 'U5')
#'
#' # Fixed graph has extra arcs 1 -> 5, 2 -> 5
#' fixed_DAG = DAG_to_restricted(DAG)
#'
#'
#' # DAG with an interfering v-structures node 3
#' DAG = create_DAG(5)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
#' DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
#'
#' # Fixed graph has extra arc 1 -> 5
#' fixed_DAG = DAG_to_restricted(DAG)
#'
#' @export
DAG_to_restricted <- function(DAG) {
  DAG_copy = DAG

  # TODO: This should be a while loop, since fixing active cycles and int_vs can
  # create more of them

  # Remove active cycles
  active_cycles = active_cycle_check(DAG)
  if (length(active_cycles) > 0) {
    # Point arcs from all nodes to the v-structure
    for (ac in active_cycles) {
      vstruc = ac[[1]]
      rest = ac[which(ac != vstruc)]
      for (node in rest) {
        DAG_copy = bnlearn::set.arc(DAG_copy, node, vstruc)
      }
    }
  }

  # Remove interfering v-structures
  res = find_B_sets(DAG)
  if (res$has_interfering_vstrucs) {
    for (node in res$nodes_with_inter_vs) {
      B_sets = res$B_sets[[node]]

      # B_sets consists of empty B-set, full B-set, at least two interfering
      # B-sets
      # Empty and full B-set never provide problems, it remains to check the others
      N_rows_B_sets = dim(B_sets)[1]
      for (i in 2:(N_rows_B_sets - 1)) {
        for (j in (i + 1):(N_rows_B_sets - 1)) {
          # B_sets[i] not in B_sets[j]
          Bset_i = B_sets[i, ]
          Bset_j = B_sets[j, ]
          increasing = all(Bset_i <= Bset_j)

          if (! increasing)
          {
            # Add arcs from all nodes in Bset_i not in Bset_j to the bq of Bset_j
            nodes_in_i_not_in_j = names(which(Bset_i > Bset_j))
            bq_j = rownames(B_sets)[j]
            for (node in nodes_in_i_not_in_j){
              DAG_copy = bnlearn::set.arc(DAG_copy, node, bq_j)
            }
          }
        }
      }
    }
  }

  return(DAG_copy)
}
