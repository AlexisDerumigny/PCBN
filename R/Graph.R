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
#' @return TRUE if the sets are d-separated and FALSE if not
#'
dsep_set <- function(DAG, X, Y, Z=c()){
  dseparated = TRUE
  for (x in X){
    for (y in Y){
      if (length(Z)==0){
        if (!bnlearn::dsep(DAG,x,y)){
          dseparated = FALSE
        }
      } else{
        if (!bnlearn::dsep(DAG,x,y,Z)){
          dseparated = FALSE
        }
      }
    }
  }
  return(dseparated)
}

#' Checks if graph has interfering v-structures
#'
#' @param DAG Directed Acyclic Graph
#'
#' @returns TRUE if graph contains interfering vs and FALSE if not
#'
interfering_vstrucs_check <- function(DAG){
  return(find_B_sets(DAG)$interfering_vstrucs)
}

#' Find all the B-sets of a given DAG
#'
#' @param DAG A bnlearn graph object
#' @param v node at which we want to find the B-sets
#'
#' @returns \code{find_B_sets} list with three elements \itemize{
#'    \item \code{B_sets} list of B-sets for each node
#'
#'    \item \code{interfering_vstrucs} a boolean specifying if the graph
#'    contains interfering v-structures or not
#'
#'    \item \code{nodes_with_inter_vs} a list containing nodes forming the
#'    interfering v-structures.
#' }
#' \code{find_B_sets_v} returns a boolean matrix with \code{(2 + length(children))}
#' columns and \code{length(parents)} rows.
#'
#' @export
#'
find_B_sets <- function(DAG)
{
  interfering_vstrucs = FALSE
  nodes_with_inter_vs = list()
  node.names = bnlearn::nodes(DAG)
  adj.mat = bnlearn::amat(DAG)
  B_set_list = list()

  for (v in node.names) {
    B_set = find_B_sets_v(DAG = DAG, v = v)
    B_set_list[[v]] = B_set

    # TODO: check code for interfering v-structures

    parents = DAG$nodes[[v]]$parent
    if (length(parents) > 0) {
      if (!increasing_B_set_check(B_set)) {
        interfering_vstrucs = TRUE
        nodes_with_inter_vs = append(nodes_with_inter_vs, v)
      }
    }
  }
  return( list( B_sets = B_set_list,
                interfering_vstrucs = interfering_vstrucs,
                nodes_with_inter_vs = nodes_with_inter_vs
  ) )
}

#' Find the B sets for a given node v
#'
#' @rdname find_B_sets
#' @export
find_B_sets_v <- function(DAG, v)
{
  parents = DAG$nodes[[v]]$parent
  nparents = length(parents)
  children = DAG$nodes[[v]]$children

  if (nparents == 0) {
    all_B_sets = matrix(nrow = 2 + length(children), ncol = 0)
    return (all_B_sets)
  }

  if (length(children) == 0) {
    all_B_sets = rbind(
      rep(FALSE, nparents) ,
      rep(TRUE, nparents)
    )
  } else if (length(children) == 1) {
    all_B_sets = rbind(
      rep(FALSE, nparents) ,
      parents %in% DAG$nodes[[children]]$parents,
      rep(TRUE, nparents)
    )
    return (all_B_sets)
  } else { # We know now that length(children) > 1

    # We put all the B-sets in a matrix
    all_B_sets = apply(X = children, MARGIN = 1, FUN = function(w){
      # This returns a vector of booleans of the same size as `parents`
      parents %in% DAG$nodes[[w]]$parents
    })

    # We add the trivial B-sets
    all_B_sets = rbind(
      rep(FALSE, nparents) ,
      all_B_sets,
      rep(TRUE, nparents)
    )
  }

  # # We remove duplicates
  # all_B_sets = unique(all_B_sets)

  rownames(all_B_sets) <- c("Empty B-set", children, "Full B-set")

  if (length(children) <= 1){
    return (all_B_sets)
  }

  # We sort the B-sets by size
  B_sets_sizes = rowSums(all_B_sets)
  order_B_sets_sizes = order(B_sets_sizes, decreasing = FALSE)
  all_B_sets = all_B_sets[order_B_sets_sizes , ]

  return (all_B_sets)
}

#' Checks if the B-sets for a particular node form an increasing sequence.
#'
#' @param B_set list containig nodes in a B-set.
#'
#' @returns TRUE if the list forms an ordered sequence, FALSE if not.
increasing_B_set_check <- function(B_set){
  increasing = TRUE
  if (length(B_set)>1){
    for (i in 1:(length(B_set)-1)){
      for (j in (i+1):length(B_set)){
        # B_set[i] not in B_set[j]
        if (!(sets::as.set(B_set[[i]]) < sets::as.set(B_set[[j]]))){
          increasing = FALSE
        }
      }
    }
  }
  return(increasing)
}


#' Checks if a graph contains active cycles
#'
#' @param DAG Directed Acyclic Graph
#'
#' @returns a list containing a boolean specifying if DAG contains active cycles,
#' number of active cycles, and list of the active cycles.
#'
active_cycle_check <- function(DAG){
  node.names = bnlearn::nodes(DAG)
  adj.mat = bnlearn::amat(DAG)
  active_cycles = FALSE
  active_cycle_list = list()

  for (v in node.names){
    parents = DAG$nodes[[v]]$parent
    children = DAG$nodes[[v]]$children
    if (length(parents)>1){
      for (i_parent in 1:length(parents)) {
        w = parents[i_parent]
        parents_up_to_w = if (i_parent == 1) {c()
        } else {parents[1:(i_parent - 1)]}

        if (length(parents_up_to_w)>0){
          for (z in parents_up_to_w){
            # Parents must be non-adjacent
            if (adj.mat[w,z]==0 & adj.mat[z,w]==0){
              # Parents in an active cycle are joined by a trail
              # with no chords with no converging connection
              # consisting of nodes not adjacent to v (the undirected cycle has no chords)

              # This means that no node in pa(v) or ch(v) can be on the trail

              # 1: Turn DAG into undirected graph
              DAG_igraph = igraph::as.undirected(bnlearn::as.igraph(DAG))
              # 2: Remove all nodes in pa(v)\{w,z}, ch(v) and v
              DAG_igraph = DAG_igraph - v - parents[which((parents != w) &
                                                            (parents != z))] - children
              # Find all undirected paths between w and z
              paths = igraph::all_simple_paths(DAG_igraph, w, z)
              # 3: Check for converging connections and chords
              if (length(paths)>0){
                for (i in 1:length(paths)){
                  if (path_check(DAG,paths[i])){ # Checks for v-strucs and chords
                    L = length(active_cycle_list)
                    path_vec = c(v,names(paths[[i]])) # v + path = active cycle
                    active_cycle_list[[L+1]] = path_vec
                    active_cycles = TRUE
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  res = list(active_cycles = active_cycles,
             N = length(active_cycle_list),
             active_cycle_list = active_cycle_list)
  return(res)
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
DAG_to_restricted <- function(DAG) {
  DAG_copy = DAG

  # Remove active cycles
  res = active_cycle_check(DAG)
  if (res$active_cycles) {
    acs = res$active_cycle_list
    N_ac = res$N
    # Point arcs from all nodes to the v-structure
    for (ac in acs) {
      vstruc = ac[[1]]
      rest = ac[which(ac != vstruc)]
      for (node in rest) {
        DAG_copy = bnlearn::set.arc(DAG_copy, node, vstruc)
      }
    }
  }

  # Remove interfering v-structures
  res = find_B_sets(DAG)
  if (res$interfering_vstrucs) {
    for (node in res$nodes_with_inter_vs) {
      B_sets = res$B_sets[[node]]
      for (i in 1:(length(B_sets) - 1)) {
        for (j in (i + 1):length(B_sets)) {
          # B_sets[i] not in B_sets[j]
          if (!(sets::as.set(B_sets[[i]]) < sets::as.set(B_sets[[j]])))
          {
            # Find all nodes b_q corresponding to B_sets[[j]]
            # A node b_q must be a child of all nodes in B_sets[[j]] and a child of node
            bqs = DAG$nodes[[node]]$children
            for (parent in B_sets[[j]]){
              bqs = intersect(DAG$nodes[[parent]]$children, bqs)
            }

            # Point arcs from each node in B_sets[[i]] to bqs
            for (a in B_sets[[i]]){
              for (b in bqs){
                DAG_copy = bnlearn::set.arc(DAG_copy, a, b)
              }
            }
          }
        }
      }
    }
  }

  return(DAG_copy)
}
