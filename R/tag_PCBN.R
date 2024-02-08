
#' Creates copula tags
#'
#' @examples
#'
#' DAG = create_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#'
#' create_copula_tag(DAG = DAG, order_hash = order_hash,
#'                   w = "U1", v = "U2", cond_set=NULL)
#'
#' create_copula_tag(DAG = DAG, order_hash = order_hash,
#'                   w = "U2", v = "U3", cond_set="U1")
#'
#' @export
#'
create_copula_tag <- function(DAG, order_hash, w, v, cond_set=NULL){
  tag_w_cond_set = create_margin_tag(DAG, order_hash, w, cond_set)
  tag_v_cond_set = create_margin_tag(DAG, order_hash, v, cond_set)
  tag = paste0(tag_w_cond_set, "_and_", tag_v_cond_set)
  return(tag)
}

#' Creates a tag for a (conditional) margin
#'
#' For a conditional margin u_{v|K} the key cannot only be of the form "v|K".
#' We also require that the string represents the recursion of h-functions
#' with which the margin has been computed.
#'
#' @examples
#'
#' DAG = create_DAG(3)
#' DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
#' DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
#'
#' order_hash = r2r::hashmap()
#' order_hash[['U3']] = c("U1", "U2")
#'
#' create_margin_tag(DAG = DAG, order_hash = order_hash,
#'                   v = "U2", cond_set=NULL)
#'
#' create_margin_tag(DAG = DAG, order_hash = order_hash,
#'                   v = "U3", cond_set="U1")
#'
#' @export
#'
create_margin_tag <- function(DAG, order_hash, v, cond_set){
  g = igraph::make_empty_graph()
  cond_set = remove_CondInd(DAG, v, cond_set)
  top_node = paste_margin(v, cond_set)

  if (length(cond_set)==0){
    return(v)
  }

  # Find specified copula
  for (w in cond_set){
    cond_set_minus_w = cond_set[!cond_set==w]
    if (is_cond_copula_specified(DAG, order_hash, w, v, cond_set_minus_w)){
      break
    }
  }

  # Initialize first branches of the tree
  nodew = paste_margin(w, cond_set_minus_w)
  nodev = paste_margin(v, cond_set_minus_w)

  g = g + igraph::vertices(top_node, nodew, nodev)
  g = g |> igraph::add_edges(c(top_node, nodew, top_node, nodev))

  tree_size = length(cond_set_minus_w)
  if (tree_size>0){
    for (j in 1:(tree_size)){
      # Each conditional margin x is again computed with two margins y and z
      # therefore we add them as nodes and point arcs from x to y and z
      g = extend_margin_tree(DAG, order_hash, top_node, g)
    }
  }
  tag = tree_to_string(g, top_node)
  return(tag)
}


# Builds a tree representing the recursions needed to compute a conditional margin
# To each leaf node (a margin) we add the two margins needed to compute it
extend_margin_tree <- function(DAG, order_hash, top_node, g){
  leaf_nodes = get_leaf_nodes(g, top_node)
  for (leaf in leaf_nodes){
    untag = unpaste_margin(leaf)
    v = untag$v
    cond_set = untag$cond_set

    # Find specified copula
    for (w in cond_set){
      cond_set_minus_w = cond_set[!cond_set==w]
      if (is_cond_copula_specified(DAG, order_hash, w, v, cond_set_minus_w)){
        break
      }
    }

    nodew = paste_margin(w, cond_set_minus_w)
    nodev = paste_margin(v, cond_set_minus_w)

    # Nodes may already be added, in particular unconditional margins
    if (sum(igraph::as_ids(igraph::V(g))==nodew)==0){
      g = g + igraph::vertices(nodew)
    }
    if (sum(igraph::as_ids(igraph::V(g))==nodev)==0){
      g = g + igraph::vertices(nodev)
    }
    g = g + igraph::edges(leaf, nodew, leaf, nodev)

  }
  return(g)
}

# Function which finds leaf nodes (I copied this from internet)
get_leaf_nodes <- function(graph, node="B"){
  path <- igraph::ego(graph, order=length(igraph::V(graph)), nodes=node, mode="out")
  nms <- names(path[[1]])
  nms[igraph::ego_size(graph, order=1, nodes=nms, mode="out", mindist=1) == 0]
}

# A function which turns a tree of margins into a unique string
tree_to_string <- function(g, top_node){
  tree_string = top_node

  continue = TRUE
  # Iteratively loop over the levels of the tree starting from the highest level;
  # the top node
  level_nodes = c(top_node)
  while(continue){
    new_level = c()
    for (node in level_nodes){
      children = sort(igraph::as_ids(igraph::neighbors(g, node, mode="out")))
      tree_string = paste0(tree_string,"-", children[1], "-", children[2])
      new_level = append(new_level, children)
    }
    level_nodes = new_level

    # Check if have reached the last level
    if (length(igraph::as_ids(igraph::neighbors(g, level_nodes[1], mode="out")))==0){
      continue = FALSE
    }
  }
  return(tree_string)
}


# Creates a simple string node|cond_set
paste_margin <- function(node, cond_set)
{
  if (length(cond_set) == 0){
    return (node)
  }
  tag = paste0("l_", node, " | ", paste(cond_set, collapse=','), "_r")

  return(tag)
}

# Turns string node|cond_set back to its vector notation
unpaste_margin <- function(tag){
  v = strsplit(tag,split="|", fixed = TRUE)[[1]][1]
  cond_set = strsplit(tag,split="|", fixed = TRUE)[[1]][2]
  cond_set = strsplit(cond_set,split=",", fixed = TRUE)[[1]]
  return(list(v = v, cond_set=cond_set))
}

