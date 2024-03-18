
###################################################
######### Creating  random PCBNs  #####################
##################################################

# Creates a PCBN with a random graph and copula assignment (to be fixed)
pick_random_PCBN <- function(N.nodes, N.arcs, familyset=c(1,3,4,5,6)){
  DAG = random_good_graph(N.nodes,N.arcs)
  order_hash = pick_random_ordering(DAG)
  node.names = bnlearn::nodes(DAG)

  # For each arc we assign a random copula
  copula_mat = pick_random_copulas(bnlearn::amat(DAG),familyset)

  return(new_PCBN(list(DAG = DAG, order_hash = order_hash, copula_mat = copula_mat)))
}


# Takes in a DAG and chooses a random order that will not require integration
pick_random_ordering <- function(DAG)
{
  order_hash = r2r::hashmap()
  node.names = bnlearn::node.ordering(DAG)

  for (v in node.names){
    parents = DAG$nodes[[v]]$parents
    if (length(parents) == 0){
      order_hash[[v]] = c()
    } else if (length(parents) == 1){
      order_hash[[v]] = parents
    } else{
      order_v = c()
      B_sets_v = find_B_sets_v(DAG = DAG, v = v)
      B_sets_v = unique(B_sets_v)
      delta_B_sets = B_sets_cut_increments(B_sets_v)

      for (i_delta_B_sets in 1:length(delta_B_sets))
      {
        delta_B_set = delta_B_sets[[i_delta_B_sets]]
        for (i in 1:length(delta_B_set))
        {
          B_minus_O = setdiff(delta_B_set, order_v)
          Poss.Cand = possible_candidates(DAG, v, order_v, order_hash, B_minus_O)

          addition = sample(Poss.Cand, 1)
          order_v = append(order_v, addition)
        }
        order_hash[[v]] = order_v
      }
    }
  }
  return(order_hash)
}


# Creates a graph with no active cycles or interfering vstrucs
random_good_graph <- function(N.nodes, N.arcs){
  DAG = create_DAG(N.nodes)
  node.names = bnlearn::nodes(DAG)

  additions = which(bnlearn::amat(DAG)+diag(N.nodes)+t(bnlearn::amat(DAG))==0, arr.ind = TRUE)
  L = length(additions[,1])
  iter = 0
  while ((L>0) & (iter<= N.arcs)) {
    k = sample(1:L, 1)
    op = list()
    op['from'] = node.names[additions[k,][1]]
    op['to'] = node.names[additions[k,][2]]
    op['operation'] = 'set'
    DAG_help = apply.operation(DAG, op)


    if (bnlearn::acyclic(DAG_help)){
      if (!(has_interfering_vstrucs(DAG_help))){
        if (length(active_cycles(DAG_help, early.stopping = TRUE)) == 0){
          DAG = DAG_help
        }
      }
    }
    additions = additions[-k,]
    L = dim(additions)[1]
    iter = iter + 1
  }
  return(DAG)
}


# Assigns a random copula to each edge
pick_random_copulas <- function(adj_mat, familyset){
  fam = adj_mat
  tau = adj_mat
  for (i in 1:nrow(adj_mat)){
    for (j in 1:ncol(adj_mat)){
      if (adj_mat[i,j]==1){
        fam[i,j] = sample(familyset,1)
        tau[i,j] = stats::runif(1,0.3,0.9)
      }
    }
  }
  return(list(fam=fam, tau=tau))
}

