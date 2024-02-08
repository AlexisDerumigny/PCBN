
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
pick_random_ordering <- function(DAG){
  assign("order_hash", r2r::hashmap(), envir = .GlobalEnv)
  node.names = bnlearn::node.ordering(DAG)
  B_sets = find_B_sets(DAG)$B_sets

  for (v in node.names){
    parents = DAG$nodes[[v]]$parents
    if (length(parents)==0){
      order_hash[[v]] = c()
    } else if (length(parents)==1){
      order_hash[[v]] = c(parents[1])
    } else{
      order_v = c()
      B_sets_v = B_sets[[v]]
      n.parents = length(parents)

      while (length(order_v)< n.parents){
        B_minus_O = find_B_minus_O(B_sets_v, order_v)

        if (length(order_v) == 0){
          Poss.Cand = B_minus_O
        } else{
          Poss.Cand = possible_candidates(DAG, v, order_v, order_hash, B_minus_O)
        }
        addition = sample(Poss.Cand, 1)
        order_v = append(order_v, addition)
      }
      order_hash[[v]] = order_v
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
      if (!(interfering_vstrucs_check(DAG_help))){
        if (!(active_cycle_check(DAG_help)$active_cycles)){
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

