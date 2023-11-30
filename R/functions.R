# Creates copula tags
create_copula_tag <- function(DAG, order_hash, w, v, cond_set=NULL){
  tag_w_cond_set = create_margin_tag(DAG, order_hash, w, cond_set)
  tag_v_cond_set = create_margin_tag(DAG, order_hash, v, cond_set)
  tag = paste0(tag_w_cond_set, tag_v_cond_set)
  return(tag)
}


# For a conditional margin u_{v|K} the key cannot only be of the form "v|K"
# We also require that the string represents the recursion of h-functions with which
# The margin has been computed
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
    if (copula_is_specified(DAG, order_hash, w, v, cond_set_minus_w)){
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
      if (copula_is_specified(DAG, order_hash, w, v, cond_set_minus_w)){
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
  # Iteratively loop over the levels of the tree starting from the highest level; the top node
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


######################################################
###########     Hill climbing   ###################
###################################################

#' Hill climbing algorithm for restricted PCBNs
#'
#' @param data data frame
#' @param start starting Directed Acyclic Graph
#'
#' @returns DAG which locally maximizes BIC based score function
hill.climbing.PCBN <- function(data, start, familyset, debug=FALSE){
  assign("copula_hash", r2r::hashmap(), envir = .GlobalEnv)
  assign("margin_hash", r2r::hashmap(), envir = .GlobalEnv)

  nodes = names(data)
  n.nodes = length(nodes)
  adj.mat = bnlearn::amat(start)
  nparents = colSums(adj.mat)
  iter = 1
  DAG = start
  fitted = fit_all_orders(data, DAG)
  reference = fitted$best_fit$metrics$BIC

  if (debug) {
    cat("----------------------------------------------------------------\n")
    cat("* starting from the following network:\n")
    print(DAG)
    cat("* current score:", reference, "\n")
  }

  repeat{
    if (iter>1){
      cat("----------------------------------------------------------------\n")
      cat("* current network:\n")
      print(DAG)
      cat("* current score:", reference, "\n")
    }

    allowed.operations = allowed.operations.general(DAG)
    # Compute data frame with the score delta of each operation
    df = operation_score_deltas(data, DAG, familyset, reference, allowed.operations)

    ## Select best operation based on column order
    bestop = df[which.max(df$score.delta),]

    ## Select best operation at random (function slice_max from tidyverse library)
    # maxima = slice_max(df, order_by = score.delta)
    # bestop = sample_n(maxima, 1)

    if (bestop$improve){
      if (debug){
        cat("----------------------------------------------------------------\n")
        cat("* possible operations:\n")
        print(df)
        cat("*best operation:\n")
        print(bestop)
      }

      # There is a function for this below
      if (bestop$operation == 'set'){
        DAG = bnlearn::set.arc(DAG, bestop$from, bestop$to)
      }
      if (bestop$operation == 'drop'){
        DAG = bnlearn::drop.arc(DAG, bestop$from, bestop$to)
      }
      if (bestop$operation == 'reverse'){
        DAG = bnlearn::reverse.arc(DAG, bestop$from, bestop$to)
      }
      reference = reference + bestop$score.delta
    }
    else{
      break;
    }
    iter = iter + 1
  }
  return(DAG)
}

# Computes the score delta for all allowed operations
operation_score_deltas = function(data, DAG, familyset, reference, allowed.operations){
  nodes = names(data)
  n.nodes = length(nodes)
  adj.mat = bnlearn::amat(DAG)

  # Create dataframe to store all operations
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(df) <- c("from", "to", "operation", "score.delta", "improve")

  # Loop over all allowed operations
  for (i in 1:nrow(allowed.operations)){
    op = allowed.operations[i,]
    DAG_new = apply.operation(DAG, op)

    # Fit all possible orders
    fitted = fit_all_orders(data, DAG_new, familyset, reuse_hash = TRUE)
    score.delta = fitted$best_fit$metrics$BIC - reference


    if (score.delta>0){improve = TRUE}
    else{improve=FALSE}
    df = rbind(df, data.frame(list(from=op$from, to=op$to, operation=op$operation, score.delta=score.delta, improve=improve)))
  }
  return(df)
}

# Applies operation op to DAG
apply.operation <- function(DAG, op){
  if (op$operation == 'set'){
    DAG_new = bnlearn::set.arc(DAG, op$from, op$to, check.cycles = FALSE)
  }
  if (op$operation == 'drop'){
    DAG_new = bnlearn::drop.arc(DAG, op$from, op$to)
  }
  if (op$operation == 'reverse'){
    DAG_new = bnlearn::reverse.arc(DAG, op$from, op$to, check.cycles = FALSE)
  }
  return(DAG_new)
}

# Finds all operations resulting in a restricted DAG
allowed.operations.general <- function(DAG){
  nodes = bnlearn::nodes(DAG)
  n.nodes = length(nodes)
  adj.mat = bnlearn::amat(DAG)

  # Create data frame to store all operations
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(df) <- c("from", "to", "operation")

  # Loop over all edges
  for (i in 1:nrow(adj.mat)){
    for (j in 1:nrow(adj.mat)){
      if (i!=j){
        from = nodes[i]
        to = nodes[j]

        # Addition
        if (adj.mat[i,j]==0 & adj.mat[j,i]==0){
          # If you set check.cycles to TRUE it stops your code I think
          DAG_new = bnlearn::set.arc(DAG, from, to, check.cycles = FALSE)
          if (bnlearn::acyclic(DAG_new)){
            if (!(interfering_vstrucs_check(DAG_new))){
              if (!(active_cycle_check(DAG_new)$active_cycles)){
                df = rbind(df, list(from = from, to = to, operation = "set"))
              }
            }
          }
        }
        # Removal
        if (adj.mat[i,j]==1){
          DAG_new = bnlearn::drop.arc(DAG, from, to)
          if (bnlearn::acyclic(DAG_new)){
            if (!(interfering_vstrucs_check(DAG_new))){
              if (!(active_cycle_check(DAG_new)$active_cycles)){
                df = rbind(df, list(from = from, to = to, operation = "drop"))
              }
            }
          }
        }
        # Reversal
        if (adj.mat[i,j]==1){
          DAG_new = bnlearn::reverse.arc(DAG, from, to, check.cycles = FALSE)
          if (bnlearn::acyclic(DAG_new)){
            if (!(interfering_vstrucs_check(DAG_new))){
              if (!(active_cycle_check(DAG_new)$active_cycles)){
                df = rbind(df, list(from = from, to = to, operation = "reverse"))
              }
            }
          }
        }
      }
    }
  }
  return(df)
}














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
        tau[i,j] = runif(1,0.3,0.9)
      }
    }
  }
  return(list(fam=fam, tau=tau))
}




####################################################
########## Computes density/CDF/log-lik of PCBN ####
####################################################

# Computes the marginal part of the log-likelihood
logLik_margins_PCBN <- function(data_normal, margins){
  log_lik = 0
  nodes = colnames(data_normal)

  for (v in nodes) {
    mean_v = margins[[v]]$mean
    sd_v = margins[[v]]$sd
    log_lik = log_lik + sum(log(stats::dnorm(data_normal[[v]], mean = mean_v, sd = sd_v)))
  }
  return(log_lik)
}

# Computes the copula part of the log-likelihood
logLik_copulas_PCBN <- function(data_uniform, PCBN){
  # Unpack PCBN object
  DAG = PCBN$DAG
  order_hash = PCBN$order_hash
  copula_mat = PCBN$copula_mat

  log_lik = 0

  well_ordering = bnlearn::node.ordering(DAG)
  # For every node v for every parent w, we need to compute the density of arc c_{ wv|pa(v; w)\{w} }
  for (v in well_ordering) {
    parents = order_hash[[v]]
    if (length(parents) > 0) {
      for (w in parents) {
        fam = copula_mat$fam[w, v]
        tau = copula_mat$tau[w, v]
        par = VineCopula::BiCopTau2Par(fam, tau)

        # Compute the required margins
        lower = parents[0:(which(parents == w) - 1)]
        v_given_lower = compute_sample_margin(PCBN, data_uniform, v, lower)
        w_given_lower = compute_sample_margin(PCBN, data_uniform, w, lower)

        log_lik_arc_w_to_v = sum(log(VineCopula::BiCopPDF(w_given_lower, v_given_lower, family = fam, par = par)))

        log_lik = log_lik + log_lik_arc_w_to_v
      }
    }
  }
  return(log_lik)
}

# Compute full log-likelihood of PCBN with margins
logLik_PCBN <- function(data, PCBN, margins){
  data_uniform = to_uniform_scale(data)
  log_lik = logLik_margins_PCBN(data, margins) + logLik_copulas_PCBN(data_uniform, PCBN)
  return(log_lik)
}


###################################################
########  Performance metrics #####################
##################################################

# Computes several performance metrics of the GBN:
# logLik, AIC, BIC, KL divergence and CvM distance
metrics_GBN <- function(data, PCBN, margins, GBN_fit, N_monte_carlo){
  metrics =list()

  log_lik = logLik(GBN_fit, data)
  AIC = -2*BIC(GBN_fit, data)
  BIC = -2*AIC(GBN_fit, data)


  # Generate data for KL and CvM
  data_help_uniform = PCBN:::sample_PCBN(PCBN, N_monte_carlo)
  data_help = PCBN:::to_normal_scale(data_help_uniform)

  KL = PCBN:::KL_divergence_GBN(data_help, PCBN, margins, GBN_fit)

  return(list(logLik = log_lik, AIC = AIC, BIC = BIC, KL = KL))
}


# Computes several performance metrics of the PCBN:
# logLik, AIC, BIC, KL divergence and CvM distance
metrics_PCBN <- function(data, PCBN, margins, PCBN_fit, margins_fit, N_monte_carlo){
  metrics =list()

  log_lik = PCBN_fit$metrics$logLik
  AIC = -PCBN_fit$metrics$AIC
  BIC = -PCBN_fit$metrics$BIC

  # To compute the logLik, AIC and BIC we need to add the information from the margins
  log_lik_margins = 0
  nodes = colnames(data)
  for (v in nodes) {
    mean_v = margins[[v]]$mean
    sd_v = margins[[v]]$sd
    log_lik_margins = log_lik_margins + sum( log( dnorm(data[[v]], mean = mean_v, sd = sd_v) ) )
  }

  log_lik = log_lik + log_lik_margins
  AIC = AIC + (2*(2 * length(margins_fit)) - 2*log_lik_margins)
  BIC = BIC + (log(length(data))*(2 * length(margins_fit)) - 2*log_lik_margins)


  # Generate data for KL and CvM
  data_help_uniform = PCBN:::sample_PCBN(PCBN, N_monte_carlo)
  data_help = PCBN:::to_normal_scale(data_help_uniform)

  KL = PCBN:::KL_divergence_PCBN(data_help, PCBN, margins, PCBN_fit, margins_fit)

  return(list(logLik = log_lik, AIC = AIC, BIC = BIC, KL = KL))
}


### Computing the marginal part of the PDF of a PCBN given normal data
PDF_margins_PCBN <- function(data_normal, margins){
  likelihood = 1
  nodes = colnames(data_normal)

  for (v in nodes) {
    mean_v = margins[[v]]$mean
    sd_v = margins[[v]]$sd
    likelihood = likelihood * prod(stats::dnorm(data_normal[[v]], mean = mean_v, sd = sd_v))
  }
  return(likelihood)
}

### Computing the copula part of the PDF of a PCBN given uniform data
PDF_copulas_PCBN <- function(data_uniform, PCBN){
  # Unpack PCBN object
  DAG = PCBN$DAG
  order_hash = PCBN$order_hash
  copula_mat = PCBN$copula_mat

  density = 1

  well_ordering = bnlearn::node.ordering(DAG)
  # For every node v for every parent w, we need to compute the density of arc c_{ wv|pa(v; w)\{w} }
  for (v in well_ordering) {
    parents = order_hash[[v]]
    if (length(parents) > 0) {
      for (w in parents) {
        fam = copula_mat$fam[w, v]
        tau = copula_mat$tau[w, v]
        par = VineCopula::BiCopTau2Par(fam, tau)

        # Compute the required margins
        lower = parents[0:(which(parents == w) - 1)]
        v_given_lower = compute_sample_margin(PCBN, data_uniform, v, lower)
        w_given_lower = compute_sample_margin(PCBN, data_uniform, w, lower)

        density_arc_w_to_v = prod(VineCopula::BiCopPDF(w_given_lower, v_given_lower, family = fam, par = par))

        density = density * density_arc_w_to_v
      }
    }
  }
  return(density)
}


### Computes the Kullback-Leibler divergence of a estimated PCBN to the true PCBN
KL_divergence_PCBN <-
  function(data, PCBN, margins, PCBN_fit, margins_fit) {
    KL = 0
    L = length(data[[1]])
    KL_vec = rep(0, L)
    data_uniform = PCBN:::to_uniform_scale(data)

    log_lik_true = PCBN:::logLik_PCBN(data, PCBN, margins)
    log_lik_est = PCBN:::logLik_PCBN(data, PCBN_fit, margins_fit)
    KL = (log_lik_true - log_lik_est)/L
    return(KL)
  }

### Computes the Kullback-Leibler divergence of a estimated GBN to the true PCBN
KL_divergence_GBN <- function(data, PCBN, margins, GBN) {
  KL = 0
  L = length(data[[1]])
  KL_vec = rep(0, L)
  data_uniform = PCBN:::to_uniform_scale(data)

  log_lik_true = PCBN:::logLik_PCBN(data, PCBN, margins)
  log_lik_est = logLik(GBN, data)
  KL = (log_lik_true - log_lik_est)/L
  return(KL)
}

# Returns the distance between two graphs
distance_DAGs = function(DAG1, DAG2){
  distance = 0
  node_set =  nodes(DAG1)
  amat1 = bnlearn::amat(DAG1)
  amat2 = bnlearn::amat(DAG2)

  ### Compute difference in the skeleton differences in the skeleton

  # Turn both into undirected graphs
  g1 = as.undirected(as.igraph(DAG1))
  g2 = as.undirected(as.igraph(DAG2))

  # Sum all arcs and substract the intersection
  int <- graph.intersection(g1,g2)
  distance = distance + ecount(g1)+ecount(g2)-2*ecount(int)


  ### Compute correct v-structures
  for (v in node_set){
    parents1 = DAG1$nodes[[v]]$parents
    parents2 = DAG2$nodes[[v]]$parents

    # One must be part of a v-structure
    if ((length(parents1)>0) | (length(parents2)>0)){
      distance = distance + length(sym_diff(parents1, parents2))
    }
  }
  return(distance)
}




####################################################
##### For data with standard normal margins  #######
####################################################

### Fits normal margins and returns a list with the mean and sd for each node
fit_normal_margins <- function(data){
  nodes = colnames(data)
  fit_list = list()
  for (node in nodes){
    fit = MASS::fitdistr(data[[node]],"normal")
    estimates = list( mean = fit$estimate[[1]], sd = fit$estimate[[2]])
    fit_list[[node]] =  estimates
  }
  return(fit_list)
}


### Scale data-frame with uniform margins to standard normal margins
to_normal_scale <- function(data){
  data_new = data
  for (i in 1:length(data)){
    data_new[[i]] = stats::qnorm(data[[i]], 0, 1)
  }
  return(data_new)
}

### Scale data-frame with uniform margins to standard normal margins
to_uniform_scale <- function(data){
  data_new = data
  for (i in 1:length(data)){
    data_new[[i]] = VineCopula::pobs(data[[i]])
  }
  return(data_new)
}






###################################
###### not used? ##################
##################################



# Creates a simple string node|cond_set
paste_margin <- function(node, cond_set){
  if (length(cond_set) > 0){
    tag = paste0(node,"|",paste(cond_set, collapse=','))
  }
  else{
    tag = paste0(node)
  }
  return(tag)
}

# Turns string node|cond_set back to its vector notation
unpaste_margin <- function(tag){
  v = strsplit(tag,split="|", fixed = TRUE)[[1]][1]
  cond_set = strsplit(tag,split="|", fixed = TRUE)[[1]][2]
  cond_set = strsplit(cond_set,split=",", fixed = TRUE)[[1]]
  return(list(v = v, cond_set=cond_set))
}

# Takes in two PCBNs and returns TRUE if they have the same ordering
check_correct_order <- function(PCBN1, PCBN2) {
  order1 = PCBN1$order_hash
  order2 = PCBN2$order_hash
  DAG1 = PCBN1$DAG
  DAG2 = PCBN2$DAG

  # TODO: Graphs should have distance equal to 0

  # For now assume they do
  node.names = bnlearn::nodes(DAG1)
  # Loop over all v-structures
  for (v in node.names){
    if (length(DAG1$nodes[[v]]$parents) > 1){
      # The orders should be the same everywhere
      if (!min(order1[[v]] == order2[[v]])) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}



# Iteratively plots a list of active cycles
plot_active_cycles = function(DAG, active_cycle_list){
  if (length(active_cycle_list)==0){
    stop("No active cycles")
    break
  }


  no_list = c("N", "n", "No", "NO")
  adj.mat = bnlearn::amat(DAG)
  L = length(active_cycle_list)
  for (i in 1:L){
    if (i>1){
      more.plots <- readline(prompt="Plot next active cycle? (Y/N): ")
      if (more.plots %in% no_list){
        break
      }
    }

    cat("Plotting active cycle ", i, "of", L, "\n")
    active_cycle = active_cycle_list[[i]]

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
    graphviz.plot(DAG , highlight = list(arcs = df, col = "red", lwd = 3))
  }
}

