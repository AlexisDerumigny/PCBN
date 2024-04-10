test_that("Hill climbing works for a small example", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U2", "U1")
  order_hash[['U4']] = c("U2", "U3")

  fam = matrix(c(0, 0, 1, 0,
                 0, 0, 1, 1,
                 0, 0, 0, 1,
                 0, 0, 0, 0), byrow = TRUE, ncol = 4)
  tau = 0.8 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  mydata = PCBN_sim(my_PCBN, N = 5000)

  result = hill.climbing.PCBN(data = mydata, familyset = 1, verbose = 0,
                              score_metric = "BIC")
  result$best_fit$copula_mat$fam
  result$best_fit$copula_mat$tau
})
