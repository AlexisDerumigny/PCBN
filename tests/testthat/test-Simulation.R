test_that("compute_sample_margin works", {

  DAG = create_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")

  fam = matrix(c(0, 1, 1,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)

  rownames(fam) <- c("U1", "U2", "U3")
  colnames(fam) <- c("U1", "U2", "U3")
  tau = 0.2 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  # Initialize data frame
  N = 100
  nodes = bnlearn::nodes(my_PCBN$DAG)
  data = data.frame(matrix(ncol = length(nodes), nrow = N))
  colnames(data) <- nodes

  data[, "U1"] = runif(N)
  data[, "U2"] = runif(N)
  u_1_given2 = compute_sample_margin(object = my_PCBN, data = data,
                                     v = "U1", cond_set = c("U2"))

  expect_identical(data[, "U1"], u_1_given2)

})
