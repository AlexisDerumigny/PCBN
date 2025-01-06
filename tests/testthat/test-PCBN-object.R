test_that("new_PCBN respect contraints on copula_mat", {
  DAG = create_empty_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")

  expect_error(new_PCBN(DAG = DAG, order_hash = order_hash,
                        copula_mat = NULL))

  expect_error(new_PCBN(DAG = DAG, order_hash = order_hash,
                        copula_mat = list() ) )

  expect_error(new_PCBN(DAG = DAG, order_hash = order_hash,
                        copula_mat = list(tau = 0.5, fam = 1) ) )

  expect_error(new_PCBN(DAG = DAG, order_hash = order_hash,
                        copula_mat = list(tau = matrix(1:9, ncol = 9),
                                          fam = matrix(1:9, ncol = 9)) ) )

  # Missing family/tau
  fam = matrix(c(0, 0, 0,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)
  tau = 0.2 * fam

  expect_error(new_PCBN(DAG = DAG, order_hash = order_hash,
                        copula_mat = list(tau = tau, fam = fam) ) )


  # Extra family/tau
  fam = matrix(c(0, 1, 1,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)
  tau = 0.2 * fam

  expect_error(new_PCBN(DAG = DAG, order_hash = order_hash,
                        copula_mat = list(tau = tau, fam = fam) ) )
})


test_that("new_PCBN controls names on copula_mat", {
  DAG = create_empty_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")

  fam = matrix(c(0, 0, 1,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)
  tau = 0.2 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  rownames(fam) <- c("U1", "U2", "U3")
  colnames(fam) <- c("U1", "U2", "U3")

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  rownames(fam) <- c("U1", "U2", "aaaaa")

  expect_error(new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam)))
})
