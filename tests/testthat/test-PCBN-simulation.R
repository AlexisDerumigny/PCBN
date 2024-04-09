test_that("PCBN_sim does not sample from a non-restricted PCBN", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U4']] = c("U2", "U3")

  fam = matrix(c(0, 1, 1, 1,
                 0, 0, 1, 1,
                 0, 0, 0, 1,
                 0, 0, 0, 0), byrow = TRUE, ncol = 4)

  tau = 0.2 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  # 1 active cycle, so no simulation is possible
  expect_error({ mydata = PCBN_sim(my_PCBN, N = 5) },
               class = "UnRestrictedPCBNError")
})

test_that("PCBN_sim does not sample if the ordering do not abide by the Bsets", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U2", "U1")
  order_hash[['U4']] = c("U1", "U3")

  fam = matrix(c(0, 1, 1, 1,
                 0, 0, 1, 1,
                 0, 0, 0, 1,
                 0, 0, 0, 0), byrow = TRUE, ncol = 4)

  tau = 0.2 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  # Order for U3 does not abide by the B-sets, so no simulation is possible
  expect_error({ mydata = PCBN_sim(my_PCBN, N = 5) },
               class = "ParentalOrderingsBsetsError")

  my_PCBN$order_hash[['U3']] = c("U1", "U2")
  # Now this works
  mydata = PCBN_sim(my_PCBN, N = 5)
})

test_that("compute_sample_margin works", {

  DAG = create_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")
  complete_and_check_orders(DAG, order_hash)

  fam = matrix(c(0, 1, 1,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)
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


test_that("PCBN_sim applies proper recursion of h-functions for small example", {

  # Initialize PCBN
  DAG = create_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")
  complete_and_check_orders(DAG, order_hash)

  fam = matrix(c(0, 1, 1,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)
  tau = 0.5 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  N = 10

  # Sample data using package
  set.seed(51)
  PCBN_sim_data = PCBN_sim(object = my_PCBN, N = N)

  # Sample data manually
  par_13  = VineCopula::BiCopTau2Par(family = fam[1,3], tau = tau[1,3])
  par_23  = VineCopula::BiCopTau2Par(family = fam[2,3], tau = tau[2,3])

  set.seed(51)
  U1 = stats::runif(N, 0, 1)
  U2 = stats::runif(N, 0, 1)
  U2_given_1 = U2

  marginal = stats::runif(N, 0, 1)
  U3_given_1 = VineCopula::BiCopHinv1(U1, marginal, family = fam[1,3], par = par_13)
  U3_given_12 = VineCopula::BiCopHinv1(U2_given_1, U3_given_1, family = fam[2,3], par = par_23)

  # Answers must agree
  expect_equal(PCBN_sim_data$U1, U1)
  expect_equal(PCBN_sim_data$U2, U2)
  expect_equal(PCBN_sim_data$U3, U3_given_12)

})



