test_that("Hill climbing works for a small example", {

  DAG = create_empty_DAG(4)
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


test_that("allowed_operations_fromDAG returns the right matrix in a small example", {
  DAG = create_empty_DAG(4)

  operations = allowed_operations_fromDAG(DAG)

  expect_identical(length(dim(operations)), 2L)

  # 12 operations are possible
  expect_identical(NROW(operations), 4L * 3L)
  expect_identical(NCOL(operations), 3L)
})


test_that("operation_do does not modify the original DAG", {
  DAG = create_empty_DAG(4)

  operations = allowed_operations_fromDAG(DAG)

  op = operations[1,]
  DAG2 = operation_do(DAG, op)

  # This should not modify the original DAG
  expect_identical(nrow(DAG$arcs), 0L)

  # But should modify the first DAG
  expect_identical(nrow(DAG2$arcs), 1L)
})


test_that("operation_do and operation_undo are invert", {
  DAG1 = create_empty_DAG(4)

  operations = allowed_operations_fromDAG(DAG1)

  op = operations[1,]
  DAG2 = operation_do(DAG1, op)

  DAG3 = operation_undo(DAG2, op)

  # This should return the original DAG
  expect_identical(DAG1, DAG3)
})
