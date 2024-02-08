
test_that("remove_CondInd works", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  # plot(DAG)

  expect_identical(
    remove_CondInd(DAG = DAG, node = "U1", cond_set = c("U2")),
    character(0)
  )

  expect_identical(
    remove_CondInd(DAG = DAG, node = "U3", cond_set = c("U1", "U2")),
    c("U1", "U2")
  )

  expect_identical(
    remove_CondInd(DAG = DAG, node = "U4", cond_set = c("U1")),
    c("U1")
  )

  expect_identical(
    remove_CondInd(DAG = DAG, node = "U4", cond_set = c("U1", "U3")),
    c("U1", "U3")
  )

  expect_identical(
    remove_CondInd(DAG = DAG, node = "U4", cond_set = c("U1", "U2", "U3")),
    c("U2", "U3")
  )


})
