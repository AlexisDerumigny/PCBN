
test_that("remove_CondInd works", {

  DAG = create_empty_DAG(7)
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

test_that("remove_CondInd works on a complicated example", {

  DAG = create_empty_DAG(9)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U6', 'U8')
  DAG = bnlearn::set.arc(DAG, 'U7', 'U9')
  DAG = bnlearn::set.arc(DAG, 'U8', 'U9')

  # plot(bnlearn::as.igraph(DAG))

  expect_identical(
    remove_CondInd(DAG = DAG, node = "U9",
                   cond_set = c("U1", "U2", "U3", "U4", "U5", "U6", "U7")),
    c("U6", "U7")
  )
  expect_identical(
    remove_CondInd(DAG = DAG, node = "U9",
                   cond_set = rev(c("U1", "U2", "U3", "U4", "U5", "U6", "U7"))),
    c("U7", "U6")
  )
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  expect_identical(
    remove_CondInd(DAG = DAG, node = "U9",
                   cond_set = rev(c("U1", "U2", "U3", "U4", "U5", "U6", "U7"))),
    c("U7", "U6")
  )
  expect_identical(
    remove_CondInd(DAG = DAG, node = "U9",
                   cond_set = c("U3", "U5", "U2", "U4", "U7")),
    c("U4", "U7")
  )

})


test_that("active_cycles works" , {

  DAG = create_empty_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  expect_equal(length(active_cycles(DAG)), 0)
  # no active cycle

  DAG = create_empty_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  expect_equal(length(active_cycles(DAG)), 1)
  # 1 active cycle

  DAG = create_empty_DAG(5)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')

  expect_equal(length(active_cycles(DAG)), 2)
  # 2 active cycles

  expect_equal(active_cycles(DAG, early.stopping = TRUE) |>
                 length(),
               1)
  # The first active cycle is returned.
})

test_that("path_check works", {

  DAG = create_empty_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  expect_true(path_check(DAG, c('U1', 'U2', 'U3', 'U4')))

  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  expect_false(path_check(DAG, c('U1', 'U2', 'U3', 'U4')))
  # has a chord

  DAG = create_empty_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U3')
  expect_false(path_check(DAG, c('U1', 'U2', 'U3', 'U4')))
  expect_false(path_check(DAG, c('U2', 'U3', 'U4')))
  # has a converging connection

})
