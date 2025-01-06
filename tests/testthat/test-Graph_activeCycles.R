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

test_that("path_hasConvergingConnections works", {

  DAG = create_empty_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  expect_false(path_hasConvergingConnections(DAG, c('U1', 'U2', 'U3', 'U4')))
  # has no converging connections

  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  expect_false(path_hasConvergingConnections(DAG, c('U1', 'U2', 'U3', 'U4')))
  # has a chord but no converging connections

  DAG = create_empty_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U3')
  expect_true(path_hasConvergingConnections(DAG, c('U1', 'U2', 'U3', 'U4')))
  expect_true(path_hasConvergingConnections(DAG, c('U2', 'U3', 'U4')))
  # has a converging connection
})

test_that("path_hasChords works", {

  DAG = create_empty_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  expect_false(path_hasChords(DAG, c('U1', 'U2', 'U3', 'U4')))
  # has no chords

  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  expect_true(path_hasChords(DAG, c('U1', 'U2', 'U3', 'U4')))
  # has a chord

  DAG = create_empty_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U3')
  expect_false(path_hasChords(DAG, c('U1', 'U2', 'U3', 'U4')))
  expect_false(path_hasChords(DAG, c('U2', 'U3', 'U4')))
  # has a converging connection but no chords
})

