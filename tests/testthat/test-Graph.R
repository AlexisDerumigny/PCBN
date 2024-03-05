
test_that("remove_CondInd works", {

  DAG = create_DAG(7)
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

  DAG = create_DAG(9)
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


test_that("DAG_to_restricted works"){

  # DAG with active cycle
  DAG = create_DAG(5)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  fixed_DAG = DAG_to_restricted(DAG)

  # Fixed graph should have extra arcs 1 -> 5, 2 -> 5
  expected_DAG = DAG
  expected_DAG = bnlearn::set.arc(expected_DAG, 'U2', 'U5')
  expected_DAG = bnlearn::set.arc(expected_DAG, 'U1', 'U5')

  false_positives = bnlearn::compare(expected_DAG, fixed_DAG)$fp
  false_negatives = bnlearn::compare(expected_DAG, fixed_DAG)$fn

  expect_identical(false_positives + false_negatives, 0)



  # DAG with an interfering v-structures node 3
  DAG = create_DAG(5)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')

  fixed_DAG = DAG_to_restricted(DAG)

  # Fixed graph should have extra arcs 1 -> 5
  expected_DAG = DAG
  expected_DAG = bnlearn::set.arc(expected_DAG, 'U1', 'U5')

  false_positives = bnlearn::compare(expected_DAG, fixed_DAG)$fp
  false_negatives = bnlearn::compare(expected_DAG, fixed_DAG)$fn

  expect_identical(false_positives + false_negatives, 0)
}


test_that("DAG_to_restricted works with complicated graph"){

  DAG = create_DAG(8)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U7')

  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U6')

  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U8')

  DAG = bnlearn::set.arc(DAG, 'U4', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U8')

  DAG = bnlearn::set.arc(DAG, 'U6', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U6', 'U8')

  find_B_sets(DAG)$nodes_with_inter_vs

  fixed_DAG = DAG_to_restricted(DAG)

  # TODO: currently test is failed since fixing the inter-v's at node
  # U3 introduces new inter-v's at node 6
  expect_identical(find_B_sets(fixed_DAG)$has_interfering_vstrucs, FALSE)
}
















