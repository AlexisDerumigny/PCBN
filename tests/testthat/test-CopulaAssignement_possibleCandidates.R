test_that("incoming_arc works", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")

  # Node of interest
  v = "U4"

  # If we start by 1, then the arc 1 -> 3 cannot be used as an incoming arc
  # (it is actually an outgoing arc)
  result = incoming_arc(DAG = DAG, w = "U3", v = v, order_v = c("U1"),
                        order_hash = order_hash)
  expect_identical(result, NULL)


  # If we add U1 then U3, note that U2 can finally be added by incoming arc U2 -> U3
  # Since u_{2 | Ovk} = u_{2 | 13} can be computed using the copula c_{23 | 1}
  result = incoming_arc(DAG = DAG, w = "U2", v = v, order_v = c("U1", "U3"),
                        order_hash = order_hash)
  expect_identical(result, "U3")

  # If we start by U3, we can only add U1 as an incoming arc, not U2.
  # This is because U1 is first in the order of parents of U3
  # So U1 | U3 can be computed but not U2 | U3.
  result = incoming_arc(DAG = DAG, w = "U1", v = v, order_v = c("U3"),
                        order_hash = order_hash)
  expect_identical(result, "U3")

  result = incoming_arc(DAG = DAG, w = "U2", v = v, order_v = c("U3"),
                        order_hash = order_hash)
  expect_identical(result, NULL)

  # We can however add U2 after U1
  result = incoming_arc(DAG = DAG, w = "U2", v = v, order_v = c("U3", "U1"),
                        order_hash = order_hash)
  expect_identical(result, "U3")

})

test_that("incoming_arc works", {

  DAG = create_DAG(5)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")
  order_hash[['U4']] = c("U1", "U3", "U2")

  # Node of interest
  v = "U5"

  # We can add U2 by U3
  result = incoming_arc(DAG = DAG, w = "U2", v = v, order_v = c("U1", "U3"),
                        order_hash = order_hash)
  expect_identical(result, "U3")

  # We add U4 by outgoing arc, and can add U2 by U4
  # we cannot add U2 by U3 anymore due to the d-separation contraint.
  result = incoming_arc(DAG = DAG, w = "U2", v = v, order_v = c("U1", "U3", "U4"),
                        order_hash = order_hash)
  expect_identical(result, "U4")

})


test_that("incoming_arc works", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")

  # Node of interest
  v = "U4"

  # If we start by 1, then the arc 1 -> 3 can be used as an outgoing arc
  result = outgoing_arc(DAG = DAG, w = "U3", v = v, order_v = c("U1"),
                        order_hash = order_hash)
  expect_identical(result, "U1")


  # If we add U1 then U3, note that U2 can finally be added by an incoming arc
  # not by an outgoing arc
  result = outgoing_arc(DAG = DAG, w = "U2", v = v, order_v = c("U1", "U3"),
                        order_hash = order_hash)
  expect_identical(result, NULL)
})

