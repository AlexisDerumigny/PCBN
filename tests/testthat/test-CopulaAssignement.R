
test_that("extend_orders works for the edge cases: 0 parent, 1 parent", {

  DAG = create_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  all_orders_1 = extend_orders(DAG = DAG, all_orders = list(), node = "U1")

  expect_equal(length(all_orders_1), 0)

  all_orders_2 = extend_orders(DAG = DAG, all_orders = list(), node = "U2")

  expect_equal(length(all_orders_2), 1)

  all_orders_3 = extend_orders(DAG = DAG, all_orders = list(), node = "U3")

  expect_equal(length(all_orders_3), 1)
})


test_that("find_all_orders_v works", {

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

  all_orders_5 = find_all_orders_v(DAG, v = "U5", order_hash = order_hash)

  # There are 8 possible orders for U5
  expect_equal(length(all_orders_5), 8)
})


test_that("complete_and_check_orders works for 4 dimensional example", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()

  expect_error({complete_and_check_orders(DAG, order_hash)})
  # Error because the order of the parents on "U3" should be specified.

  order_hash[['U3']] = c("U1", "U2")
  complete_and_check_orders(DAG, order_hash)

  expect_identical(r2r::keys(order_hash) |> unlist(),
                   c("U3", "U4"))
  # U4 should have been added since it has only parent.

  expect_identical(order_hash[["U4"]], c("U3"))
  # The parent of "U4" should be "U3"
})



