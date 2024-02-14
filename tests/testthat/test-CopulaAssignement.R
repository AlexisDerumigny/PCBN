test_that("is_cond_copula_specified works for 4 dimensional example", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")
  complete_and_check_orders(DAG, order_hash)

  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U1", v = "U3", cond = c()) )
  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U3", v = "U1", cond = c()) )
  # returns TRUE because the copula c_{1,3} is known

  expect_false(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                        w = "U2", v = "U3", cond = c()) )
  expect_false(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                        w = "U3", v = "U2", cond = c()) )
  # returns FALSE because the copula c_{2,3} is not known

  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U2", v = "U3", cond = c("U1")))
  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U3", v = "U2", cond = c("U1")))
  # returns TRUE because the copula c_{2,3 | 1} is known

  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U3", v = "U4", cond = c()) )
  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U4", v = "U3", cond = c()) )
  # returns TRUE because the copula c_{3,4} is known
})


test_that("is_cond_copula_specified works for another 4 dimensional example", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U4']] = c("U2", "U1", "U3")
  complete_and_check_orders(DAG, order_hash)

  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U2", v = "U4",
                                       cond = NULL))

  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U4", v = "U2",
                                       cond = NULL))
})


test_that("find_cond_copula_specified works for 4 dimensional example", {


  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U2", "U1")
  order_hash[['U4']] = c("U2", "U3")
  complete_and_check_orders(DAG, order_hash)

  expect_identical(
    find_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                               v = "U3", cond = c("U2"))[["w"]] ,
    "U2"
  )
  # because c_{2,3} is known

  expect_identical(
    find_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                               v = "U3", cond = c("U1", "U2"))[["w"]] ,
    "U1"
  )
  # because c_{1,3 | 2} is known

  expect_identical(
    find_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                               v = "U4", cond = c("U2"))[["w"]] ,
    "U2"
  )
  # because c_{2,4} is known

  expect_identical(
    find_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                               v = "U4", cond = c("U2", "U3"))[["w"]] ,
    "U3"
  )
  # because c_{3,4 | 2} is known
})

test_that("is_cond_copula_specified works for another 4 dimensional example", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U4']] = c("U2", "U1", "U3")
  complete_and_check_orders(DAG, order_hash)

  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U2", v = "U4",
                                       cond = NULL))

  expect_true(is_cond_copula_specified(DAG = DAG, order_hash = order_hash,
                                       w = "U4", v = "U2",
                                       cond = NULL))
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



