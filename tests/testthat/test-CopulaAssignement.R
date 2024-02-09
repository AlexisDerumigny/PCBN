test_that("is_cond_copula_specified works for 3 dimensional example", {

  DAG = create_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")

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
})


test_that("find_cond_copula_specified works for 1 dimensional example", {


  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U2", "U1")
  order_hash[['U4']] = c("U2", "U3")

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


