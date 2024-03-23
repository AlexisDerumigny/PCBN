test_that("has_interfering_vstrucs works", {

  DAG = create_DAG(7)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  DAG = bnlearn::set.arc(DAG, 'U1', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U7')

  # There is one interfering v-structure
  expect_true(has_interfering_vstrucs(DAG))

  DAG = bnlearn::set.arc(DAG, 'U1', 'U7')
  # Now no interfering v-structure
  expect_false(has_interfering_vstrucs(DAG))
})

test_that("B_sets_are_increasing works", {

  B_sets = matrix(c(FALSE, FALSE, FALSE, FALSE,
                    TRUE , FALSE, FALSE, FALSE,
                    TRUE , TRUE , FALSE, FALSE,
                    TRUE , TRUE ,  TRUE,  TRUE),
                  nrow = 4, byrow = TRUE)

  expect_true(B_sets_are_increasing(B_sets))

  B_sets = matrix(c(FALSE, FALSE, FALSE, FALSE,
                    TRUE , FALSE, TRUE, FALSE,
                    TRUE , TRUE , FALSE, FALSE,
                    TRUE , TRUE ,  TRUE,  TRUE),
                  nrow = 4, byrow = TRUE)

  expect_false(B_sets_are_increasing(B_sets))
})


test_that("find_B_sets works", {

  DAG = create_DAG(6)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  DAG = bnlearn::set.arc(DAG, 'U1', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U6')

  B_sets = find_B_sets_v(DAG, v = 'U5')

  expect_equal(nrow(B_sets), 3)
  expect_equal(ncol(B_sets), 4)
  expect_identical(colnames(B_sets), c("U1", "U2", "U3", "U4"))
  expect_identical(rownames(B_sets), c("Empty B-set", "U6", "Full B-set"))
})


test_that("find_B_sets works on a more complicated example", {

  DAG = create_DAG(7)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  DAG = bnlearn::set.arc(DAG, 'U1', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U7')

  B_sets = find_B_sets_v(DAG, v = 'U5')

  expect_equal(nrow(B_sets), 4)
  expect_equal(ncol(B_sets), 4)
  expect_identical(colnames(B_sets), c("U1", "U2", "U3", "U4"))
  expect_identical(rownames(B_sets), c("Empty B-set", "U6", "U7", "Full B-set"))
})


test_that("find_B_sets_v works in an example with 1 parent and 2 children", {

  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')

  B_sets = find_B_sets_v(DAG, v = 'U2')

  # We test that B_sets does not drop, i.e. that it stays an object of the class
  # "matrix", even if it has only 1 column.
  expect_identical(class(B_sets), c("matrix", "array"))
  expect_equal(nrow(B_sets), 4)
  expect_equal(ncol(B_sets), 1)
})


test_that("B_sets_are_increasing", {
  B_sets = matrix(c(FALSE, FALSE, FALSE, FALSE,
                    TRUE , FALSE, FALSE, FALSE,
                    TRUE , TRUE , FALSE, FALSE,
                    TRUE , TRUE , TRUE , TRUE),
                  nrow = 4, byrow = TRUE)

  expect_true(B_sets_are_increasing(B_sets))

  B_sets = matrix(c(FALSE, FALSE, FALSE, FALSE,
                    TRUE , FALSE, TRUE , FALSE,
                    TRUE , TRUE , FALSE, FALSE,
                    TRUE , TRUE , TRUE , TRUE),
                  nrow = 4, byrow = TRUE)

  expect_false(B_sets_are_increasing(B_sets))
})

test_that("find_interfering_v_from_B_sets works", {

  DAG = create_DAG(7)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  DAG = bnlearn::set.arc(DAG, 'U1', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U7')

  B_sets = find_B_sets_v(DAG, v = 'U5')
  interf_v = find_interfering_v_from_B_sets(B_sets)

  # `interf_v` is a data.frame
  expect_s3_class(interf_v, "data.frame")

  # There is one interfering v-structure
  expect_equal(nrow(interf_v), 1)

  expect_equal(object = interf_v[1,] |> unlist(),
               expected = c(A = "U6", B = "U7",
                            parents.A.but.not.parents.B = "U1",
                            parents.B.but.not.parents.A = "U2"
               )
  )

})


test_that("find_interfering_v_from_B_sets works with 3 interfering v-structures", {

  DAG = create_DAG(8)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  DAG = bnlearn::set.arc(DAG, 'U1', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U8')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U8')

  B_sets = find_B_sets_v(DAG, v = 'U5')
  unique_B_sets = B_sets_make_unique(B_sets)
  interf_v = find_interfering_v_from_B_sets(B_sets)

  # There is 3 interfering v-structures
  expect_equal(nrow(interf_v), 3)

})


test_that("find_interfering_v_from_B_sets works with 3 interfering v-structures", {

  DAG = create_DAG(8)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  DAG = bnlearn::set.arc(DAG, 'U1', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U8')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U6')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U7')
  DAG = bnlearn::set.arc(DAG, 'U5', 'U8')

  B_sets = find_B_sets_v(DAG, v = 'U5')
  unique_B_sets = B_sets_make_unique(B_sets)
  interf_v = find_interfering_v_from_B_sets(B_sets)

  # There is 3 interfering v-structures
  expect_equal(nrow(interf_v), 1)

})

test_that("B_sets_cut_increments works with the output of B_sets_make_unique", {

  DAG = create_DAG(5)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')

  B_sets = find_B_sets_v(DAG, v = 'U3') |>
    B_sets_make_unique()

  B_sets_incr = B_sets_cut_increments(B_sets = B_sets)

  # There are 2 B-sets increments: first "U1" and then "U2"
  expect_identical(B_sets_incr, list("U1", "U2"))
})

