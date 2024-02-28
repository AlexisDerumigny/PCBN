
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

