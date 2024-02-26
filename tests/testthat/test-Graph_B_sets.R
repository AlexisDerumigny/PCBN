
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

