
# Setup code ===================================================================

DAG = create_empty_DAG(4)
DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

order_hash = r2r::hashmap()
order_hash[['U4']] = c("U2", "U1", "U3")
complete_and_check_orders(DAG, order_hash)

fam = matrix(c(0, 1, 1, 1,
               0, 0, 0, 1,
               0, 0, 0, 1,
               0, 0, 0, 0), byrow = TRUE, ncol = 4)

tau = 0.2 * fam

my_PCBN = new_PCBN(
  DAG, order_hash,
  copula_mat = list(tau = tau, fam = fam))

N = 100
mydata = PCBN_sim(my_PCBN, N = N)

# Silent the output by default
verbose = 0

# Tests ========================================================================


test_that("BiCopCondFit estimates the copulas well",  {
  e = default_envir()

  C_12 = BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U2",
                      cond_set = c(), familyset = 1, order_hash = order_hash,
                      e = e, verbose = verbose, method = "mle")

  C_12_direct = VineCopula::BiCopSelect(mydata[, "U1"], mydata[, "U2"], familyset = 1)

  expect_equal(C_12$tau, C_12_direct$tau)

  C_13 = BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U3",
                      cond_set = c(), familyset = 1, order_hash = order_hash,
                      e = e, verbose = verbose, method = "mle")

  C_13_direct = VineCopula::BiCopSelect(mydata[, "U1"], mydata[, "U3"], familyset = 1)

  expect_equal(C_13$tau, C_13_direct$tau)

  BiCopCondFit(data = mydata, DAG = DAG, v = "U2", w = "U4",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U4",
               cond_set = c("U2"), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")

  BiCopCondFit(data = mydata, DAG = DAG, v = "U3", w = "U4",
               cond_set = c("U1", "U2"), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")


  C_12_again = BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U2",
                            cond_set = c(), familyset = 1, order_hash = order_hash,
                            e = e, verbose = verbose, method = "mle")

  expect_equal(C_12_again$tau, C_12_direct$tau)

})


test_that("BiCopCondFit makes all the keys as specified",  {
  e = default_envir()

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U2",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")

  # We get all the keys in a textual form
  all_keys_keychain = e$keychain |>
    r2r::keys() |>
    lapply(FUN = print_key_keychain) |>
    unlist() |>
    sort()

  expect_identical(all_keys_keychain,
                   c("U1", "U1 | U2", "U1, U2", "U2", "U2 | U1") )

  expect_identical(e$keychain[[list(margin = "U1", cond = "U2")]]$name,
                   "U1 | U2")
  expect_identical(e$keychain[[list(margin = "U2", cond = "U1")]]$name,
                   "U2 | U1")
})


test_that("BiCopCondFit completes the hashmaps as needed",  {
  e = default_envir()

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U2",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")
  # ls(e)
  # length(r2r::keys(e$keychain))
  # r2r::keys(e$keychain)[[1]] |> print_key_keychain()
  # r2r::keys(e$keychain)[[2]]
  # r2r::keys(e$keychain)[[3]]
  # r2r::keys(e$keychain)[[4]]
  # r2r::keys(e$keychain)[[5]]

  expect_equal(e$margin_hash[[ list(name = "U1", margin = c("U1"),
                                    cond = character(0), copula = NULL) ]] |>
                 length() ,
               N)

  # e$keychain[[r2r::keys(e$keychain)[[2]]]]
  # e$keychain[[r2r::keys(e$keychain)[[3]]]]
  # e$keychain[[r2r::keys(e$keychain)[[4]]]]
  # e$keychain[[r2r::keys(e$keychain)[[5]]]]

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U3",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")

  # Both copulas are different, so the keys should also be different
  expect_true(! identical(e$keychain[[list(margins = c("U1", "U2"),
                                           cond = character(0))]] ,

                          e$keychain[[list(margins = c("U1", "U3"),
                                           cond = character(0))]] ) )

  expect_equal(e$margin_hash[[ e$keychain[[ list(margin = c("U3"),
                                                 cond = character(0)) ]] ]] |>
                 length(),
               N)


  BiCopCondFit(data = mydata, DAG = DAG, v = "U2", w = "U4",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U4",
               cond_set = c("U2"), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")

  BiCopCondFit(data = mydata, DAG = DAG, v = "U3", w = "U4",
               cond_set = c("U1", "U2"), familyset = 1, order_hash = order_hash,
               e = e, verbose = verbose, method = "mle")

  all_keys_keychain = e$keychain |>
    r2r::keys() |>
    lapply(FUN = print_key_keychain) |>
    unlist() |>
    sort()
})




test_that("ComputeCondMargin works", {
  e = default_envir()

  U1 = ComputeCondMargin(data = mydata, DAG = my_PCBN, v = "U1", cond_set = NULL,
                         familyset = 1, order_hash = order_hash,
                         e = e, verbose = verbose, method = "mle")

  expect_identical(U1, mydata[, "U1"])

  key_U1 = list(margin = "U1", cond = character(0))
  expect_true(r2r::has_key(x = e$keychain, key_U1))

  expect_identical(object = e$keychain[[key_U1]],
                   expected = list(name = "U1", margin = "U1",
                                   cond = character(0), copula = NULL) )
})


test_that("fit_copulas gives the same results as what is in the keychain", {
  e = default_envir()

  result = fit_copulas(data = mydata, DAG = DAG, order_hash = order_hash,
                       familyset = 1, e = e, verbose = verbose, method = "mle")

  expect_equal(
    e$copula_hash[[
      e$keychain[[list(margins = c("U1", "U2"), cond = character(0))]]
    ]]$tau ,

    result$copula_mat$tau["U1", "U2"]
  )
  expect_equal(
    e$copula_hash[[
      e$keychain[[list(margins = c("U1", "U3"), cond = character(0))]]
    ]]$tau ,

    result$copula_mat$tau["U1", "U3"]
  )
  expect_equal(
    e$copula_hash[[
      e$keychain[[list(margins = c("U2", "U4"), cond = character(0))]]
    ]]$tau ,

    result$copula_mat$tau["U2", "U4"]
  )
  expect_equal(
    e$copula_hash[[
      e$keychain[[list(margins = c("U1", "U4"), cond = c("U2"))]]
    ]]$tau ,

    result$copula_mat$tau["U1", "U4"]
  )
  expect_equal(
    e$copula_hash[[
      e$keychain[[list(margins = c("U3", "U4"), cond = c("U1", "U2"))]]
    ]]$tau ,

    result$copula_mat$tau["U3", "U4"]
  )
})


test_that("fit_copulas respects the adjacency matrix", {

  DAG = create_empty_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U2")

  fam = matrix(c(0, 0, 0,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)

  tau = 0.8 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  mydata = PCBN_sim(my_PCBN, N = 5)
  e = default_envir()

  result = fit_copulas(data = mydata, DAG = DAG, order_hash = order_hash,
                       familyset = 1, e = e, verbose = verbose, method = "mle")

  expect_equal(result$copula_mat$fam, my_PCBN$copula_mat$fam)
})


test_that("fit_copulas can take the family object from familyMatrix", {

  DAG = create_empty_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U2")

  fam = matrix(c(0, 0, 0,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)

  tau = 0.8 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  mydata = PCBN_sim(my_PCBN, N = 5)

  result1 = fit_copulas(data = mydata, DAG = DAG, order_hash = order_hash,
                        familyset = 1, e = default_envir(), verbose = verbose)

  result3 = fit_copulas(data = mydata, DAG = DAG, order_hash = order_hash,
                        familyset = 3, e = default_envir(), verbose = verbose)

  colnames(fam) <- rownames(fam) <- paste0("U", 1:3)

  result_matrix = fit_copulas(
    data = mydata, DAG = DAG, order_hash = order_hash,
    familyset = 3,  # should be ignore since familyMatrix is given
    familyMatrix = fam, e = default_envir(), verbose = verbose)

  expect_true (identical(result1$copula_mat, result_matrix$copula_mat) )
  expect_false(identical(result3$copula_mat, result_matrix$copula_mat) )
})


test_that("fit_copulas works for an example of dimension 5", {

  # Initialize PCBN
  DAG = create_empty_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")
  complete_and_check_orders(DAG, order_hash)

  fam = matrix(c(0, 0, 1,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)
  tau = fam * 0.8

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  N = 5000
  mydata = PCBN_sim(object = my_PCBN, N = N)

  parent_given_lower = compute_sample_margin(object = my_PCBN,
                                             data = mydata,
                                             v = "U2",
                                             cond_set = "U1",
                                             check_PCBN = FALSE)

  par_ = VineCopula::BiCopTau2Par(family = 1, tau = 0.8)
  U3_correct =
    VineCopula::BiCopHinv1(u1 = mydata[, "U1"],
                           u2 =  VineCopula::BiCopHinv1(u1 = mydata[, "U2"],
                                                        u2 = runif(N),
                                                        family = 1,
                                                        par = par_),
                           family = 1,
                           par = par_)

  plot(mydata[, "U1"], U3_correct)
  cor(mydata[, "U1"], U3_correct, method = "kendall")
  plot(mydata[, "U1"], mydata[, "U3"])


  e = default_envir()

  result = fit_copulas(data = mydata, DAG = DAG,
                       order_hash = order_hash, familyset = 1,
                       e = e, verbose = verbose)

  expect_equal(object = result$copula_mat$tau |> unname(),
               expected = tau, tolerance = 0.05)

})



test_that("fit_copulas works for an example of dimension 5", {

  # Initialize PCBN
  DAG = create_empty_DAG(5)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U5')
  DAG = bnlearn::set.arc(DAG, 'U4', 'U5')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")
  order_hash[['U5']] = c("U4", "U3", "U1", "U2")
  complete_and_check_orders(DAG, order_hash)

  fam = matrix(c(0, 0, 1, 0, 1,
                 0, 0, 1, 0, 1,
                 0, 0, 0, 1, 1,
                 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 0), byrow = TRUE, ncol = 5)
  tau = fam * 0.8

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  N = 5000
  mydata = PCBN_sim(object = my_PCBN, N = N)

  e = default_envir()

  expect_no_error({
    result = fit_copulas(data = mydata, DAG = DAG,
                         order_hash = order_hash, familyset = 1,
                         e = e, verbose = verbose)
  })

  result$copula_mat$tau

  expect_equal(
    e$copula_hash[[
      e$keychain[[list(margins = c("U1", "U3"), cond = character(0))]]
    ]]$tau ,

    result$copula_mat$tau["U1", "U3"]
  )

  expect_equal(
    BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U3",
                 cond_set = c(), familyset = 1, order_hash = order_hash,
                 e = e, verbose = verbose, method = "mle")$tau ,

    BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U3",
                 cond_set = c(), familyset = 1, order_hash = order_hash,
                 e = default_envir(), verbose = verbose, method = "mle")$tau
  )
})


test_that("fit_all_orders works", {

  DAG = create_empty_DAG(3)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U3')

  order_hash = r2r::hashmap()
  order_hash[['U3']] = c("U1", "U2")

  fam = matrix(c(0, 0, 1,
                 0, 0, 1,
                 0, 0, 0), byrow = TRUE, ncol = 3)

  tau = 0.2 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  mydata = PCBN_sim(my_PCBN, N = 5)
  e = default_envir()

  expect_error(
    fit_all_orders(data = mydata, DAG = DAG,
                   familyset = 1, e = e, score_metric = "aaaaa")
  )

})


