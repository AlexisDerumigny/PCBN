test_that("BiCopCondFit estimates the copulas well", {

  # Setup code =============================================================
  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U4']] = c("U2", "U1", "U3")
  complete_and_check_orders(DAG, order_hash)

  fam = matrix(c(0, 1, 1, 1,
                 0, 0, 1, 1,
                 0, 0, 0, 1,
                 0, 0, 0, 0), byrow = TRUE, ncol = 4)

  tau = 0.2 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  mydata = sample_PCBN(my_PCBN, N = 100)

  e = default_envir()

  # Estimation code ===========================================================

  C_12 = BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U2",
                      cond_set = c(), familyset = 1, order_hash = order_hash,
                      e = e)

  C_12_direct = VineCopula::BiCopSelect(mydata[, "U1"], mydata[, "U2"], familyset = 1)

  expect_equal(C_12$tau, C_12_direct$tau)

  C_13 = BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U3",
                      cond_set = c(), familyset = 1, order_hash = order_hash,
                      e = e)

  C_13_direct = VineCopula::BiCopSelect(mydata[, "U1"], mydata[, "U3"], familyset = 1)

  expect_equal(C_13$tau, C_13_direct$tau)

  BiCopCondFit(data = mydata, DAG = DAG, v = "U2", w = "U4",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e)

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U4",
               cond_set = c("U2"), familyset = 1, order_hash = order_hash,
               e = e)

  BiCopCondFit(data = mydata, DAG = DAG, v = "U3", w = "U4",
               cond_set = c("U1", "U2"), familyset = 1, order_hash = order_hash,
               e = e)


  C_12_again = BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U2",
                            cond_set = c(), familyset = 1, order_hash = order_hash,
                            e = e)

  expect_equal(C_12_again$tau, C_12_direct$tau)

})


test_that("BiCopCondFit completes the hashmaps as needed", {
  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U4']] = c("U2", "U1", "U3")
  complete_and_check_orders(DAG, order_hash)

  fam = matrix(c(0, 1, 1, 1,
                 0, 0, 1, 1,
                 0, 0, 0, 1,
                 0, 0, 0, 0), byrow = TRUE, ncol = 4)

  tau = 0.2 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  N = 100
  mydata = sample_PCBN(my_PCBN, N = N)

  e = default_envir()

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U2",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e)
  # ls(e)
  # length(r2r::keys(e$keychain))
  # r2r::keys(e$keychain)[[1]] |> print_key_keychain()
  # r2r::keys(e$keychain)[[2]]
  # r2r::keys(e$keychain)[[3]]
  # r2r::keys(e$keychain)[[4]]
  # r2r::keys(e$keychain)[[5]]

  all_keys_keychain = e$keychain |>
    r2r::keys() |>
    lapply(FUN = print_key_keychain) |>
    unlist() |>
    sort()

  expect_identical(all_keys_keychain,
                   c("U1", "U1 | U2", "U1, U2", "U2", "U2 | U1") )

  expect_equal(e$margin_hash[[ list(margin = c("U1"), cond = character(0)) ]] |>
                 length() ,
               N)

  # e$keychain[[r2r::keys(e$keychain)[[2]]]]
  # e$keychain[[r2r::keys(e$keychain)[[3]]]]
  # e$keychain[[r2r::keys(e$keychain)[[4]]]]
  # e$keychain[[r2r::keys(e$keychain)[[5]]]]

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U3",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e)

  # Both copulas are different, so the keys should also be different
  expect_true(! identical(e$keychain[[list(margins = c("U1", "U2"), cond = character(0))]] ,
                          e$keychain[[list(margins = c("U1", "U3"), cond = character(0))]] ) )

  e$margin_hash[[ e$keychain[[r2r::keys(e$keychain)[[3]]]] ]]


  BiCopCondFit(data = mydata, DAG = DAG, v = "U2", w = "U4",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e)

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U4",
               cond_set = c("U2"), familyset = 1, order_hash = order_hash,
               e = e)

  BiCopCondFit(data = mydata, DAG = DAG, v = "U3", w = "U4",
               cond_set = c("U1", "U2"), familyset = 1, order_hash = order_hash,
               e = e)

  all_keys_keychain = e$keychain |>
    r2r::keys() |>
    lapply(FUN = print_key_keychain) |>
    unlist() |>
    sort()
})




test_that("ComputeCondMargin works", {
  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U4']] = c("U2", "U1", "U3")
  complete_and_check_orders(DAG, order_hash)

  fam = matrix(c(0, 1, 1, 1,
                 0, 0, 1, 1,
                 0, 0, 0, 1,
                 0, 0, 0, 0), byrow = TRUE, ncol = 4)

  tau = 0.2 * fam

  my_PCBN = new_PCBN(
    DAG, order_hash,
    copula_mat = list(tau = tau, fam = fam))

  mydata = sample_PCBN(my_PCBN, N = 100)

  e = default_envir()

  U1 = ComputeCondMargin(data = mydata, DAG = my_PCBN, v = "U1", cond_set = NULL,
                         familyset = 1, order_hash = order_hash,
                         e = e, verbose = 0)

  expect_identical(U1, mydata[, "U1"])

  key_U1 = list(margin = "U1", cond = character(0))
  expect_true(r2r::has_key(x = e$keychain, key_U1))

  expect_identical(key_U1, e$keychain[[key_U1]])
})

