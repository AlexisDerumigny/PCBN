test_that("BiCopCondFit works", {
  DAG = create_DAG(4)
  DAG = bnlearn::set.arc(DAG, 'U1', 'U2')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U3')
  DAG = bnlearn::set.arc(DAG, 'U1', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U2', 'U4')
  DAG = bnlearn::set.arc(DAG, 'U3', 'U4')

  order_hash = r2r::hashmap()
  order_hash[['U4']] = c("U2", "U1", "U3")

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

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U2",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e)
  ls(e)

  length(r2r::keys(e$keychain))
  r2r::keys(e$keychain)[[1]]
  r2r::keys(e$keychain)[[2]]
  r2r::keys(e$keychain)[[3]]
  r2r::keys(e$keychain)[[4]]
  r2r::keys(e$keychain)[[5]]

  e$keychain[[r2r::keys(e$keychain)[[1]]]]
  e$keychain[[r2r::keys(e$keychain)[[2]]]]
  e$keychain[[r2r::keys(e$keychain)[[3]]]]
  e$keychain[[r2r::keys(e$keychain)[[4]]]]
  e$keychain[[r2r::keys(e$keychain)[[5]]]]

  e$keychain[[r2r::keys(e$keychain)[[1]]]]

  BiCopCondFit(data = mydata, DAG = DAG, v = "U1", w = "U3",
               cond_set = c(), familyset = 1, order_hash = order_hash,
               e = e)

  length(r2r::keys(e$keychain))
  for (i in 1:9){
    print(r2r::keys(e$keychain)[[i]])
  }

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

})
