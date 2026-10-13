test_that("set_arc adds arcs and rejects invalid edits with either backend", {
  backends = "bnlearn"
  if (.caugi_available()) {
    backends = c(backends, "caugi")
  }
  for (backend in backends) {
    graph = create_empty_DAG(3, graph_backend = backend)
    graph = set_arc(graph, "U1", "U2")
    expect_identical(.graph_backend(graph), backend)
    expect_identical(.graph_parents(graph, "U2"), "U1")
    expect_identical(.graph_nodes(graph), c("U1", "U2", "U3"))
    expect_equal(.graph_adjacency(set_arc(graph, "U1", "U2")),
                 .graph_adjacency(graph))
    graph = set_arc(graph, "U2", "U3")
    expect_error(set_arc(graph, "U3", "U1"))
    expect_error(set_arc(graph, "missing", "U2"))
    expect_error(set_arc(graph, "U1", "missing"))
    expect_error(set_arc(graph, "U1", "U1"))
  }
  expect_error(set_arc(list(), "U1", "U2"), "bnlearn or caugi graph")
  expect_error(set_arc(create_empty_DAG(2), "U1", "U2",
                       graph_backend = "other"), "graph_backend")
})

test_that("set_arc converts graphs when a backend is specified", {
  skip_if_not_installed("caugi", minimum_version = "1.3.0")
  graph = set_arc(create_empty_DAG(3), "U1", "U2", graph_backend = "caugi")
  expect_true(caugi::is_caugi(graph))
  graph = set_arc(graph, "U2", "U3", graph_backend = "bnlearn")
  expect_s3_class(graph, "bn")
  expect_identical(bnlearn::parents(graph, "U2"), "U1")
  expect_identical(bnlearn::parents(graph, "U3"), "U2")
})
