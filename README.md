Development of the PCBN package
================

# Files that have been checked and are fully working

## 1. PCBN object and structure

- File `PCBN-object.R`

  - `new_PCBN()`: working, but may be changed / developed further. So
    the API is still not stable.

  - `plot.PCBN()` and `print.PCBN()` are working but are somehow basic
    and can be improved.

- File `PCBN-tag.R` : contains some utilities about the tag/hashmap
  system for storing estimated copulas and estimated conditional
  margins.

  The functions are: `paste_margin()`, `unpaste_margin()`,
  `print_key_keychain()`, `make_and_store_keyCopula()` and
  `make_and_store_keyMargin()`.

## 2. PCBN and graphical models

- File `PCBN-find-copula-specified.R`

  - `is_cond_copula_specified()`, `find_cond_copula_specified()` are
    already checked and working.

- File `CopulaAssignement_possibleCandidates.R`

  - `possible_candidates()`, `possible_candidate_incoming_arc()`,
    `possible_candidate_outgoing_arc()` are already checked and working.

- File `Parents_ordering.R`

  - `find_all_orders()`, `extend_orders()`, `find_all_orders_v()`,
    `complete_and_check_orders()`, `is_order_abiding_Bsets`,
    `is_order_abiding_Bsets_v` are already checked and working.

## 3. Simulation and estimation

- File `PCBN-simulation.R`

  - `PCBN_sim()`, `compute_sample_margin()`

  Both functions are finished, but will return an error if the PCBN does
  not satisfy the restrictions. This is the case for example if the
  order of the parents are not chosen in the correct way (even though
  the graph may not have any active cycle nor interfering v-structure).

  They also can be further optimized by storing already computed margins
  in a hashmap.

- File `PCBN-logLik.R`

  - `logLik.PCBN()` and `PCBN_PDF()` are working.

  **TODO:** add methods `AIC` and `BIC`.

- File `Estimation.R`

  - `BiCopCondFit()` and `ComputeCondMargin()` have been rewritten.

  A vignette *“How to use the estimation procedures”* has been written
  to explain how these functions work.

  - `fit_copulas()` and `fit_all_orders()` are working.

## 4. Graph utilities

- File `Graph_B_sets.R`
  - `has_interfering_vstrucs()`, `find_B_sets()`, `find_B_sets_v()`,
    `find_interfering_v_from_B_sets()`, `B_sets_make_unique()`,
    `B_sets_cut_increments()` are already checked and working.

  A vignette *“B-sets and interfering v-structures”* has been written to
  explain how these functions work.

# Plan for future work

## Files to be checked soon

- File `Hill_climbing.R`

## Files to be checked and rearranged

- File `Graph.R`

  - `remove_CondInd()`, `dsep_set()`, `active_cycles()`, `path_check()`
    are already checked and working.

- File `PCBN-restrictions.R`

  - `is_restrictedDAG()` is already checked and working.

- File `Utilities.R`

## Other files to be checked later

- File `Performance_metrics.R`

- File `Random_PCBN.R`
