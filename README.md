Development of the PCBN package
================

# Files that have been checked and are fully working

- File `PCBN-object.R`

  - `new_PCBN()`: working, but may be changed / developed further. So
    the API is still not stable.

  - `plot.PCBN()` and `print.PCBN()` are working but are somehow basic
    and can be improved.

- File `PCBN-logLik.R`

  - `logLik.PCBN()` and `PCBN_PDF` are working. TODO: add methods `AIC`
    and `BIC`.

- File `Simulation.R`

  - `sample_PCBN()`
  - `compute_sample_margin()`

  Both functions are finished, but will return an error if the PCBN does
  not satisfy the restrictions. This is the case for example if the
  order of the parents are not chosen in the correct way (even though
  the graph may not have any active cycle nor interfering v-structure).

- File `PCBN-tag.R` : contains some utilities about the tag/hashmap
  system for storing estimated copulas and estimated conditional
  margins.

  The functions are: `paste_margin()`, `unpaste_margin()`,
  `print_key_keychain()`, `make_and_store_keyCopula()` and
  `make_and_store_keyMargin()`.

- File `Graph_B_sets.R`

  - `has_interfering_vstrucs()`
  - `find_B_sets()`, `find_B_sets_v()`
  - `find_interfering_v()`
  - `B_sets_make_unique()`
  - `B_sets_cut_increments()`

  are already checked and working.

- File `CopulaAssignement_possibleCandidates.R`

  - `possible_candidates()`
  - `possible_candidate_incoming_arc()`
  - `possible_candidate_outgoing_arc()`

  are already checked and working.

# Plan for future work

- File `Estimation.R` : rewriting in progress.

  - `BiCopCondFit()` and `ComputeCondMargin()` have been rewritten.

  A vignette *“How to use the estimation procedures”* has been written
  to explain how these functions work.

  - **TODO:** `fit_copulas()` and `fit_all_orders()` will be checked
    soon.

## Files to be checked soon

- File `Hill_climbing.R`

## Files to be checked and rearranged

- File `CopulaAssignement.R`

  - `is_cond_copula_specified()`
  - `find_cond_copula_specified()`
  - `complete_and_check_orders()`

  are already checked and working.

- File `Graph.R`

  - `remove_CondInd()` is already checked and working.

- File `Utilities.R`

## Other files to be checked later

- File `Performance_metrics.R`

- File `Random_PCBN.R`
