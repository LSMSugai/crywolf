#' Stability of the crywolf interface
#'
#' What crywolf promises not to change between versions, and which
#' functions are still experimental.
#'
#' @section What is frozen:
#' From crywolf 2.0.0, for the stable functions listed below:
#'
#' * **Arguments**: their names, their order and their default values,
#'   and those of the methods of crywolf's classes. A minor version may add
#'   an argument after the existing named ones (and before a trailing
#'   `...`), with a default that keeps the earlier behaviour.
#' * **Classes**: the crywolf classes their results inherit from
#'   (`fp_data`, `fp_occu`, `fp_priors`, `fp_sim`, `fp_calibration`,
#'   `fp_sbc`, `fp_ppc`, `fp_prior_predictive`, `fp_effects`, `fp_sites`,
#'   `fp_threshold_sweep`, `fp_sheepdog`, `summary.fp_occu`,
#'   `summary.fp_calibration`, `summary.fp_sbc`), and the class of the
#'   design warnings of [fp_data()] (`crywolf_design_warning`), which code
#'   may catch. Classes of other packages (the `CmdStanModel` of
#'   [fp_compile()], the `psis_loo` of `loo()`) are theirs to keep.
#' * **Results**: the names documented in [crywolf-output] (parameters,
#'   derived quantities, columns of summaries and of the other results
#'   described there) and the components and columns described in the
#'   Value section of each stable function, for example the tables of
#'   [fp_read_validation()], [fp_read_sheepdog()] and [fp_aggregate()]. New
#'   columns or components may be added. Columns that a reader passes on
#'   from its input unchanged (the selection-table columns of
#'   [fp_read_birdnet()]; the `scores` and `sites` tables and the `schema`
#'   of [fp_read_sheepdog()]) follow BirdNET and 'sheepdog' and are not
#'   covered.
#'
#' These change only in a major version, announced in the NEWS of the
#' previous minor version. Not covered: arguments passed through
#' `...` to 'cmdstanr', the fields of `stan_data` and `design`, raw Stan
#' names, the text of messages and warnings, and the numerical results
#' themselves, which may move within Monte Carlo error when the sampler or
#' its settings change (for example `threads_per_chain`, whose default
#' became `"auto"` in 2.0.0). A test of the package (run on its continuous
#' integration, not on CRAN) records the arguments of every stable function
#' and method, so that a change cannot pass unnoticed.
#'
#' @section Stable functions:
#' Data and reading: [fp_data()], [fp_data_wide()], [fp_read_birdnet()],
#' [fp_read_validation()], [fp_read_sheepdog()], [fp_aggregate()],
#' [fp_threshold_sweep()].
#'
#' Fitting and checking: [fp_occu()], [fp_priors()], [fp_gap_prior()],
#' [fp_compile()], [fp_stan_file()], [fp_prior_predictive()], [fp_ppc()],
#' [fp_loo_weights()], the methods in [fp_occu-methods] (`print()`,
#' `summary()`, `coef()`, `predict()`, `loo()`, `as_draws_df()`) and
#' [simulate.fp_occu()].
#'
#' Results and simulation: [fp_effects()], [fp_sites()], [simulate_fp()],
#' [fp_calibrate()], [fp_verification_curve()], [fp_calibrate_fit()],
#' [fp_sbc()].
#'
#' The `print()`, `summary()` and `plot()` methods of the classes above
#' are stable too.
#'
#' @section Experimental functions:
#' [fp_occu_mc()], [fp_draws_data()] and [fp_to_unmarked()], with the
#' `fp_occu_mc` class and its methods, are new in 2.0.0. They propagate a
#' classifier calibration from 'sheepdog' into the model or pass data to
#' 'unmarked', and are still being evaluated. Their arguments, classes and
#' results may change in a minor version; the NEWS will say so. Each says
#' so in its help. [fp_occu_dyn()], the dynamic (multi-season) model, with
#' the `fp_occu_dyn` class and its methods, is experimental in the same
#' way; the `season` and `season_covs` arguments of [fp_data()] that feed
#' it are stable. [fp_read_sheepdog()], by contrast, reads a versioned
#' export format (`sheepdog-occupancy`, schema 1) and is stable.
#'
#' @name crywolf-api
#' @keywords documentation
NULL
