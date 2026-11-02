# crywolf (development version)

## Dynamic model (experimental)

* New `fp_occu_dyn()`, experimental: the false-positive occupancy model
  over several seasons, with initial occupancy (`psi1`), colonization
  (`col`) and extinction (`ext`) formulas on site and season covariates.
  Within each season the observation model is that of `fp_occu()`, with
  all its options. The states are summed out by the forward algorithm.
  It has `print()`, `summary()` and `as_draws_df()` methods; its output
  names (`psi1_*`, `col_*`, `ext_*`) are documented in `crywolf-output`.
* With a `fp_occu_dyn()` fit, `fp_sites()` returns one row per site and
  season, with the occupancy expected from the covariates (`psi`) and the
  probability that the site-season is occupied given all the data of the
  site (`occupied`, by forward-backward over seasons). `fp_ppc()` checks
  detections per site-season, with the latent states drawn jointly over
  seasons by forward filtering, backward sampling, and `simulate()` draws
  new states and detections with the structure of the data.
* `fp_occu_dyn()` gives derived quantities per season: expected
  occupancy over sites (`psi_mean[t]`), realized occupancy with the states
  drawn from their joint posterior (`psi_fs[t]`), turnover (`turnover[t]`)
  and, with intercept-only colonization and extinction, the equilibrium
  occupancy (`psi_eq`). `summary()` lists them with a `season` column.
* New `simulate_fp_dyn()`, experimental: simulates data for
  `fp_occu_dyn()`, with colonization on a season covariate and initial
  occupancy and extinction on a site covariate.
* New vignette `dynamic`: `fp_occu_dyn()` on the `curiango` data with the
  nights of each site as seasons (use from one night to the next),
  compared with the model of independent site-nights by cross-validation
  on blocks of nights.
* `fp_data()` gains `season_covs`, covariates of each season or of each
  site and season.
* `fp_data()` fills in seasons that are dates (`Date`, or text as
  `YYYY-MM-DD` or `YYYYMMDD`, like `night` and `recording_date` of a
  'sheepdog' export): the dates missing between the first and the last are
  added as seasons without data, so that every transition spans the same
  time (one night for nightly seasons). It warns when dates or numbers are
  unevenly spaced and cannot be filled in.
* `fp_data()` gains `season`, the column of the season (primary period)
  of each automated and verified occasion, for `fp_occu_dyn()`.
  Seasons follow the levels of a factor or the sorted values, so `night`
  or `recording_date` of a 'sheepdog' export work as they are. `print()`
  shows the number of seasons and of sites with data in each, and
  `fp_data()` warns about seasons without data at any site and about
  seasons without verified or confirmed data. `fp_occu()` ignores the
  column.
* The Stan program keeps the likelihood of each site's observations given
  that it is occupied and given that it is not (its emissions) in
  `fp_emissions.stanfunctions`, which `fp_occu.stan` includes and
  `fp_occu_dyn.stan` shares. `fp_occu()` gives the same draws as before.

## Integration with 'sheepdog' (experimental)

* `fp_draws_data()` reads exports of schema 1.3, whose calibration has
  continuous covariates (`schema$calibration$covariates`, such as a
  standardised temperature of the night). The probability that a
  prediction is true then depends on its occasion, so sheepdog draws the
  histories when it exports and writes no calibration draws:
  `fp_draws_data()` reads those histories as they are, takes `rule`,
  `p_min` and `seed` from the export when they are left at their defaults,
  and stops with an error when asked for other values or more draws.
  crywolf 2.1.0 stops with an error on these exports, rather than ignoring
  the covariates.
* The `print()` method of `fp_read_sheepdog()` shows the formula of the
  calibration.

# crywolf 2.1.0

## Integration with 'sheepdog' (experimental)

* `fp_read_sheepdog()` reads exports of the `sheepdog-occupancy` schema
  1.1 (any 1.x; 2.x is still refused). Its occasion covariates
  (`schema$occasion_covariates`) are extra columns of `occasions`, and the
  keys of the calibration are read as text.
* `fp_draws_data()` matches every key of the calibration
  (`schema$calibration$keys` beyond site and species, such as month,
  season or recorder) by the column of that name in the occasions of the
  export. An occasion without a value of a key is treated as not surveyed,
  like one without calibration draws. Exports of schema 1.0 calibrated by
  month still take the month from the recording date.
* `fp_draws_data()` reads the detection histories that sheepdog drew
  (`occasion_draws.parquet`, schema 1.2) when they were drawn for the same
  species, rule, seed or `p_min`, and occasions; they are the histories it
  would draw itself.

# crywolf 2.0.0

The first major version: the interface of the stable functions is frozen,
the Stan model is compiled when the package is installed (the route to
CRAN), and fits use several threads per chain by default. Three
experimental functions start the integration with 'sheepdog'.

## Stable interface

* The interface of the stable functions is frozen (`help("crywolf-api")`),
  extending the freeze of the names in `help("crywolf-output")` made in
  1.9.0: the names, order and defaults of the arguments of these functions
  and of their methods, the crywolf classes of their results (and of the
  design warnings of `fp_data()`), and the components and columns
  described in their help. These change only in a major version, announced
  in the NEWS of the previous minor version; a minor version may add
  arguments and columns. A test records the arguments of every stable
  function and method. `fp_occu_mc()`, `fp_draws_data()` and
  `fp_to_unmarked()` are experimental and say so in their help.

## Breaking changes

* `threads_per_chain` of `fp_occu()` defaults to `"auto"`, as announced in
  1.9.0: each chain uses its share of the physical cores. Fits are faster,
  but with more than one thread a fit is no longer reproduced exactly by
  its `seed` (the posterior is the same). Set `threads_per_chain = 1` to
  keep exact reproduction, for example to rerun an analysis made with an
  earlier version. The once-per-session notice of 1.9.0 and the advice
  after slow one-thread fits are gone. `fp_calibrate()`,
  `fp_calibrate_fit()` and `fp_sbc()` keep one thread per chain unless
  `fit_args` says otherwise, since their fits run side by side in parallel
  workers.

## Installation

* The Stan model is compiled when the package is installed from source
  (through 'instantiate'), with threading, if CmdStan 2.40 or later is
  found; otherwise, as for binary packages, `fp_compile()` compiles it into
  the user cache on first use, as before. 'cmdstanr' moves from Imports to
  Suggests (it is not on CRAN) and 'instantiate' joins Imports; the Stan
  program moves from `inst/stan/` to `src/stan/`, and `fp_stan_file()`
  finds it in either case. Draws are unchanged.

## New features (experimental)

* `fp_occu_mc()` fits the model to several versions of the same data, such
  as detection histories drawn from the posterior of a classifier
  calibration (`sheepdog::history_draws()`), and pools the posteriors in
  equal shares (Monte Carlo propagation; a "cut" model, in which the
  calibration does not learn from occupancy). `summary()` splits each
  variance into parts within and between fits, and reports convergence fit
  by fit. With a single version and one thread per chain it reproduces
  `fp_occu()` draw for draw. Fits that did not converge (largest R-hat
  above `max_rhat`, 1.05 by default) are left out of the pool with a
  warning that names them.
* `fp_draws_data()` builds those versions from an export of sheepdog with
  calibration draws (`calibration_draws.parquet`): in each draw a
  prediction is judged true with its calibrated probability (`rule =
  "sample"`, the default) or when that probability reaches `p_min` (`rule =
  "threshold"`). The sites of the model can be other than the recorders
  (for instance, site-nights). A calibration by month (`"month"` among the
  keys of the export) is matched by the month of recording; an occasion with
  a prediction that has no calibration is left unsurveyed (`NA`), not counted
  as a non-detection. Needs the 'arrow' package (Suggests).
* `fp_to_unmarked()` turns an `fp_data` object into the
  `unmarkedFrameOccuFP` of `unmarked::occuFP()`, the reverse of
  `fp_data_wide()`: verified occasions first, a `METH` covariate that tells
  the methods apart, and optional observation covariates of the automated
  occasions. With `fp_read_sheepdog()` and `fp_aggregate()`, an export of
  sheepdog can be fitted by maximum likelihood as well.

# crywolf 1.9.0

This release freezes the API on the way to crywolf 2.0.0, and adds a site
random effect on true detection.

## API

* The names documented in `help("crywolf-output")` (readable parameter
  names, the columns of the results described there, the fields of the
  fitted object, now including `random_det`) are frozen from this version:
  they change only in a major version, announced one minor version
  earlier. New options may add names without changing existing ones.
  Results not described on that page are not covered yet.

## Upcoming changes

* The default of `threads_per_chain` in `fp_occu()` changes from 1 to
  `"auto"` in crywolf 2.0.0. Fits will be faster, but a fit will no longer
  be reproduced exactly by its `seed` (the posterior is the same), so
  rerunning an old analysis with the same seed will give slightly
  different draws. A fit that does not set `threads_per_chain` now says so
  once per session; set `threads_per_chain = 1` to keep exact reproduction.
  `?fp_occu` explains the three choices (`"auto"`, 1, or a number).
* `fp_read_validation()` will stop reading the legacy format of the sampler
  (`validacao_por_clipe.csv`, column names in Portuguese) in crywolf 3.0.0.
  It is read as before throughout crywolf 2.x, with no warning at run time.

## New features

* `fp_read_sheepdog()` reads the folder written by
  `sheepdog::export_occupancy_data()` (schema `sheepdog-occupancy`,
  version 1): occasions, predictions, reviewed clips with their sampling
  design, and recordings annotated by ear, in the forms of `fp_aggregate()`
  and `fp_data()`. It drops clips drawn with `sampling_mode = "top"` with a
  warning (`keep_top = TRUE` keeps them) and refuses other schemas and
  versions. Needs the 'jsonlite' package (Suggests).
* `fp_occu(random_det = TRUE)` adds a site random intercept to the logit of
  the true-detection probability of the automated method, for detection
  that varies among occupied sites in ways the covariates in `det` do not
  capture. It is off by default, reported as `sd_det_site`, with the prior
  `fp_priors(det_site_sd = 1)`, and can be combined with `random_fp`; the
  combination mixes slowly, and `fp_occu()` warns about it. With the option
  off the draws are identical to those of the previous version with the same
  seed. `loo()` integrates the site effects out of each site's
  log-likelihood also with `random_det`, with a two-dimensional
  Gauss-Hermite grid of 30 nodes per effect when both site effects are on.
  `simulate()`, `fp_ppc()`, `fp_calibrate_fit()`, `fp_sbc()` and
  `fp_prior_predictive()` support it. The new article
  `vignette("site_effects")` shows in simulations (100 replicates of four
  scenarios) what each site effect recovers and costs, a simulation-based
  calibration, and the `curiango` data, where leave-one-site-out comparison
  prefers `random_det` together with the run term of `vignette("markov")`.
* `vignette("threshold")` reports refits at thresholds 0.2 and 0.4 with
  `fp = ~ temp + rh`.

## Changes

* Under `R CMD check` (when `_R_CHECK_LIMIT_CORES_` is set, as on CRAN), a
  fit uses at most two cores: `threads_per_chain` is capped at 2 and
  `parallel_chains` at `2 %/% threads_per_chain`. Fewer parallel chains do
  not change the draws, only the wall time: the chains run one after
  another.
* With `random_fp = TRUE`, `loo()` now integrates each site's effect on the
  false-positive rate out of its log-likelihood, by Gauss-Hermite quadrature
  (new arguments `integrate = TRUE` and `nodes = NULL`, which means 100
  nodes with one site effect), before the leave-one-out computation.
  A left-out site's effect is not known from its own data, so this is the
  right predictive density; conditional on the fitted effects, many sites
  had Pareto k above 0.7 (in the example of `vignette("fp_variation")`, 0
  integrated against 26 to 29 conditional). `integrate = FALSE` gives the
  previous values. The exact refits of `reloo = TRUE` integrate the
  left-out site's effect the same way, and no longer draw random numbers.
* Tests and vignettes skip their model fits when the CmdStan found is older
  than 2.40, instead of failing in `fp_occu()`.

# crywolf 1.7.0

## Requirements

* CmdStan 2.40 or later is now required (it was 2.34). The data of the
  likelihood now travel in tuples, which `reduce_sum()` accepts only from
  CmdStan 2.40. Update with `cmdstanr::install_cmdstan()`. The posterior is
  unchanged: at the same parameter values the log-likelihood is exactly
  equal to that of crywolf 1.6.0 in seven designs, and with the same seed
  the draws are identical in six (`dev/profile/validate-tuples*.R`).

## Documentation

* `help("crywolf-output")` gathers readable parameter names, their scales
  and availability, posterior summaries, diagnostics, predictions, site
  estimates, checks, calibration, simulation, and raw Stan mappings.

## New features

* `fp_threshold_sweep()` can derive MaxF1 and MaxPrecision thresholds from
  reviewed classifier clips. Ties favour recall, then the smaller threshold;
  selected thresholds and validation metrics are retained as attributes.
  Numeric thresholds remain supported and duplicate thresholds are fitted once.

## Bug fixes

* `print()` of an `fp_data` object counted the sites known to be occupied
  (from verified detections and confirmed true observations together) on
  the line of verified sites, so data with only confirmations showed "0
  verified sites ... 24 with a confirmed detection". Verified sites,
  confirmed occasions, and sites known to be occupied now have their own
  lines.
* With `parallel = TRUE`, the workers of 'furrr' use the L'Ecuyer-CMRG
  generator, so `simulate_fp(seed = s)` inside `fp_calibrate()`,
  `fp_verification_curve()`, `fp_calibrate_fit()`, and `fp_sbc()` simulated
  other datasets than a sequential run with the same seeds. Every function
  that takes a `seed` now sets R's default generators along with it, so the
  two modes give identical results. In a session with the default
  generators nothing changes: `simulate_fp()` still reproduces the datasets
  of the fpcalib compendium exactly.
* Two fits made within the same minute after the same `set.seed()` could
  share CmdStan's output file names, which 'cmdstanr' draws from R's random
  numbers; discarding the first fit then deleted the files of the second,
  and reading its draws failed with "File does not exist" or "Supplied CSV
  file is corrupt". `fp_occu()` now gives every fit its own file names,
  without using R's random numbers.

## Tests

* The precomputed results of the vignettes carry the MD5 of the Stan
  program that made them, and a vignette shows a note when its results come
  from another version of the program. `dev/precompute-all.R` refreshes them
  all; `dev/restamp-precomputed.R` marks them current without recomputing,
  only after the posterior regression test passes. The MD5 is of the whole
  file, so even an edit to a comment in the Stan program shows the note
  until the results are restamped.
* A posterior regression test checks, in four designs, that the
  log-likelihood at saved parameter values is unchanged and that posterior
  means stay within Monte Carlo error of saved references
  (`tests/testthat/test-regression.R`; references made by
  `dev/make-regression-refs.R`).

## New articles

* "From BirdNET output to an occupancy model": BirdNET selection tables and
  'sheepdog' verdicts to the data of `fp_data()`.
* "Choosing a classifier threshold": `fp_threshold_sweep()` on the
  `curiango` data.
* "Reporting a crywolf analysis": the items of the Bayesian Analysis
  Reporting Guidelines, with a template for methods and results.
* "Detections that come in runs": the Markov term on the `curiango` data.
* "False-positive rates that vary among sites and occasions": `random_fp`
  and occasion-level covariates on simulated data.
* "Planning a study": decisions to make before recording starts.
* "Performance": computing time over the versions, threads, and what did
  not help.

Results that need long fits are precomputed by `dev/precompute-markov.R`,
`dev/precompute-threshold.R`, and `dev/precompute-fp-variation.R`.

## Installation

* 'sheepdog', suggested for reading clip reviews, is now installed from
  https://diogoprov.r-universe.dev (listed in `Additional_repositories`)
  instead of from GitHub through `Remotes`. The README says how.

# crywolf 1.6.0

## Speed

* `fp_occu(threads_per_chain = "auto")` uses the physical cores divided
  among the parallel chains. On the `curiango` data, 150 iterations of one
  chain took 55 s with one thread, 16 s with four, and 10 s with eight. The default stays at
  one thread, because with more the order of the partial sums varies from
  run to run and a fit is no longer reproduced exactly by its `seed` (the
  posterior is the same). A fit whose chains took over a minute on one core
  now says so and suggests `"auto"`.
* `fp_calibrate()` documents that two chains per fit give an effective
  sample size of about 2000 and, with the chains of a fit run one after
  another as in a parallel calibration, half the cost of four.
* Measured and not adopted (`dev/profile/rejected-2026-10-21.md`): a dense
  mass matrix, sampling the false-positive intercept with the gap, a
  centred recorder effect, and `reduce_sum_static()`. Considered and not
  pursued, from the same measurements: initialization by Pathfinder
  (warm-up is 40% of the time) and a QR decomposition (no posterior
  correlation among covariate effects above 0.5). The confirmation design
  keeps its current likelihood path, for lack of a large dataset to
  measure it on.
* With `random_fp = TRUE`, 3 of 20 fits of the unstructured regime had an
  R-hat above 1.01 with the default settings; a warm-up of 2000 (2 of 20)
  or `adapt_delta = 0.99` (3 of 20) took longer for no clear gain, so the
  defaults stay.

## Data

* `fp_read_validation()` reads `validation_by_clip.csv` from 'sheepdog'
  2.0.0 and later (English column names and verdicts: `site`,
  `source_audio`, `clip_start`, `verdict = "positive"`, ...) as well as the
  legacy `validacao_por_clipe.csv` of the sampler up to 1.2.0. The format is
  detected from the columns, and both give the same output, which now
  keeps `occasion_id` when the file has it.

# crywolf 1.5.0

## Model options

All three are off by default; with them off, the draws are identical to
those of crywolf 1.4.0.

* `fp_occu(ver = ~ ...)` gives the verified method covariates of its own
  (from the verified table or the sites), reported as `verified_*`.
* False-positive formulas can use covariates measured on each occasion
  (columns of the automated table), so the rate varies among the occasions
  of a site.
* `fp_occu(markov = "time")` adds the detection on the previous occasion of
  the same site as a covariate of true detection (`det_y_prev`), for
  species recorded in runs. The previous occasion must be one time step
  earlier: after a gap (a missing occasion) a new run starts. `simulate()`
  simulates the runs in order, and `fp_ppc()` checks them one step ahead.

`dev/validate-model-options.R` checks the identity with 1.4.0 in three
designs, compares the R likelihood used by `fp_ppc()`, `fp_sites()`, and
`simulate()` with the Stan `log_lik` in six combinations of the options,
and checks in simulations that each option's effect is recovered; its
output is in `dev/validate-model-options.txt`.

## Checking and comparing models

* `fp_sbc()`: simulation-based calibration of the computation (Talts et al.
  2018; Modrák et al. 2025), with the log-likelihood of the data as a test
  quantity besides the parameters, and `summary()` and `plot()`.
* `loo(reloo = TRUE)` refits the model without each site whose Pareto k
  exceeds a threshold and uses the exact leave-one-out value for it. A site
  that is the only one with its recorder, or with a level of a factor,
  keeps its approximate value, with a warning.
* `fp_loo_weights()`: stacking or pseudo-BMA weights of several fits.
* `fp_gap_prior()` turns the precision measured on a random sample of
  reviewed detections into a prior for the gap.
* `predict(process = "fp")` predicts for the automated occasions by default
  when the false-positive rate has occasion covariates.
* `fp_threshold_sweep()` fits the same model at several classifier
  thresholds and plots the estimates against the threshold.

# crywolf 1.4.0

## Plots and site-level estimates

* `plot()` for fits draws the posterior intervals of the coefficients, one
  panel per process, or sampler diagnostics (`type = "trace"`, `"rank"`,
  `"pairs"`, via 'bayesplot', with divergent transitions marked).
* `fp_effects()` gives the occupancy probability or the false-positive rate
  along one covariate, with the others at their means, and a `transform`
  argument for covariates standardized before fitting; `plot()` draws it.
* `fp_sites()` returns, for every site, the occupancy probability, the
  posterior probability that it is occupied given its data, and its
  false-positive rate; `plot()` ranks the sites or, with coordinates, maps
  them.
* `plot()` for `fp_ppc()` results: replicated against observed discrepancy,
  or the distribution of detections per site.
* `fp_prior_predictive()` simulates data from the priors with the design of
  the real data and compares simple summaries with the observed ones.
* `plot()` for calibrations can show verification effort as the expected
  number of verified sites or occasions, and mark a target coverage.

## Reading classifier output

* `fp_read_birdnet()` reads BirdNET selection tables (Raven format).
* `fp_read_validation()` reads the clip verdicts of the
  birdnet-validation-sampler.
* `fp_aggregate()` turns detections into automated occasions at a
  confidence threshold, and verdicts into counts for the
  observation-confirmation design, warning when a checked occasion has
  detections that were not reviewed.
* `fp_data()` warns about sites without automated occasions and about
  automated data that are all detections or all non-detections.

## Other changes

* Fits keep every draw in memory, so a fit saved with `saveRDS()` works in
  later sessions, after CmdStan's temporary files are gone.

# crywolf 1.3.0

## Faster likelihood, same posterior

The Stan program evaluates the likelihood up to three times faster, with
the log density, its gradient and `log_lik` unchanged up to floating-point
rounding (checked in nine designs; scripts and results in `dev/profile/`
of the repository). Executables compiled by earlier versions are rebuilt
once, and CmdStan 2.34 or later is now required.

* The unoccupied branch of each site is computed from its totals of
  automated occasions and detections, because the false-positive rate is
  constant within a site.
* When `det` has covariates, the occupied branch uses one
  `binomial_logit_glm_lpmf()` call per site, which evaluates the linear
  predictor, the log-likelihood and its gradient in one pass instead of
  creating an autodiff variable per occasion.
* With `random_fp = FALSE` the site effects on the false-positive rate enter
  as data rather than as parameters fixed at zero.

The gain grows with the number of occasion-level rows per site. Time per
gradient evaluation with one thread, same parameter values:

| Data | Rows | 1.2.0 | 1.3.0 | Gain |
|---|---|---|---|---|
| 300 sites x 20 occasions, `det = ~ 1` (rows collapsed) | 300 | 0.064 ms | 0.056 ms | 1.1x |
| 300 sites x 20 occasions, `det = ~ bio + aci`, recorder effect | 6,000 | 0.58 ms | 0.23 ms | 2.5x |
| `curiango`, `det = ~ bio + aci`, recorder effect | 97,806 | 8.5 ms | 2.6 ms | 3.3x |

In a short run of the sampler on the `curiango` data (four chains), each
leapfrog step took 2.3 times less with one thread and 2.6 times less with
four threads per chain.

# crywolf 1.2.0

* `fp_occu(random_fp = TRUE)` adds a site random intercept to the logit of
  the false-positive rate, for rates that vary among sites in ways the
  covariates do not capture; its prior scale is `fp_priors(fp_site_sd = )`.
  It is off by default: it slows the fits and mixes less well, and it does
  not replace covariates when false positives follow an environmental
  gradient (see `?fp_occu`). With the option off, the draws are identical
  to those of crywolf 1.1.0.
* `fp_calibrate_fit()` checks the calibration of candidate models on data
  simulated from a fitted model, which keep the design of the real data
  (sites, covariates, occasions, recorders, and verification).
* `simulate()` method for `fp_occu` fits, from the posterior means or from
  random posterior draws.
* The summary of a calibration prints every occupancy slope, not only `occ_x`.

# crywolf 1.1.0

* The likelihood is summed over sites with Stan's `reduce_sum()`, and the
  model is compiled with threading, so `fp_occu(threads_per_chain = )` can
  split each chain among several cores. On the curiango data (97,806
  automated occasions) a fit with four chains of four threads each took
  16.5 minutes instead of 70, and each gradient evaluation was 2.7 times faster than with the
  previous model in a short benchmark on the same computer. The log
  density and its gradient are unchanged up to floating-point rounding, and
  the automated occasions use `y * lp - n * log1p_exp(lp)`, an equivalent
  form of the binomial log-likelihood with half the transcendental
  functions. Executables compiled by earlier versions are rebuilt once.

# crywolf 1.0.0

First release, accompanying the simulation study of Provete et al.
(fpcalib research compendium, v0.1.0, doi:10.5281/zenodo.23379100).

## Models

* `fp_occu()` fits the multiple-methods false-positive occupancy model of
  Miller et al. (2011) in Stan, with any number of occupancy, detection, and
  false-positive covariates, an optional recorder random intercept, and
  configurable priors (`fp_priors()`). `fp = NULL` gives the standard
  occupancy model.
* Observation-confirmation design of Chambert et al. (2015): argument
  `confirmed` of `fp_data()` and `frac_confirmed` of `simulate_fp()`.
* Methods `print`, `summary`, `coef`, `predict`, `loo`, and `as_draws_df`.
* `fp_ppc()`: posterior predictive check with a Bayesian p-value.

## Simulation and calibration

* `simulate_fp()` generates homogeneous, unstructured, and structured
  false-positive regimes; with the same seed it returns the same datasets
  as the simulator of the fpcalib compendium.
* `fp_calibrate()` and `fp_verification_curve()` report the coverage and
  bias of credible intervals for a design, with a `plot()` method and a
  progress bar that shows the time remaining (needs 'progressr'). They fit
  the recorder random intercept only when the simulated design has one,
  which avoids divergent transitions.

## Data and documentation

* `curiango` dataset: acoustic monitoring of the Common Pauraque in the
  Pantanal (Piatti et al. 2026).
* Articles: introduction, simulation, calibration, a case study with the
  `curiango` data, coming from `unmarked::occuFP()`, and the Stan model.
