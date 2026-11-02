#' Fit a Bayesian false-positive occupancy model
#'
#' Fits the multiple-methods false-positive occupancy model of Miller et al.
#' (2011) in Stan. Detections from an automated method may be false positives,
#' whereas detections from a verified subset are certain. The occupancy
#' probability depends on site covariates, true detection on observation
#' (or site) covariates and an optional recorder random intercept, and the
#' false-positive rate is either constant or a function of site covariates.
#'
#' @param data An `fp_data` object from [fp_data()], [fp_data_wide()], or
#'   [simulate_fp()].
#' @param occ One-sided formula for the logit of occupancy, using columns of
#'   `data$sites`.
#' @param det One-sided formula for the logit of the true-detection probability
#'   of the automated method. It may use columns of `data$auto` and of
#'   `data$sites`. It must keep the intercept.
#' @param fp One-sided formula for the logit of the false-positive rate, using
#'   columns of `data$sites` or, for covariates measured on each occasion
#'   (background noise, for instance), of `data$auto`; the rate then varies
#'   among the occasions of a site. `~ 1` (the default) fits a constant rate; `~ w`
#'   lets the rate change along `w`. Use `NULL` to fit a standard occupancy
#'   model without false positives. The intercept is fixed by the model (see
#'   Details) and cannot be removed.
#' @param random_recorder Logical. Add a recorder random intercept to true
#'   detection? Defaults to `TRUE` when `data` has a recorder column. Data
#'   from [simulate_fp()] always have one; with its default `sd_recorder = 0`
#'   set this to `FALSE`, because a random effect whose true variance is zero
#'   causes divergent transitions.
#' @param ver One-sided formula for the logit of the detection probability of
#'   the verified method, using columns of `data$verified` or `data$sites`
#'   (for example, who listened, or the length of the clip). The default,
#'   `~ 1`, keeps it constant.
#' @param markov Optional name of a numeric, date, or date-time column of
#'   `data$auto` that places the occasions of each site in time (an occasion
#'   number or a recording time). With it,
#'   true detection depends on whether the previous occasion of the same site
#'   was a detection, through the covariate `y_prev` added to `det` (reported
#'   as `det_y_prev`): a first-order Markov term for runs of detections, such
#'   as a bird calling for several consecutive minutes. See Details.
#' @param random_fp Logical. Add a site random intercept to the logit of the
#'   false-positive rate, for rates that vary among sites in ways the
#'   covariates in `fp` do not capture. `FALSE` by default; see the section
#'   "Site random effect on the false-positive rate" before turning it on.
#' @param random_det Logical. Add a site random intercept to the logit of the
#'   true-detection probability of the automated method, for detection that
#'   varies among occupied sites in ways the covariates in `det` do not
#'   capture (distance to the recorder, abundance, local noise). `FALSE` by
#'   default; it can be combined with `random_fp`. See the section "Site
#'   random effect on true detection".
#' @param priors Prior settings from [fp_priors()].
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed Passed to the
#'   `$sample()` method of 'cmdstanr'.
#' @param refresh How often each chain reports its progress, in iterations.
#'   The default, `NULL`, reports every 10% of the iterations in interactive
#'   sessions and nothing otherwise; `0` silences it. Warm-up iterations are
#'   usually slower than sampling ones, so the first half of the run takes
#'   more than half of the time.
#' @param parallel_chains Number of chains run in parallel.
#' @param threads_per_chain Number of threads each chain uses to compute the
#'   likelihood, which is summed over sites in parallel with Stan's
#'   `reduce_sum()`. One of:
#'   * `"auto"`: the physical cores divided among the parallel chains of
#'     this fit (not among other processes, such as calibration workers).
#'     The fastest choice for most data, but a fit is not reproduced
#'     exactly by its `seed`: with more than one thread the partial sums are
#'     added in an order that varies from run to run, so the draws drift
#'     apart after a while. The posterior is the same.
#'   * `1`: each chain runs on one core, and the same `seed` gives the same
#'     draws. Use it when exact reproduction matters (published analyses,
#'     tests, saved results).
#'   * a larger number: that many threads per chain, with the same caveat
#'     as `"auto"`. Keep `parallel_chains * threads_per_chain` at or below
#'     the number of cores.
#'
#'   Threads speed up fits with many sites or occasions (on the `curiango`
#'   data, 4 threads made each chain 3.4 times faster). The default is
#'   `"auto"` since crywolf 2.0.0 (it was 1 before). Under `R CMD check`
#'   (when `_R_CHECK_LIMIT_CORES_` is set, as on CRAN) a fit uses at most
#'   two cores, counting `parallel_chains * threads_per_chain`.
#' @param ... Further arguments passed to `$sample()`.
#'
#' @details
#' **Identifiability.** The intercept of the false-positive rate equals the
#' intercept of true detection minus a positive gap, \eqn{\gamma_0 = \alpha_1 -
#' \Delta}, so that at the reference values of the covariates (all zero) a
#' true detection is more likely than a false positive. This removes the label
#' switching between the occupied and unoccupied states. Centre the
#' covariates, so that the reference is an average site; [simulate_fp()]
#' already returns them standardized.
#'
#' **Verified data.** The verified detections anchor the false-positive rate.
#' With none, the model still runs, but the rate is informed only by the priors
#' and the gap, and the fit should not be trusted.
#'
#' **Observation-confirmation design.** When `data` has confirmed occasions
#' (argument `confirmed` of [fp_data()]), the model follows the
#' observation-confirmation design of Chambert et al. (2015): an expert checks
#' every observation of some automated occasions and finds only false ones,
#' only true ones, both, or none. False observations can then occur at
#' occupied sites too, and `det` and `fp` model \eqn{s_1} and \eqn{s_0}, the
#' probabilities of at least one true and at least one false observation on an
#' occasion; a site-level detection at an occupied site has probability
#' \eqn{p_{11} = 1 - (1 - s_1)(1 - s_0)}, and `p11_ref` reports it at the
#' reference. The confirmed occasions identify \eqn{s_0} directly, so verified
#' data are optional; they can be combined with this design. The gap still
#' applies, so that \eqn{s_1 > s_0} at the reference.
#'
#' **Site random effect on the false-positive rate.** With
#' `random_fp = TRUE` the logit of the false-positive rate of site \eqn{i}
#' gains a term \eqn{e_i \sim \mathrm{Normal}(0, \sigma_{fp})}, reported as
#' `sd_fp_site`, with a half-Normal prior (`fp_priors(fp_site_sd = )`).
#' In simulations with 300 sites, 20 automated occasions, and 30% of sites
#' verified, it restored the coverage of the mean false-positive rate
#' (`p10_bar`) when the rate varied among sites at random (from 0.00 and 0.55
#' with a constant rate to 0.95) and cost nothing when the rate was
#' homogeneous. Use it with care:
#' * It adds one parameter per site, and fits took 6 to 9 times longer.
#' * It mixes less well: in those simulations the largest R-hat exceeded 1.01
#'   in a quarter to three quarters of the fits. Check the diagnostics, and
#'   run longer chains (`iter_warmup`, `iter_sampling`) when they are poor.
#' * Each site has a single latent state, so a site effect on false
#'   positives is told apart from occupancy mainly by the verified data; with
#'   few occasions or little verification it may be weakly identified.
#' * It does not fix false positives that follow an environmental gradient.
#'   Those bias the occupancy effects, and the remedy is to put the
#'   occupancy covariates, and any measured correlate of the confusion, in
#'   `fp` (see below). Prefer covariates, and add the random effect for the
#'   variation that remains.
#'
#' Predictions with `predict(process = "fp")` are for a site with
#' \eqn{e_i = 0}, a typical site.
#'
#' **Site random effect on true detection.** With `random_det = TRUE` the
#' logit of the true-detection probability of the automated method at site
#' \eqn{i} gains a term \eqn{d_i \sim \mathrm{Normal}(0, \sigma_{det})},
#' reported as `sd_det_site`, with a half-Normal prior
#' (`fp_priors(det_site_sd = )`). It is for detection that varies among
#' occupied sites for reasons the covariates in `det` do not capture, such as
#' the distance from the recorder to the calling animals, their abundance, or
#' local noise. The detection probability of the verified method has no site
#' effect. `p11_ref` is the true-detection probability of a typical site
#' (\eqn{d_i = 0}).
#'
#' The two site effects answer different questions, and one cannot stand in
#' for the other: in simulations where only detection varied among sites, a
#' model with `random_fp = TRUE` alone attributed that variation to false
#' positives and estimated a large `sd_fp_site` although the true rate was
#' constant. Both can be fitted together. They are then identified by
#' different data (true detection only at occupied sites, false positives
#' only at unoccupied ones, with the verified data telling the two apart),
#' but they mix slowly, and `fp_occu()` warns about it: check R-hat and the
#' effective sample sizes, and run longer chains when they are poor.
#' `vignette("site_effects")` shows both effects in simulations and on the
#' `curiango` data.
#'
#' **Detections that come in runs.** Occasions close in time are not
#' independent when the species, once present, is active for a while: a
#' calling bird is recorded in consecutive minutes, which makes detections
#' more clustered among sites than independent occasions allow (a posterior
#' predictive check flags it). `markov` adds the detection on the previous
#' occasion as a covariate of true detection, so that \eqn{p_{11}} can be
#' higher right after a detection. It acts on the occupied state only:
#' false positives stay independent among occasions. `y_prev` is the
#' detection on the occasion one time step earlier, where the step is the
#' shortest interval between consecutive occasions of a site in the data
#' (one occasion, or 15 minutes for recordings every 15 minutes). The first
#' occasion of each site, and an occasion that follows a gap (for instance,
#' a missing occasion), has `y_prev` = 0: its probability of detection is
#' the one after no detection, so a run is not carried across a gap. Two
#' occasions of a site at the same time are an error. [fp_ppc()] checks one
#' step ahead, replicating each occasion given the observed previous one;
#' [simulate.fp_occu()] simulates the runs sequentially.
#'
#' **When to model the false-positive rate.** In the simulations behind this
#' package, a constant rate gave calibrated intervals for the occupancy slope
#' when false positives were homogeneous; random variation of the rate among
#' sites made the estimate of the rate unreliable but largely spared the slope.
#' When the rate tracked an environmental gradient correlated with the
#' occupancy covariate and the confusion was strong, the slope was biased and
#' verifying more sites did not fix it. Placing the occupancy covariate in `fp`
#' recovered part of the calibration, and placing the gradient itself
#' recovered all of it. Use [fp_calibrate()] to check your own design.
#'
#' **Parameter names.** Coefficients are labelled by process and column of the
#' design matrix: `occ_*` (occupancy), `det_*` (true detection), `fp_*` (false
#' positive, where `fp_(Intercept)` is \eqn{\gamma_0}), plus `fp_gap`
#' (\eqn{\Delta}), `sd_recorder`, `p_verified` (detection probability of the
#' verified method), `p11_ref` and `p10_ref` (true- and false-positive
#' probabilities at the reference), and `p10_bar` (false-positive probability
#' averaged over sites or automated occasions, depending on the FP covariates).
#' See [crywolf-output] for all names, scales, and availability conditions.
#'
#' @return An object of class `fp_occu`, with methods for [base::print()],
#'   [base::summary()], [stats::coef()], [stats::predict()], [loo::loo()], and
#'   [posterior::as_draws_df()]. The 'cmdstanr' fit is in element `fit`.
#' @references
#' Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
#' positive detections in species occurrence data under different study
#' designs. *Ecology*, 96, 332–339. \doi{10.1890/14-1507.1}
#'
#' Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
#' Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
#' two types of observational error occur: non-detection and species
#' misidentification. *Ecology*, 92, 1422–1428. \doi{10.1890/10-1396.1}
#' @seealso [crywolf-output], [fp_data()], [fp_priors()], [fp_calibrate()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
#'                    conf_cor = 0.7, seed = 42)
#'
#' # constant false-positive rate (the simulated recorders have no effect,
#' # so the recorder term is left out)
#' m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
#' # false-positive rate that tracks the occupancy covariate
#' m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
#'
#' coef(m0)["occ_x"]; coef(m1)["occ_x"]   # truth: -0.8
#' loo::loo_compare(loo(m0), loo(m1))
#' }
fp_occu <- function(data, occ = ~ 1, det = ~ 1, fp = ~ 1,
                    random_recorder = !is.null(data$recorder),
                    ver = ~ 1, markov = NULL, random_fp = FALSE,
                    random_det = FALSE, priors = fp_priors(),
                    chains = 4, parallel_chains = chains, threads_per_chain = "auto",
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.95, seed = NULL, refresh = NULL, ...) {
  .check_fp_data(data)
  no_ver <- is.null(data$verified) || nrow(data$verified) == 0
  no_conf <- is.null(data$confirmed) || all(is.na(data$auto[[data$confirmed[1]]]))
  if (no_ver && no_conf) {
    warning("No verified data: the false-positive rate is identified only by ",
            "the priors. Treat the estimates with caution.", call. = FALSE)
  }
  cores <- .sampler_cores(threads_per_chain, parallel_chains)
  threads_per_chain <- cores$threads_per_chain
  parallel_chains <- cores$parallel_chains
  sd_list <- .build_stan_data(data, occ, det, fp, random_recorder, priors,
                              random_fp = random_fp, ver = ver, markov = markov,
                              random_det = random_det)
  if (isTRUE(random_fp) && isTRUE(random_det)) {
    warning("Site effects on both the false-positive rate and true detection ",
            "mix slowly: check R-hat and the effective sample sizes, and run ",
            "longer chains (`iter_warmup`, `iter_sampling`) if they are poor.",
            call. = FALSE)
  }
  if (is.null(refresh)) {
    refresh <- if (interactive()) max(1, round((iter_warmup + iter_sampling) / 10)) else 0
  }
  mod <- fp_compile()
  # cmdstanr names output files by the minute plus a suffix drawn from R's
  # random numbers: after the same set.seed(), two fits in one minute share
  # file names, and discarding the first deletes the files of the second.
  # A name of our own unless the user gives one in `...`.
  run_sampler <- function(..., output_basename = .output_basename()) {
    mod$sample(..., output_basename = output_basename)
  }
  fit <- run_sampler(data = sd_list$stan_data, chains = chains,
                     parallel_chains = parallel_chains,
                     threads_per_chain = threads_per_chain,
                     iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                     adapt_delta = adapt_delta, seed = seed, refresh = refresh,
                     ...)
  # read every draw and the sampler diagnostics into memory, so that a fit
  # saved with saveRDS() still works after CmdStan's temporary files are gone
  fit$draws()
  fit$sampler_diagnostics()
  draws <- .pretty_draws(fit, sd_list)
  diag <- .diagnostics(fit, draws)

  out <- structure(
    list(fit = fit, draws = draws, stan_data = sd_list$stan_data,
         formulas = list(occ = occ, det = det, fp = fp, ver = ver),
         markov = markov,
         design = sd_list$design, data = data,
         random_recorder = isTRUE(sd_list$stan_data$estimate_re == 1),
         random_fp = isTRUE(sd_list$stan_data$estimate_fp_re == 1),
         random_det = isTRUE(sd_list$stan_data$estimate_det_re == 1),
         priors = priors,
         sampler = list(chains = chains, parallel_chains = parallel_chains,
                        threads_per_chain = threads_per_chain,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = adapt_delta), diagnostics = diag, call = match.call()),
    class = "fp_occu"
  )
  .warn_diagnostics(diag)
  out
}

# threads per chain ("auto" or a positive integer) and parallel chains,
# within the two cores that R CMD check allows
.sampler_cores <- function(threads_per_chain, parallel_chains) {
  if (identical(threads_per_chain, "auto")) {
    threads_per_chain <- .auto_threads(parallel_chains)
  }
  threads_per_chain <- suppressWarnings(as.integer(threads_per_chain))
  if (length(threads_per_chain) != 1 || is.na(threads_per_chain) ||
      threads_per_chain < 1) {
    stop("`threads_per_chain` must be a positive integer or \"auto\".", call. = FALSE)
  }
  # CRAN allows two cores under R CMD check: count chains x threads
  if (.check_limits_cores()) {
    threads_per_chain <- min(threads_per_chain, 2L)
    parallel_chains <- max(1L, min(parallel_chains, 2L %/% threads_per_chain))
  }
  list(threads_per_chain = threads_per_chain, parallel_chains = parallel_chains)
}

# a name for CmdStan's output files that no other fit in any session shares,
# made without R's random numbers so that the user's stream is untouched
.output_basename <- function() {
  .fp_env$n_fits <- if (is.null(.fp_env$n_fits)) 1L else .fp_env$n_fits + 1L
  paste0("fp_occu-", gsub("[^0-9]", "", format(Sys.time(), "%Y%m%d%H%M%OS6")),
         "-", Sys.getpid(), "-", .fp_env$n_fits)
}

# physical cores divided among the parallel chains, at most two in total
# under R CMD check
.auto_threads <- function(parallel_chains) {
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores)) cores <- 1L
  if (.check_limits_cores()) cores <- min(cores, 2L)
  max(1L, as.integer(cores %/% max(1L, parallel_chains)))
}

# TRUE when R CMD check limits the cores (set by --as-cran and on CRAN)
.check_limits_cores <- function() {
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  nzchar(chk) && chk != "false"
}

#' Prior settings for fp_occu()
#'
#' Normal priors on the logit scale. The defaults are weakly informative for
#' standardized covariates.
#'
#' @param intercept_sd Standard deviation of the Normal(0, sd) priors on the
#'   intercepts of occupancy, true detection, and verified detection.
#' @param slope_sd Standard deviation of the Normal(0, sd) priors on all
#'   slopes, including those of the false-positive rate.
#' @param gap_mean,gap_sd Mean and standard deviation of the Normal prior,
#'   truncated at zero, on the gap between the true-detection and
#'   false-positive intercepts. The default, Normal(2, 1.5), puts the
#'   false-positive rate at the reference roughly an order of magnitude below
#'   true detection, while leaving room for the data.
#' @param recorder_sd Scale of the half-Normal(0, sd) prior on the standard
#'   deviation of the recorder random intercept.
#' @param fp_site_sd Scale of the half-Normal(0, sd) prior on the standard
#'   deviation of the site random intercept on the false-positive rate
#'   (`fp_occu(random_fp = TRUE)`).
#' @param det_site_sd Scale of the half-Normal(0, sd) prior on the standard
#'   deviation of the site random intercept on true detection
#'   (`fp_occu(random_det = TRUE)`).
#'
#' @return A list of class `fp_priors`.
#' @export
#' @examples
#' fp_priors()
#' fp_priors(slope_sd = 1)
fp_priors <- function(intercept_sd = 1.5, slope_sd = 2.5,
                      gap_mean = 2, gap_sd = 1.5, recorder_sd = 2,
                      fp_site_sd = 1, det_site_sd = 1) {
  stopifnot(intercept_sd > 0, slope_sd > 0, gap_sd > 0, recorder_sd > 0,
            fp_site_sd > 0, det_site_sd > 0)
  structure(list(intercept_sd = intercept_sd, slope_sd = slope_sd,
                 gap_mean = gap_mean, gap_sd = gap_sd, recorder_sd = recorder_sd,
                 fp_site_sd = fp_site_sd, det_site_sd = det_site_sd),
            class = "fp_priors")
}

#' @export
print.fp_priors <- function(x, ...) {
  cat("<fp_priors> (logit scale)\n")
  cat(sprintf("  intercepts      ~ Normal(0, %g)\n", x$intercept_sd))
  cat(sprintf("  slopes          ~ Normal(0, %g)\n", x$slope_sd))
  cat(sprintf("  fp_gap          ~ Normal(%g, %g), truncated at 0\n", x$gap_mean, x$gap_sd))
  cat(sprintf("  sd_recorder     ~ half-Normal(0, %g)\n", x$recorder_sd))
  cat(sprintf("  sd_fp_site      ~ half-Normal(0, %g)\n", x$fp_site_sd))
  cat(sprintf("  sd_det_site     ~ half-Normal(0, %g)\n", x$det_site_sd))
  invisible(x)
}

# internal: data list for Stan ------------------------------------------------

.design <- function(form, df, what) {
  if (!inherits(form, "formula") || length(form) != 2) {
    stop("`", what, "` must be a one-sided formula, e.g. ~ x.", call. = FALSE)
  }
  vars <- all.vars(form)
  miss <- setdiff(vars, names(df))
  if (length(miss)) {
    stop("Variable(s) in `", what, "` not found: ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  }
  mf <- stats::model.frame(form, df, na.action = stats::na.pass)
  if (anyNA(mf)) {
    stop("Missing values in the covariates of `", what, "`.", call. = FALSE)
  }
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, mf)
  list(X = X, terms = tt, xlevels = stats::.getXlevels(tt, mf),
       contrasts = attr(X, "contrasts"),
       intercept = attr(tt, "intercept") == 1)
}

.build_stan_data <- function(data, occ, det, fp, random_recorder, priors,
                             random_fp = FALSE, ver = ~ 1, markov = NULL,
                             random_det = FALSE) {
  s <- data$site
  sites <- data$sites
  auto <- data$auto
  nsite <- nrow(sites)

  d_occ <- .design(occ, sites, "occ")
  if (!d_occ$intercept) stop("`occ` must keep the intercept.", call. = FALSE)

  site2 <- match(auto[[s]], sites[[s]])
  # the likelihood sums each site's occasions over contiguous rows; with a
  # Markov term they are also in time order within each site
  if (!is.null(markov)) {
    .check_df(auto, "data$auto", markov)
    ord <- order(site2, auto[[markov]])
  } else {
    ord <- order(site2)
  }
  auto <- auto[ord, , drop = FALSE]
  site2 <- site2[ord]
  has_prev <- NULL
  if (!is.null(markov)) {
    has_prev <- .markov_adjacent(site2, auto[[markov]], markov)
    # detection on the previous occasion of the same site; 0 on the first
    # occasion and after a gap in time
    prev <- c(0L, as.integer(auto[[data$y]])[-nrow(auto)])
    prev[!has_prev] <- 0L
    auto$y_prev <- prev
    det <- stats::update(det, ~ . + y_prev)
  }
  # det and fp may use site covariates: bring them into the automated table
  site_only <- setdiff(names(sites), names(auto))
  auto_full <- cbind(auto, sites[site2, site_only, drop = FALSE])
  d_det <- .design(det, auto_full, "det")
  if (!d_det$intercept) stop("`det` must keep the intercept.", call. = FALSE)

  estimate_fp <- !is.null(fp)
  # false-positive covariates measured on each occasion make the rate vary
  # within sites; otherwise it is a property of the site
  fp_by_row <- estimate_fp && !all(all.vars(fp) %in% names(sites))
  if (estimate_fp) {
    d_fp <- .design(fp, if (fp_by_row) auto_full else sites, "fp")
    if (!d_fp$intercept) {
      stop("`fp` must keep the intercept; it is tied to the detection ",
           "intercept by the gap (see ?fp_occu).", call. = FALSE)
    }
    Xfp_all <- d_fp$X[, -1, drop = FALSE]
  } else {
    d_fp <- NULL
    Xfp_all <- matrix(0, nsite, 0)
  }
  Xfp <- if (fp_by_row) matrix(0, nsite, 0) else Xfp_all

  use_fp_re <- isTRUE(random_fp)
  if (use_fp_re && !estimate_fp) {
    stop("`random_fp = TRUE` needs a false-positive process; do not use ",
         "`fp = NULL`.", call. = FALSE)
  }
  use_re <- isTRUE(random_recorder)
  if (use_re && is.null(data$recorder)) {
    stop("`random_recorder = TRUE` but `data` has no recorder column.", call. = FALSE)
  }
  if (use_re) {
    rec_f <- factor(auto[[data$recorder]])
    rec2 <- as.integer(rec_f)
    rec_levels <- levels(rec_f)
  } else {
    rec2 <- rep(1L, nrow(auto))
    rec_levels <- character(0)
  }

  ver_data <- data$verified
  if (is.null(ver_data) || nrow(ver_data) == 0) {
    y1 <- integer(0); site1 <- integer(0)
    X1s <- matrix(0, 0, 0)
    d_ver <- NULL
  } else {
    site1 <- match(ver_data[[s]], sites[[s]])
    # sorted by site, so that the rows of each site are contiguous
    o1 <- order(site1)
    ver_data <- ver_data[o1, , drop = FALSE]
    site1 <- site1[o1]
    y1 <- as.integer(ver_data[[data$y]])
    ver_full <- cbind(ver_data, sites[site1, setdiff(names(sites), names(ver_data)),
                                      drop = FALSE])
    d_ver <- .design(ver, ver_full, "ver")
    if (!d_ver$intercept) stop("`ver` must keep the intercept.", call. = FALSE)
    X1s <- d_ver$X[, -1, drop = FALSE]
  }
  if (is.null(d_ver) && length(all.vars(ver))) {
    stop("`ver` has covariates but the data have no verified occasions.", call. = FALSE)
  }

  # observation confirmation: checked occasions leave the automated rows and
  # enter with their confirmation state
  conf <- data$confirmed
  is3 <- if (is.null(conf)) rep(FALSE, nrow(auto)) else !is.na(auto[[conf[1]]])
  obs_confirm <- any(is3)
  if (obs_confirm && !estimate_fp) {
    stop("The observation-confirmation design needs a false-positive ",
         "process; do not use `fp = NULL`.", call. = FALSE)
  }
  # occasion-level structure, kept for simulating from a fit
  Xfp_occ <- if (fp_by_row) Xfp_all else matrix(0, nrow(auto), 0)
  occasions <- list(order = ord, site = site2, rec = rec2, Z = d_det$X,
                    confirmed = is3, Xfp = Xfp_occ,
                    verified_order = if (length(y1)) o1 else integer(0),
                    X1s = X1s, markov = !is.null(markov),
                    has_prev = has_prev)
  v3 <- if (obs_confirm) {
    as.integer(2L * (auto[[conf[1]]][is3] > 0) + (auto[[conf[2]]][is3] > 0))
  } else integer(0)
  site3 <- site2[is3]
  rec3 <- rec2[is3]
  Z3 <- d_det$X[is3, , drop = FALSE]
  Xfp3 <- Xfp_occ[is3, , drop = FALSE]
  auto <- auto[!is3, , drop = FALSE]
  site2 <- site2[!is3]
  rec2 <- rec2[!is3]
  X_det <- d_det$X[!is3, , drop = FALSE]
  Xfp_rows <- Xfp_occ[!is3, , drop = FALSE]

  # collapse occasions with identical covariates within a site into
  # binomial rows (exact, and much faster when detection is site-level)
  y_raw <- as.integer(auto[[data$y]])
  key <- do.call(paste, c(list(site2, rec2), as.data.frame(X_det),
                          as.data.frame(Xfp_rows), list(sep = "\r")))
  grp <- match(key, unique(key))
  first <- !duplicated(grp)
  n2 <- tabulate(grp)
  y2 <- as.integer(rowsum(y_raw, grp, reorder = FALSE))
  site2 <- site2[first]
  rec2 <- rec2[first]
  Z <- X_det[first, , drop = FALSE]
  Xfp2 <- Xfp_rows[first, , drop = FALSE]
  len2 <- tabulate(site2, nbins = nsite)
  start2 <- cumsum(c(1L, len2))[seq_len(nsite)]

  stan_data <- list(
    nsite = nsite, Kbeta = ncol(d_occ$X), X = unname(d_occ$X),
    N2 = length(n2), n2 = as.integer(n2), y2 = y2, site2 = as.integer(site2),
    start2 = as.integer(start2), len2 = as.integer(len2),
    Kalpha = ncol(Z), Z = unname(Z),
    nrec = max(1L, length(rec_levels)), rec2 = as.integer(rec2),
    estimate_re = as.integer(use_re),
    Kfp = ncol(Xfp_all), fp_by_row = as.integer(fp_by_row), Xfp = unname(Xfp),
    Xfp2 = unname(Xfp2), estimate_fp = as.integer(estimate_fp),
    N1 = length(y1), y1 = y1, site1 = as.integer(site1),
    Kver_s = ncol(X1s), X1s = unname(if (length(y1)) X1s else matrix(0, 0, 0)),
    known_present = as.integer(data$known_present),
    obs_confirm = as.integer(obs_confirm), N3 = length(v3), v3 = v3,
    site3 = as.integer(site3), Z3 = unname(Z3), rec3 = as.integer(rec3),
    Xfp3 = unname(Xfp3),
    prior_sd_int = priors$intercept_sd, prior_sd_slope = priors$slope_sd,
    gap_mean = priors$gap_mean, gap_sd = priors$gap_sd,
    recorder_sd = priors$recorder_sd,
    estimate_fp_re = as.integer(use_fp_re),
    # fits saved before this prior existed carry no value for it
    fp_site_sd = if (is.null(priors$fp_site_sd)) 1 else priors$fp_site_sd,
    estimate_det_re = as.integer(isTRUE(random_det)),
    det_site_sd = if (is.null(priors$det_site_sd)) 1 else priors$det_site_sd
  )
  design <- list(occ = d_occ, det = d_det, fp = d_fp, ver = d_ver,
                 occ_names = colnames(d_occ$X), det_names = colnames(d_det$X),
                 fp_names = if (estimate_fp) colnames(d_fp$X) else character(0),
                 ver_names = if (is.null(d_ver)) "(Intercept)" else colnames(d_ver$X),
                 recorders = rec_levels)
  list(stan_data = stan_data, design = design, occasions = occasions)
}

# internal: readable parameter names ------------------------------------------

.name_map <- function(design, estimate_fp, use_re, has_verified = TRUE,
                      use_fp_re = FALSE, use_det_re = FALSE) {
  m <- c(
    stats::setNames(paste0("occ_", design$occ_names),
                    sprintf("beta[%d]", seq_along(design$occ_names))),
    stats::setNames(paste0("det_", design$det_names),
                    sprintf("alpha[%d]", seq_along(design$det_names))),
    if (has_verified) c("delta0" = "verified_(Intercept)", "r11" = "p_verified"),
    if (has_verified && length(design$ver_names) > 1) {
      stats::setNames(paste0("verified_", design$ver_names[-1]),
                      sprintf("delta_s[%d]", seq_along(design$ver_names[-1])))
    },
    c("p11_ref" = "p11_ref")
  )
  if (estimate_fp) {
    slopes <- design$fp_names[-1]
    m <- c(m, c("gamma0" = "fp_(Intercept)"))
    if (length(slopes)) {
      m <- c(m, stats::setNames(paste0("fp_", slopes),
                                sprintf("gamma_fp[%d]", seq_along(slopes))))
    }
    m <- c(m, c("fp_gap" = "fp_gap", "p10_ref" = "p10_ref", "p10_bar" = "p10_bar"))
  }
  if (use_re) m <- c(m, c("sigma_u" = "sd_recorder"))
  if (use_fp_re) m <- c(m, c("sigma_fp[1]" = "sd_fp_site"))
  if (use_det_re) m <- c(m, c("sigma_det[1]" = "sd_det_site"))
  m
}

.pretty_draws <- function(fit, sd_list) {
  sd <- sd_list$stan_data
  map <- .name_map(sd_list$design, sd$estimate_fp == 1, sd$estimate_re == 1,
                   has_verified = sd$N1 > 0,
                   use_fp_re = isTRUE(sd$estimate_fp_re == 1),
                   use_det_re = isTRUE(sd$estimate_det_re == 1))
  dr <- posterior::as_draws_df(fit$draws(variables = unique(sub("\\[.*", "", names(map)))))
  dr <- posterior::subset_draws(dr, variable = names(map))
  posterior::variables(dr) <- unname(map[posterior::variables(dr)])
  dr
}

.diagnostics <- function(fit, draws) {
  ds <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
  summ <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
  list(
    divergences = if (is.null(ds)) NA_integer_ else sum(ds$num_divergent),
    max_treedepth = if (is.null(ds)) NA_integer_ else sum(ds$num_max_treedepth),
    max_rhat = max(summ$rhat, na.rm = TRUE),
    min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(summ$ess_tail, na.rm = TRUE)
  )
}

.warn_diagnostics <- function(diag) {
  msgs <- c(
    if (isTRUE(diag$divergences > 0))
      sprintf("%d divergent transition(s); consider a higher `adapt_delta`.",
              diag$divergences),
    if (isTRUE(diag$max_rhat > 1.01))
      sprintf("Largest R-hat is %.3f (> 1.01); run longer chains.", diag$max_rhat),
    if (isTRUE(diag$min_ess_bulk < 400))
      sprintf("Smallest bulk ESS is %.0f (< 400).", diag$min_ess_bulk)
  )
  if (length(msgs)) warning(paste(msgs, collapse = "\n"), call. = FALSE)
  invisible(NULL)
}

# Is each occasion (sorted by site and time) one time step after the previous
# row of the same site? The step is the smallest interval between consecutive
# occasions in the data, so missing occasions (e.g. rows with NA dropped by
# fp_data()) break the chain instead of joining non-adjacent occasions.
.markov_adjacent <- function(site, time, name) {
  if (!is.numeric(time) && !inherits(time, c("POSIXt", "Date"))) {
    stop("The `markov` column '", name, "' must be numeric, a date, or a ",
         "date-time.", call. = FALSE)
  }
  t <- as.numeric(time)
  if (anyNA(t)) stop("The `markov` column '", name, "' has missing values.",
                     call. = FALSE)
  n <- length(t)
  same <- c(FALSE, site[-1] == site[-n])
  dt <- c(NA, diff(t))
  if (any(same & dt == 0)) {
    stop("Two occasions of the same site share the same `markov` value; ",
         "their order in time is undefined.", call. = FALSE)
  }
  if (!any(same)) return(same)
  step <- min(dt[same])
  same & abs(dt - step) <= 1e-8 * step
}
