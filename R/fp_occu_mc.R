#' Fit the model to draws of the detection histories
#'
#' Fits [fp_occu()] to each of several versions of the same data, which
#' differ only in which occasions count as detections, and pools the
#' posterior draws in equal shares. The versions are typically detection
#' histories drawn from the posterior of a calibration of the classifier,
#' as `sheepdog::history_draws()` produces, so that the uncertainty about
#' which detections are true reaches the occupancy estimates (Monte Carlo
#' propagation).
#'
#' @param data_list A list of `fp_data` objects with the same sites, one per
#'   draw of the detection histories.
#' @param ... Arguments passed to [fp_occu()] for every fit (formulas,
#'   priors, sampler settings).
#' @param seed Optional integer: fit `i` uses seed `seed + i - 1`. The
#'   fits are reproduced draw for draw only with one thread per chain
#'   (`threads_per_chain = 1`, the setting when `parallel = TRUE`).
#' @param parallel Fit the versions in parallel with 'furrr', on the
#'   workers set by [future::plan()]. Unless set in `...`,
#'   `parallel_chains` and `threads_per_chain` are then 1, so that each
#'   worker uses one core.
#' @param max_rhat Fits whose largest R-hat exceeds this value are left out
#'   of the pooled draws, with a warning that names them. `Inf` keeps every
#'   fit.
#' @param keep_fits Keep the [fp_occu()] fit of every version in the
#'   result (`fits`)? They hold the backend objects and can be large. Their
#'   draws are in memory, but with `parallel = TRUE` the CmdStan output
#'   files of a worker may be gone once it ends.
#'
#' @details
#' **What the pooled posterior is.** Each fit conditions on one version of
#' the detection histories; pooling the fits in equal shares averages the
#' posterior over the versions. When the versions are draws from the
#' posterior of a calibration, this approximates the posterior that accounts
#' for the uncertainty of the calibration, but the calibration itself is not
#' updated by the occupancy model (a "cut" model): information about the
#' state of a site does not flow back to how likely its detections are to
#' be true.
#'
#' **Between and within.** `summary()` splits the posterior variance of
#' each quantity into the part within fits (uncertainty given one version of
#' the histories) and the part between fits (uncertainty about the
#' histories), by the law of total variance: the squares of `sd_within`
#' and `sd_between` add up to the variance of the pooled draws. A large
#' share between fits means that the result depends on which detections
#' are taken as true. The between part includes the Monte Carlo error of
#' each fit's mean, which is small next to it unless the fits are short.
#'
#' **Failed fits.** A version whose fit stops with an error is left out
#' with a warning; the versions should weigh equally, so the function also
#' warns when the fits have different numbers of draws.
#'
#' **Convergence.** R-hat over the pooled draws would mix differences
#' between fits with lack of convergence, so it is not reported for the
#' pool. Each fit is checked on its own: a fit whose largest R-hat exceeds
#' `max_rhat` (1.05 by default) has not converged and is left out, with a
#' warning; report how many were left out. This is looser than the 1.01
#' that [fp_occu()] warns about, because mixing error within one fit mostly
#' widens the within-fit part of the pooled variance; fits between the two
#' values are kept, and their warnings show in `diagnostics`. Divergent
#' transitions do not exclude a fit. Leaving fits out conditions the
#' propagation on the versions that converged: if non-convergence depends
#' on the version (for instance, histories with very few detections), the
#' pooled posterior is biased and the between-fit variance understated.
#' Rather than dropping them, refit those versions with longer chains or a
#' higher `adapt_delta`. The summary reports the largest R-hat and the
#' smallest effective sample size over the fits kept, and `diagnostics`
#' lists every fit, with `kept` marking those in the pool.
#'
#' **Number of versions.** The cost is that of one fit times the number of
#' versions. The between-fit variance is estimated from the versions, so
#' it needs several of them; check that the estimates of interest are
#' stable as versions are added.
#'
#' @return An object of class `fp_occu_mc`: `draws` (the pooled draws, in
#'   the readable names of [fp_occu()], the chains of fit `i` following
#'   those of fit `i - 1`), `summaries` (the posterior summary of each fit),
#'   `diagnostics` (one row per fit that ran: number of draws, divergences,
#'   largest R-hat, smallest bulk and tail effective sample sizes, and
#'   `kept`), `fits` (with `keep_fits = TRUE`, the fits kept), `n_fits` (the
#'   fits in the pool), and `call`.
#' @seealso [fp_occu()], [fp_read_sheepdog()]
#' @section Experimental:
#' This function is experimental (new in 2.0.0, part of the integration with
#' 'sheepdog'): its arguments, its class and methods, and its results may
#' change in a minor version.
#' See [crywolf-api].
#'
#' @export
#' @examples
#' \dontrun{
#' # three versions of a simulated dataset, with detections switched off at
#' # random as if they had been judged false
#' sim <- simulate_fp(nsite = 100, nocc = 8, seed = 1)
#' versions <- purrr::map(1:3, function(i) {
#'   d <- sim$data
#'   set.seed(i)
#'   hit <- d$auto$y == 1
#'   d$auto$y[hit] <- stats::rbinom(sum(hit), 1, 0.9)
#'   d
#' })
#' mc <- fp_occu_mc(versions, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
#'                  threads_per_chain = 1, seed = 1)
#' summary(mc)
#' }
fp_occu_mc <- function(data_list, ..., seed = NULL, parallel = FALSE,
                       max_rhat = 1.05, keep_fits = FALSE) {
  if (inherits(data_list, "fp_data") || !is.list(data_list) || !length(data_list) ||
      !all(purrr::map_lgl(data_list, inherits, "fp_data"))) {
    stop("`data_list` must be a list of `fp_data` objects.", call. = FALSE)
  }
  .check_fp_data(data_list[[1]])
  dots <- list(...)
  if (!is.numeric(max_rhat) || length(max_rhat) != 1 || is.na(max_rhat) || max_rhat < 1) {
    stop("`max_rhat` must be a single number of at least 1 (Inf keeps every fit).",
         call. = FALSE)
  }
  sites <- function(d) d$sites
  if (!all(purrr::map_lgl(data_list, function(d) identical(sites(d), sites(data_list[[1]]))))) {
    stop("The versions in `data_list` must have the same sites table (sites, ",
         "order and covariates); they should differ only in the detections.",
         call. = FALSE)
  }
  # one core per fit when the fits run in parallel, unless set by the user
  # (not modifyList(), which would drop an argument set to NULL, as fp = NULL)
  if (parallel) {
    one <- list(parallel_chains = 1, threads_per_chain = 1)
    dots <- c(dots, one[setdiff(names(one), names(dots))])
  }
  fp_compile()
  jobs <- tibble::tibble(i = seq_along(data_list), d = data_list)
  run_one <- function(i, d) {
    s <- if (is.null(seed)) NULL else as.integer(seed) + i - 1L
    m <- tryCatch(do.call(fp_occu, c(list(d, seed = s), dots)),
                  error = function(e) conditionMessage(e))
    if (is.character(m)) return(list(ok = FALSE, msg = m))
    list(ok = TRUE, draws = m$draws, summary = summary(m)$table,
         diagnostics = data.frame(fit = i, n_draws = posterior::ndraws(m$draws),
                                  divergences = m$diagnostics$divergences,
                                  max_rhat = m$diagnostics$max_rhat,
                                  min_ess_bulk = m$diagnostics$min_ess_bulk,
                                  min_ess_tail = m$diagnostics$min_ess_tail),
         fit = if (keep_fits) m)
  }
  res <- .run_jobs(jobs, run_one, parallel)
  ok <- purrr::map_lgl(res, "ok")
  if (!any(ok)) {
    stop("Every fit failed; the first error: ", res[[1]]$msg, call. = FALSE)
  }
  if (!all(ok)) {
    warning(sum(!ok), " fit(s) failed and were left out (version(s) ",
            paste(which(!ok), collapse = ", "), "); the first error: ",
            res[[which(!ok)[1]]]$msg, call. = FALSE)
  }
  res <- res[ok]
  diagnostics <- purrr::list_rbind(purrr::map(res, "diagnostics"))
  converged <- !((diagnostics$max_rhat > max_rhat) %in% TRUE)
  diagnostics$kept <- converged
  if (!any(converged)) {
    stop("No fit converged (largest R-hat above ", max_rhat, " in every fit); ",
         "run longer chains.", call. = FALSE)
  }
  if (!all(converged)) {
    warning(sum(!converged), " fit(s) did not converge (largest R-hat above ",
            max_rhat, ") and were left out: version(s) ",
            paste(diagnostics$fit[!converged], collapse = ", "), ". Refit them with ",
            "longer chains rather than dropping them: leaving versions out can bias ",
            "the pooled posterior.", call. = FALSE)
  }
  res <- res[converged]
  n_draws <- purrr::map_int(res, function(r) r$diagnostics$n_draws)
  if (length(unique(n_draws)) > 1) {
    warning("The fits have different numbers of draws (", paste(range(n_draws), collapse = " to "),
            "), for instance after a failed chain, so the versions do not weigh ",
            "equally in the pooled draws.", call. = FALSE)
  }
  draws <- purrr::reduce(purrr::map(res, "draws"), posterior::bind_draws, along = "chain")
  # summaries named by version, as `fit` in diagnostics
  summaries <- stats::setNames(purrr::map(res, "summary"),
                               diagnostics$fit[diagnostics$kept])
  structure(list(draws = draws, summaries = summaries,
                 diagnostics = diagnostics,
                 fits = if (keep_fits) purrr::map(res, "fit"),
                 n_fits = length(res), call = match.call()),
            class = "fp_occu_mc")
}

#' @rdname fp_occu_mc
#' @param object,x An `fp_occu_mc` object.
#' @param prob Width of the equal-tailed intervals.
#' @param digits Digits to print.
#' @return `summary()` returns, for each quantity, the pooled posterior mean,
#'   standard deviation and interval, the standard deviations within and
#'   between fits (`sd_within`, `sd_between`), the share of the variance
#'   between fits (`share_between`), and the largest R-hat and smallest bulk
#'   effective sample size over the fits.
#' @export
summary.fp_occu_mc <- function(object, prob = 0.95, ...) {
  a <- (1 - prob) / 2
  pooled <- posterior::summarise_draws(
    object$draws, mean = mean, sd = stats::sd,
    lower = function(v) unname(stats::quantile(v, a)),
    upper = function(v) unname(stats::quantile(v, 1 - a)))
  per_fit <- purrr::list_rbind(purrr::imap(object$summaries, function(t, i) {
    # i is the version number (the names of `summaries`), or the position in
    # objects made before they were named
    data.frame(fit = i, variable = t$variable, mean = t$mean, sd = t$sd,
               rhat = t$rhat, ess_bulk = t$ess_bulk)
  }))
  parts <- purrr::map(split(per_fit, factor(per_fit$variable, levels = pooled$variable)),
                      function(g) {
    # variance of the fit means over the fits (divisor M), so that the two
    # parts add up to the variance of the pooled draws
    between <- mean((g$mean - mean(g$mean))^2)
    within <- mean(g$sd^2)
    data.frame(sd_within = sqrt(within), sd_between = sqrt(between),
               share_between = if (between + within > 0) between / (between + within) else 0,
               max_rhat = max(g$rhat, na.rm = TRUE),
               min_ess_bulk = min(g$ess_bulk, na.rm = TRUE))
  }) |> purrr::list_rbind()
  structure(list(table = tibble::as_tibble(cbind(as.data.frame(pooled), parts)),
                 prob = prob, n_fits = object$n_fits, diagnostics = object$diagnostics,
                 call = object$call),
            class = "summary.fp_occu_mc")
}

#' @rdname fp_occu_mc
#' @export
print.summary.fp_occu_mc <- function(x, digits = 3, ...) {
  cat(sprintf("Posterior pooled over %d fits (%g%% equal-tailed intervals)\n\n",
              x$n_fits, 100 * x$prob))
  print(as.data.frame(x$table), digits = digits, row.names = FALSE)
  d <- x$diagnostics
  # fits made before `kept` existed
  if (is.null(d$kept)) d$kept <- TRUE
  k <- d[d$kept, , drop = FALSE]
  cat(sprintf(paste0("\n  Per fit: %d divergent transitions in all, largest R-hat %.3f, ",
                     "smallest bulk ESS %.0f\n"),
              sum(k$divergences), max(k$max_rhat), min(k$min_ess_bulk)))
  if (any(!d$kept)) {
    cat(sprintf("  %d fit(s) left out for lack of convergence: version(s) %s\n",
                sum(!d$kept), paste(d$fit[!d$kept], collapse = ", ")))
  }
  invisible(x)
}

#' @rdname fp_occu_mc
#' @export
print.fp_occu_mc <- function(x, digits = 3, ...) {
  print(summary(x), digits = digits)
  invisible(x)
}

#' @rdname fp_occu_mc
#' @export
as_draws_df.fp_occu_mc <- function(x, ...) {
  x$draws
}
