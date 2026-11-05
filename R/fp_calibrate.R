#' Simulation-based calibration check for a study design
#'
#' Simulates many datasets from a design with [simulate_fp()], fits one or more
#' false-positive models to each with [fp_occu()], and scores how often the
#' credible intervals contain the true values and how biased the posterior
#' means are. Use it before collecting data, to ask whether the correction
#' will be calibrated for the number of sites, occasions, and verified sites
#' you can afford, and whether a constant false-positive rate is enough under
#' the confusion you suspect.
#'
#' @param design A data frame (or list) whose columns are arguments of
#'   [simulate_fp()]; each row is one cell of the design. Arguments not given
#'   take the defaults of [simulate_fp()].
#' @param models Named list of formulas for the false-positive rate, one model
#'   per element, as in the `fp` argument of [fp_occu()]. The default fits only
#'   the constant-rate model.
#' @param occ,det Formulas for occupancy and true detection, shared by all
#'   models. Simulated data have the occupancy covariate `x` and the detection
#'   covariates `bio` and `aci`.
#' @param reps Number of replicate datasets per cell. Each costs one Stan fit
#'   per model. Fifty replicates give a Monte Carlo standard error of about
#'   0.03 for a coverage of 0.95.
#' @param prob Width of the credible intervals scored.
#' @param base_seed Integer. Replicate `j` (counting over cells and replicates)
#'   uses seed `base_seed + j`, so results are reproducible and every model
#'   sees exactly the same datasets.
#' @param parallel Logical. Fit replicates in parallel with 'furrr'. Set the
#'   plan first, e.g. `future::plan(future::multisession, workers = 8)`, and
#'   keep `parallel_chains = 1` in `fit_args`. In interactive sessions with
#'   'progressr' installed, a progress bar shows the estimated time remaining,
#'   also when fitting in parallel.
#' @param fit_args List of further arguments for [fp_occu()]. Unless set
#'   here, `random_recorder` is `TRUE` only when the design simulates a
#'   recorder effect (`sd_recorder > 0`): fitting one whose true standard
#'   deviation is zero produces divergent transitions. With the chains of
#'   each fit run one after another, the cost is proportional to `chains`:
#'   `chains = 2` takes half the time of four and still gives an effective
#'   sample size of about 2000 per parameter in designs like the default
#'   one, ample for the 95% intervals that are scored. Each fit runs on one
#'   thread per chain (`threads_per_chain = 1`, unlike [fp_occu()], whose
#'   default is `"auto"`): `"auto"` divides the cores among the chains of one
#'   fit, not among parallel workers, and one thread keeps every replicate
#'   reproducible by its seed.
#'
#' @details
#' The scored quantities are the occupancy intercept and slope
#' (`occ_(Intercept)`, `occ_x`), the false-positive rate averaged over sites
#' (`p10_bar`, against the realized mean rate), the detection probability of
#' the verified method (`p_verified`), and, for models that include the
#' confusing gradient `w`, its effect on the false-positive rate (`fp_w`).
#'
#' @return A tibble of class `fp_calibration`, one row per replicate, model,
#'   and scored quantity, with the truth, the posterior mean and interval,
#'   whether the interval covered the truth, and sampler diagnostics. Use
#'   [base::summary()] for coverage, its Monte Carlo standard error, bias, and RMSE
#'   per cell and model.
#' @seealso [fp_verification_curve()], [simulate_fp()]; for dynamic designs,
#'   which this function does not simulate, [simulate_fp_dyn()]
#' @export
#' @examples
#' \dontrun{
#' design <- data.frame(p10 = 0.05, conf_effect = c(0, 1.5), conf_cor = 0.7)
#' cal <- fp_calibrate(design,
#'                     models = list(constant = ~ 1, fp_on_x = ~ x),
#'                     occ = ~ x, reps = 50)
#' summary(cal)
#' }
fp_calibrate <- function(design, models = list(constant = ~ 1),
                         occ = ~ x, det = ~ 1, reps = 50, prob = 0.95,
                         base_seed = 1L, parallel = FALSE,
                         fit_args = list(chains = 4, parallel_chains = 1,
                                         iter_warmup = 1000, iter_sampling = 1000,
                                         adapt_delta = 0.95)) {
  design <- as.data.frame(design)
  bad <- setdiff(names(design), setdiff(names(formals(simulate_fp)), "seed"))
  if (length(bad)) {
    stop("Unknown simulate_fp() argument(s) in `design`: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    stop("`models` must be a named list of formulas.", call. = FALSE)
  }
  design$cell <- seq_len(nrow(design))
  jobs <- expand.grid(rep = seq_len(reps), cell = design$cell)
  jobs$seed <- as.integer(base_seed) + seq_len(nrow(jobs))

  fp_compile()   # compile once, before any worker starts

  run_one <- function(cell, rep, seed) {
    args <- as.list(design[design$cell == cell, setdiff(names(design), "cell"),
                           drop = FALSE])
    sim <- do.call(simulate_fp, c(args, list(seed = seed)))
    purrr::imap(models, function(fp, model) {
      .score_fit(sim, fp = fp, occ = occ, det = det, model = model,
                 prob = prob, seed = seed, fit_args = fit_args) |>
        cbind(cell = cell, rep = rep, seed = seed)
    }) |>
      purrr::list_rbind()
  }

  res <- .run_jobs(jobs, run_one, parallel)
  out <- purrr::list_rbind(res)
  out <- merge(out, design, by = "cell", sort = FALSE)
  out <- tibble::as_tibble(out)
  structure(out, class = c("fp_calibration", class(out)), prob = prob)
}

#' Coverage as a function of verification effort
#'
#' Runs [fp_calibrate()] over a grid of verification fractions and prints the
#' coverage of the occupancy slope against the fraction of sites verified. If
#' the constant-rate model does not approach nominal coverage as verification
#' grows, the false positives are probably structured, and the remedy is to
#' model the false-positive rate (for instance on the occupancy covariate),
#' not to verify more.
#'
#' @param ... Arguments of [simulate_fp()] that define the scenario (held
#'   fixed), e.g. `p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7`.
#' @param frac_verified Fractions of sites verified to sweep.
#' @param models,occ,det,reps,base_seed,parallel,fit_args As in
#'   [fp_calibrate()].
#'
#' @return The `fp_calibration` object (invisibly), after printing the curve.
#' @export
#' @examples
#' \dontrun{
#' fp_verification_curve(p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7,
#'                       models = list(constant = ~ 1, fp_on_x = ~ x),
#'                       reps = 50)
#' }
fp_verification_curve <- function(..., frac_verified = c(0.1, 0.3, 0.5, 0.7, 0.9),
                                  models = list(constant = ~ 1),
                                  occ = ~ x, det = ~ 1, reps = 50,
                                  base_seed = 1L, parallel = FALSE,
                                  fit_args = list(chains = 4, parallel_chains = 1,
                                                  iter_warmup = 1000,
                                                  iter_sampling = 1000,
                                                  adapt_delta = 0.95)) {
  fixed <- list(...)
  design <- do.call(expand.grid, c(list(frac_verified = frac_verified),
                                   purrr::map(fixed, unique),
                                   list(stringsAsFactors = FALSE)))
  cal <- fp_calibrate(design, models = models, occ = occ, det = det,
                      reps = reps, base_seed = base_seed, parallel = parallel,
                      fit_args = fit_args)
  s <- summary(cal)
  curve <- s[s$parameter == "occ_x",
             c("model", "frac_verified", "coverage", "mc_se", "bias")]
  curve <- curve[order(curve$model, curve$frac_verified), ]
  cat("\nCoverage of the occupancy slope against verification effort\n")
  print(as.data.frame(curve), digits = 2, row.names = FALSE)
  const <- curve[curve$model == names(models)[1], ]
  if (nrow(const) > 1 && max(const$coverage) < 0.85) {
    cat("\n  The first model stays below nominal coverage at every verification",
        "\n  level. Verifying more will not fix it; model the false-positive rate.\n")
  }
  invisible(cal)
}

#' @export
summary.fp_calibration <- function(object, ...) {
  prob <- attr(object, "prob")
  keys <- c("cell", "model", "parameter")
  grp <- split(as.data.frame(object), as.data.frame(object)[keys], drop = TRUE)
  design_cols <- setdiff(names(object),
                         c(keys, "rep", "seed", "truth", "mean", "lower",
                           "upper", "covered", "divergences", "max_rhat",
                           "ok", "seconds"))
  out <- purrr::map(grp, function(g) {
    ok <- g[g$ok, , drop = FALSE]
    n <- nrow(ok)
    cov <- if (n) mean(ok$covered) else NA_real_
    err <- ok$mean - ok$truth
    data.frame(
      g[1, c(keys, design_cols), drop = FALSE],
      n = n, failed = sum(!g$ok),
      coverage = cov,
      mc_se = if (n) sqrt(cov * (1 - cov) / n) else NA_real_,
      bias = if (n) mean(err) else NA_real_,
      rmse = if (n) sqrt(mean(err^2)) else NA_real_,
      prop_divergent = if (n) mean(ok$divergences > 0) else NA_real_,
      row.names = NULL
    )
  }) |>
    purrr::list_rbind()
  out <- out[order(out$cell, out$model, out$parameter), ]
  structure(tibble::as_tibble(out), prob = prob,
            class = c("summary.fp_calibration", "tbl_df", "tbl", "data.frame"))
}

#' @export
print.fp_calibration <- function(x, ...) {
  cat(sprintf("<fp_calibration> %d fits: %d cell(s) x %d model(s) x %d replicate(s)\n",
              nrow(unique(as.data.frame(x)[c("cell", "model", "rep")])),
              length(unique(x$cell)), length(unique(x$model)),
              length(unique(x$rep))))
  cat("Use summary() for coverage and bias per cell and model.\n")
  invisible(x)
}

#' @export
print.summary.fp_calibration <- function(x, ...) {
  cat(sprintf("Calibration summary (%g%% intervals)\n\n", 100 * attr(x, "prob")))
  is_slope <- grepl("^occ_", x$parameter) & x$parameter != "occ_(Intercept)"
  slope <- as.data.frame(x)[is_slope, , drop = FALSE]
  if (nrow(slope)) {
    several <- length(unique(slope$parameter)) > 1
    cat(if (several) "Occupancy slopes:\n" else "Occupancy slope:\n")
    slope$verdict <- vapply(slope$coverage, .verdict, character(1))
    keep <- if (several) names(slope) else setdiff(names(slope), "parameter")
    print(slope[, keep], digits = 2, row.names = FALSE)
    cat("\n")
  }
  cat("All quantities: inspect the returned tibble.\n")
  invisible(x)
}

# internal ---------------------------------------------------------------------

.verdict <- function(cov, nominal = 0.95, tol = 0.05) {
  if (is.na(cov)) return("NA")
  if (cov >= nominal - tol) "ok" else if (cov >= 0.80) "low" else "collapsed"
}

.score_fit <- function(sim, fp, occ, det, model, prob, seed, fit_args) {
  tr <- sim$truth
  # fit a recorder effect only when the data were simulated with one: with
  # sd_recorder = 0 its variance collapses to zero and the sampler diverges
  fit_args <- utils::modifyList(list(show_messages = FALSE,
                                     show_exceptions = FALSE, refresh = 0,
                                     threads_per_chain = 1,
                                     random_recorder = tr$sd_recorder > 0),
                                fit_args)
  targets <- c(`occ_(Intercept)` = tr$b0, occ_x = tr$b1)
  if (tr$n_verified_sites > 0) targets <- c(targets, p_verified = tr$r11)
  if (!is.null(fp)) {
    targets <- c(targets, p10_bar = tr$p10_marginal)
    if ("w" %in% all.vars(fp)) targets <- c(targets, fp_w = tr$conf_effect)
  }
  if (!"x" %in% all.vars(occ)) targets <- targets[names(targets) != "occ_x"]
  .fit_and_score(sim$data, targets, fp = fp, occ = occ, det = det,
                 model = model, prob = prob, seed = seed, fit_args = fit_args)
}

# fits one model and scores the posterior intervals of `targets`, a named
# vector of true values of output parameters
.fit_and_score <- function(data, targets, fp, occ, det, model, prob, seed,
                           fit_args) {
  t0 <- Sys.time()
  fit <- tryCatch(
    suppressWarnings(do.call(fp_occu, c(list(data = data, occ = occ,
                                             det = det, fp = fp, seed = seed),
                                        fit_args))),
    error = function(e) NULL
  )
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(fit)) {
    return(data.frame(model = model, parameter = names(targets),
                      truth = unname(targets), mean = NA_real_,
                      lower = NA_real_, upper = NA_real_, covered = NA,
                      divergences = NA_integer_, max_rhat = NA_real_,
                      ok = FALSE, seconds = secs))
  }
  a <- (1 - prob) / 2
  purrr::imap(targets, function(truth, par) {
    v <- fit$draws[[par]]
    lo <- unname(stats::quantile(v, a)); hi <- unname(stats::quantile(v, 1 - a))
    data.frame(model = model, parameter = par, truth = truth, mean = mean(v),
               lower = lo, upper = hi, covered = truth >= lo & truth <= hi,
               divergences = fit$diagnostics$divergences,
               max_rhat = fit$diagnostics$max_rhat, ok = TRUE, seconds = secs)
  }) |>
    purrr::list_rbind()
}

#' Plot the results of a calibration check
#'
#' Plots the coverage (or bias) of one scored quantity against the column of
#' the design that varies, one colour per model. Coverage is shown with an
#' interval of two Monte Carlo standard errors and a dashed line at the nominal
#' level; bias, with a dashed line at zero.
#'
#' @param x An `fp_calibration` object from [fp_calibrate()] or
#'   [fp_verification_curve()].
#' @param parameter Scored quantity to plot, by default the occupancy slope
#'   (`"occ_x"`). See [fp_calibrate()] for the others.
#' @param type `"coverage"` (default) or `"bias"`.
#' @param by Name of the design column for the horizontal axis. By default, the
#'   first design column that varies among cells; further varying columns
#'   become panels. With no varying column, the models are on the axis.
#' @param verified How to show the verification effort when `by` is
#'   `"frac_verified"`: as the fraction of sites (default), as the expected
#'   number of verified sites (`"sites"`, `frac_verified * nsite`), or as the
#'   expected number of verified occasions, e.g. clips listened to
#'   (`"occasions"`, times `nver`), which is closer to the cost of the
#'   design.
#' @param target Optional coverage to mark with a second line, for instance
#'   the lowest coverage you are willing to accept (0.90).
#' @param ... Unused.
#'
#' @return A 'ggplot2' object, which can be modified with the usual `+`
#'   syntax.
#' @export
#' @examples
#' \dontrun{
#' cal <- fp_verification_curve(p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7,
#'                              models = list(constant = ~ 1, fp_on_x = ~ x),
#'                              reps = 50)
#' plot(cal)
#' plot(cal, type = "bias")
#' }
plot.fp_calibration <- function(x, parameter = "occ_x",
                                 type = c("coverage", "bias"), by = NULL,
                                 verified = c("fraction", "sites", "occasions"),
                                 target = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() needs the 'ggplot2' package.", call. = FALSE)
  }
  type <- match.arg(type)
  s <- as.data.frame(summary(x))
  if (!parameter %in% s$parameter) {
    stop("`parameter` must be one of: ",
         paste(unique(s$parameter), collapse = ", "), ".", call. = FALSE)
  }
  s <- s[s$parameter == parameter, , drop = FALSE]
  varying <- .varying_design(s)
  if (is.null(by)) {
    by <- if (length(varying)) varying[1] else "model"
  } else if (!by %in% names(s)) {
    stop("`by` must be a column of the design.", call. = FALSE)
  }
  panels <- setdiff(varying, by)
  verified <- match.arg(verified)
  xlab <- by
  if (verified != "fraction" && by != "frac_verified") {
    stop("`verified` applies only when the axis is `frac_verified`.", call. = FALSE)
  }
  if (by == "frac_verified" && verified != "fraction") {
    # design values, or the defaults of simulate_fp() when not in the design
    design_value <- function(nm) {
      if (nm %in% names(s)) s[[nm]] else eval(formals(simulate_fp)[[nm]])
    }
    s$frac_verified <- s$frac_verified * design_value("nsite") *
      if (verified == "occasions") design_value("nver") else 1
    xlab <- if (verified == "sites") "Expected number of verified sites" else
      "Expected number of verified occasions"
  }

  df <- data.frame(
    x = s[[by]], model = s$model,
    y = if (type == "coverage") s$coverage else s$bias,
    lo = if (type == "coverage") pmax(0, s$coverage - 2 * s$mc_se) else NA_real_,
    hi = if (type == "coverage") pmin(1, s$coverage + 2 * s$mc_se) else NA_real_,
    panel = if (length(panels)) {
      do.call(paste, c(purrr::map(panels, function(p) paste(p, "=", s[[p]])),
                       list(sep = ", ")))
    } else ""
  )
  dodge <- ggplot2::position_dodge(width = if (is.numeric(df$x))
    0.03 * diff(range(df$x, 1)) else 0.3)
  ref <- if (type == "coverage") attr(x, "prob") else 0
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = model, group = model)) +
    ggplot2::geom_hline(yintercept = ref, linetype = "dashed", colour = "grey40")
  if (is.numeric(df$x)) p <- p + ggplot2::geom_line(position = dodge, alpha = 0.6)
  if (!is.null(target) && type == "coverage") {
    p <- p + ggplot2::geom_hline(yintercept = target, linetype = "dotted",
                                 colour = "#FF3EA5")
  }
  p <- if (type == "coverage") {
    p + ggplot2::geom_pointrange(ggplot2::aes(ymin = lo, ymax = hi),
                                 position = dodge) +
      ggplot2::coord_cartesian(ylim = c(0, 1))
  } else {
    p + ggplot2::geom_point(position = dodge, size = 2)
  }
  if (length(panels)) p <- p + ggplot2::facet_wrap(~ panel)
  p + ggplot2::labs(
    x = xlab, colour = "Model",
    y = if (type == "coverage")
      sprintf("Coverage of the %g%% interval of %s", 100 * attr(x, "prob"), parameter)
    else sprintf("Bias of %s", parameter)
  ) +
    ggplot2::theme_minimal()
}

# design columns of a calibration summary that vary among cells
.varying_design <- function(s) {
  fixed <- c("cell", "model", "parameter", "n", "failed", "coverage", "mc_se",
             "bias", "rmse", "prop_divergent")
  cols <- setdiff(names(s), fixed)
  cols[purrr::map_lgl(cols, function(cl) length(unique(s[[cl]])) > 1)]
}

utils::globalVariables(c("y", "lo", "hi", "model"))

# runs one replicate per row of `jobs`, sequentially or with 'furrr', with a
# progress bar that shows the time remaining in interactive sessions
.run_jobs <- function(jobs, run_one, parallel) {
  if (parallel && !requireNamespace("furrr", quietly = TRUE)) {
    stop("`parallel = TRUE` needs the 'furrr' package.", call. = FALSE)
  }
  run_all <- function(tick) {
    job <- function(...) {
      out <- run_one(...)
      tick()
      out
    }
    if (parallel) {
      furrr::future_pmap(jobs, job, .options = furrr::furrr_options(seed = TRUE))
    } else {
      purrr::pmap(jobs, job)
    }
  }
  if (interactive() && requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p <- progressr::progressor(steps = nrow(jobs))
      run_all(function() p())
    }, handlers = progressr::handler_cli(
      format = "{cli::pb_current}/{cli::pb_total} replicates {cli::pb_bar} {cli::pb_percent} | ETA {cli::pb_eta}"))
  } else {
    run_all(function() NULL)
  }
}
