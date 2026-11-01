# Simulation and calibration from a fitted model ------------------------------

#' Simulate data from a fitted model
#'
#' Generates new detections from an `fp_occu` fit, keeping everything else
#' of the original data: the sites and their covariates, the automated
#' occasions with their observation covariates and recorders, which sites
#' and occasions were verified, and which occasions were checked under the
#' observation-confirmation design. Only the latent occupancy states and the
#' detections are drawn anew, from the fitted model.
#'
#' @param object An `fp_occu` object.
#' @param nsim Number of datasets.
#' @param seed Optional integer seed.
#' @param truth Parameter values to simulate from: `"mean"` (default) uses
#'   the posterior means for every dataset; `"posterior"` uses a different
#'   random posterior draw for each dataset, which carries the uncertainty
#'   of the fit into the simulations.
#' @param ... Unused.
#'
#' @details
#' The recorder effects, when the model has them, are kept at their values in
#' the chosen draw (or their posterior means), so the simulated data have the
#' same recorders with the same effects. Site effects on the false-positive
#' rate (`fp_occu(random_fp = TRUE)`) and on true detection (`random_det =
#' TRUE`) are kept at their values in the draw with `truth = "posterior"`;
#' with `truth = "mean"`, whose site effects would be shrunk towards zero,
#' new ones are drawn from a normal distribution with the posterior mean of
#' `sd_fp_site` (or `sd_det_site`) as standard deviation. Under the
#' multiple-methods design, an automated occasion at an occupied site is a detection with
#' probability \eqn{p_{11}} and at an unoccupied site with probability
#' \eqn{p_{10}}; verified occasions are detections with probability
#' \eqn{r_{11}} at occupied sites only. Under the observation-confirmation
#' design, true and false observations are drawn independently, as in
#' [simulate_fp()].
#'
#' @return A list of `nsim` elements, each a list with `data` (an `fp_data`
#'   object with the structure of the original data) and `truth` (the
#'   generating values on the scale of the output of [fp_occu()], the latent
#'   states `z`, the false-positive rate of each site `p10_site`, and its mean
#'   `p10_bar`).
#' @seealso [fp_calibrate_fit()], [simulate_fp()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
#' m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
#' new <- simulate(m, nsim = 2, seed = 1)
#' new[[1]]$data
#' new[[1]]$truth$p10_bar
#' }
simulate.fp_occu <- function(object, nsim = 1, seed = NULL,
                             truth = c("mean", "posterior"), ...) {
  truth <- match.arg(truth)
  gen <- .fit_generator(object)
  if (!is.null(seed)) .set_seed(seed)
  purrr::map(seq_len(nsim), function(i) {
    .simulate_from(gen, .pick_pars(gen$draws, truth))
  })
}

#' Calibration of a fitted model on data like the real ones
#'
#' Asks whether the credible intervals of a model would cover the truth in
#' data like the ones it was fitted to. The function simulates datasets from
#' the fit with [simulate.fp_occu()], so that they keep the design of the
#' real data (sites, covariates, occasions, verification, recorders) and
#' follow the fitted model, refits each candidate model to every dataset,
#' and scores the coverage and bias of their intervals, as [fp_calibrate()]
#' does for designed simulations.
#'
#' @param fit An `fp_occu` object: the generating model.
#' @param models Named list of false-positive formulas (or `NULL` for no
#'   false positives) to refit. The default refits the generating model.
#'   The occupancy and detection formulas, the recorder effect, and the
#'   priors are those of `fit`.
#' @param reps Number of simulated datasets.
#' @param truth `"mean"` (default) or `"posterior"`; see
#'   [simulate.fp_occu()].
#' @param prob Width of the credible intervals.
#' @param base_seed Seeds are `base_seed + 1, ..., base_seed + reps`.
#' @param parallel,fit_args As in [fp_calibrate()].
#'
#' @details
#' The scored quantities are the occupancy coefficients, the detection
#' coefficients, the intercept of the verified method, the mean
#' false-positive rate over sites (`p10_bar`, when both the generating and
#' the refitted model have false positives), and the slopes of the
#' false-positive rate in the refitted model. A slope that is not in the
#' generating model has a true value of zero.
#'
#' The most useful comparison pits the model you fitted against a simpler
#' one. For example, if a model with `fp = ~ temp + rh` estimated clear
#' effects of the covariates on the false-positive rate, refitting
#' `fp = ~ 1` to data simulated from it shows how much a constant rate would
#' bias the occupancy effects at your sample size and verification effort.
#' This answers the question for the fitted structure only: a dependence of
#' the false-positive rate on something you did not measure is not part of
#' the generating model.
#'
#' Each replicate refits every model, so the cost is that of `reps` times
#' the number of models fits of the size of the real data. Use `parallel`
#' for large datasets; each fit runs on one thread per chain unless
#' `fit_args` sets `threads_per_chain`.
#'
#' @return An object of class `fp_calibration`, as from [fp_calibrate()],
#'   with a single design cell.
#' @seealso [fp_calibrate()], [simulate.fp_occu()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
#'                    conf_cor = 0.7, seed = 1)
#' m_x <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
#' cal <- fp_calibrate_fit(m_x, models = list(constant = ~ 1, fp_on_x = ~ x),
#'                         reps = 20)
#' summary(cal)
#' plot(cal, parameter = "occ_x")
#' }
fp_calibrate_fit <- function(fit, models = list(fitted = fit$formulas$fp),
                             reps = 50, truth = c("mean", "posterior"),
                             prob = 0.95, base_seed = 1L, parallel = FALSE,
                             fit_args = list(chains = 4, parallel_chains = 1,
                                             iter_warmup = 1000,
                                             iter_sampling = 1000,
                                             adapt_delta = 0.95)) {
  if (!inherits(fit, "fp_occu")) {
    stop("`fit` must be an `fp_occu` object.", call. = FALSE)
  }
  truth <- match.arg(truth)
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    stop("`models` must be a named list of formulas.", call. = FALSE)
  }
  gen <- .fit_generator(fit)
  fit_args <- utils::modifyList(
    list(show_messages = FALSE, show_exceptions = FALSE, refresh = 0,
         threads_per_chain = 1, random_recorder = fit$random_recorder,
         random_fp = isTRUE(fit$random_fp), random_det = isTRUE(fit$random_det),
         priors = fit$priors,
         ver = if (is.null(fit$formulas$ver)) ~ 1 else fit$formulas$ver,
         markov = fit$markov),
    fit_args)
  jobs <- data.frame(rep = seq_len(reps), cell = 1L,
                     seed = as.integer(base_seed) + seq_len(reps))

  fp_compile()   # compile once, before any worker starts

  run_one <- function(rep, cell, seed) {
    .set_seed(seed)
    sim <- .simulate_from(gen, .pick_pars(gen$draws, truth))
    purrr::imap(models, function(fp, model) {
      .fit_and_score(sim$data, .fit_targets(sim$truth, fp, gen),
                     fp = fp, occ = gen$formulas$occ, det = gen$formulas$det,
                     model = model, prob = prob, seed = seed,
                     fit_args = fit_args) |>
        cbind(cell = cell, rep = rep, seed = seed)
    }) |>
      purrr::list_rbind()
  }
  out <- purrr::list_rbind(.run_jobs(jobs, run_one, parallel))
  out <- tibble::as_tibble(out)
  structure(out, class = c("fp_calibration", class(out)), prob = prob)
}

# internal ---------------------------------------------------------------------

# everything needed to simulate from a fit, without the CmdStan objects (so
# that it can be sent to parallel workers)
.fit_generator <- function(fit) {
  f <- fit$formulas
  sl <- .build_stan_data(fit$data, f$occ, f$det, f$fp, fit$random_recorder,
                         fit$priors, random_fp = isTRUE(fit$random_fp),
                         random_det = isTRUE(fit$random_det),
                         ver = if (is.null(f$ver)) ~ 1 else f$ver,
                         markov = fit$markov)
  list(data = fit$data, formulas = f, stan_data = sl$stan_data,
       design = sl$design, occasions = sl$occasions,
       draws = .ppc_pars(fit$fit, sl$stan_data, NULL))
}

# one set of parameter values: posterior means, or one random draw
.pick_pars <- function(draws, truth) {
  i <- if (truth == "mean") NULL else sample.int(length(draws$delta0), 1)
  pick <- function(m) {
    if (is.null(m)) return(NULL)
    if (is.null(dim(m))) return(if (is.null(i)) mean(m) else unname(m[i]))
    unname(if (is.null(i)) colMeans(m) else m[i, ])
  }
  p <- purrr::map(draws, pick)
  # posterior means of site effects are shrunk towards zero: with
  # truth = "mean", the simulation draws new site effects from their
  # distribution instead
  if (truth == "mean" && !is.null(p$sigma_fp)) p$e_fp <- NULL
  if (truth == "mean" && !is.null(p$sigma_det)) p$d_det <- NULL
  p
}

# with z, the latent states are those given (the dynamic model draws them
# over seasons), not drawn from the occupancy probability
.simulate_from <- function(gen, p, z = NULL) {
  sd <- gen$stan_data
  oc <- gen$occasions
  d <- gen$data
  nsite <- sd$nsite
  has_fp <- sd$estimate_fp == 1
  by_row <- isTRUE(sd$fp_by_row == 1)

  if (is.null(z)) z <- stats::rbinom(nsite, 1, stats::plogis(drop(sd$X %*% p$beta)))
  lp10 <- if (has_fp) {
    v <- rep(p$gamma0, nsite)
    if (sd$Kfp > 0 && !by_row) v <- v + drop(sd$Xfp %*% p$gamma_fp)
    if (!is.null(p$e_fp)) {
      v <- v + p$e_fp
    } else if (!is.null(p$sigma_fp)) {
      v <- v + stats::rnorm(nsite, 0, p$sigma_fp)
    }
    v
  } else rep(-Inf, nsite)

  # automated occasions, in the sorted order of the design
  lp11 <- drop(oc$Z %*% p$alpha)
  if (sd$estimate_re == 1) lp11 <- lp11 + p$u[oc$rec]
  if (!is.null(p$d_det)) {
    lp11 <- lp11 + p$d_det[oc$site]
  } else if (!is.null(p$sigma_det)) {
    lp11 <- lp11 + stats::rnorm(nsite, 0, p$sigma_det)[oc$site]
  }
  zi <- z[oc$site]
  lp10_occ <- lp10[oc$site]
  if (by_row) lp10_occ <- lp10_occ + drop(oc$Xfp %*% p$gamma_fp)
  s0 <- stats::plogis(lp10_occ)
  n <- length(oc$site)
  auto <- d$auto[oc$order, , drop = FALSE]
  confirm <- sd$obs_confirm == 1
  if (isTRUE(oc$markov)) {
    # runs of detections: simulate the occasions of each site in time order,
    # true detection depending on the simulated previous occasion
    lag_col <- which(colnames(oc$Z) == "y_prev")
    theta <- p$alpha[lag_col]
    base <- lp11 - theta * oc$Z[, lag_col]
    # position of each occasion within its site (rows are sorted by site)
    pos <- seq_len(n) - match(oc$site, oc$site) + 1L
    sims <- purrr::reduce(seq_len(max(pos)), function(acc, k) {
      idx <- which(pos == k)
      # as in the fit: no run carried over from before a gap in time
      prev <- if (k == 1) 0 else acc$y[idx - 1] * oc$has_prev[idx]
      s1 <- stats::plogis(base[idx] + theta * prev)
      if (confirm) {
        acc$true[idx] <- stats::rbinom(length(idx), 1, zi[idx] * s1)
        acc$false[idx] <- stats::rbinom(length(idx), 1, s0[idx])
        acc$y[idx] <- pmax(acc$true[idx], acc$false[idx])
      } else {
        acc$y[idx] <- stats::rbinom(length(idx), 1, ifelse(zi[idx] == 1, s1, s0[idx]))
      }
      acc
    }, .init = list(y = integer(n), true = integer(n), false = integer(n)))
    y <- sims$y
    true_obs <- sims$true
    false_obs <- sims$false
  } else if (confirm) {
    # true and false observations occur independently on every occasion
    true_obs <- stats::rbinom(n, 1, zi * stats::plogis(lp11))
    false_obs <- stats::rbinom(n, 1, s0)
    y <- pmax(true_obs, false_obs)
  } else {
    y <- stats::rbinom(n, 1, ifelse(zi == 1, stats::plogis(lp11), s0))
  }
  if (confirm) {
    conf <- d$confirmed
    auto[[conf[1]]] <- ifelse(oc$confirmed, true_obs, NA_integer_)
    auto[[conf[2]]] <- ifelse(oc$confirmed, false_obs, NA_integer_)
  }
  auto[[d$y]] <- as.integer(y)

  verified <- d$verified
  if (!is.null(verified) && nrow(verified) > 0) {
    vs <- match(verified[[d$site]], d$sites[[d$site]])
    # verified detection of each occasion, back in the order of the data
    eta <- rep(p$delta0, nrow(verified))
    if (isTRUE(sd$Kver_s > 0)) {
      eta[oc$verified_order] <- p$delta0 + drop(oc$X1s %*% p$delta_s)
    }
    verified[[d$y]] <- stats::rbinom(nrow(verified), 1, z[vs] * stats::plogis(eta))
  }

  # warnings about the design are for the user's data, not for every draw
  data <- withCallingHandlers(
    fp_data(auto = auto, sites = d$sites, verified = verified,
            site = d$site, y = d$y, recorder = d$recorder,
            confirmed = d$confirmed),
    crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
  truth <- c(
    stats::setNames(p$beta, paste0("occ_", gen$design$occ_names)),
    stats::setNames(p$alpha, paste0("det_", gen$design$det_names)),
    `verified_(Intercept)` = p$delta0,
    if (isTRUE(sd$Kver_s > 0))
      stats::setNames(p$delta_s, paste0("verified_", gen$design$ver_names[-1])),
    if (has_fp) c(`fp_(Intercept)` = p$gamma0),
    if (has_fp && sd$Kfp > 0)
      stats::setNames(p$gamma_fp, paste0("fp_", gen$design$fp_names[-1]))
  )
  # false-positive rate of each site: with occasion covariates, its mean over
  # the site's automated occasions (confirmed ones excluded), as in the fit
  auto_rows <- if (confirm) !oc$confirmed else rep(TRUE, n)
  p10_site <- if (has_fp) {
    if (by_row) {
      purrr::map_dbl(split(s0[auto_rows], factor(oc$site[auto_rows],
                                                 levels = seq_len(nsite))), mean)
    } else stats::plogis(lp10)
  }
  list(data = data,
       truth = list(pars = truth, z = z, p10_site = unname(p10_site),
                    p10_bar = if (!has_fp) 0 else if (by_row) mean(s0[auto_rows])
                              else mean(p10_site)))
}

# true values of the quantities scored for one refitted model
.fit_targets <- function(truth, fp, gen) {
  tp <- truth$pars
  keep <- grepl("^(occ_|det_)", names(tp))
  targets <- tp[keep]
  if (gen$stan_data$N1 > 0) {
    targets <- c(targets, tp[grepl("^verified_", names(tp))])
  }
  if (!is.null(fp)) {
    if (gen$stan_data$estimate_fp == 1) {
      targets <- c(targets, p10_bar = truth$p10_bar)
    }
    slopes <- colnames(.design(fp, gen$data$sites, "fp")$X)[-1]
    if (length(slopes)) {
      nm <- paste0("fp_", slopes)
      targets <- c(targets, stats::setNames(
        purrr::map_dbl(nm, function(v) if (v %in% names(tp)) tp[[v]] else 0), nm))
    }
  }
  targets
}
