# Plots of fitted models -------------------------------------------------------

#' Plot a fitted model
#'
#' Posterior intervals of the coefficients, or sampler diagnostics.
#'
#' @param x An `fp_occu` object.
#' @param type `"intervals"` (default) draws the posterior median with 50%
#'   and `prob` intervals of every coefficient on the logit scale, one panel
#'   per process, and the standard deviations of the random effects.
#'   `"trace"`, `"rank"`, and `"pairs"` draw sampler diagnostics with
#'   'bayesplot': trace plots and rank histograms by chain, and pairs plots
#'   with divergent transitions marked.
#' @param pars Output names of the parameters for the diagnostics. The
#'   default takes the occupancy slopes and the parameters that are hardest
#'   to sample in this model: the false-positive intercept, the gap, and the
#'   standard deviations of the random effects.
#' @param prob Width of the outer interval for `type = "intervals"`.
#' @param ... Passed to the 'bayesplot' function for the diagnostics.
#'
#' @return A 'ggplot2' object (for `"pairs"`, a 'bayesplot' grid).
#' @seealso [fp_effects()], [fp_sites()], [plot.fp_ppc()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
#' m <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
#' plot(m)
#' plot(m, type = "trace")
#' plot(m, type = "pairs", pars = c("fp_(Intercept)", "fp_gap"))
#' }
plot.fp_occu <- function(x, type = c("intervals", "trace", "rank", "pairs"),
                         pars = NULL, prob = 0.95, ...) {
  .require_pkg("ggplot2")
  type <- match.arg(type)
  if (type == "intervals") return(.plot_intervals(x, prob))
  .require_pkg("bayesplot")
  vars <- posterior::variables(x$draws)
  if (is.null(pars)) {
    pars <- intersect(c(grep("^occ_", vars, value = TRUE), "fp_(Intercept)",
                        "fp_gap", "sd_recorder", "sd_fp_site", "sd_det_site"), vars)
    pars <- setdiff(pars, "occ_(Intercept)")
    if (!length(pars)) {
      stop("The model has none of the default parameters; name them in `pars`.",
           call. = FALSE)
    }
  }
  bad <- setdiff(pars, vars)
  if (length(bad)) {
    stop("Unknown parameter(s): ", paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  np <- .nuts_params(x$fit)
  switch(type,
    trace = bayesplot::mcmc_trace(x$draws, pars = pars, np = np, ...),
    rank = bayesplot::mcmc_rank_overlay(x$draws, pars = pars, ...),
    pairs = bayesplot::mcmc_pairs(x$draws, pars = pars, np = np, ...)
  )
}

.plot_intervals <- function(x, prob) {
  d <- posterior::as_draws_matrix(x$draws)
  keep <- grepl("^(occ_|det_|fp_|verified_|sd_)", colnames(d)) &
    colnames(d) != "fp_gap"
  d <- d[, keep, drop = FALSE]
  a <- (1 - prob) / 2
  q <- function(p) .col_quantile(d, p)
  process <- sub("_.*", "", colnames(d))
  labels <- c(occ = "Occupancy", det = "True detection", fp = "False positives",
              verified = "Verified detection", sd = "Random effects (SD)")
  df <- data.frame(
    parameter = factor(colnames(d), rev(colnames(d))),
    process = factor(labels[process], unique(labels[process])),
    median = q(0.5), l50 = q(0.25), u50 = q(0.75), lo = q(a), hi = q(1 - a)
  )
  ggplot2::ggplot(df, ggplot2::aes(median, parameter)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi), colour = "#241042") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = l50, xmax = u50), linewidth = 1.6,
                            colour = "#FF3EA5") +
    ggplot2::geom_point(size = 2, colour = "#241042") +
    ggplot2::facet_wrap(~ process, scales = "free", ncol = 1) +
    ggplot2::labs(x = sprintf("Posterior median, 50%% and %g%% intervals (logit scale)",
                              100 * prob), y = NULL) +
    ggplot2::theme_minimal()
}

# sampler diagnostics in the long format of bayesplot::nuts_params()
.nuts_params <- function(fit) {
  sd <- posterior::as_draws_df(fit$sampler_diagnostics())
  cols <- intersect(c("accept_stat__", "stepsize__", "treedepth__",
                      "n_leapfrog__", "divergent__", "energy__"), names(sd))
  purrr::map(cols, function(cl) {
    data.frame(Chain = sd$.chain, Iteration = sd$.iteration,
               Parameter = cl, Value = sd[[cl]])
  }) |>
    purrr::list_rbind() |>
    transform(Parameter = factor(Parameter))
}

#' Marginal effect of a covariate
#'
#' Posterior of the occupancy probability (or the false-positive rate) along
#' one covariate, with the other covariates of that process held at their
#' means (numeric) or first level (factors) over the sites.
#'
#' @param object An `fp_occu` object.
#' @param covariate Name of a covariate of the process.
#' @param process `"occ"` (default) or `"fp"`.
#' @param values Values of the covariate, on the scale used in the model.
#'   The default spans its observed range in 100 steps, or takes every level
#'   of a factor.
#' @param transform Optional function that maps model-scale values to the
#'   scale shown on the x axis, for covariates standardized or transformed
#'   before fitting; for a covariate standardized with [base::scale()],
#'   `function(z) z * sd(raw) + mean(raw)`.
#' @param prob Width of the credible interval.
#'
#' @return A data frame of class `fp_effects` with the covariate value
#'   (`value`, model scale), the value shown on the plot (`x`), and the
#'   posterior mean and interval of the probability; [plot()] draws it.
#' @seealso [predict.fp_occu()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
#' m <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
#' plot(fp_effects(m, "x"))
#' plot(fp_effects(m, "x", process = "fp"))
#' }
fp_effects <- function(object, covariate, process = c("occ", "fp"),
                       values = NULL, transform = NULL, prob = 0.95) {
  process <- match.arg(process)
  f <- object$formulas[[process]]
  if (is.null(f)) stop("This model has no false-positive process.", call. = FALSE)
  vars <- all.vars(f)
  if (!covariate %in% vars) {
    stop("`", covariate, "` is not in the ", process, " formula.", call. = FALSE)
  }
  # covariates measured on each occasion come from the automated table
  sites <- object$data$sites
  if (!all(vars %in% names(sites))) {
    a <- object$data$auto
    sites <- cbind(a, sites[match(a[[object$data$site]], sites[[object$data$site]]),
                            setdiff(names(sites), names(a)), drop = FALSE])
  }
  if (is.null(values)) {
    v <- sites[[covariate]]
    values <- if (is.numeric(v)) {
      seq(min(v, na.rm = TRUE), max(v, na.rm = TRUE), length.out = 100)
    } else {
      factor(levels(factor(v)), levels = levels(factor(v)))
    }
  }
  typical <- purrr::map(sites[vars], function(v) {
    if (is.numeric(v)) mean(v, na.rm = TRUE) else factor(levels(factor(v))[1], levels(factor(v)))
  })
  nd <- as.data.frame(typical)[rep(1, length(values)), , drop = FALSE]
  nd[[covariate]] <- values
  p <- predict.fp_occu(object, newdata = nd, process = process, prob = prob)
  out <- data.frame(covariate = covariate, process = process, value = values,
                    x = if (is.null(transform)) values else transform(values),
                    estimate = p$estimate, lower = p$lower, upper = p$upper)
  structure(out, class = c("fp_effects", "data.frame"), prob = prob)
}

#' @rdname fp_effects
#' @param x An `fp_effects` object.
#' @param xlab Label of the x axis.
#' @param ... Unused.
#' @export
plot.fp_effects <- function(x, xlab = x$covariate[1], ...) {
  .require_pkg("ggplot2")
  ylab <- if (x$process[1] == "occ") "Occupancy probability" else "False-positive rate"
  if (!is.numeric(x$x)) {
    return(
      ggplot2::ggplot(x, ggplot2::aes(x, estimate)) +
        ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper),
                                 colour = "#241042") +
        ggplot2::labs(x = xlab, y = ylab) +
        ggplot2::theme_minimal()
    )
  }
  ggplot2::ggplot(x, ggplot2::aes(x, estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "#2DE2E6", alpha = 0.35) +
    ggplot2::geom_line(colour = "#241042", linewidth = 0.9) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_minimal()
}

#' Site-level estimates
#'
#' For every site: the occupancy probability given its covariates, the
#' probability that it is occupied given its covariates and its data
#' (detections and verification), and its false-positive rate. The second
#' is the quantity to map or rank sites by: a site with automated detections
#' but a high false-positive rate can still be unlikely to be occupied.
#'
#' @param object An `fp_occu` object, or an `fp_occu_dyn` object
#'   (see Details).
#' @param prob Width of the credible intervals.
#' @param ndraws Number of posterior draws to use (all by default).
#' @param seed Optional seed for the choice of draws.
#'
#' @details
#' For a dynamic model ([fp_occu_dyn()]) the result has one row per site and
#' season, with the site and season columns and the season covariates.
#' `psi` is then the probability of occupancy in that season given the
#' covariates alone (initial occupancy carried forward by colonization and
#' extinction), and `occupied` the probability that the site-season is
#' occupied given all the data of the site, in every season (computed by the
#' forward-backward algorithm for each posterior draw). `detections`,
#' `verified` and `p10` refer to the site-season.
#'
#' @return A data frame of class `fp_sites`, one row per site (per
#'   site-season for a dynamic model), with the
#'   site columns of the data, `psi` (occupancy probability, with `psi_lower`
#'   and `psi_upper`), `occupied` (posterior probability that the site is
#'   occupied given its data), `occupied_lower` and `occupied_upper`
#'   (quantiles over the posterior draws of that conditional probability;
#'   they show how much it depends on the uncertain parameters, not the
#'   uncertainty of the occupancy state itself, which `occupied` already
#'   gives), `p10` (posterior
#'   mean false-positive rate, when the model has false positives),
#'   `detections` (number of automated detections), and `verified` (whether
#'   a verified detection proves occupancy). [plot()] draws it.
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
#' m <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
#' st <- fp_sites(m)
#' head(st)
#' plot(st)
#' }
fp_sites <- function(object, prob = 0.95, ndraws = NULL, seed = NULL) {
  if (inherits(object, "fp_occu_dyn")) return(.fp_sites_dyn(object, prob, ndraws, seed))
  sd <- object$stan_data
  if (!is.null(seed)) .set_seed(seed)
  pars <- .ppc_pars(object$fit, sd, ndraws)
  r <- .ppc_core(sd, pars, "automated", "freeman-tukey")
  psi <- stats::plogis(pars$beta %*% t(sd$X))
  a <- (1 - prob) / 2
  q <- .col_quantile
  new <- c("psi", "psi_lower", "psi_upper", "occupied", "occupied_lower",
           "occupied_upper", "detections", "verified", "p10")
  clash <- intersect(new, names(object$data$sites))
  if (length(clash)) {
    stop("The sites table already has column(s) ", paste(clash, collapse = ", "),
         "; rename them before calling fp_sites().", call. = FALSE)
  }
  out <- data.frame(
    object$data$sites,
    psi = colMeans(psi), psi_lower = q(psi, a), psi_upper = q(psi, 1 - a),
    occupied = colMeans(r$pz), occupied_lower = q(r$pz, a),
    occupied_upper = q(r$pz, 1 - a),
    detections = tabulate(rep(sd$site2, sd$y2), nbins = sd$nsite),
    verified = sd$known_present == 1,
    row.names = NULL
  )
  if (sd$estimate_fp == 1) out$p10 <- colMeans(.p10_site_draws(sd, pars))
  structure(out, class = c("fp_sites", "data.frame"), prob = prob)
}

#' @rdname fp_sites
#' @param x An `fp_sites` object.
#' @param coords Optional names of two columns of `x` with the coordinates
#'   of the sites (x then y, e.g. `c("lon", "lat")`). With them the plot is
#'   a map of the probability that each site is occupied; without them,
#'   the sites are ranked by it.
#' @param ... Unused.
#' @export
plot.fp_sites <- function(x, coords = NULL, ...) {
  .require_pkg("ggplot2")
  x$evidence <- ifelse(x$verified, "verified detection",
                       ifelse(x$detections > 0, "automated detections only",
                              "no detection"))
  if (!is.null(coords)) {
    if (length(coords) != 2 || !all(coords %in% names(x))) {
      stop("`coords` must name two columns of `x`.", call. = FALSE)
    }
    return(
      ggplot2::ggplot(x, ggplot2::aes(.data[[coords[1]]], .data[[coords[2]]],
                                      colour = occupied, shape = evidence)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::scale_colour_viridis_c(limits = c(0, 1)) +
        ggplot2::labs(colour = "P(occupied)", shape = NULL) +
        ggplot2::coord_equal() +
        ggplot2::theme_minimal()
    )
  }
  x <- x[order(x$occupied), ]
  x$rank <- seq_len(nrow(x))
  ggplot2::ggplot(x, ggplot2::aes(rank, occupied, colour = evidence)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = occupied_lower, ymax = occupied_upper),
                            alpha = 0.3) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_colour_manual(values = c("verified detection" = "#241042",
                                            "automated detections only" = "#FF3EA5",
                                            "no detection" = "#2DE2E6")) +
    ggplot2::labs(x = "Sites, ranked", y = "Posterior probability of occupancy",
                  colour = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot a posterior predictive check
#'
#' @param x An `fp_ppc` object.
#' @param type `"discrepancy"` (default) plots the replicated against the
#'   observed discrepancy for every draw, with the line of equality; the
#'   Bayesian p-value is the share of points above it. `"histogram"` shows
#'   the same as the distribution of their difference, replicated minus
#'   observed, with the p-value as the share above zero. `"counts"` compares
#'   the distribution of detections per site with the replicates, with
#'   [bayesplot::ppc_bars()].
#' @param ... Passed to [bayesplot::ppc_bars()].
#'
#' @return A 'ggplot2' object.
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
#' m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
#' p <- fp_ppc(m, ndraws = 500, seed = 1)
#' plot(p)
#' plot(p, type = "counts")
#' }
plot.fp_ppc <- function(x, type = c("discrepancy", "histogram", "counts"), ...) {
  .require_pkg("ggplot2")
  type <- match.arg(type)
  if (type == "counts") {
    .require_pkg("bayesplot")
    return(bayesplot::ppc_bars(x$y, x$yrep, ...))
  }
  subtitle <- sprintf("%s, Bayesian p-value %.3f", x$fit_stat, x$p_value)
  if (type == "histogram") {
    return(
      ggplot2::ggplot(data.frame(d = x$fit_rep - x$fit_obs), ggplot2::aes(d)) +
        ggplot2::geom_histogram(bins = 40, fill = "#2DE2E6", colour = "white") +
        ggplot2::geom_vline(xintercept = 0, colour = "#FF3EA5", linewidth = 1) +
        ggplot2::labs(x = "Replicated minus observed discrepancy", y = "Draws",
                      subtitle = subtitle) +
        ggplot2::theme_minimal()
    )
  }
  df <- data.frame(observed = x$fit_obs, replicated = x$fit_rep)
  ggplot2::ggplot(df, ggplot2::aes(observed, replicated)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(alpha = 0.3, colour = "#241042") +
    ggplot2::labs(x = "Discrepancy of the observed data",
                  y = "Discrepancy of the replicated data",
                  subtitle = subtitle) +
    ggplot2::theme_minimal()
}

# quantile of every column of a draws matrix
.col_quantile <- function(m, p) {
  stats::setNames(purrr::map_dbl(seq_len(ncol(m)), function(j) stats::quantile(m[, j], p)),
                  colnames(m))
}

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("This function needs the '", pkg, "' package.", call. = FALSE)
  }
  invisible(TRUE)
}

utils::globalVariables(c("median", "parameter", "lo", "hi", "l50", "u50",
                         "estimate", "lower", "upper", "occupied",
                         "occupied_lower", "occupied_upper", "evidence",
                         "observed", "replicated", "Parameter", ".data", "d"))
