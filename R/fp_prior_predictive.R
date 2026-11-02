#' Prior predictive check
#'
#' Draws parameter values from the priors, simulates data with the design of
#' `data` (its sites, covariates, occasions, recorders, and verification), and
#' compares a few summaries of the simulated data with the observed ones. If
#' the observed values sit far in the tails of what the priors produce, the
#' priors (or the model) are at odds with the data before any fitting.
#'
#' @param data An `fp_data` object.
#' @param occ,det,fp,ver,markov,random_recorder,random_fp,random_det,priors As in
#'   [fp_occu()].
#' @param ndraws Number of prior draws, one simulated dataset each.
#' @param seed Optional integer seed.
#'
#' @details
#' Each prior draw takes the intercepts and slopes from their Normal priors,
#' the gap from its truncated Normal prior, and the standard deviations of
#' the random effects from their half-Normal priors, with new recorder and
#' site effects drawn from them. Covariates are taken as they are in `data`,
#' so the priors are judged on the scale the model will use. The summaries
#' are the share of sites with at least one automated detection (naive
#' occupancy), the share of automated occasions with a detection, the share
#' of verified occasions with a detection, the mean occupancy probability,
#' and the mean false-positive rate over sites.
#'
#' @return An object of class `fp_prior_predictive`: a data frame with one
#'   row per draw and one column per summary, with the observed values as
#'   the attribute `observed`. [print()] shows the prior predictive
#'   intervals and where the observed values fall; [plot()] draws them.
#' @seealso [fp_priors()], [simulate.fp_occu()]
#' @export
#' @examples
#' sim <- simulate_fp(nsite = 100, nocc = 10, seed = 1)
#' pp <- fp_prior_predictive(sim$data, occ = ~ x, fp = ~ 1,
#'                           random_recorder = FALSE, ndraws = 50, seed = 1)
#' pp
fp_prior_predictive <- function(data, occ = ~ 1, det = ~ 1, fp = ~ 1,
                                ver = ~ 1, markov = NULL,
                                random_recorder = !is.null(data$recorder),
                                random_fp = FALSE, random_det = FALSE,
                                priors = fp_priors(),
                                ndraws = 200, seed = NULL) {
  .check_fp_data(data)
  sl <- .build_stan_data(data, occ, det, fp, random_recorder, priors,
                         random_fp = random_fp, ver = ver, markov = markov,
                         random_det = random_det)
  gen <- list(data = data, formulas = list(occ = occ, det = det, fp = fp, ver = ver),
              stan_data = sl$stan_data, design = sl$design,
              occasions = sl$occasions)
  if (!is.null(seed)) .set_seed(seed)
  out <- purrr::map(seq_len(ndraws), function(i) {
    p <- .prior_draw(sl$stan_data, priors)
    sim <- .simulate_from(gen, p)
    data.frame(t(.data_summaries(sim$data)),
               mean_psi = mean(stats::plogis(drop(sl$stan_data$X %*% p$beta))),
               p10_bar = sim$truth$p10_bar)
  }) |> purrr::list_rbind()
  structure(out, class = c("fp_prior_predictive", "data.frame"),
            observed = .data_summaries(data))
}

# one draw of every parameter from the priors
.prior_draw <- function(sd, priors) {
  slopes <- function(k) stats::rnorm(k, 0, priors$slope_sd)
  int <- function() stats::rnorm(1, 0, priors$intercept_sd)
  # Normal(gap_mean, gap_sd) truncated at zero, by inversion
  lo <- stats::pnorm(0, priors$gap_mean, priors$gap_sd)
  gap <- stats::qnorm(stats::runif(1, lo, 1), priors$gap_mean, priors$gap_sd)
  alpha <- c(int(), slopes(sd$Kalpha - 1))
  # drawn in the same order as in crywolf 1.4.0, so that seeds reproduce
  beta <- c(int(), slopes(sd$Kbeta - 1))
  delta0 <- int()
  delta_s <- slopes(if (is.null(sd$Kver_s)) 0 else sd$Kver_s)
  sigma_u <- if (sd$estimate_re == 1) abs(stats::rnorm(1, 0, priors$recorder_sd))
  list(
    beta = beta,
    alpha = alpha,
    delta0 = delta0,
    delta_s = delta_s,
    sigma_u = sigma_u,
    u = if (sd$estimate_re == 1) stats::rnorm(sd$nrec, 0, sigma_u),
    gamma0 = alpha[1] - gap,
    gamma_fp = slopes(sd$Kfp),
    sigma_fp = if (isTRUE(sd$estimate_fp_re == 1)) {
      abs(stats::rnorm(1, 0, priors$fp_site_sd))
    },
    sigma_det = if (isTRUE(sd$estimate_det_re == 1)) {
      abs(stats::rnorm(1, 0, priors$det_site_sd))
    }
  )
}

# summaries compared between the simulated and the observed data
.data_summaries <- function(d) {
  a <- d$auto
  # every site in the denominator, also those without automated occasions
  detected <- d$sites[[d$site]] %in% a[[d$site]][a[[d$y]] == 1]
  v <- d$verified
  c(naive_occupancy = mean(detected),
    automated_rate = mean(a[[d$y]]),
    verified_rate = if (is.null(v) || nrow(v) == 0) NA_real_ else mean(v[[d$y]]))
}

#' @export
print.fp_prior_predictive <- function(x, digits = 3, ...) {
  obs <- attr(x, "observed")
  cat(sprintf("Prior predictive check (%d draws)\n\n", nrow(x)))
  rows <- purrr::map(names(x), function(v) {
    q <- stats::quantile(x[[v]], c(0.05, 0.5, 0.95), na.rm = TRUE)
    o <- if (v %in% names(obs)) obs[[v]] else NA_real_
    data.frame(summary = v, q05 = q[[1]], median = q[[2]], q95 = q[[3]],
               observed = o,
               share_below_observed = if (is.na(o)) NA_real_ else mean(x[[v]] < o))
  }) |> purrr::list_rbind()
  print(rows, digits = digits, row.names = FALSE)
  cat("\nObserved values near 0 or 1 in `share_below_observed` sit in the tails",
      "of the prior predictive distribution.\n")
  invisible(x)
}

#' @rdname fp_prior_predictive
#' @param x An `fp_prior_predictive` object.
#' @param ... Unused.
#' @export
plot.fp_prior_predictive <- function(x, ...) {
  .require_pkg("ggplot2")
  obs <- attr(x, "observed")
  long <- purrr::map(names(x), function(v) {
    data.frame(summary = v, value = x[[v]])
  }) |> purrr::list_rbind()
  long <- long[!is.na(long$value), ]
  ref <- data.frame(summary = names(obs), value = unname(obs))
  ref <- ref[!is.na(ref$value) & ref$summary %in% long$summary, ]
  ggplot2::ggplot(long, ggplot2::aes(value)) +
    ggplot2::geom_histogram(bins = 30, fill = "#2DE2E6", colour = "white") +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = value),
                        colour = "#FF3EA5", linewidth = 1) +
    ggplot2::facet_wrap(~ summary, scales = "free") +
    ggplot2::labs(x = "Value in data simulated from the priors", y = "Draws",
                  caption = "Line: observed value") +
    ggplot2::theme_minimal()
}

utils::globalVariables("value")
