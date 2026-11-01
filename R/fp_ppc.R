#' Posterior predictive check for a false-positive occupancy model
#'
#' Replicates the detection data from the posterior of an [fp_occu()] fit and
#' compares observed and replicated data with a discrepancy statistic, giving a
#' Bayesian p-value. Detections are summed per site, so the check asks whether
#' the model reproduces how many detections each site gets.
#'
#' @param object An `fp_occu` object, or an `fp_occu_dyn` object
#'   (experimental), for which detections are summed per site and season
#'   and the latent states are drawn jointly over seasons by forward
#'   filtering, backward sampling.
#' @param method Which data to check: `"automated"` (default; every site) or
#'   `"verified"` (verified sites only).
#' @param fit_stat Discrepancy statistic: `"freeman-tukey"` (default),
#'   \eqn{\sum_i (\sqrt{y_i} - \sqrt{E_i})^2}, or `"chi-squared"`,
#'   \eqn{\sum_i (y_i - E_i)^2 / E_i}, where \eqn{y_i} is the number of
#'   detections at site \eqn{i} and \eqn{E_i} its expected value given the
#'   parameters and the latent state.
#' @param ndraws Number of posterior draws to use, sampled without replacement.
#'   `NULL` uses all of them.
#' @param seed Optional integer seed for the subsampling and replication.
#'
#' @details
#' The latent occupancy state is marginalized in the Stan model, so for each
#' posterior draw \eqn{z_i} is first drawn from its conditional posterior given
#' the data and the parameters,
#' \deqn{\Pr(z_i = 1 \mid \mathbf{y}_i, \theta) = \frac{\psi_i\,
#' f(\mathbf{y}_i \mid z_i = 1, \theta)}{\psi_i\, f(\mathbf{y}_i \mid z_i = 1,
#' \theta) + (1 - \psi_i)\, f(\mathbf{y}_i \mid z_i = 0, \theta)},}
#' which is 1 at sites with a verified detection. Given \eqn{z_i}, automated
#' detections are replicated with probability \eqn{z_i p_{11} + (1 - z_i)
#' p_{10}} and verified detections with probability \eqn{z_i r_{11}}. The
#' discrepancy is computed for the observed and the replicated data with the
#' same \eqn{z} and parameters, and the Bayesian p-value is the proportion of
#' draws in which the replicated discrepancy exceeds the observed one. Values
#' near 0 or 1 indicate lack of fit; values near 0.5 indicate that the
#' observed data are typical of the fitted model.
#'
#' Under the observation-confirmation design, the conditional probability of
#' \eqn{z_i} also uses the confirmed occasions, and `"automated"` checks the
#' occasions that were not confirmed.
#'
#' A p-value near 0.5 does not show that the false-positive process is
#' correctly specified: a check of detection counts can have little power
#' against a false-positive rate that tracks an unmeasured gradient. Use
#' [fp_calibrate()] to ask whether the intervals of a given model are
#' calibrated for a design.
#'
#' The elements `y` and `yrep` follow the conventions of the 'bayesplot'
#' package, so its `ppc_*()` functions can be used directly, for example
#' `bayesplot::ppc_bars(ppc$y, ppc$yrep)`.
#'
#' @return An object of class `fp_ppc`, a list with `y` (observed detections
#'   per site), `yrep` (matrix of replicated detections, draws in rows and
#'   sites in columns), `fit_obs` and `fit_rep` (discrepancy of the observed and
#'   replicated data, one value per draw), `p_value`, `site` (the site
#'   identifiers of the columns), `method`, and `fit_stat`; for a dynamic
#'   model, also `season`, the season of each column.
#' @seealso [fp_occu()], [fp_calibrate()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
#' m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
#' ppc <- fp_ppc(m, ndraws = 500, seed = 1)
#' ppc
#' bayesplot::ppc_bars(ppc$y, ppc$yrep)
#' }
fp_ppc <- function(object, method = c("automated", "verified"),
                   fit_stat = c("freeman-tukey", "chi-squared"),
                   ndraws = 1000, seed = NULL) {
  dyn <- inherits(object, "fp_occu_dyn")
  if (!inherits(object, "fp_occu") && !dyn) {
    stop("`object` must be an `fp_occu` or `fp_occu_dyn` fit.", call. = FALSE)
  }
  method <- match.arg(method)
  fit_stat <- match.arg(fit_stat)
  sd <- object$stan_data
  if (method == "verified" && sd$N1 == 0) {
    stop("The data have no verified occasions.", call. = FALSE)
  }
  if (!is.null(seed)) .set_seed(seed)
  if (dyn) {
    # site-seasons, with the states drawn over seasons
    fb <- .dyn_ffbs(sd, .dyn_pars(object$fit, sd, ndraws))
    res <- .ppc_core(fb$em$sd_u, fb$em$pars_u, method, fit_stat, z = fb$z)
    units <- .dyn_units(object$data)$sites
    site_ids <- units[[object$data$site]]
  } else {
    pars <- .ppc_pars(object$fit, sd, ndraws)
    res <- .ppc_core(sd, pars, method, fit_stat)
    site_ids <- object$data$sites[[object$data$site]]
  }
  structure(
    list(y = res$y, yrep = res$yrep, fit_obs = res$fit_obs,
         fit_rep = res$fit_rep, p_value = mean(res$fit_rep > res$fit_obs),
         site = site_ids[res$keep],
         season = if (dyn) units[[object$data$season]][res$keep],
         method = method, fit_stat = fit_stat),
    class = "fp_ppc"
  )
}

#' @export
print.fp_ppc <- function(x, digits = 3, ...) {
  cat("Posterior predictive check (", x$method, " detections per site)\n\n",
      sep = "")
  cat(sprintf("  statistic:          %s\n", x$fit_stat))
  cat(sprintf("  draws, sites:       %d, %d\n", nrow(x$yrep), ncol(x$yrep)))
  cat(sprintf("  mean discrepancy:   observed %s, replicated %s\n",
              format(mean(x$fit_obs), digits = digits),
              format(mean(x$fit_rep), digits = digits)))
  cat(sprintf("  Bayesian p-value:   %s\n", format(x$p_value, digits = digits)))
  if (x$p_value < 0.1 || x$p_value > 0.9) {
    cat("\n  The p-value is far from 0.5: the model may not fit these data.\n")
  }
  invisible(x)
}

# internal ---------------------------------------------------------------------

# posterior draws of the parameters the replication needs, as matrices with
# one row per draw
.ppc_pars <- function(fit, sd, ndraws) {
  vars <- c("beta", "alpha", "delta0",
            if (isTRUE(sd$Kver_s > 0)) "delta_s",
            if (sd$estimate_re == 1) "u",
            if (sd$estimate_fp == 1) "gamma0",
            if (sd$estimate_fp == 1 && sd$Kfp > 0) "gamma_fp",
            if (isTRUE(sd$estimate_fp_re == 1)) c("e_raw", "sigma_fp"),
            if (isTRUE(sd$estimate_det_re == 1)) c("d_raw", "sigma_det"))
  dm <- posterior::as_draws_matrix(fit$draws(variables = vars))
  total <- nrow(dm)
  if (!is.null(ndraws) && ndraws < total) {
    dm <- dm[sort(sample.int(total, ndraws)), , drop = FALSE]
  }
  dm <- unclass(dm)
  cols <- function(v) {
    dm[, grep(paste0("^", v, "(\\[|$)"), colnames(dm)), drop = FALSE]
  }
  S <- nrow(dm)
  list(
    beta = cols("beta"), alpha = cols("alpha"), delta0 = cols("delta0")[, 1],
    delta_s = if (isTRUE(sd$Kver_s > 0)) cols("delta_s") else matrix(0, S, 0),
    u = if (sd$estimate_re == 1) cols("u") else NULL,
    gamma0 = if (sd$estimate_fp == 1) cols("gamma0")[, 1] else NULL,
    gamma_fp = if (sd$estimate_fp == 1 && sd$Kfp > 0) cols("gamma_fp") else
      matrix(0, S, 0),
    # site effects on the false-positive rate, S x nsite
    e_fp = if (isTRUE(sd$estimate_fp_re == 1))
      cols("e_raw") * cols("sigma_fp")[, 1] else NULL,
    sigma_fp = if (isTRUE(sd$estimate_fp_re == 1)) cols("sigma_fp")[, 1] else NULL,
    # site effects on true detection, S x nsite
    d_det = if (isTRUE(sd$estimate_det_re == 1))
      cols("d_raw") * cols("sigma_det")[, 1] else NULL,
    sigma_det = if (isTRUE(sd$estimate_det_re == 1)) cols("sigma_det")[, 1] else NULL
  )
}

# sum the columns of an S x N matrix by group, into an S x nsite matrix
.site_sum <- function(M, g, nsite) {
  out <- matrix(0, nrow(M), nsite)
  if (length(g) == 0) return(out)
  r <- rowsum(t(M), g)
  out[, as.integer(rownames(r))] <- t(r)
  out
}

# with ll_only = TRUE, returns only the log-likelihood, without drawing
# random numbers (for loo()); with emissions_only = TRUE, only the
# log-likelihood of each site's observations given that it is occupied (occ)
# and given that it is not (unocc), without the occupancy probability (for
# the units of the dynamic model); with z (S draws x nsite), the latent
# states are those given, not drawn from their conditional posterior (the
# states of the dynamic model, from .dyn_ffbs()), and pz and ll are NULL
.ppc_core <- function(sd, pars, method, fit_stat, ll_only = FALSE,
                      emissions_only = FALSE, z = NULL) {
  S <- nrow(pars$beta)
  nsite <- sd$nsite
  rows <- function(v) rep(v, each = S)   # recycle a per-column vector over draws

  # linear predictors
  use_psi <- !emissions_only && is.null(z)
  lpsi <- if (use_psi) pars$beta %*% t(sd$X)             # S x nsite
  lp11 <- pars$alpha %*% t(sd$Z)                         # S x N2
  if (sd$estimate_re == 1) lp11 <- lp11 + pars$u[, sd$rec2, drop = FALSE]
  if (!is.null(pars$d_det)) lp11 <- lp11 + pars$d_det[, sd$site2, drop = FALSE]
  fp <- sd$estimate_fp == 1
  by_row <- isTRUE(sd$fp_by_row == 1)
  if (fp) {
    lp10 <- .lp10_draws(sd, pars)                        # S x nsite
    lp10_row <- lp10[, sd$site2, drop = FALSE]           # S x N2
    if (by_row) lp10_row <- lp10_row + pars$gamma_fp %*% t(sd$Xfp2)
  }
  y2 <- sd$y2
  f2 <- sd$n2 - sd$y2
  conf <- isTRUE(sd$obs_confirm == 1)

  # conditional posterior of z, as in loglik_sites() of the Stan program
  if (conf) {
    # observation confirmation: log(1 - p11) = log(1 - s1) + log(1 - s0)
    lq11 <- stats::plogis(-lp11, log.p = TRUE) +
      stats::plogis(-lp10_row, log.p = TRUE)
    l11 <- rows(y2) * log(-expm1(lq11)) + rows(f2) * lq11
  } else {
    l11 <- rows(y2) * stats::plogis(lp11, log.p = TRUE) +
      rows(f2) * stats::plogis(-lp11, log.p = TRUE)
  }
  occ <- (if (use_psi) stats::plogis(lpsi, log.p = TRUE) else 0) +
    .site_sum(l11, sd$site2, nsite)
  if (sd$N1 > 0) {
    eta1 <- .eta_verified(sd, pars)                      # S x N1
    l1 <- rows(sd$y1) * stats::plogis(eta1, log.p = TRUE) +
      rows(1 - sd$y1) * stats::plogis(-eta1, log.p = TRUE)
    occ <- occ + .site_sum(l1, sd$site1, nsite)
  }
  unocc <- if (use_psi) stats::plogis(-lpsi, log.p = TRUE) else matrix(0, S, nsite)
  if (fp) {
    l10 <- rows(y2) * stats::plogis(lp10_row, log.p = TRUE) +
      rows(f2) * stats::plogis(-lp10_row, log.p = TRUE)
    unocc <- unocc + .site_sum(l10, sd$site2, nsite)
  } else {
    any_det <- tabulate(sd$site2[y2 > 0], nbins = nsite) > 0
    unocc[, any_det] <- -Inf
  }
  if (conf) {
    v3 <- sd$v3
    lp3 <- pars$alpha %*% t(sd$Z3)
    if (sd$estimate_re == 1) lp3 <- lp3 + pars$u[, sd$rec3, drop = FALSE]
    if (!is.null(pars$d_det)) lp3 <- lp3 + pars$d_det[, sd$site3, drop = FALSE]
    lp10_3 <- lp10[, sd$site3, drop = FALSE]
    if (by_row) lp10_3 <- lp10_3 + pars$gamma_fp %*% t(sd$Xfp3)
    occ3 <- rows(v3 >= 2) * stats::plogis(lp3, log.p = TRUE) +
      rows(v3 < 2) * stats::plogis(-lp3, log.p = TRUE) +
      rows(v3 %% 2 == 1) * stats::plogis(lp10_3, log.p = TRUE) +
      rows(v3 %% 2 == 0) * stats::plogis(-lp10_3, log.p = TRUE)
    unocc3 <- rows(v3 == 1) * stats::plogis(lp10_3, log.p = TRUE) +
      rows(v3 == 0) * stats::plogis(-lp10_3, log.p = TRUE)
    occ <- occ + .site_sum(occ3, sd$site3, nsite)
    unocc <- unocc + .site_sum(unocc3, sd$site3, nsite)
    unocc[, unique(sd$site3[v3 >= 2])] <- -Inf
  }
  if (emissions_only) return(list(occ = occ, unocc = unocc))
  ll <- pz <- NULL
  if (is.null(z)) {
    # marginal log-likelihood per site (the log_lik of the Stan program)
    m <- pmax(occ, unocc)
    ll <- m + log(exp(occ - m) + exp(unocc - m))
    ll[, sd$known_present == 1] <- occ[, sd$known_present == 1]
    if (ll_only) return(list(ll = ll))
    pz <- stats::plogis(occ - unocc)
    pz[, sd$known_present == 1] <- 1
    z <- matrix(stats::rbinom(S * nsite, 1, pz), S, nsite)
  }

  if (method == "automated") {
    zr <- z[, sd$site2, drop = FALSE]
    p11_row <- if (conf) -expm1(lq11) else stats::plogis(lp11)
    p_row <- zr * p11_row +
      if (fp) (1 - zr) * stats::plogis(lp10_row) else 0
    n_row <- rows(sd$n2)
    yrep_row <- matrix(stats::rbinom(length(p_row), n_row, p_row), S)
    keep <- sd$len2 > 0
    y <- tabulate(rep(sd$site2, y2), nbins = nsite)
    E <- .site_sum(n_row * p_row, sd$site2, nsite)
    yrep <- .site_sum(yrep_row, sd$site2, nsite)
  } else {
    keep <- tabulate(sd$site1, nbins = nsite) > 0
    y <- tabulate(sd$site1[sd$y1 == 1], nbins = nsite)
    p_row <- z[, sd$site1, drop = FALSE] * stats::plogis(eta1)   # S x N1
    E <- .site_sum(p_row, sd$site1, nsite)
    yrep <- .site_sum(matrix(stats::rbinom(length(p_row), 1, p_row), S),
                      sd$site1, nsite)
  }
  y <- y[keep]
  E <- E[, keep, drop = FALSE]
  yrep <- yrep[, keep, drop = FALSE]
  Y <- matrix(y, S, length(y), byrow = TRUE)
  disc <- switch(fit_stat,
    "freeman-tukey" = function(obs) rowSums((sqrt(obs) - sqrt(E))^2),
    "chi-squared"   = function(obs) rowSums((obs - E)^2 / (E + 1e-6))
  )
  list(y = y, yrep = yrep, fit_obs = disc(Y), fit_rep = disc(yrep),
       z = z, pz = pz, ll = ll, keep = which(keep))
}

# logit of the false-positive rate of each site, S draws x nsite
.lp10_draws <- function(sd, pars) {
  lp10 <- matrix(pars$gamma0, length(pars$gamma0), sd$nsite)
  if (sd$Kfp > 0 && !isTRUE(sd$fp_by_row == 1)) lp10 <- lp10 + pars$gamma_fp %*% t(sd$Xfp)
  if (!is.null(pars$e_fp)) lp10 <- lp10 + pars$e_fp
  lp10
}

# false-positive rate of each site, S x nsite: with occasion-level covariates,
# its mean over the site's automated occasions
.p10_site_draws <- function(sd, pars) {
  if (!isTRUE(sd$fp_by_row == 1)) return(stats::plogis(.lp10_draws(sd, pars)))
  lp <- .lp10_draws(sd, pars)[, sd$site2, drop = FALSE] + pars$gamma_fp %*% t(sd$Xfp2)
  w <- matrix(sd$n2, nrow(lp), ncol(lp), byrow = TRUE)
  .site_sum(w * stats::plogis(lp), sd$site2, sd$nsite) /
    matrix(pmax(1, tabulate(rep(sd$site2, sd$n2), nbins = sd$nsite)),
           nrow(lp), sd$nsite, byrow = TRUE)
}

# logit of the verified detection probability of each verified occasion, S x N1
.eta_verified <- function(sd, pars) {
  eta <- matrix(pars$delta0, length(pars$delta0), sd$N1)
  if (isTRUE(sd$Kver_s > 0)) eta <- eta + pars$delta_s %*% t(sd$X1s)
  eta
}
