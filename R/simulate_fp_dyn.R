#' Simulate data for a dynamic false-positive occupancy model
#'
#' Generates one dataset from known parameters for [fp_occu_dyn()]: the
#' occupancy state of each site changes between seasons by colonization and
#' extinction, and within each season it is observed as in [simulate_fp()],
#' by an automated method with false positives and, at a random subset of
#' site-seasons, a verified method without them. **This function is
#' experimental** (see [crywolf-api]).
#'
#' @param nsite Number of sites.
#' @param nseason Number of seasons.
#' @param nocc Number of automated occasions per site and season.
#' @param psi1 Occupancy probability in the first season, at the mean of the
#'   site covariate `x`.
#' @param b_psi1 Slope of initial occupancy on `x` (logit scale).
#' @param col Colonization probability at the mean of the season covariate
#'   `rain`.
#' @param b_col Slope of colonization on `rain`.
#' @param ext Extinction probability at the mean of `x`.
#' @param b_ext Slope of extinction on `x`.
#' @param p11 True-detection probability of the automated method.
#' @param p10 False-positive probability of the automated method.
#' @param r11 Detection probability of the verified method.
#' @param frac_verified Probability that a site-season is verified.
#' @param nver Number of verified occasions per verified site-season.
#' @param frac_confirmed Probability that an automated occasion is checked
#'   under the observation-confirmation design, as in [simulate_fp()]. 0
#'   (default) for none.
#' @param seed Optional integer seed.
#'
#' @details
#' The site covariate `x` and the site-season covariate `rain` are
#' standardized normal draws; `rain` is returned in `season_covs` of the
#' data. A transition into season \eqn{t} uses the covariates of season
#' \eqn{t}, as in [fp_occu_dyn()]. With `frac_confirmed > 0`, true and false
#' observations occur independently on every occasion, and checked
#' occasions get the columns `n_true` and `n_false` of `auto`.
#'
#' @return A list of class `fp_sim_dyn` with elements `data` (an `fp_data`
#'   object with a `season` column and the season covariate `rain`) and
#'   `truth` (the generating values, their logits named as the output of
#'   [fp_occu_dyn()] in `pars`, and the latent states `z`, a matrix with
#'   sites in rows and seasons in columns).
#' @seealso [fp_occu_dyn()], [simulate_fp()]
#' @export
#' @examples
#' s <- simulate_fp_dyn(nsite = 50, nseason = 4, seed = 1)
#' s
#' colMeans(s$truth$z)   # realized occupancy per season
simulate_fp_dyn <- function(nsite = 100, nseason = 5, nocc = 4,
                            psi1 = 0.5, b_psi1 = 0, col = 0.2, b_col = 0,
                            ext = 0.2, b_ext = 0, p11 = 0.6, p10 = 0.05,
                            r11 = 0.8, frac_verified = 0.2, nver = 2,
                            frac_confirmed = 0, seed = NULL) {
  stopifnot(nsite >= 2, nseason >= 1, nocc >= 1, nver >= 1,
            psi1 > 0, psi1 < 1, col > 0, col < 1, ext > 0, ext < 1,
            p11 > 0, p11 < 1, p10 > 0, p10 < 1, r11 > 0, r11 <= 1,
            frac_verified >= 0, frac_verified <= 1,
            frac_confirmed >= 0, frac_confirmed <= 1)
  if (!is.null(seed)) .set_seed(seed)
  x <- as.numeric(scale(stats::rnorm(nsite)))
  rain <- matrix(stats::rnorm(nsite * nseason), nsite, nseason)

  # occupancy states, sites in rows and seasons in columns
  z <- matrix(0L, nsite, nseason)
  z[, 1] <- stats::rbinom(nsite, 1, stats::plogis(stats::qlogis(psi1) + b_psi1 * x))
  for (t in seq_len(nseason)[-1]) {
    p_col <- stats::plogis(stats::qlogis(col) + b_col * rain[, t])
    p_ext <- stats::plogis(stats::qlogis(ext) + b_ext * x)
    z[, t] <- stats::rbinom(nsite, 1, ifelse(z[, t - 1] == 1, 1 - p_ext, p_col))
  }

  # automated occasions
  auto <- expand.grid(site = seq_len(nsite), season = seq_len(nseason),
                      occasion = seq_len(nocc))
  zi <- z[cbind(auto$site, auto$season)]
  if (frac_confirmed > 0) {
    true_obs <- stats::rbinom(nrow(auto), 1, zi * p11)
    false_obs <- stats::rbinom(nrow(auto), 1, p10)
    auto$y <- pmax(true_obs, false_obs)
    checked <- stats::runif(nrow(auto)) < frac_confirmed
    auto$n_true <- ifelse(checked, true_obs, NA_integer_)
    auto$n_false <- ifelse(checked, false_obs, NA_integer_)
  } else {
    auto$y <- stats::rbinom(nrow(auto), 1, ifelse(zi == 1, p11, p10))
  }

  # verified occasions at a random subset of site-seasons
  units <- expand.grid(site = seq_len(nsite), season = seq_len(nseason))
  units <- units[stats::runif(nrow(units)) < frac_verified, , drop = FALSE]
  verified <- units[rep(seq_len(nrow(units)), each = nver), , drop = FALSE]
  verified$occasion <- rep(seq_len(nver), nrow(units))
  verified$y <- stats::rbinom(nrow(verified), 1,
                              r11 * z[cbind(verified$site, verified$season)])
  rownames(verified) <- NULL

  season_covs <- data.frame(site = rep(seq_len(nsite), nseason),
                            season = rep(seq_len(nseason), each = nsite),
                            rain = as.vector(rain))
  data <- fp_data(auto = auto, sites = data.frame(site = seq_len(nsite), x = x),
                  verified = verified, season = "season", season_covs = season_covs,
                  confirmed = if (frac_confirmed > 0) c("n_true", "n_false"))
  pars <- c(`psi1_(Intercept)` = stats::qlogis(psi1), psi1_x = b_psi1,
            `col_(Intercept)` = stats::qlogis(col), col_rain = b_col,
            `ext_(Intercept)` = stats::qlogis(ext), ext_x = b_ext,
            `det_(Intercept)` = stats::qlogis(p11),
            `fp_(Intercept)` = stats::qlogis(p10),
            `verified_(Intercept)` = stats::qlogis(r11))
  truth <- list(psi1 = psi1, col = col, ext = ext, p11 = p11, p10 = p10,
                r11 = r11, b_psi1 = b_psi1, b_col = b_col, b_ext = b_ext,
                frac_verified = frac_verified, frac_confirmed = frac_confirmed,
                pars = pars, z = z)
  structure(list(data = data, truth = truth), class = "fp_sim_dyn")
}

#' @export
print.fp_sim_dyn <- function(x, ...) {
  t <- x$truth
  cat(sprintf("<fp_sim_dyn> %d sites over %d seasons\n", nrow(t$z), ncol(t$z)))
  cat(sprintf(paste0("  truth: psi1 = %.2f, col = %.2f, ext = %.2f, p11 = %.2f, ",
                     "p10 = %.3f, r11 = %.2f\n"),
              t$psi1, t$col, t$ext, t$p11, t$p10, t$r11))
  cat("  realized occupancy per season:",
      paste(sprintf("%.2f", colMeans(t$z)), collapse = " "), "\n")
  print(x$data)
  invisible(x)
}
