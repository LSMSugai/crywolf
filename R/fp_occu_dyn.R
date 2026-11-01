#' Dynamic (multi-season) false-positive occupancy model
#'
#' Fits a false-positive occupancy model over several seasons (primary
#' periods), in which the occupancy state of each site changes between
#' seasons by colonization and extinction. Within each season the
#' observation model is that of [fp_occu()]: the multiple-methods and
#' observation-confirmation designs, verified detections, recorder and site
#' effects, false-positive covariates and runs of detections. **This
#' function is experimental** (see [crywolf-api]): its arguments, results
#' and output names may change in a minor version.
#'
#' The occupancy state of site \eqn{i} in the first season is Bernoulli with
#' probability \eqn{\psi_{1,i}}. In each later season \eqn{t}, an occupied
#' site stays occupied with probability \eqn{1 - \epsilon_{it}} (extinction
#' \eqn{\epsilon}) and an unoccupied one becomes occupied with probability
#' \eqn{\gamma_{it}} (colonization), on the logit scale linear in the
#' covariates of `col` and `ext`. The covariates of a transition are those
#' of the season it enters (row \eqn{t} of the site-season table for the
#' transition from \eqn{t - 1} to \eqn{t}). The states are summed out by
#' the forward algorithm, so the log-likelihood (`log_lik`) is per site,
#' over all its seasons.
#'
#' A site without data in a season contributes to that season only through
#' the transitions. With seasons that are dates, [fp_data()] fills in the
#' missing dates, so colonization and extinction are probabilities per step
#' (per night, for nightly seasons). A single false positive looks like a
#' colonization followed by an extinction: without verified or confirmed
#' data in the seasons, the two processes are hard to separate from false
#' positives.
#'
#' @param data An `fp_data` object with a season column; see the `season`
#'   and `season_covs` arguments of [fp_data()].
#' @param psi1 One-sided formula for the occupancy probability in the first
#'   season. It may use site covariates and covariates of the first season.
#' @param col,ext One-sided formulas for the colonization and extinction
#'   probabilities. They may use site covariates and season covariates
#'   (`fp_data(season_covs = )`). All three formulas keep their intercept,
#'   with the priors of [fp_priors()] for intercepts and slopes.
#' @param det,fp,random_recorder,ver,markov,random_fp,random_det,priors As in
#'   [fp_occu()]; formulas may also use season covariates. `markov` runs of
#'   detections do not cross seasons. Site effects (`random_fp`,
#'   `random_det`) belong to sites and are shared by their seasons.
#' @param chains,parallel_chains,threads_per_chain,iter_warmup,iter_sampling,adapt_delta,seed,refresh,...
#'   Passed to the sampler, as in [fp_occu()].
#'
#' @return An object of class `fp_occu_dyn`, a list with `fit`, `draws`
#'   (readable names: see [crywolf-output]), `stan_data`, `formulas`
#'   (`psi1`, `col`, `ext`, `det`, `fp`, `ver`), `markov`, `design`, `data`,
#'   `random_recorder`, `random_fp`, `random_det`, `priors`, `sampler`,
#'   `diagnostics` and `call`. Its methods are `print()`, `summary()` and
#'   `as_draws_df()`.
#' @seealso [fp_data()], [fp_occu()], [crywolf-output]
#' @export
#' @examples
#' \dontrun{
#' # three seasons of two occasions each, from a single-season simulation
#' sim <- simulate_fp(nsite = 150, nocc = 6, frac_verified = 0.3, seed = 1)
#' d <- sim$data
#' d$auto$season <- ceiling(d$auto$occasion / 2)
#' d$verified$season <- ceiling(d$verified$occasion * 3 / max(d$verified$occasion))
#' dd <- fp_data(d$auto, d$sites, verified = d$verified, season = "season")
#' m <- fp_occu_dyn(dd, psi1 = ~ x, col = ~ 1, ext = ~ 1, fp = ~ 1,
#'                  random_recorder = FALSE)
#' m
#' }
fp_occu_dyn <- function(data, psi1 = ~ 1, col = ~ 1, ext = ~ 1, det = ~ 1, fp = ~ 1,
                        random_recorder = !is.null(data$recorder),
                        ver = ~ 1, markov = NULL, random_fp = FALSE,
                        random_det = FALSE, priors = fp_priors(),
                        chains = 4, parallel_chains = chains, threads_per_chain = "auto",
                        iter_warmup = 1000, iter_sampling = 1000,
                        adapt_delta = 0.95, seed = NULL, refresh = NULL, ...) {
  if (!inherits(data, "fp_data") || is.null(data$season)) {
    stop("`data` must be an `fp_data` object with a season column; see ",
         "`fp_data(season = )`.", call. = FALSE)
  }
  no_ver <- is.null(data$verified) || nrow(data$verified) == 0
  no_conf <- is.null(data$confirmed) || all(is.na(data$auto[[data$confirmed[1]]]))
  if (no_ver && no_conf) {
    warning("No verified data: the false-positive rate is identified only by ",
            "the priors. Treat the estimates with caution.", call. = FALSE)
  }
  cores <- .sampler_cores(threads_per_chain, parallel_chains)
  threads_per_chain <- cores$threads_per_chain
  parallel_chains <- cores$parallel_chains
  sd_list <- .build_dyn_stan_data(data, psi1, col, ext, det, fp, random_recorder, priors,
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
  mod <- .compile_model("fp_occu_dyn")
  run_sampler <- function(..., output_basename = .output_basename()) {
    mod$sample(..., output_basename = output_basename)
  }
  fit <- run_sampler(data = sd_list$stan_data, chains = chains,
                     parallel_chains = parallel_chains,
                     threads_per_chain = threads_per_chain,
                     iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                     adapt_delta = adapt_delta, seed = seed, refresh = refresh, ...)
  fit$draws()
  fit$sampler_diagnostics()
  draws <- .pretty_draws_dyn(fit, sd_list)
  diag <- .diagnostics(fit, draws)
  out <- structure(
    list(fit = fit, draws = draws, stan_data = sd_list$stan_data,
         formulas = list(psi1 = psi1, col = col, ext = ext, det = det, fp = fp, ver = ver),
         markov = markov, design = sd_list$design, data = data,
         random_recorder = isTRUE(sd_list$stan_data$estimate_re == 1),
         random_fp = isTRUE(sd_list$stan_data$estimate_fp_re == 1),
         random_det = isTRUE(sd_list$stan_data$estimate_det_re == 1),
         priors = priors,
         sampler = list(chains = chains, parallel_chains = parallel_chains,
                        threads_per_chain = threads_per_chain,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = adapt_delta),
         diagnostics = diag, call = match.call()),
    class = "fp_occu_dyn"
  )
  .warn_diagnostics(diag)
  out
}

#' @rdname fp_occu_dyn
#' @param x,object An `fp_occu_dyn` object.
#' @param digits Number of digits to print.
#' @export
print.fp_occu_dyn <- function(x, digits = 3, ...) {
  sd <- x$stan_data
  f <- x$formulas
  cat("Dynamic Bayesian false-positive occupancy model (",
      if (isTRUE(sd$obs_confirm == 1)) "observation-confirmation" else "multiple-methods",
      " design)\n\n", sep = "")
  cat("  initial occupancy:", .fmt_formula(f$psi1), "\n")
  cat("  colonization:     ", .fmt_formula(f$col), "\n")
  cat("  extinction:       ", .fmt_formula(f$ext), "\n")
  cat("  true detection:   ", .fmt_formula(f$det),
      if (x$random_recorder) "+ (1 | recorder)" else "",
      if (isTRUE(x$random_det)) "+ (1 | site)" else "", "\n")
  cat("  false positive:   ",
      if (is.null(f$fp)) "none (standard occupancy model)" else .fmt_formula(f$fp),
      if (isTRUE(x$random_fp)) "+ (1 | site)" else "", "\n")
  if (!is.null(f$ver) && length(all.vars(f$ver))) {
    cat("  verified:         ", .fmt_formula(f$ver), "\n")
  }
  if (!is.null(x$markov)) {
    cat("  runs:              previous detection (y_prev) within seasons, ordered by",
        x$markov, "\n")
  }
  cat(sprintf("  data:              %d sites over %d seasons, %d automated and %d verified occasions\n\n",
              sd$nsite, sd$nseason, sum(sd$n2) + sd$N3, sd$N1))
  print(as.data.frame(summary(x)$table), digits = digits, row.names = FALSE)
  .print_diag(x$diagnostics)
  invisible(x)
}

#' @rdname fp_occu_dyn
#' @param prob Width of the equal-tailed credible interval.
#' @return `summary()` returns an object of class `summary.fp_occu`, as for
#'   [fp_occu()].
#' @export
summary.fp_occu_dyn <- function(object, prob = 0.95, ...) {
  summary.fp_occu(object, prob = prob)
}

#' @rdname fp_occu_dyn
#' @export
as_draws_df.fp_occu_dyn <- function(x, ...) {
  x$draws
}

# The data of fp_occu_dyn.stan. Each site-season is a unit, numbered
# (site - 1) * nseason + season; the units form an fp_data whose "sites" are
# the units, so that .build_stan_data() builds the observation model of each
# unit as fp_occu() builds that of each site. The designs of initial
# occupancy, colonization and extinction come from the table of units, with
# the season covariates: the first season for psi1, the later ones for col
# and ext (rows of the first season are zero and unused).
.build_dyn_stan_data <- function(data, psi1, col, ext, det, fp, random_recorder, priors,
                                 random_fp = FALSE, ver = ~ 1, markov = NULL,
                                 random_det = FALSE) {
  units <- .dyn_units(data)
  sl <- .build_stan_data(units, ~ 1, det, fp, random_recorder, priors,
                         random_fp = random_fp, ver = ver, markov = markov,
                         random_det = random_det)
  nsite <- nrow(data$sites)
  nseason <- length(data$seasons)
  t_unit <- rep(seq_len(nseason), nsite)
  d_psi1 <- .design(psi1, units$sites[t_unit == 1, , drop = FALSE], "psi1")
  trans <- function(form, what) {
    if (nseason == 1) {
      if (length(all.vars(form))) {
        stop("`", what, "` has covariates but the data have one season.", call. = FALSE)
      }
      X <- matrix(1, nsite, 1, dimnames = list(NULL, "(Intercept)"))
      return(list(d = list(X = X, intercept = TRUE), X = X))
    }
    d <- .design(form, units$sites[t_unit > 1, , drop = FALSE], what)
    X <- matrix(0, length(t_unit), ncol(d$X), dimnames = list(NULL, colnames(d$X)))
    X[t_unit > 1, ] <- d$X
    list(d = d, X = X)
  }
  d_col <- trans(col, "col")
  d_ext <- trans(ext, "ext")
  for (d in list(psi1 = d_psi1, col = d_col$d, ext = d_ext$d)) {
    if (!d$intercept) {
      stop("`psi1`, `col` and `ext` must keep the intercept.", call. = FALSE)
    }
  }
  sd <- sl$stan_data
  sd$nunit <- sd$nsite
  sd$nsite <- nsite
  sd$nseason <- nseason
  sd$Kbeta <- ncol(d_psi1$X)
  sd$X <- unname(d_psi1$X)
  sd$Kcol <- ncol(d_col$X)
  sd$Xcol <- unname(d_col$X)
  sd$Kext <- ncol(d_ext$X)
  sd$Xext <- unname(d_ext$X)
  sl$stan_data <- sd
  sl$design$occ <- d_psi1
  sl$design$occ_names <- colnames(d_psi1$X)
  sl$design$col <- d_col$d
  sl$design$ext <- d_ext$d
  sl$design$col_names <- colnames(d_col$X)
  sl$design$ext_names <- colnames(d_ext$X)
  sl
}

# readable draws: the names of fp_occu() with psi1_ for occ_, plus col_ and
# ext_
.pretty_draws_dyn <- function(fit, sd_list) {
  sd <- sd_list$stan_data
  des <- sd_list$design
  map <- .name_map(des, sd$estimate_fp == 1, sd$estimate_re == 1,
                   has_verified = sd$N1 > 0,
                   use_fp_re = isTRUE(sd$estimate_fp_re == 1),
                   use_det_re = isTRUE(sd$estimate_det_re == 1))
  map <- sub("^occ_", "psi1_", map)
  map <- c(map,
           stats::setNames(paste0("col_", des$col_names),
                           sprintf("beta_col[%d]", seq_along(des$col_names))),
           stats::setNames(paste0("ext_", des$ext_names),
                           sprintf("beta_ext[%d]", seq_along(des$ext_names))))
  dr <- posterior::as_draws_df(fit$draws(variables = unique(sub("\\[.*", "", names(map)))))
  dr <- posterior::subset_draws(dr, variable = names(map))
  posterior::variables(dr) <- unname(map[posterior::variables(dr)])
  dr
}

# the site-season units of a dynamic data set, as an fp_data with one "site"
# per unit (column .unit); the table of units carries the site covariates,
# the season and the season covariates
.dyn_units <- function(data) {
  s <- data$site
  sites <- data$sites
  nt <- length(data$seasons)
  unit_of <- function(df) {
    (match(df[[s]], sites[[s]]) - 1L) * nt + .season_index(df, data$season, data$seasons)
  }
  auto <- data$auto
  auto$.unit <- unit_of(auto)
  verified <- data$verified
  if (!is.null(verified)) verified$.unit <- unit_of(verified)
  usites <- sites[rep(seq_len(nrow(sites)), each = nt), , drop = FALSE]
  rownames(usites) <- NULL
  usites$.unit <- seq_len(nrow(usites))
  usites[[data$season]] <- rep(data$seasons, nrow(sites))
  sc <- data$season_covs
  if (!is.null(sc)) {
    covs <- setdiff(names(sc), c(s, data$season))
    clash <- intersect(covs, names(sites))
    if (length(clash)) {
      stop("Season covariates share names with site covariates: ",
           paste(clash, collapse = ", "), ".", call. = FALSE)
    }
    t_sc <- .season_index(sc, data$season, data$seasons)
    if (s %in% names(sc)) {
      row_sc <- (match(sc[[s]], sites[[s]]) - 1L) * nt + t_sc
      idx <- match(usites$.unit, row_sc)
    } else {
      idx <- match(rep(seq_len(nt), nrow(sites)), t_sc)
    }
    for (v in covs) usites[[v]] <- sc[[v]][idx]
  }
  # units with a verified detection or a confirmed true observation
  known_present <- rep(0L, nrow(usites))
  if (!is.null(verified)) known_present[verified$.unit[verified[[data$y]] == 1]] <- 1L
  if (!is.null(data$confirmed)) {
    known_present[auto$.unit[which(auto[[data$confirmed[1]]] > 0)]] <- 1L
  }
  data$auto <- auto
  data$verified <- verified
  data$sites <- usites
  data$site <- ".unit"
  data$known_present <- known_present
  data
}

# Marginal log-likelihood of each site, S draws x nsite, as the log_lik of
# fp_occu_dyn.stan: the emissions of the units from .ppc_core(), then the
# forward algorithm over seasons
.dyn_loglik <- function(sd, pars) {
  nsite <- sd$nsite
  nt <- sd$nseason
  unit_site <- rep(seq_len(nsite), each = nt)
  pars_u <- pars
  if (!is.null(pars$e_fp)) pars_u$e_fp <- pars$e_fp[, unit_site, drop = FALSE]
  if (!is.null(pars$d_det)) pars_u$d_det <- pars$d_det[, unit_site, drop = FALSE]
  sd_u <- sd
  sd_u$nsite <- sd$nunit
  em <- .ppc_core(sd_u, pars_u, "automated", "freeman-tukey", emissions_only = TRUE)
  em$unocc[, sd$known_present == 1] <- -Inf
  lse <- function(a, b) {
    m <- pmax(a, b)
    out <- m + log(exp(a - m) + exp(b - m))
    out[m == -Inf] <- -Inf
    out
  }
  lpsi <- pars$beta %*% t(sd$X)
  lcol <- pars$beta_col %*% t(sd$Xcol)
  lext <- pars$beta_ext %*% t(sd$Xext)
  k <- (seq_len(nsite) - 1L) * nt + 1L
  l1 <- stats::plogis(lpsi, log.p = TRUE) + em$occ[, k, drop = FALSE]
  l0 <- stats::plogis(-lpsi, log.p = TRUE) + em$unocc[, k, drop = FALSE]
  for (t in seq_len(nt)[-1]) {
    k <- k + 1L
    n1 <- lse(l1 + stats::plogis(-lext[, k, drop = FALSE], log.p = TRUE),
              l0 + stats::plogis(lcol[, k, drop = FALSE], log.p = TRUE)) +
      em$occ[, k, drop = FALSE]
    n0 <- lse(l1 + stats::plogis(lext[, k, drop = FALSE], log.p = TRUE),
              l0 + stats::plogis(-lcol[, k, drop = FALSE], log.p = TRUE)) +
      em$unocc[, k, drop = FALSE]
    l1 <- n1
    l0 <- n0
  }
  lse(l1, l0)
}

# posterior draws for .dyn_loglik(): those of .ppc_pars() plus the
# coefficients of colonization and extinction
.dyn_pars <- function(fit, sd) {
  pars <- .ppc_pars(fit, sd, NULL)
  dm <- unclass(posterior::as_draws_matrix(fit$draws(variables = c("beta_col", "beta_ext"))))
  pars$beta_col <- dm[, grep("^beta_col\\[", colnames(dm)), drop = FALSE]
  pars$beta_ext <- dm[, grep("^beta_ext\\[", colnames(dm)), drop = FALSE]
  pars
}
