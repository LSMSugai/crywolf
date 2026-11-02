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
#' The fit also gives, per season, the expected and the realized occupancy
#' over sites (`psi_mean`, `psi_fs`) and the turnover, and, with constant
#' transitions, the equilibrium occupancy (`psi_eq`); see [crywolf-output].
#' `summary()` lists them with their season; `print()` shows only the
#' parameters.
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
#'   `diagnostics` and `call`. Its methods are `print()`, `summary()`,
#'   `as_draws_df()` and `simulate()`; [fp_sites()] gives the occupancy of
#'   each site-season and [fp_ppc()] checks the fit.
#' @seealso [fp_data()], [fp_occu()], [simulate_fp_dyn()], [crywolf-output];
#'   `vignette("dynamic")` fits it to the `curiango` data, with the nights
#'   of each site as seasons.
#' @export
#' @examples
#' \dontrun{
#' s <- simulate_fp_dyn(nsite = 150, nseason = 4, b_col = 1, b_ext = -1, seed = 1)
#' m <- fp_occu_dyn(s$data, psi1 = ~ x, col = ~ rain, ext = ~ x,
#'                  random_recorder = FALSE)
#' m
#' head(fp_sites(m))   # one row per site-season
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
  # convergence is judged on the parameters, not on the derived quantities
  # (psi_fs is constant when every state is known)
  vars <- posterior::variables(draws)
  diag <- .diagnostics(fit, posterior::subset_draws(draws, variable = vars[!.is_derived(vars)]))
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
  tab <- summary(x)$table
  print(as.data.frame(tab[tab$season == "", names(tab) != "season"]), digits = digits,
        row.names = FALSE)
  cat("\n  Occupancy per season (psi_mean, psi_fs), turnover",
      if ("psi_eq" %in% tab$variable) "and equilibrium occupancy (psi_eq)",
      "in summary().\n")
  .print_diag(x$diagnostics)
  invisible(x)
}

#' @rdname fp_occu_dyn
#' @param prob Width of the equal-tailed credible interval.
#' @return `summary()` returns an object of class `summary.fp_occu`, as for
#'   [fp_occu()], whose table has also a `season` column: the season of the
#'   derived quantities `psi_mean[t]`, `psi_fs[t]` and `turnover[t]`, and
#'   empty for the other rows (see [crywolf-output]).
#' @export
summary.fp_occu_dyn <- function(object, prob = 0.95, ...) {
  s <- summary.fp_occu(object, prob = prob)
  t <- suppressWarnings(as.integer(sub("^[a-z_]+\\[(\\d+)\\]$", "\\1", s$table$variable)))
  t[!.is_derived(s$table$variable)] <- NA
  season <- as.character(object$data$seasons)[t]
  season[is.na(season)] <- ""
  s$table <- tibble::add_column(s$table, season = season, .after = "variable")
  s
}

# derived quantities of fp_occu_dyn.stan, by readable name
.is_derived <- function(v) grepl("^(psi_mean|psi_fs|turnover)\\[|^psi_eq$", v)

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
  # derived quantities; turnover[t] is numbered by the season t it enters
  nt <- sd$nseason
  map <- c(map,
           stats::setNames(sprintf("psi_mean[%d]", seq_len(nt)), sprintf("psi_mean[%d]", seq_len(nt))),
           stats::setNames(sprintf("psi_fs[%d]", seq_len(nt)), sprintf("psi_fs[%d]", seq_len(nt))),
           if (nt > 1) stats::setNames(sprintf("turnover[%d]", 2:nt),
                                       sprintf("turnover[%d]", seq_len(nt - 1))),
           if (nt > 1 && length(des$col_names) == 1 && length(des$ext_names) == 1)
             c(`psi_eq[1]` = "psi_eq"))
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

# Emissions of the units, S draws x nunit: the log-likelihood of each
# site-season's observations given that it is occupied (occ) and given that
# it is not (unocc, -Inf at units known to be occupied), from .ppc_core()
# with the site effects expanded to units. Also returns the unit-level
# stan_data and parameters, for replicating observations.
.dyn_emissions <- function(sd, pars) {
  unit_site <- rep(seq_len(sd$nsite), each = sd$nseason)
  pars_u <- pars
  if (!is.null(pars$e_fp)) pars_u$e_fp <- pars$e_fp[, unit_site, drop = FALSE]
  if (!is.null(pars$d_det)) pars_u$d_det <- pars$d_det[, unit_site, drop = FALSE]
  sd_u <- sd
  sd_u$nsite <- sd$nunit
  em <- .ppc_core(sd_u, pars_u, "automated", "freeman-tukey", emissions_only = TRUE)
  em$unocc[, sd$known_present == 1] <- -Inf
  c(em, list(sd_u = sd_u, pars_u = pars_u))
}

# log(exp(a) + exp(b)), elementwise, with -Inf for both -Inf
.lse <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[m == -Inf] <- -Inf
  out
}

# Forward algorithm over seasons, from the emissions: lists, one element per
# season, of the S x nsite matrices log P(z_t = 1, y_1..t) (a1) and
# log P(z_t = 0, y_1..t) (a0), and the logits of colonization (lcol) and
# extinction (lext) of the units
.dyn_forward <- function(sd, pars, em) {
  nsite <- sd$nsite
  nt <- sd$nseason
  lpsi <- pars$beta %*% t(sd$X)
  lcol <- pars$beta_col %*% t(sd$Xcol)
  lext <- pars$beta_ext %*% t(sd$Xext)
  k <- (seq_len(nsite) - 1L) * nt + 1L
  a1 <- a0 <- vector("list", nt)
  a1[[1]] <- stats::plogis(lpsi, log.p = TRUE) + em$occ[, k, drop = FALSE]
  a0[[1]] <- stats::plogis(-lpsi, log.p = TRUE) + em$unocc[, k, drop = FALSE]
  for (t in seq_len(nt)[-1]) {
    k <- k + 1L
    a1[[t]] <- .lse(a1[[t - 1]] + stats::plogis(-lext[, k, drop = FALSE], log.p = TRUE),
                    a0[[t - 1]] + stats::plogis(lcol[, k, drop = FALSE], log.p = TRUE)) +
      em$occ[, k, drop = FALSE]
    a0[[t]] <- .lse(a1[[t - 1]] + stats::plogis(lext[, k, drop = FALSE], log.p = TRUE),
                    a0[[t - 1]] + stats::plogis(-lcol[, k, drop = FALSE], log.p = TRUE)) +
      em$unocc[, k, drop = FALSE]
  }
  list(a1 = a1, a0 = a0, lcol = lcol, lext = lext)
}

# Marginal log-likelihood of each site, S draws x nsite, as the log_lik of
# fp_occu_dyn.stan
.dyn_loglik <- function(sd, pars) {
  fw <- .dyn_forward(sd, pars, .dyn_emissions(sd, pars))
  nt <- sd$nseason
  .lse(fw$a1[[nt]], fw$a0[[nt]])
}

# Forward filtering, backward sampling. Returns, S draws x nunit (units
# numbered (site - 1) * nseason + season), the probability that each
# site-season is occupied given all the data and the parameters of the draw
# (pz, by the forward-backward algorithm) and one draw of the states from
# their joint conditional posterior (z), plus the emissions.
.dyn_ffbs <- function(sd, pars) {
  em <- .dyn_emissions(sd, pars)
  fw <- .dyn_forward(sd, pars, em)
  nsite <- sd$nsite
  nt <- sd$nseason
  S <- nrow(pars$beta)
  unit <- function(t) (seq_len(nsite) - 1L) * nt + t
  pz <- z <- matrix(NA_real_, S, sd$nunit)
  # backward: log P(y_t+1..T | z_t) for z_t = 1 (b1) and 0 (b0)
  b1 <- b0 <- matrix(0, S, nsite)
  zt <- NULL
  for (t in rev(seq_len(nt))) {
    k <- unit(t)
    if (t < nt) {
      k1 <- unit(t + 1)
      lcol <- fw$lcol[, k1, drop = FALSE]
      lext <- fw$lext[, k1, drop = FALSE]
      e1 <- em$occ[, k1, drop = FALSE] + b1
      e0 <- em$unocc[, k1, drop = FALSE] + b0
      nb1 <- .lse(stats::plogis(-lext, log.p = TRUE) + e1, stats::plogis(lext, log.p = TRUE) + e0)
      nb0 <- .lse(stats::plogis(lcol, log.p = TRUE) + e1, stats::plogis(-lcol, log.p = TRUE) + e0)
      b1 <- nb1
      b0 <- nb0
      # state of season t given the sampled state of t + 1
      logit <- ifelse(zt == 1,
                      (fw$a1[[t]] + stats::plogis(-lext, log.p = TRUE)) -
                        (fw$a0[[t]] + stats::plogis(lcol, log.p = TRUE)),
                      (fw$a1[[t]] + stats::plogis(lext, log.p = TRUE)) -
                        (fw$a0[[t]] + stats::plogis(-lcol, log.p = TRUE)))
    } else {
      logit <- fw$a1[[t]] - fw$a0[[t]]
    }
    zt <- matrix(stats::rbinom(S * nsite, 1, stats::plogis(logit)), S, nsite)
    z[, k] <- zt
    pz[, k] <- stats::plogis((fw$a1[[t]] + b1) - (fw$a0[[t]] + b0))
  }
  list(pz = pz, z = z, em = em)
}

# posterior draws for the dynamic model: those of .ppc_pars() plus the
# coefficients of colonization and extinction, all from the same draws
.dyn_pars <- function(fit, sd, ndraws = NULL) {
  pars <- .ppc_pars(fit, sd, NULL)
  dm <- unclass(posterior::as_draws_matrix(fit$draws(variables = c("beta_col", "beta_ext"))))
  pars$beta_col <- dm[, grep("^beta_col\\[", colnames(dm)), drop = FALSE]
  pars$beta_ext <- dm[, grep("^beta_ext\\[", colnames(dm)), drop = FALSE]
  total <- nrow(pars$beta)
  if (!is.null(ndraws) && ndraws < total) {
    keep <- sort(sample.int(total, ndraws))
    pars <- purrr::map(pars, function(m) {
      if (is.null(m)) NULL else if (is.null(dim(m))) m[keep] else m[keep, , drop = FALSE]
    })
  }
  pars
}

# fp_sites() for the dynamic model: one row per site-season
.fp_sites_dyn <- function(object, prob, ndraws, seed) {
  sd <- object$stan_data
  if (!is.null(seed)) .set_seed(seed)
  pars <- .dyn_pars(object$fit, sd, ndraws)
  fb <- .dyn_ffbs(sd, pars)
  # occupancy from the covariates alone: initial occupancy carried forward
  # by colonization and extinction
  nt <- sd$nseason
  col <- stats::plogis(pars$beta_col %*% t(sd$Xcol))
  ext <- stats::plogis(pars$beta_ext %*% t(sd$Xext))
  psi <- matrix(NA_real_, nrow(pars$beta), sd$nunit)
  k <- (seq_len(sd$nsite) - 1L) * nt + 1L
  psi[, k] <- stats::plogis(pars$beta %*% t(sd$X))
  for (t in seq_len(nt)[-1]) {
    psi[, k + 1L] <- psi[, k] * (1 - ext[, k + 1L]) + (1 - psi[, k]) * col[, k + 1L]
    k <- k + 1L
  }
  units <- .dyn_units(object$data)$sites
  units$.unit <- NULL
  new <- c("psi", "psi_lower", "psi_upper", "occupied", "occupied_lower",
           "occupied_upper", "detections", "verified", "p10")
  clash <- intersect(new, names(units))
  if (length(clash)) {
    stop("The sites table already has column(s) ", paste(clash, collapse = ", "),
         "; rename them before calling fp_sites().", call. = FALSE)
  }
  a <- (1 - prob) / 2
  q <- .col_quantile
  out <- data.frame(
    units,
    psi = colMeans(psi), psi_lower = q(psi, a), psi_upper = q(psi, 1 - a),
    occupied = colMeans(fb$pz), occupied_lower = q(fb$pz, a),
    occupied_upper = q(fb$pz, 1 - a),
    detections = tabulate(rep(sd$site2, sd$y2), nbins = sd$nunit),
    verified = sd$known_present == 1,
    row.names = NULL
  )
  if (sd$estimate_fp == 1) out$p10 <- colMeans(.p10_site_draws(fb$em$sd_u, fb$em$pars_u))
  structure(out, class = c("fp_sites", "data.frame"), prob = prob)
}

#' Simulate data from a fitted dynamic model
#'
#' Generates new detections from an [fp_occu_dyn()] fit, keeping the
#' structure of the original data (sites, seasons, occasions and their
#' covariates, which site-seasons were verified and which occasions were
#' checked), as [simulate.fp_occu()] does for a single season. The latent
#' states are drawn anew over seasons, from initial occupancy, colonization
#' and extinction, and the detections given them. **Experimental**, like
#' [fp_occu_dyn()].
#'
#' @param object An `fp_occu_dyn` object.
#' @inheritParams simulate.fp_occu
#' @return A list of `nsim` elements, each a list with `data` (an `fp_data`
#'   object with the structure of the original data) and `truth` (the
#'   generating values named as the output of [fp_occu_dyn()] in `pars`, the
#'   latent states `z` as a matrix with sites in rows and seasons in
#'   columns, the false-positive rate of each site-season `p10_site`, and
#'   its mean `p10_bar`).
#' @seealso [fp_occu_dyn()], [simulate_fp_dyn()]
#' @export
#' @examples
#' \dontrun{
#' s <- simulate_fp_dyn(nsite = 100, nseason = 4, seed = 1)
#' m <- fp_occu_dyn(s$data, col = ~ rain, ext = ~ x, random_recorder = FALSE)
#' new <- simulate(m, nsim = 2, seed = 1)
#' new[[1]]$truth$pars
#' }
simulate.fp_occu_dyn <- function(object, nsim = 1, seed = NULL,
                                 truth = c("mean", "posterior"), ...) {
  truth <- match.arg(truth)
  f <- object$formulas
  sl <- .build_dyn_stan_data(object$data, f$psi1, f$col, f$ext, f$det, f$fp,
                             object$random_recorder, object$priors,
                             random_fp = isTRUE(object$random_fp),
                             random_det = isTRUE(object$random_det),
                             ver = if (is.null(f$ver)) ~ 1 else f$ver,
                             markov = object$markov)
  sd <- sl$stan_data
  sd_u <- sd
  sd_u$nsite <- sd$nunit
  units <- .dyn_units(object$data)
  gen <- list(data = units, stan_data = sd_u, design = sl$design,
              occasions = sl$occasions)
  draws <- .dyn_pars(object$fit, sd)
  if (!is.null(seed)) .set_seed(seed)
  purrr::map(seq_len(nsim), function(i) {
    .simulate_dyn_from(gen, sd, .pick_pars(draws, truth), object$data)
  })
}

.simulate_dyn_from <- function(gen, sd, p, data) {
  nsite <- sd$nsite
  nt <- sd$nseason
  unit_site <- rep(seq_len(nsite), each = nt)
  # site effects belong to sites: new ones (truth = "mean") or those of the
  # draw, expanded to the site-seasons
  site_effect <- function(e, sigma) {
    if (!is.null(e)) return(e[unit_site])
    if (!is.null(sigma)) return(stats::rnorm(nsite, 0, sigma)[unit_site])
    NULL
  }
  p$e_fp <- site_effect(p$e_fp, p$sigma_fp)
  p$d_det <- site_effect(p$d_det, p$sigma_det)
  col <- stats::plogis(drop(sd$Xcol %*% p$beta_col))
  ext <- stats::plogis(drop(sd$Xext %*% p$beta_ext))
  z <- matrix(0L, nsite, nt)
  z[, 1] <- stats::rbinom(nsite, 1, stats::plogis(drop(sd$X %*% p$beta)))
  for (t in seq_len(nt)[-1]) {
    k <- (seq_len(nsite) - 1L) * nt + t
    z[, t] <- stats::rbinom(nsite, 1, ifelse(z[, t - 1] == 1, 1 - ext[k], col[k]))
  }
  res <- .simulate_from(gen, p, z = as.vector(t(z)))
  strip <- function(df) {
    if (is.null(df)) return(NULL)
    df$.unit <- NULL
    df
  }
  new <- withCallingHandlers(
    fp_data(auto = strip(res$data$auto), sites = data$sites,
            verified = strip(res$data$verified), site = data$site, y = data$y,
            recorder = data$recorder, confirmed = data$confirmed,
            season = data$season, season_covs = data$season_covs),
    crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
  # in the order of the output of fp_occu_dyn(): psi1, col, ext, then the rest
  pars <- res$truth$pars
  names(pars) <- sub("^occ_", "psi1_", names(pars))
  first <- grepl("^psi1_", names(pars))
  pars <- c(pars[first],
            stats::setNames(p$beta_col, paste0("col_", gen$design$col_names)),
            stats::setNames(p$beta_ext, paste0("ext_", gen$design$ext_names)),
            pars[!first])
  list(data = new,
       truth = list(pars = pars, z = z, p10_site = res$truth$p10_site,
                    p10_bar = res$truth$p10_bar))
}
