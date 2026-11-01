# Dynamic (multi-season) false-positive occupancy model -------------------------
#
# Not exported yet: the formulas of colonization and extinction, and the
# names of their parameters, come in the next step (dev/roadmap-v2.md). For
# now colonization and extinction have intercepts only, and the result keeps
# the raw Stan names. The observation model and its arguments are those of
# fp_occu(); `data` needs a season column (fp_data(season = )).

fp_occu_dyn <- function(data, occ = ~ 1, det = ~ 1, fp = ~ 1,
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
  sd_list <- .build_dyn_stan_data(data, occ, det, fp, random_recorder, priors,
                                  random_fp = random_fp, ver = ver, markov = markov,
                                  random_det = random_det)
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
  structure(
    list(fit = fit, stan_data = sd_list$stan_data, design = sd_list$design,
         data = data, priors = priors, call = match.call()),
    class = "fp_occu_dyn"
  )
}

# The data of fp_occu_dyn.stan. Each site-season is a unit, numbered
# (site - 1) * nseason + season; the units form an fp_data whose "sites" are
# the units, so that .build_stan_data() builds the observation model of each
# unit as fp_occu() builds that of each site. The initial occupancy design
# is that of the sites.
.build_dyn_stan_data <- function(data, occ, det, fp, random_recorder, priors,
                                 random_fp = FALSE, ver = ~ 1, markov = NULL,
                                 random_det = FALSE) {
  units <- .dyn_units(data)
  sl <- .build_stan_data(units, occ, det, fp, random_recorder, priors,
                         random_fp = random_fp, ver = ver, markov = markov,
                         random_det = random_det)
  nsite <- nrow(data$sites)
  nseason <- length(data$seasons)
  sd <- sl$stan_data
  d_occ <- .design(occ, data$sites, "occ")
  sd$nunit <- sd$nsite
  sd$nsite <- nsite
  sd$nseason <- nseason
  sd$X <- unname(d_occ$X)
  sd$Kcol <- 1L
  sd$Xcol <- matrix(1, sd$nunit, 1)
  sd$Kext <- 1L
  sd$Xext <- matrix(1, sd$nunit, 1)
  sl$stan_data <- sd
  sl$design$occ <- d_occ
  sl
}

# the site-season units of a dynamic data set, as an fp_data with one "site"
# per unit (column .unit)
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
