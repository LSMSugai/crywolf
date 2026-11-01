# Dynamic model (dev/roadmap-v2.md): the season column of fp_data(), the
# likelihood of fp_occu_dyn.stan, checked against fp_occu() with one season
# and against its R version with several, and the latent states.

test_that("fp_data() orders seasons and reports them", {
  sites <- data.frame(site = 1:3)
  auto <- data.frame(site = c(1, 1, 2, 3), y = c(1, 0, 0, 1),
                     night = c("2024-10-02", "2024-10-01", "2024-10-01", "2024-10-03"))
  ver <- data.frame(site = 1, y = 1, night = "2024-10-03")
  # two seasons lack verified data; that warning is tested below
  quiet <- function(x) suppressWarnings(x, classes = "crywolf_design_warning")
  d <- quiet(fp_data(auto, sites, verified = ver, season = "night"))
  expect_equal(d$seasons, c("2024-10-01", "2024-10-02", "2024-10-03"))
  expect_output(print(d), "seasons: +3, 1 day\\(s\\) apart \\(sites with data per season: 1 to 2, median 2\\)")
  # without season the object is as before
  expect_null(quiet(fp_data(auto, sites, verified = ver))$season)

  f <- transform(auto, night = factor(night, levels = c("2024-10-03", "2024-10-01",
                                                        "2024-10-02")))
  expect_equal(quiet(fp_data(auto = f, sites, verified = ver, season = "night"))$seasons,
               levels(f$night))
  # numbers are matched by value, not by their printed form
  num <- quiet(fp_data(transform(auto, night = c(0.3, 0.1 + 0.2, 0.3, 0.3)), sites,
                       season = "night"))
  expect_equal(crywolf:::.season_index(num$auto, "night", num$seasons), c(1, 2, 1, 1))
  expect_error(fp_data(transform(auto, night = NA), sites, season = "night"),
               "missing values")
  expect_error(fp_data(auto, sites, verified = transform(ver, night = as.Date(night)),
                       season = "night"), "same class")
  expect_error(fp_data(auto, sites, verified = ver[, 1:2], season = "night"),
               "lacks column")
})

test_that("fp_data() warns about seasons without data or verification", {
  sites <- data.frame(site = 1:2)
  auto <- data.frame(site = c(1, 2, 1, 2), y = c(1, 0, 0, 1), t = c(1, 1, 2, 2))
  ver <- data.frame(site = 1, y = 1, t = 1)
  expect_warning(fp_data(auto, sites, verified = ver, season = "t"),
                 "1 season\\(s\\) have no verified or confirmed data \\(e.g. 2\\)",
                 class = "crywolf_design_warning")
  auto$t <- factor(auto$t, levels = 1:3)
  w <- character(0)
  withCallingHandlers(fp_data(auto, sites, verified = rbind(ver, transform(ver, t = 2)),
                              season = "t"),
                      warning = function(cnd) {
                        w <<- c(w, conditionMessage(cnd))
                        invokeRestart("muffleWarning")
                      })
  expect_match(w, "1 season\\(s\\) have no data at any site \\(e.g. 3\\)", all = FALSE)
})

# the value of an expression and the messages of its warnings
with_warnings <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(cnd) {
    msgs <<- c(msgs, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

test_that("fp_data() fills in missing dates and warns about uneven seasons", {
  sites <- data.frame(site = 1:2)
  auto <- data.frame(site = c(1, 2, 1, 2), y = c(1, 0, 0, 1),
                     night = c("20241001", "20241002", "20241004", "20241004"))
  r <- with_warnings(fp_data(auto, sites, season = "night"))
  d <- r$value
  expect_equal(d$seasons, c("20241001", "20241002", "20241003", "20241004"))
  expect_equal(d$season_step, 1)
  expect_match(r$warnings, "1 season\\(s\\) have no data at any site \\(e.g. 20241003\\)",
               all = FALSE)
  expect_output(print(d), "seasons: +4, 1 day\\(s\\) apart")
  # the step is the shortest interval: Date seasons every other day stay so
  dd <- with_warnings(fp_data(transform(auto, night = as.Date("2024-10-01") +
                                          c(0, 2, 6, 6)), sites, season = "night"))$value
  expect_equal(dd$seasons, as.Date(c("2024-10-01", "2024-10-03", "2024-10-05",
                                     "2024-10-07")))
  expect_equal(dd$season_step, 2)
  # dates not on a common step, and numbers, are taken as they come
  r3 <- with_warnings(fp_data(transform(auto, night = c("2024-10-01", "2024-10-03",
                                                        "2024-10-06", "2024-10-06")),
                              sites, season = "night"))
  expect_match(r3$warnings, "not evenly spaced", all = FALSE)
  expect_length(r3$value$seasons, 3)
  expect_null(r3$value$season_step)
  r4 <- with_warnings(fp_data(transform(auto, night = c(1, 2, 4, 4)), sites,
                              season = "night"))
  expect_match(r4$warnings, "unequally spaced", all = FALSE)
  expect_equal(r4$value$seasons, c(1, 2, 4))
  # evenly spaced fractions do not warn, and factor levels are not filled in
  r5 <- with_warnings(fp_data(transform(auto, night = c(0.1, 0.2, 0.3, 0.3)), sites,
                              season = "night"))
  expect_false(any(grepl("unequally", r5$warnings)))
  f <- factor(auto$night, levels = c("20241004", "20241001", "20241002"))
  r6 <- with_warnings(fp_data(transform(auto, night = f), sites, season = "night"))
  expect_equal(r6$value$seasons, levels(f))
  expect_null(r6$value$season_step)
})

test_that("season covariates are validated and reach the units", {
  sites <- data.frame(site = c("a", "b"), x = c(-1, 1))
  auto <- data.frame(site = c("a", "b", "a", "b"), y = c(1, 0, 0, 1), t = c(1, 1, 2, 2))
  ver <- data.frame(site = "a", y = 1, t = 1:2)
  sc <- data.frame(t = 1:2, rain = c(5, 7))
  quiet <- function(x) suppressWarnings(x, classes = "crywolf_design_warning")
  d <- quiet(fp_data(auto, sites, verified = ver, season = "t", season_covs = sc))
  expect_output(print(d), "season covariates: +rain")
  u <- crywolf:::.dyn_units(d)
  expect_equal(u$sites$rain, c(5, 7, 5, 7))
  sd <- crywolf:::.build_dyn_stan_data(d, ~ x, ~ rain, ~ x, ~ 1, ~ 1, FALSE,
                                       fp_priors())$stan_data
  # rows of the first season are unused by colonization and extinction
  expect_equal(sd$Xcol, cbind(c(0, 1, 0, 1), c(0, 7, 0, 7)))
  expect_equal(sd$Xext, cbind(c(0, 1, 0, 1), c(0, -1, 0, 1)))
  expect_equal(sd$X, cbind(1, c(-1, 1)), ignore_attr = TRUE)
  site_sc <- data.frame(site = c("a", "b", "a"), t = c(1, 1, 2), rain = 1:3)
  d2 <- quiet(fp_data(auto, sites, verified = ver, season = "t", season_covs = site_sc))
  expect_equal(crywolf:::.dyn_units(d2)$sites$rain, c(1, 3, 2, NA))
  # (b, 2) has no rain: an error only for a formula that uses it
  expect_error(crywolf:::.build_dyn_stan_data(d2, ~ 1, ~ rain, ~ 1, ~ 1, ~ 1, FALSE,
                                              fp_priors()), "Missing values")
  expect_error(fp_data(auto, sites, season = "t", season_covs = rbind(sc, sc)),
               "one row per t")
  expect_error(fp_data(auto, sites, season = "t", season_covs = data.frame(t = 3, r = 1)),
               "not seasons of the data")
  expect_error(fp_data(auto, sites, season_covs = sc), "needs `season`")
  d3 <- quiet(fp_data(auto, sites, season = "t", season_covs = transform(sc, x = 1)))
  expect_error(crywolf:::.dyn_units(d3), "share names")
  expect_error(crywolf:::.build_dyn_stan_data(d, ~ 0 + x, ~ 1, ~ 1, ~ 1, ~ 1, FALSE,
                                              fp_priors()), "keep the intercept")
})

test_that("site-season units carry the data and known presences of their season", {
  sites <- data.frame(site = c("a", "b"), x = c(-1, 1))
  auto <- data.frame(site = c("a", "a", "b"), y = c(1, 0, 1), t = c(1, 2, 2),
                     n_true = c(NA, NA, 1), n_false = c(NA, NA, 0))
  ver <- data.frame(site = "a", y = 1, t = 2)
  d <- suppressWarnings(fp_data(auto, sites, verified = ver, season = "t",
                                confirmed = c("n_true", "n_false")))
  u <- crywolf:::.dyn_units(d)
  expect_equal(u$auto$.unit, c(1, 2, 4))
  expect_equal(u$verified$.unit, 2)
  expect_equal(u$sites$x, c(-1, -1, 1, 1))
  # unit 2 (a, season 2) has a verified detection, unit 4 (b, 2) a confirmed one
  expect_equal(u$known_present, c(0L, 1L, 0L, 1L))
})

test_that("the R forward algorithm matches a hand computation", {
  # site 1 has no data; site 2 one automated detection in season 2; site 3 a
  # verified detection in season 1 and nothing in season 2
  sites <- data.frame(site = 1:3)
  auto <- data.frame(site = 2, y = 1, t = 2)
  ver <- data.frame(site = 3, y = 1, t = 1)
  d <- suppressWarnings(fp_data(auto, sites, verified = ver, season = "t"))
  sd <- crywolf:::.build_dyn_stan_data(d, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, FALSE, fp_priors())$stan_data
  one <- function(p) matrix(stats::qlogis(p), 1, 1)
  pars <- list(beta = one(0.6), alpha = one(0.7), delta0 = stats::qlogis(0.8),
               delta_s = matrix(0, 1, 0), gamma0 = stats::qlogis(0.1),
               gamma_fp = matrix(0, 1, 0), beta_col = one(0.3), beta_ext = one(0.2))
  psi2 <- 0.6 * (1 - 0.2) + 0.4 * 0.3
  expect_equal(c(crywolf:::.dyn_loglik(sd, pars)),
               c(0, log(psi2 * 0.7 + (1 - psi2) * 0.1), log(0.6 * 0.8)),
               tolerance = 1e-12)
})

# data in which the occasions of each site fall into seasons; a fraction
# `drop` of the site-seasons lose all their automated data; with rain, a
# covariate of each site and season
dyn_data <- function(sim, nseason, drop = 0.1, seed = 1, rain = FALSE) {
  set.seed(seed)
  d <- sim$data
  nocc <- max(d$auto$occasion)
  d$auto$season <- ceiling(d$auto$occasion * nseason / nocc)
  unit <- paste(d$auto$site, d$auto$season)
  gone <- unique(unit)[stats::runif(length(unique(unit))) < drop]
  d$auto <- d$auto[!unit %in% gone, , drop = FALSE]
  if (!is.null(d$verified)) {
    d$verified$season <- sample.int(nseason, nrow(d$verified), replace = TRUE)
  }
  sc <- if (rain) {
    expand.grid(site = d$sites$site, season = seq_len(nseason)) |>
      transform(rain = stats::rnorm(nrow(d$sites) * nseason))
  }
  suppressWarnings(fp_data(d$auto, d$sites, verified = d$verified,
                           recorder = d$recorder, confirmed = d$confirmed,
                           season = "season", season_covs = sc))
}

test_that("with one season the dynamic likelihood equals that of fp_occu()", {
  skip_if_no_cmdstan()
  s_miller <- simulate_fp(nsite = 60, nocc = 6, a_bio = 0.5, sd_recorder = 0.5, seed = 3)
  s_chambert <- simulate_fp(nsite = 60, nocc = 6, frac_confirmed = 0.4,
                            frac_verified = 0.2, seed = 6)
  one_season <- function(sim) dyn_data(sim, 1, drop = 0)
  gap <- function(data, ...) {
    m <- suppressWarnings(fp_occu(data, occ = ~ x, ..., chains = 1, iter_warmup = 100,
                                  iter_sampling = 20, threads_per_chain = 1, seed = 3,
                                  refresh = 0, show_messages = FALSE, sig_figs = 18))
    sl <- crywolf:::.build_dyn_stan_data(data, ~ x, ~ 1, ~ 1, ..., priors = fp_priors())
    # the parameters of fp_occu.stan, then those of colonization and extinction
    pars <- m$fit$draws(format = "draws_array")
    pars <- posterior::subset_draws(pars, variable = posterior::variables(pars)[
      !grepl("^(lp__|u|gamma0|log_lik|p11_ref|p10_ref|p10_bar|r11)(\\[|$)",
             posterior::variables(pars))])
    set.seed(2)
    pars <- posterior::bind_draws(pars, posterior::draws_array(
      `beta_col[1]` = array(stats::rnorm(20), c(20, 1)),
      `beta_ext[1]` = array(stats::rnorm(20), c(20, 1))), along = "variable")
    gq <- crywolf:::.compile_model("fp_occu_dyn")$generate_quantities(
      fitted_params = pars, data = sl$stan_data, seed = 1, sig_figs = 18,
      threads_per_chain = 1)
    a <- unclass(posterior::as_draws_matrix(m$fit$draws("log_lik")))
    b <- unclass(posterior::as_draws_matrix(gq$draws("log_lik")))
    max(abs(a - b) / pmax(1, abs(a)))
  }
  gaps <- c(
    miller = gap(one_season(s_miller), det = ~ bio, fp = ~ w, random_recorder = TRUE),
    chambert_random = gap(one_season(s_chambert), det = ~ 1, fp = ~ 1,
                          random_recorder = FALSE, random_fp = TRUE, random_det = TRUE),
    markov = gap(one_season(s_miller), det = ~ 1, fp = ~ 1, random_recorder = FALSE,
                 markov = "occasion")
  )
  purrr::iwalk(gaps, function(g, name) {
    expect_lt(g, 1e-12, label = paste(name, "relative difference dynamic vs fp_occu"))
  })
})

test_that("the R likelihood matches the Stan log_lik over several seasons", {
  skip_if_no_cmdstan()
  s_miller <- simulate_fp(nsite = 60, nocc = 6, a_bio = 0.5, sd_recorder = 0.5, seed = 3)
  s_chambert <- simulate_fp(nsite = 60, nocc = 6, frac_confirmed = 0.4,
                            frac_verified = 0.3, seed = 6)
  gap <- function(data, ...) {
    m <- suppressWarnings(fp_occu_dyn(data, psi1 = ~ x, ..., chains = 1,
                                      iter_warmup = 150, iter_sampling = 20,
                                      threads_per_chain = 1, seed = 3, refresh = 0,
                                      show_messages = FALSE, sig_figs = 15))
    ll_r <- crywolf:::.dyn_loglik(m$stan_data, crywolf:::.dyn_pars(m$fit, m$stan_data))
    ll <- unclass(posterior::as_draws_matrix(m$fit$draws("log_lik")))
    max(abs(ll_r - ll) / pmax(1, abs(ll)))
  }
  gaps <- c(
    miller = gap(dyn_data(s_miller, 3), col = ~ x, ext = ~ w, det = ~ bio, fp = ~ w,
                 random_recorder = TRUE),
    chambert_random = gap(dyn_data(s_chambert, 3, rain = TRUE), col = ~ rain,
                          ext = ~ rain + x, fp = ~ 1, random_recorder = FALSE,
                          random_fp = TRUE, random_det = TRUE),
    markov_noise = gap(dyn_data(s_chambert, 2), fp = ~ aci, random_recorder = FALSE,
                       markov = "occasion")
  )
  purrr::iwalk(gaps, function(g, name) {
    expect_lt(g, 1e-8, label = paste(name, "relative difference R vs Stan"))
  })
})

test_that("fp_occu_dyn() names its parameters and prints", {
  skip_if_no_cmdstan()
  sim <- simulate_fp(nsite = 60, nocc = 6, frac_verified = 0.3, seed = 6)
  m <- suppressWarnings(fp_occu_dyn(dyn_data(sim, 3, rain = TRUE), psi1 = ~ x,
                                    col = ~ rain, ext = ~ 1, fp = ~ 1,
                                    random_recorder = FALSE, chains = 1,
                                    iter_warmup = 100, iter_sampling = 50,
                                    threads_per_chain = 1, seed = 1, refresh = 0,
                                    show_messages = FALSE))
  expect_s3_class(m, "fp_occu_dyn")
  expect_setequal(posterior::variables(as_draws_df(m)),
                  c("psi1_(Intercept)", "psi1_x", "col_(Intercept)", "col_rain",
                    "ext_(Intercept)", "det_(Intercept)", "fp_(Intercept)", "fp_gap",
                    "verified_(Intercept)", "p_verified", "p11_ref", "p10_ref",
                    "p10_bar"))
  # short chains: posterior warns that the ESS is capped
  expect_s3_class(suppressWarnings(summary(m)), "summary.fp_occu")
  out <- capture.output(suppressWarnings(print(m)))
  expect_match(out, "colonization: +~rain", all = FALSE)
  expect_match(out, "60 sites over 3 seasons", all = FALSE)
})

# small data for the latent states: site 1 has no data in season 2, site 3 a
# verified detection in season 2, so its state there is known
ffbs_case <- function() {
  set.seed(4)
  sites <- data.frame(site = 1:4, x = c(-1, 0, 0.5, 1))
  auto <- expand.grid(site = 1:4, t = 1:3, occasion = 1:3)
  auto <- auto[!(auto$site == 1 & auto$t == 2), ]
  auto$y <- stats::rbinom(nrow(auto), 1, 0.4)
  ver <- data.frame(site = 3, y = 1, t = 2)
  d <- suppressWarnings(fp_data(auto, sites, verified = ver, season = "t"))
  sd <- crywolf:::.build_dyn_stan_data(d, ~ 1, ~ x, ~ 1, ~ 1, ~ 1, FALSE,
                                       fp_priors())$stan_data
  pars <- list(beta = matrix(c(0.2, -0.5), 2), alpha = matrix(c(0.4, 1), 2),
               delta0 = c(1, 1.5), delta_s = matrix(0, 2, 0),
               gamma0 = c(-2, -1.5), gamma_fp = matrix(0, 2, 0),
               beta_col = cbind(c(-1, -0.5), c(0.8, -0.3)),
               beta_ext = matrix(c(-1.2, -0.4), 2))
  list(sd = sd, pars = pars)
}

# joint posterior of the state sequences of each site, by enumerating them:
# a list over draws of nsite x 2^nseason matrices, sequences in the rows of seqs
enumerate_states <- function(sd, pars, seqs) {
  em <- crywolf:::.dyn_emissions(sd, pars)
  nt <- sd$nseason
  purrr::map(seq_len(nrow(pars$beta)), function(s) {
    t(vapply(seq_len(sd$nsite), function(i) {
      k <- (i - 1) * nt + seq_len(nt)
      col <- stats::plogis(drop(sd$Xcol[k, , drop = FALSE] %*% pars$beta_col[s, ]))
      ext <- stats::plogis(drop(sd$Xext[k, , drop = FALSE] %*% pars$beta_ext[s, ]))
      psi <- stats::plogis(sum(sd$X[i, ] * pars$beta[s, ]))
      lp <- apply(seqs, 1, function(z) {
        prior <- log(ifelse(z[1] == 1, psi, 1 - psi))
        for (t in seq_len(nt)[-1]) {
          p1 <- if (z[t - 1] == 1) 1 - ext[t] else col[t]
          prior <- prior + log(ifelse(z[t] == 1, p1, 1 - p1))
        }
        prior + sum(ifelse(z == 1, em$occ[s, k], em$unocc[s, k]))
      })
      exp(lp - max(lp)) / sum(exp(lp - max(lp)))
    }, numeric(nrow(seqs))))
  })
}

test_that("forward-backward and FFBS match the enumeration of state sequences", {
  cs <- ffbs_case()
  sd <- cs$sd
  seqs <- as.matrix(expand.grid(rep(list(0:1), sd$nseason)))
  post <- enumerate_states(sd, cs$pars, seqs)
  fb <- crywolf:::.dyn_ffbs(sd, cs$pars)
  for (s in 1:2) {
    exact <- post[[s]] %*% seqs   # nsite x nseason
    expect_equal(fb$pz[s, ], c(t(exact)), tolerance = 1e-10)
  }
  # the known state has probability 1
  expect_equal(fb$pz[, (3 - 1) * sd$nseason + 2], c(1, 1))
  # sampled sequences: draw 1 repeated many times
  n <- 20000
  rep1 <- purrr::map(cs$pars, function(m) {
    if (is.null(dim(m))) rep(m[1], n) else m[rep(1, n), , drop = FALSE]
  })
  set.seed(1)
  z <- crywolf:::.dyn_ffbs(sd, rep1)$z
  nt <- sd$nseason
  for (i in seq_len(sd$nsite)) {
    code <- z[, (i - 1) * nt + seq_len(nt)] %*% 2^(seq_len(nt) - 1)
    freq <- tabulate(code + 1, nbins = 2^nt) / n
    expect_lt(max(abs(freq - post[[1]][i, ])), 0.015,
              label = paste("site", i, "sampled vs exact sequence frequencies"))
  }
})

test_that("fp_sites(), fp_ppc() and simulate() work with the dynamic model", {
  skip_if_no_cmdstan()
  s <- simulate_fp_dyn(nsite = 150, nseason = 4, b_col = 1, b_ext = -1, seed = 2)
  m <- suppressWarnings(fp_occu_dyn(s$data, psi1 = ~ x, col = ~ rain, ext = ~ x,
                                    random_recorder = FALSE, chains = 2,
                                    iter_warmup = 400, iter_sampling = 400,
                                    threads_per_chain = 1, seed = 1, refresh = 0,
                                    show_messages = FALSE))
  # recovery: the generating values within the 95% intervals
  dr <- as_draws_df(m)
  for (p in names(s$truth$pars)) {
    ci <- stats::quantile(dr[[p]], c(0.025, 0.975))
    expect_true(s$truth$pars[[p]] > ci[1] && s$truth$pars[[p]] < ci[2],
                label = paste(p, "within its 95% interval"))
  }

  st <- fp_sites(m, ndraws = 200, seed = 1)
  expect_s3_class(st, "fp_sites")
  expect_equal(nrow(st), 150 * 4)
  expect_equal(st$season, rep(1:4, 150))
  expect_true(all(c("x", "rain", "psi", "occupied", "p10") %in% names(st)))
  expect_true(all(st$occupied[st$verified] > 1 - 1e-12))
  expect_gt(cor(st$occupied, c(t(s$truth$z))), 0.8)

  pp <- fp_ppc(m, ndraws = 100, seed = 1)
  expect_s3_class(pp, "fp_ppc")
  expect_length(pp$season, length(pp$site))
  expect_equal(ncol(pp$yrep), length(pp$site))
  expect_true(pp$p_value > 0 && pp$p_value < 1)

  new <- simulate(m, nsim = 2, seed = 1)
  expect_length(new, 2)
  expect_s3_class(new[[1]]$data, "fp_data")
  expect_equal(dim(new[[1]]$truth$z), c(150, 4))
  expect_setequal(names(new[[1]]$truth$pars), names(s$truth$pars))
  expect_equal(nrow(new[[1]]$data$auto), nrow(s$data$auto))
  expect_equal(new[[1]]$data$season_covs, s$data$season_covs)
})
