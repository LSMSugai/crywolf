# Dynamic model, first steps (dev/roadmap-v2.md): the season column of
# fp_data(), and the likelihood of fp_occu_dyn.stan, checked against
# fp_occu() with one season and against its R version with several.

test_that("fp_data() orders seasons and reports them", {
  sites <- data.frame(site = 1:3)
  auto <- data.frame(site = c(1, 1, 2, 3), y = c(1, 0, 0, 1),
                     night = c("2024-10-02", "2024-10-01", "2024-10-01", "2024-10-03"))
  ver <- data.frame(site = 1, y = 1, night = "2024-10-03")
  # two seasons lack verified data; that warning is tested below
  quiet <- function(x) suppressWarnings(x, classes = "crywolf_design_warning")
  d <- quiet(fp_data(auto, sites, verified = ver, season = "night"))
  expect_equal(d$seasons, c("2024-10-01", "2024-10-02", "2024-10-03"))
  expect_output(print(d), "seasons: +3 \\(sites with data per season: 1 to 2, median 2\\)")
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
  sd <- crywolf:::.build_dyn_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors())$stan_data
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
# `drop` of the site-seasons lose all their automated data
dyn_data <- function(sim, nseason, drop = 0.1, seed = 1) {
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
  suppressWarnings(fp_data(d$auto, d$sites, verified = d$verified,
                           recorder = d$recorder, confirmed = d$confirmed,
                           season = "season"))
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
    sl <- crywolf:::.build_dyn_stan_data(data, ~ x, ..., priors = fp_priors())
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
    m <- suppressWarnings(crywolf:::fp_occu_dyn(data, occ = ~ x, ..., chains = 1,
                                                iter_warmup = 150, iter_sampling = 20,
                                                threads_per_chain = 1, seed = 3,
                                                refresh = 0, show_messages = FALSE,
                                                sig_figs = 15))
    ll_r <- crywolf:::.dyn_loglik(m$stan_data, crywolf:::.dyn_pars(m$fit, m$stan_data))
    ll <- unclass(posterior::as_draws_matrix(m$fit$draws("log_lik")))
    max(abs(ll_r - ll) / pmax(1, abs(ll)))
  }
  gaps <- c(
    miller = gap(dyn_data(s_miller, 3), det = ~ bio, fp = ~ w, random_recorder = TRUE),
    chambert_random = gap(dyn_data(s_chambert, 3), fp = ~ 1, random_recorder = FALSE,
                          random_fp = TRUE, random_det = TRUE),
    markov_noise = gap(dyn_data(s_chambert, 2), fp = ~ aci, random_recorder = FALSE,
                       markov = "occasion")
  )
  purrr::iwalk(gaps, function(g, name) {
    expect_lt(g, 1e-8, label = paste(name, "relative difference R vs Stan"))
  })
})
