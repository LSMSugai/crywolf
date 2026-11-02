# two species on the same sites and occasions: the data of one simulation,
# with the detections of the second species drawn again
ms_sim <- function(seed, ...) {
  sim <- simulate_fp(nsite = 40, nocc = 6, a_bio = 0.5, sd_recorder = 0.5, seed = seed, ...)
  d <- sim$data
  set.seed(seed)
  auto <- transform(d$auto, y = stats::rbinom(nrow(d$auto), 1, 0.3))
  verified <- transform(d$verified, y = stats::rbinom(nrow(d$verified), 1, 0.4))
  b <- suppressWarnings(fp_data(auto, d$sites, verified, recorder = d$recorder))
  tag <- function(x, sp) transform(x, sp = sp)
  ms <- suppressWarnings(fp_data(rbind(tag(d$auto, "a"), tag(auto, "b")), d$sites,
                                 rbind(tag(d$verified, "a"), tag(verified, "b")),
                                 recorder = d$recorder, species = "sp"))
  list(a = d, b = b, ms = ms)
}

test_that("fp_occu_ms() checks its data and arguments", {
  s <- ms_sim(1)
  expect_error(fp_occu_ms(s$a), "species column")
  expect_error(fp_occu_ms(s$ms, markov = "occasion"), "not available")
  sl <- crywolf:::.build_ms_stan_data(s$ms, ~ x, ~ bio, ~ 1, TRUE, fp_priors())
  sd <- sl$stan_data
  expect_equal(c(sd$nsite, sd$nspecies, sd$nunit), c(40, 2, 80))
  expect_equal(sd$sp_unit, rep(1:2, 40))
  # unit (site - 1) * 2 + species; the verified detections of site 1 by species
  kp <- s$ms$known_present
  expect_equal(sd$known_present, as.vector(t(kp)))
  expect_equal(sl$design$species, c("a", "b"))
  expect_equal(sd$X[1:2, 2], rep(s$ms$sites$x[1], 2))
})

test_that("the multispecies likelihood equals that of fp_occu() species by species", {
  skip_if_no_cmdstan()
  gq <- function(model, pars, stan_data) {
    # a folder per call, as cmdstanr names files by the minute
    out <- tempfile("gq")
    dir.create(out)
    g <- crywolf:::.compile_model(model)$generate_quantities(
      fitted_params = pars, data = stan_data, seed = 1, sig_figs = 18,
      threads_per_chain = 1, output_dir = out)
    unclass(posterior::as_draws_matrix(g$draws("log_lik")))
  }
  gap <- function(s, det, fp) {
    m <- suppressWarnings(fp_occu(s$a, occ = ~ x, det = det, fp = fp, chains = 1,
                                  iter_warmup = 100, iter_sampling = 20,
                                  threads_per_chain = 1, seed = 3, refresh = 0,
                                  show_messages = FALSE, sig_figs = 18))
    pa <- posterior::as_draws_df(m$fit$draws())
    pa <- pa[, !grepl("^(lp__|u|gamma0|log_lik|p11_ref|p10_ref|p10_bar|r11)(\\[|$)",
                      names(pa))]
    # species b: its own occupancy coefficients, detection intercept and gap
    set.seed(4)
    pb <- pa
    pb[["beta[1]"]] <- stats::rnorm(20)
    pb[["beta[2]"]] <- stats::rnorm(20)
    pb[["alpha[1]"]] <- stats::rnorm(20)
    pb$fp_gap <- stats::runif(20, 0.5, 3)
    sl_a <- crywolf:::.build_stan_data(s$a, ~ x, det, fp, TRUE, fp_priors())
    sl_b <- crywolf:::.build_stan_data(s$b, ~ x, det, fp, TRUE, fp_priors())
    ll_a <- gq("fp_occu", posterior::as_draws_array(pa), sl_a$stan_data)
    ll_b <- gq("fp_occu", posterior::as_draws_array(pb), sl_b$stan_data)
    # the multispecies parameters that give each species its own values:
    # community means 0 (gap 1), standard deviations 1
    pm <- pa[, !grepl("^(beta|fp_gap)", names(pa))]
    pm[["alpha[1]"]] <- 0
    pm[c("mu_beta[1]", "mu_beta[2]", "sigma_beta[1]", "sigma_beta[2]")] <-
      list(0, 0, 1, 1)
    pm[c("zb[1,1]", "zb[2,1]", "zb[1,2]", "zb[2,2]")] <-
      list(pa[["beta[1]"]], pa[["beta[2]"]], pb[["beta[1]"]], pb[["beta[2]"]])
    pm$sigma_alpha <- 1
    pm[c("za[1]", "za[2]")] <- list(pa[["alpha[1]"]], pb[["alpha[1]"]])
    pm[c("mu_gap", "sigma_gap")] <- list(1, 1)
    pm[c("zg[1]", "zg[2]")] <- list(log(pa$fp_gap), log(pb$fp_gap))
    # CmdStan reads the parameters in the order of their declaration
    first <- c("^mu_beta", "^sigma_beta", "^zb", "^alpha", "^sigma_alpha$", "^za",
               "^mu_gap$", "^sigma_gap$", "^zg", "^u_raw", "^sigma_u$", "^gamma_fp",
               "^delta0$", "^delta_s")
    pm <- pm[, c(unlist(lapply(first, grep, names(pm), value = TRUE)),
                 ".chain", ".iteration", ".draw")]
    sl_m <- crywolf:::.build_ms_stan_data(s$ms, ~ x, det, fp, TRUE, fp_priors())
    ll_m <- gq("fp_occu_ms", posterior::as_draws_array(pm), sl_m$stan_data)
    a <- cbind(ll_a, ll_b)
    b <- ll_m[, c(seq(1, by = 2, length.out = 40), seq(2, by = 2, length.out = 40))]
    max(abs(a - b) / pmax(1, abs(a)))
  }
  s <- ms_sim(3)
  gaps <- c(site_fp = gap(s, ~ bio, ~ x), occasion_fp = gap(s, ~ bio, ~ aci))
  purrr::iwalk(gaps, function(g, name) {
    expect_lt(g, 1e-12, label = paste(name, "relative difference multispecies vs fp_occu"))
  })
})

test_that("a short multispecies fit has readable names and methods", {
  skip_if_no_cmdstan()
  s <- ms_sim(5)
  m <- suppressWarnings(fp_occu_ms(s$ms, occ = ~ x, det = ~ bio, chains = 1,
                                   iter_warmup = 150, iter_sampling = 100,
                                   threads_per_chain = 1, seed = 1, refresh = 0,
                                   show_messages = FALSE))
  expect_s3_class(m, "fp_occu_ms")
  v <- posterior::variables(m$draws)
  expect_true(all(c("occ_(Intercept)", "occ_x", "sd_occ_(Intercept)", "sd_occ_x",
                    "det_(Intercept)", "det_bio", "sd_det_(Intercept)",
                    "fp_(Intercept)", "fp_gap", "sd_log_fp_gap", "p11_ref", "p10_ref",
                    "sd_recorder", "verified_(Intercept)", "p_verified",
                    "occ_x[a]", "occ_(Intercept)[b]", "det_(Intercept)[a]",
                    "fp_gap[b]", "p10_ref[a]", "p11_ref[b]") %in% v))
  expect_false(any(grepl("^(mu_|z|beta|alpha|gap_sp)", v)))
  # every species keeps p11 above p10
  expect_true(all(m$draws[["p11_ref[a]"]] > m$draws[["p10_ref[a]"]]))
  expect_output(print(m), "2 species at 40 sites")
  expect_s3_class(summary(m), "summary.fp_occu")
  expect_identical(as_draws_df(m), m$draws)
})
