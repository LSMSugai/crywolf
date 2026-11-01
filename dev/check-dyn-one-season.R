# With one season, the log density of fp_occu_dyn.stan equals that of
# fp_occu.stan, apart from the priors of beta_col and beta_ext (criterion of
# step 2 of the dynamic model in dev/roadmap-v2.md). Needs the model methods
# of cmdstanr, which compile both programs again, so it is not part of the
# tests; tests/testthat/test-fp_occu_dyn.R compares their log_lik instead.
# Run from the package root: Rscript dev/check-dyn-one-season.R
devtools::load_all(quiet = TRUE)
s <- simulate_fp(nsite = 60, nocc = 6, frac_confirmed = 0.3, frac_verified = 0.2,
                 seed = 6)
d <- s$data
d$auto$season <- 1L
d$verified$season <- 1L
d <- suppressWarnings(fp_data(d$auto, d$sites, d$verified, recorder = "recorder",
                              confirmed = d$confirmed, season = "season"))
args <- list(d, det = ~ 1, fp = ~ 1, random_recorder = TRUE,
             priors = fp_priors(), random_fp = TRUE, random_det = TRUE)
sd0 <- do.call(crywolf:::.build_stan_data, c(args, occ = ~ x))$stan_data
sd1 <- do.call(crywolf:::.build_dyn_stan_data,
               c(args, psi1 = ~ x, col = ~ 1, ext = ~ 1))$stan_data

dir <- tempfile()
dir.create(dir)
invisible(file.copy(file.path("src", "stan", c("fp_occu.stan", "fp_occu_dyn.stan",
                                               "fp_emissions.stanfunctions")), dir))
compile <- function(f) {
  cmdstanr::cmdstan_model(file.path(dir, f), compile_model_methods = TRUE,
                          cpp_options = list(stan_threads = TRUE), quiet = TRUE)
}
run <- function(mod, sd) {
  fit <- mod$sample(data = sd, chains = 1, iter_warmup = 100, iter_sampling = 10,
                    seed = 1, refresh = 0, show_messages = FALSE, threads_per_chain = 1)
  fit$init_model_methods()
  fit
}
f0 <- run(compile("fp_occu.stan"), sd0)
f1 <- run(compile("fp_occu_dyn.stan"), sd1)

# the unconstrained parameters of fp_occu.stan, then beta_col and beta_ext
up <- unclass(f0$unconstrain_draws(draws = f0$draws(), format = "draws_matrix"))
set.seed(1)
gap <- vapply(seq_len(nrow(up)), function(i) {
  b <- stats::rnorm(2)
  # `~` drops constants, so the priors of beta_col and beta_ext add only
  # their kernels
  prior <- sum(stats::dnorm(b, 0, sd1$prior_sd_int, log = TRUE) -
                 stats::dnorm(0, 0, sd1$prior_sd_int, log = TRUE))
  f1$log_prob(c(up[i, ], b)) - f0$log_prob(up[i, ]) - prior
}, numeric(1))
print(signif(gap, 3))
stopifnot(max(abs(gap)) < 1e-12)
cat("One season: log densities agree to", signif(max(abs(gap)), 3), "\n")
