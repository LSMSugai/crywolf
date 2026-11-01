# Builds the reference values of the posterior regression test
# (tests/testthat/test-regression.R) and saves them in
# tests/testthat/fixtures/regression-refs.rds.
#
# Refresh the references only after a change that is MEANT to alter the
# posterior or the likelihood (a new prior, a bug fix), and say so in NEWS.md.
# A speed-up or refactor must pass the test with the old references.
#
# Run from the package root: Rscript dev/make-regression-refs.R
devtools::load_all(quiet = TRUE)
source(file.path("tests", "testthat", "helper-regression.R"))

refs <- purrr::imap(.regression_designs(), function(design, name) {
  message("Fitting ", name)
  fit <- .regression_fit(design)
  # 20 parameter draws (5 iterations of each chain), and their log_lik
  pars <- .regression_pars(fit)
  params <- posterior::subset_draws(fit$fit$draws(variables = pars),
                                    iteration = round(seq(1, 500, length.out = 5)))
  summ <- posterior::summarise_draws(fit$draws, "mean", "sd", "mcse_mean")
  list(params = params, log_lik = .regression_log_lik(fit$stan_data, params),
       summary = as.data.frame(summ), diagnostics = fit$diagnostics)
})
attr(refs, "made") <- list(
  date = Sys.Date(), crywolf = as.character(utils::packageVersion("crywolf")),
  cmdstan = cmdstanr::cmdstan_version(), platform = R.version$platform,
  stan_md5 = crywolf:::.stan_md5())
dir.create(file.path("tests", "testthat", "fixtures"), showWarnings = FALSE)
saveRDS(refs, file.path("tests", "testthat", "fixtures", "regression-refs.rds"))
str(attr(refs, "made"))
