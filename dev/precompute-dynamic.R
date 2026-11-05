# Pre-computes the results of vignettes/dynamic.Rmd: the curiango data with
# the nights of each site as seasons (dev/curiango-dyn-common.R) fitted with
# fp_occu_dyn() on all the data, its posterior predictive check and the
# probability that each site was in use on each night, plus the tables of
# the blocked cross-validation (dev/dyn-kfold/results.rds, from the AWS run
# of dev/dyn-kfold/run_kfold.R, which this script does not rerun). About an
# hour on 8 cores. The results keep 1,000 draws of the parameters plotted in the
# vignette, and the tables of dev/dyn-val/results.rds and dev/dyn-pika/results.rds
# that it plots. Run those first.
# Run from the package root after installing the package:
#   Rscript dev/precompute-dynamic.R
library(crywolf)
source(file.path("dev", "curiango-dyn-common.R"))
d <- suppressWarnings(fp_data(occasions, data.frame(site = sort(unique(occasions$site))),
                              verified, site = "site", recorder = "recorder",
                              season = "night", season_covs = daily))
t0 <- Sys.time()
fit <- fp_occu_dyn(d, psi1 = forms$psi1, col = forms$col, ext = forms$ext,
                   det = forms$det, fp = forms$fp, random_recorder = TRUE, seed = 1,
                   chains = 4, parallel_chains = 4, threads_per_chain = 2,
                   refresh = 0, show_messages = FALSE)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
keep <- c("col_(Intercept)", "col_temp_d", "col_rh_d", "ext_(Intercept)",
          "ext_temp_d", "ext_rh_d", "p11_ref", "p10_ref")
draws <- as.data.frame(as_draws_df(fit))[, keep]
draws <- draws[seq(1, nrow(draws), length.out = 1000), ]
ppc <- fp_ppc(fit, ndraws = 300, seed = 1)
cv <- readRDS(file.path("dev", "dyn-kfold", "results.rds"))
val <- readRDS(file.path("dev", "dyn-val", "results.rds"))
pika <- readRDS(file.path("dev", "dyn-pika", "results.rds"))
out <- list(
  summary = summary(fit)$table,
  draws = draws,
  diagnostics = fit$diagnostics,
  minutes = minutes,
  sites_per_night = table(factor(unique(occasions[c("site", "night")])$night,
                                 levels = as.character(d$seasons))),
  ppc = list(p_value = ppc$p_value, fit_obs = ppc$fit_obs, fit_rep = ppc$fit_rep),
  sites = as.data.frame(fp_sites(fit, ndraws = 1000, seed = 1)),
  cv = cv[c("total", "pairs", "by_site", "by_fold", "diagnostics")],
  cv_units = nrow(cv$units),
  val = val[c("lowp", "lowp_fits", "calib", "randet", "markov")],
  pika = pika[c("diagnostics", "summary", "draws", "occ_park", "sensitivity", "loo", "ppc")]
)
saveRDS(crywolf:::.stamp(out), file.path("vignettes", "dynamic-results.rds"), compress = "xz")
