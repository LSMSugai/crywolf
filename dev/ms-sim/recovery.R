# Short recovery check of fp_occu_ms() (PR 2 of the multispecies model):
# one simulated community, fitted once, posterior against the truth.
# The simulator here is provisional; simulate_fp_ms() comes in PR 4.
# Run from the package root: Rscript dev/ms-sim/recovery.R
devtools::load_all(quiet = TRUE)

sim_ms <- function(nsite = 80, nspecies = 12, nocc = 8, frac_verified = 0.3, nver = 3,
                   mu_beta = c(0, 0.8), sigma_beta = c(1, 0.5),
                   alpha = c(0.2, 0.5), sigma_alpha = 0.7,
                   mu_gap = 2.5, sigma_gap = 0.3, delta0 = 0.5, seed = 1) {
  set.seed(seed)
  sp <- sprintf("sp%02d", seq_len(nspecies))
  sites <- data.frame(site = seq_len(nsite), x = stats::rnorm(nsite))
  beta <- mu_beta + sigma_beta * matrix(stats::rnorm(2 * nspecies), 2)
  alpha_sp <- alpha[1] + sigma_alpha * stats::rnorm(nspecies)
  gap_sp <- mu_gap * exp(sigma_gap * stats::rnorm(nspecies))
  z <- sapply(seq_len(nspecies), function(s) {
    stats::rbinom(nsite, 1, stats::plogis(beta[1, s] + beta[2, s] * sites$x))
  })
  auto <- expand.grid(site = seq_len(nsite), occasion = seq_len(nocc), s = seq_len(nspecies))
  auto$bio <- stats::rnorm(nrow(auto))
  zi <- z[cbind(auto$site, auto$s)]
  p11 <- stats::plogis(alpha_sp[auto$s] + alpha[2] * auto$bio)
  p10 <- stats::plogis(alpha_sp[auto$s] - gap_sp[auto$s])
  auto$y <- stats::rbinom(nrow(auto), 1, ifelse(zi == 1, p11, p10))
  auto$species <- sp[auto$s]
  vs <- sample(nsite, round(frac_verified * nsite))
  ver <- expand.grid(site = vs, k = seq_len(nver), s = seq_len(nspecies))
  ver$y <- stats::rbinom(nrow(ver), 1, z[cbind(ver$site, ver$s)] * stats::plogis(delta0))
  ver$species <- sp[ver$s]
  truth <- c(`occ_(Intercept)` = mu_beta[1], occ_x = mu_beta[2],
             `sd_occ_(Intercept)` = sigma_beta[1], sd_occ_x = sigma_beta[2],
             `det_(Intercept)` = alpha[1], det_bio = alpha[2],
             `sd_det_(Intercept)` = sigma_alpha, fp_gap = mu_gap, sd_log_fp_gap = sigma_gap,
             `verified_(Intercept)` = delta0,
             stats::setNames(beta[1, ], sprintf("occ_(Intercept)[%s]", sp)),
             stats::setNames(alpha_sp, sprintf("det_(Intercept)[%s]", sp)),
             stats::setNames(stats::plogis(alpha_sp - gap_sp), sprintf("p10_ref[%s]", sp)))
  d <- fp_data(auto[c("site", "occasion", "bio", "y", "species")], sites,
               ver[c("site", "y", "species")], species = "species")
  list(data = d, truth = truth)
}

s <- sim_ms()
print(s$data)
t0 <- Sys.time()
m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, fp = ~ 1, chains = 4,
                threads_per_chain = 2, seed = 1)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
tab <- summary(m)$table
tab$truth <- s$truth[tab$variable]
tab <- tab[!is.na(tab$truth), ]
tab$covered <- tab$lower <= tab$truth & tab$truth <= tab$upper
res <- list(table = tab, diagnostics = m$diagnostics, minutes = minutes)
saveRDS(res, file.path("dev", "ms-sim", "recovery.rds"))
print(as.data.frame(tab[c("variable", "truth", "mean", "lower", "upper", "rhat", "covered")]),
      digits = 3)
cat(sprintf("\n%.1f minutes; coverage %d of %d; max R-hat %.3f; divergences %d\n",
            minutes, sum(tab$covered), nrow(tab), m$diagnostics$max_rhat,
            m$diagnostics$divergences))
