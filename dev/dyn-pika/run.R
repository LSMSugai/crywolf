# =============================================================================
# dev/dyn-pika/run.R
# Criterion 8 of dev/dyn-val/README.md: fp_occu_dyn() on a second data set,
# independent of curiango. American pika surveys of Goldman et al. (2023,
# Ecology and Evolution 13: e10019), post_process.csv from Dryad
# (doi:10.5061/dryad.ttdz08m2n, CC0), collected by the OSU Cascades Human and
# Ecosystem Resilience and Sustainability Lab with support from the National
# Park Service.
#
# Plots are sites and years are seasons, 2010 to 2021, with 2015 and 2016 as
# seasons without data. Each survey gives three automated occasions, one per
# ambiguous sign (fresh scat without an unambiguous mark, fresh haypile, call),
# and one verified occasion (a sighting or unambiguous fresh scat), following
# the classification of the article. The article models each sign with its own
# detection and false-detection rates; here sign enters det and fp as a
# covariate, and p11 and the gap are otherwise shared.
#
# Models: dyn (psi1 ~ park, col ~ 1, ext ~ 1), three priors of fp_gap, and ind
# (fp_occu() with plot-years as sites, occ ~ park + year). Both are compared by
# PSIS-LOO with whole plots left out, summing the log-likelihood of ind over the
# years of each plot. Every fit: 4 chains of 1000 + 1000 iterations, adapt_delta
# = 0.95, seed 2026.
# Run from the package root: Rscript dev/dyn-pika/run.R
# =============================================================================

library(crywolf)
dir <- file.path("dev", "dyn-pika")
raw <- read.csv(file.path(dir, "post_process.csv"), na.strings = c("NA", ""))
raw$FS_U[is.na(raw$FS_U)] <- 0
years <- as.character(2010:2021)

occasions <- function(unit) {
  sign <- function(name, y) {
    data.frame(site = unit, sign = name, y = as.integer(y), year = factor(raw$Year, levels = years))
  }
  rbind(sign("scat", raw$FS > 0 & raw$FS_U == 0), sign("haypile", raw$FH > 0),
        sign("call", raw$PC > 0))
}
verified <- function(unit) {
  data.frame(site = unit, y = as.integer(raw$PS > 0 | raw$FS_U > 0),
             year = factor(raw$Year, levels = years))
}

plots <- unique(raw[c("plot", "Park")])
d_dyn <- fp_data(occasions(raw$plot), data.frame(site = plots$plot, park = plots$Park),
                 verified(raw$plot), season = "year")
unit <- paste(raw$plot, raw$Year)
units <- unique(data.frame(site = unit, plot = raw$plot, park = raw$Park,
                           year = factor(raw$Year)))
d_ind <- fp_data(occasions(unit)[c("site", "sign", "y")], units, verified(unit)[c("site", "y")])

args <- list(det = ~ sign, fp = ~ sign, random_recorder = FALSE, chains = 4, parallel_chains = 4,
             threads_per_chain = 2, iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
             seed = 2026, refresh = 0, show_messages = FALSE)
dyn <- function(priors) {
  do.call(fp_occu_dyn, c(list(d_dyn, psi1 = ~ park, col = ~ 1, ext = ~ 1, priors = priors), args))
}
base <- fp_priors()
fits <- list(default = dyn(base),
             gap_low = dyn(fp_priors(gap_mean = base$gap_mean - 1)),
             gap_high = dyn(fp_priors(gap_mean = base$gap_mean + 1)),
             gap_wide = dyn(fp_priors(gap_sd = 2 * base$gap_sd)))
ind <- do.call(fp_occu, c(list(d_ind, occ = ~ park + year, priors = base), args))

keep <- c("col_(Intercept)", "ext_(Intercept)", "fp_(Intercept)", "fp_gap", "p11_ref", "p10_ref",
          "p10_bar", "p_verified", "psi_eq")
sens <- do.call(rbind, lapply(names(fits), function(n) {
  s <- as.data.frame(summary(fits[[n]])$table)
  cbind(prior = n, s[s$variable %in% keep, c("variable", "mean", "lower", "upper")])
}))

# log-likelihood per plot, iterations x chains x plots
ll_dyn <- posterior::as_draws_array(fits$default$fit$draws("log_lik"))
ll_ind <- posterior::as_draws_array(ind$fit$draws("log_lik"))
to_plot <- split(seq_len(dim(ll_ind)[3]), match(units$plot[match(d_ind$sites$site, units$site)],
                                                 d_dyn$sites$site))
ll_sum <- array(vapply(to_plot, function(i) apply(unclass(ll_ind)[, , i, drop = FALSE], 1:2, sum),
                       matrix(0, dim(ll_ind)[1], dim(ll_ind)[2])),
                c(dim(ll_ind)[1:2], length(to_plot)))
psis <- function(ll) loo::loo(ll, r_eff = loo::relative_eff(exp(sweep(ll, 3, apply(ll, 3, max), "-"))))
l_dyn <- psis(unclass(ll_dyn))
l_ind <- psis(ll_sum)

ppc <- list(dyn_auto = fp_ppc(fits$default, ndraws = 300, seed = 1),
            dyn_ver = fp_ppc(fits$default, method = "verified", ndraws = 300, seed = 1),
            ind_auto = fp_ppc(ind, ndraws = 300, seed = 1))
res <- list(sizes = c(plots = nrow(plots), surveys = nrow(raw)),
            diagnostics = c(lapply(fits, `[[`, "diagnostics"), list(ind = ind$diagnostics)),
            summary = as.data.frame(summary(fits$default)$table), ind_summary = as.data.frame(summary(ind)$table),
            sensitivity = sens, loo = loo::loo_compare(list(dyn = l_dyn, ind = l_ind)),
            high_k = c(dyn = sum(l_dyn$diagnostics$pareto_k > 0.7),
                       ind = sum(l_ind$diagnostics$pareto_k > 0.7)),
            ppc = lapply(ppc, function(p) list(p_value = p$p_value, fit_obs = stats::median(p$fit_obs),
                                               fit_rep = stats::median(p$fit_rep))))
saveRDS(res, file.path(dir, "results.rds"))
str(res[c("sizes", "diagnostics", "high_k", "ppc")], max.level = 3)
print(res$loo)
print(sens, digits = 3)
