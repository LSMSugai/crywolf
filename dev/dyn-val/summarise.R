# =============================================================================
# dev/dyn-val/summarise.R
# Reads the job files of run_val.R (dev/dyn-val/jobs/), applies criteria 1 to 7
# of dev/dyn-val/README.md and writes dev/dyn-val/results.rds and results.md.
# Colonization, extinction and initial occupancy are on the probability scale,
# from their logit-scale intercepts; psi_fs is pooled over seasons. A fit
# counts as converged with R-hat < 1.01 for every parameter with a true value
# and no divergent transitions. Criterion 2 pools the cells of each combination
# of p11, verification and seasons, that is, over p10 and occasions.
# Run from the package root:  Rscript dev/dyn-val/summarise.R
# =============================================================================

job_dir <- file.path("dev", "dyn-val", "jobs")
read_part <- function(part) {
  f <- list.files(job_dir, pattern = paste0("^", part, "-"), full.names = TRUE)
  do.call(rbind, lapply(f, readRDS))
}
parts <- c("lowp", "calib", "randet", "markov")
d <- stats::setNames(lapply(parts, read_part), parts)
sbc <- read_part("sbc")

core <- c("psi1", "col", "ext", "p10_bar", "psi_fs")
fit_keys <- c("rep", "cell", "model")

prepare <- function(x) {
  x$param <- sub("\\[\\d+\\]$", "", x$variable)
  inv <- x$param %in% c("psi1_(Intercept)", "col_(Intercept)", "ext_(Intercept)")
  x$param[inv] <- sub("_\\(Intercept\\)$", "", x$param[inv])
  for (v in c("mean", "lower", "upper", "truth")) x[[v]][inv] <- stats::plogis(x[[v]][inv])
  bad <- !x$ok | (!is.na(x$rhat) & x$rhat >= 1.01) | (!is.na(x$divergences) & x$divergences > 0)
  x$converged <- !ave(bad, interaction(x[fit_keys], drop = TRUE), FUN = any)
  x
}
d <- lapply(d, prepare)

# share of converged fits per cell and model
fits <- function(x, by) {
  u <- unique(x[c(fit_keys, by, "converged", "ok")])
  stats::aggregate(cbind(converged, ok) ~ ., data = u[c(by, "model", "converged", "ok")],
                   FUN = mean)
}
# bias, RMSE, coverage and width per group and parameter
metrics <- function(x, by) {
  key <- interaction(x[c(by, "param")], drop = TRUE, sep = "|")
  do.call(rbind, lapply(split(x, key), function(g) {
    data.frame(g[1, c(by, "param")], n = nrow(g), bias = mean(g$mean - g$truth),
               rmse = sqrt(mean((g$mean - g$truth)^2)),
               coverage = mean(g$lower <= g$truth & g$truth <= g$upper),
               width = mean(g$upper - g$lower), row.names = NULL)
  }))
}
# the bias threshold is for the quantities of criterion 2; criterion 7 asks only
# coverage of sd_det_site and det_y_prev; no criterion covers fp_x
only <- function(m, params = core) {
  m <- m[m$param %in% params, ]
  m$pass <- m$coverage >= 0.90 & (abs(m$bias) <= 0.05 | !m$param %in% core)
  m$pass[m$param == "fp_x"] <- NA
  m[do.call(order, m[setdiff(names(m), c("n", "bias", "rmse", "coverage", "width", "pass"))]), ]
}

lowp_cell <- c("p11", "p10", "frac_verified", "nseason", "nocc")
res <- list(
  lowp_fits = fits(d$lowp, lowp_cell),
  calib_fits = fits(d$calib, c("frac_verified", "frac_confirmed", "b_fp")),
  randet_fits = fits(d$randet, "sd_det"),
  markov_fits = fits(d$markov, "a"),
  lowp = only(metrics(d$lowp, c("p11", "frac_verified", "nseason"))),
  lowp_cells = only(metrics(d$lowp, lowp_cell)),
  calib = only(metrics(d$calib, c("frac_verified", "frac_confirmed", "b_fp", "model")),
               c(core, "fp_x")),
  randet = only(metrics(d$randet, c("sd_det", "model")), c(core, "sd_det_site")),
  markov = only(metrics(d$markov, c("a", "model")), c(core, "det_y_prev"))
)

# criterion 6: 100 possible ranks (0 to 99) in 10 bins of 10
res$sbc <- do.call(rbind, lapply(split(sbc, sbc$variable), function(x) {
  counts <- tabulate(floor(x$rank / 10) + 1, nbins = 10)
  data.frame(variable = x$variable[1], n = nrow(x), chisq_p = stats::chisq.test(counts)$p.value,
             mean_rank = mean(x$rank), row.names = NULL)
}))
res$sbc$pass <- res$sbc$chisq_p > 0.005
sf <- unique(sbc[c("rep", "cell", "ok", "max_rhat", "divergences")])
res$sbc_fits <- c(n = nrow(sf), failed = sum(!sf$ok), with_divergences = sum(sf$divergences > 0),
                  rhat_over_1.01 = sum(sf$max_rhat >= 1.01))

# criterion 1, over all cells and over the designs of criterion 2
lf <- res$lowp_fits
in2 <- lf$p11 >= 0.30 & lf$frac_verified >= 0.1
one <- function(name, f) data.frame(set = name, cells = nrow(f), below_0.95 = sum(f$converged < 0.95),
                                    lowest = min(f$converged), errors = sum(f$ok < 1))
res$criterion1 <- rbind(
  one("lowp, all", lf), one("lowp, p11 = 0.30 and verification >= 10%", lf[in2, ]),
  one("calib, fp ~ x", res$calib_fits[res$calib_fits$model == "fp_x", ]),
  one("calib, fp ~ 1 with b_fp = 0", res$calib_fits[res$calib_fits$model == "fp_const" & res$calib_fits$b_fp == 0, ]),
  one("randet, random_det", res$randet_fits[res$randet_fits$model == "random_det", ]),
  one("markov, markov", res$markov_fits[res$markov_fits$model == "markov", ]))

saveRDS(res, file.path("dev", "dyn-val", "results.rds"))

md <- function(x, digits = 3) {
  x[] <- lapply(x, function(v) if (is.numeric(v)) signif(v, digits) else v)
  c(paste("|", paste(names(x), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(trimws(r), collapse = " | "), "|")))
}
c2 <- res$lowp[res$lowp$p11 >= 0.30 & res$lowp$frac_verified >= 0.1, ]
c4 <- res$calib[res$calib$frac_confirmed > 0 & res$calib$b_fp == 0, ]
c5 <- res$calib[res$calib$b_fp == 1, ]
out <- c("# Validation of the dynamic model: results", "",
         "Written by `dev/dyn-val/summarise.R`. Criteria in `README.md`.", "",
         "## Criterion 1: convergence", "", md(res$criterion1), "",
         "Cells of `lowp` with fewer than 95% of fits converged:", "",
         md(lf[lf$converged < 0.95, ]), "",
         "## Criterion 2: p11 = 0.30, verification at least 10%", "", md(c2), "",
         "## Criterion 3: limits (p11 = 0.15 or no verification)", "",
         md(res$lowp[!(res$lowp$p11 >= 0.30 & res$lowp$frac_verified >= 0.1), ]), "",
         "## Criterion 4: observation-confirmation, b_fp = 0", "", md(c4), "",
         "## Criterion 5: false positives that vary, b_fp = 1", "", md(c5), "",
         "All cells of `calib`, converged share:", "", md(res$calib_fits), "",
         "## Criterion 6: SBC", "",
         paste0("Fits: ", res$sbc_fits[["n"]], "; failed: ", res$sbc_fits[["failed"]],
                "; with divergences: ", res$sbc_fits[["with_divergences"]],
                "; with R-hat at or above 1.01: ", res$sbc_fits[["rhat_over_1.01"]], "."), "",
         md(res$sbc), "",
         "## Criterion 7: random_det and markov", "", md(res$randet_fits), "", md(res$randet), "",
         md(res$markov_fits), "", md(res$markov), "")
writeLines(out, file.path("dev", "dyn-val", "results.md"))
message("wrote dev/dyn-val/results.rds and results.md")
