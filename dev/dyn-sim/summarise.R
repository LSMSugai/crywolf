# =============================================================================
# dev/dyn-sim/summarise.R
# Reads the job files of run_sims.R (dev/dyn-sim/jobs/) and writes
# dev/dyn-sim/results.rds (the summaries below and the SBC ranks) and the
# tables of dev/dyn-sim/results.md. Colonization, extinction and initial
# occupancy are summarised on the probability scale, from the logit-scale
# intercepts of models without covariates (ident); the realized occupancy
# psi_fs is pooled over seasons. Fits count as converged with R-hat < 1.01
# for every parameter with a true value and no divergent transitions.
# Run from the package root:  Rscript dev/dyn-sim/summarise.R
# =============================================================================

job_dir <- file.path("dev", "dyn-sim", "jobs")
read_part <- function(part) {
  f <- list.files(job_dir, pattern = paste0("^", part, "-"), full.names = TRUE)
  if (!length(f)) return(NULL)
  do.call(rbind, lapply(f, readRDS))
}
ident <- read_part("ident")
recovery <- read_part("recovery")
sbc <- read_part("sbc")

# one row per parameter (and grouping), with bias, RMSE, 95% coverage and width
metrics <- function(d, by) {
  d$param <- sub("\\[\\d+\\]$", "", d$variable)
  inv <- d$param %in% c("psi1_(Intercept)", "col_(Intercept)", "ext_(Intercept)")
  d$param[inv] <- sub("_\\(Intercept\\)$", "", d$param[inv])
  for (v in c("mean", "lower", "upper", "truth")) d[[v]][inv] <- stats::plogis(d[[v]][inv])
  key <- interaction(d[c(by, "param")], drop = TRUE, sep = "|")
  out <- lapply(split(d, key), function(x) {
    data.frame(x[1, c(by, "param")], n = nrow(x),
               bias = mean(x$mean - x$truth), rmse = sqrt(mean((x$mean - x$truth)^2)),
               coverage = mean(x$lower <= x$truth & x$truth <= x$upper),
               width = mean(x$upper - x$lower), row.names = NULL)
  })
  do.call(rbind, out)
}
fit_keys <- c("part", "rep", "cell", "model")
fit_ok <- function(d) {
  bad <- !d$ok | (!is.na(d$rhat) & d$rhat >= 1.01) | (!is.na(d$divergences) & d$divergences > 0)
  key <- interaction(d[intersect(fit_keys, names(d))], drop = TRUE)
  !ave(bad, key, FUN = any)
}

res <- list()
if (!is.null(ident)) {
  ident$converged <- fit_ok(ident)
  res$ident_fits <- stats::aggregate(converged ~ frac_verified + p10 + nseason + model,
                                     data = unique(ident[c(fit_keys, "frac_verified", "p10",
                                                           "nseason", "converged")]),
                                     FUN = mean)
  res$ident <- metrics(ident[ident$ok, ], c("frac_verified", "p10", "nseason", "model"))
}
if (!is.null(recovery)) {
  recovery$converged <- fit_ok(recovery)
  res$recovery_fits <- stats::aggregate(converged ~ nsite + nseason,
                                        data = unique(recovery[c(fit_keys, "nsite", "nseason",
                                                                 "converged")]), FUN = mean)
  res$recovery <- metrics(recovery[recovery$ok, ], c("nsite", "nseason"))
}
if (!is.null(sbc)) {
  ok <- sbc[sbc$ok, ]
  # 100 possible ranks (0 to 99) in 10 bins of 10; chi-squared test of uniformity
  res$sbc_ranks <- ok[c("rep", "variable", "truth", "rank", "max_rhat", "divergences")]
  res$sbc <- do.call(rbind, lapply(split(ok, ok$variable), function(x) {
    counts <- tabulate(floor(x$rank / 10) + 1, nbins = 10)
    data.frame(variable = x$variable[1], n = nrow(x),
               chisq_p = stats::chisq.test(counts)$p.value,
               mean_rank = mean(x$rank), row.names = NULL)
  }))
  res$sbc_fits <- c(n = length(unique(sbc$rep)), failed = length(unique(sbc$rep[!sbc$ok])),
                    with_divergences = length(unique(ok$rep[ok$divergences > 0])))
}
saveRDS(res, file.path("dev", "dyn-sim", "results.rds"))

md <- function(x, digits = 3) {
  x[] <- lapply(x, function(v) if (is.numeric(v)) signif(v, digits) else v)
  c(paste("|", paste(names(x), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
out <- c("# Simulations of the dynamic model: results", "",
         "Written by `dev/dyn-sim/summarise.R`.", "")
if (!is.null(res$ident)) {
  keep <- res$ident$param %in% c("col", "ext", "psi_fs", "p10_bar")
  out <- c(out, "## Identification (step 6)", "", "Share of converged fits:", "",
           md(res$ident_fits), "", "Colonization, extinction, realized occupancy:", "",
           md(res$ident[keep, ]), "")
}
if (!is.null(res$recovery)) {
  out <- c(out, "## Recovery with covariates (step 7)", "", md(res$recovery_fits), "",
           md(res$recovery), "")
}
if (!is.null(res$sbc)) {
  out <- c(out, "## Simulation-based calibration (step 7)", "",
           paste0("Fits: ", res$sbc_fits[["n"]], "; failed: ", res$sbc_fits[["failed"]],
                  "; with divergent transitions: ", res$sbc_fits[["with_divergences"]], "."),
           "", md(res$sbc), "")
}
writeLines(out, file.path("dev", "dyn-sim", "results.md"))
message("wrote dev/dyn-sim/results.rds and results.md")
