# Reads the job files of dev/dyn-kfold/run_kfold.R and writes results.rds
# and the tables of results.md: the held-out elpd of the dynamic model (dyn)
# and of the model of independent site-nights (ind), by site and fold, and
# their difference. Run from the package root: Rscript dev/dyn-kfold/summarise.R
dir <- file.path("dev", "dyn-kfold")
jobs <- lapply(list.files(file.path(dir, "jobs"), full.names = TRUE), readRDS)
lme <- function(x) max(x) + log(mean(exp(x - max(x))))

# one row per model, fold and site with held-out minutes
rows <- do.call(rbind, lapply(jobs, function(j) {
  data.frame(model = j$model, fold = j$fold, site = j$sites, n = j$n_heldout,
             elpd = apply(j$ll, 2, lme))
}))
rows <- rows[rows$n > 0, ]
w <- reshape(rows, idvar = c("fold", "site", "n"), timevar = "model", direction = "wide")
w <- w[stats::complete.cases(w), ]
w$diff <- w$elpd.dyn - w$elpd.ind
se <- function(x) sqrt(length(x)) * stats::sd(x)
total <- data.frame(model = c("dyn", "ind", "dyn - ind"),
                    elpd = c(sum(w$elpd.dyn), sum(w$elpd.ind), sum(w$diff)),
                    se = c(se(w$elpd.dyn), se(w$elpd.ind), se(w$diff)))
by_site <- stats::aggregate(cbind(n, elpd.dyn, elpd.ind, diff) ~ site, w, sum)
by_fold <- stats::aggregate(cbind(n, elpd.dyn, elpd.ind, diff) ~ fold, w, sum)
diag <- do.call(rbind, lapply(jobs, function(j) {
  d <- j$diagnostics
  data.frame(model = j$model, fold = j$fold, max_rhat = d$max_rhat,
             min_ess_bulk = d$min_ess_bulk, divergences = d$divergences,
             minutes = round(j$minutes, 1))
}))
diag <- diag[order(diag$model, diag$fold), ]
coefs <- do.call(rbind, lapply(jobs, function(j) {
  s <- j$summary
  data.frame(model = j$model, fold = j$fold, variable = s$variable, mean = s$mean)
}))
res <- list(units = w, total = total, by_site = by_site, by_fold = by_fold,
            diagnostics = diag, coefs = coefs)
saveRDS(res, file.path(dir, "results.rds"))

md <- function(df, digits = 1) {
  df[] <- lapply(df, function(v) if (is.numeric(v)) round(v, digits) else v)
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
writeLines(c(
  "# Blocked cross-validation of the dynamic model on curiango (step 8)",
  "",
  "Written by `summarise.R` from the job files of `run_kfold.R` (design in its header).",
  sprintf("Held-out elpd over %d site x fold units; SE from the spread of the units.", nrow(w)),
  "",
  "## Total", "", md(total), "",
  "## By site", "", md(by_site), "",
  "## By fold", "", md(by_fold), "",
  "## Convergence of the training fits", "", md(diag, 3)
), file.path(dir, "results.md"))
print(total)
