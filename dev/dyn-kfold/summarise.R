# Reads the job files of dev/dyn-kfold/run_kfold.R and writes results.rds
# and the tables of results.md: the held-out elpd of each model, by site and
# fold, and its difference from the dynamic model (dyn). Run from the package
# root: Rscript dev/dyn-kfold/summarise.R
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
models <- c("dyn", "ind", "dyn_mk", "ind_mk", "ind_rd")
cols <- paste0("elpd.", models)
se <- function(x) sqrt(length(x)) * stats::sd(x)
total <- do.call(rbind, lapply(models, function(m) {
  e <- w[[paste0("elpd.", m)]]
  d <- e - w$elpd.dyn
  data.frame(model = m, elpd = sum(e), se = se(e), diff_from_dyn = sum(d), se_diff = se(d))
}))
pairs <- do.call(rbind, lapply(list(c("dyn", "ind"), c("dyn_mk", "ind_mk"), c("dyn_mk", "dyn"),
                                    c("dyn_mk", "ind_rd"), c("ind_mk", "ind"), c("ind_rd", "ind")),
                               function(p) {
  d <- w[[paste0("elpd.", p[1])]] - w[[paste0("elpd.", p[2])]]
  data.frame(comparison = paste(p, collapse = " - "), diff = sum(d), se = se(d))
}))
by_unit <- function(by) stats::aggregate(w[c("n", cols)], w[by], sum)
by_site <- by_unit("site")
by_fold <- by_unit("fold")
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
res <- list(units = w, total = total, pairs = pairs, by_site = by_site, by_fold = by_fold,
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
  sprintf(paste("Held-out elpd over %d site x fold units with held-out minutes; SE from",
                "the spread of the units, which treats them as independent although",
                "the units of a fold share a training fit."), nrow(w)),
  "",
  "## Total", "", md(total), "",
  "## Pairs of models", "", md(pairs), "",
  "## By site", "", md(by_site), "",
  "## By fold", "", md(by_fold), "",
  "## Convergence of the training fits", "", md(diag, 3)
), file.path(dir, "results.md"))
print(total)
