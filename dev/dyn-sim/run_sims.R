# =============================================================================
# dev/dyn-sim/run_sims.R
# Steps 6 and 7 of the dynamic model (dev/roadmap-v2.md): simulations for
# fp_occu_dyn(), all from simulate_fp_dyn() (100 or 200 sites, 4 automated
# occasions per site and season, verified site-seasons with 2 occasions).
#
#   ident    (step 6) how much verification separates colonization and
#            extinction from false positives. 45 cells: verified fraction of
#            site-seasons in {0, 0.05, 0.10, 0.20, 0.40} x p10 in {0.02,
#            0.05, 0.10} x seasons in {3, 6, 10}; 100 sites, psi1 = 0.5,
#            col = ext = 0.2, p11 = 0.6, r11 = 0.8, no covariates. Two models
#            per dataset: fp ~ 1 and fp = NULL (a dynamic occupancy model
#            that ignores false positives). 50 replicates per cell.
#   recovery (step 7) recovery with covariates: psi1 ~ x (b = 0.5), col ~ rain
#            (b = 1), ext ~ x (b = -1), p10 = 0.05, 20% of site-seasons
#            verified; 4 designs, 100 or 200 sites x 5 or 10 seasons; 300
#            replicates per design.
#   sbc      (step 7) simulation-based calibration: the parameters of
#            psi1 ~ x, col ~ rain, ext ~ x, det ~ 1, fp ~ 1, ver ~ 1 drawn
#            from fp_priors(), 100 sites x 5 seasons, 20% verified; the rank
#            of each true value among 99 thinned posterior draws. 1000
#            replicates.
#
# Every fit: 4 chains of 1000 + 1000 iterations, adapt_delta = 0.95, one
# core (threads_per_chain = 1). Seeds: 61001 on (ident), 71001 on
# (recovery), 81001 on (sbc), one per task, so a task gives the same data
# and fit on a rerun.
#
# Each task is saved in dev/dyn-sim/jobs/ as soon as it ends and skipped on
# a rerun, so an interrupted run resumes where it stopped. Run from the
# package root, with crywolf installed from this checkout:
#   WORKERS=16 Rscript dev/dyn-sim/run_sims.R
# REPS scales the replicates (REPS = 1: one per cell or design, for a test);
# PARTS (default "ident,recovery,sbc") selects the parts. Do not edit this
# file while it runs: Rscript reads it as it goes.
# =============================================================================

library(crywolf)

reps <- as.integer(Sys.getenv("REPS", "0"))   # 0: the replicates above
workers <- as.integer(Sys.getenv("WORKERS", "16"))
parts <- strsplit(Sys.getenv("PARTS", "ident,recovery,sbc"), ",")[[1]]
job_dir <- file.path("dev", "dyn-sim", "jobs")
dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)
n_rep <- function(full) if (reps > 0) reps else full

ident_cells <- expand.grid(frac_verified = c(0, 0.05, 0.10, 0.20, 0.40),
                           p10 = c(0.02, 0.05, 0.10), nseason = c(3, 6, 10))
recovery_cells <- expand.grid(nsite = c(100, 200), nseason = c(5, 10))

tasks <- rbind(
  if ("ident" %in% parts) {
    t <- expand.grid(rep = seq_len(n_rep(50)), cell = seq_len(nrow(ident_cells)))
    data.frame(part = "ident", t, seed = 61000 + seq_len(nrow(t)))
  },
  if ("recovery" %in% parts) {
    t <- expand.grid(rep = seq_len(n_rep(300)), cell = seq_len(nrow(recovery_cells)))
    data.frame(part = "recovery", t, seed = 71000 + seq_len(nrow(t)))
  },
  if ("sbc" %in% parts) {
    t <- data.frame(rep = seq_len(n_rep(1000)), cell = 1L)
    data.frame(part = "sbc", t, seed = 81000 + seq_len(nrow(t)))
  }
)
# longest tasks first, so that the last ones to end are short
tasks <- tasks[order(-ifelse(tasks$part == "ident",
                             ident_cells$nseason[pmin(tasks$cell, nrow(ident_cells))],
                             ifelse(tasks$part == "recovery", 10, 5))), ]

fit_dyn <- function(d, seed, ...) {
  suppressWarnings(fp_occu_dyn(
    d, ..., random_recorder = FALSE, chains = 4, parallel_chains = 1,
    threads_per_chain = 1, iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.95, seed = seed, refresh = 0, show_messages = FALSE))
}

# summary rows of the parameters with a true value, plus psi_fs[t] against
# the realized occupancy of season t
summary_rows <- function(m, truth_pars, z) {
  s <- summary(m)$table
  real <- stats::setNames(colMeans(z), sprintf("psi_fs[%d]", seq_len(ncol(z))))
  truth <- c(truth_pars, real)
  s <- s[s$variable %in% names(truth), c("variable", "mean", "lower", "upper",
                                         "rhat", "ess_bulk")]
  s$truth <- unname(truth[s$variable])
  as.data.frame(s)
}

# one model fitted to one dataset -> rows, or one row with the error
fit_rows <- function(model, expr) {
  t0 <- Sys.time()
  res <- tryCatch(expr, error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "error")) {
    return(data.frame(model = model, ok = FALSE, msg = conditionMessage(res),
                      variable = NA_character_, mean = NA_real_, lower = NA_real_,
                      upper = NA_real_, rhat = NA_real_, ess_bulk = NA_real_,
                      truth = NA_real_, divergences = NA_integer_, seconds = secs))
  }
  data.frame(model = model, ok = TRUE, msg = NA_character_, res$rows,
             divergences = res$m$diagnostics$divergences, seconds = secs)
}

run_ident <- function(cell, seed) {
  cfg <- ident_cells[cell, ]
  s <- simulate_fp_dyn(nsite = 100, nseason = cfg$nseason, nocc = 4, psi1 = 0.5,
                       col = 0.2, ext = 0.2, p11 = 0.6, p10 = cfg$p10, r11 = 0.8,
                       frac_verified = cfg$frac_verified, nver = 2, seed = seed)
  keep <- c("psi1_(Intercept)", "col_(Intercept)", "ext_(Intercept)",
            "det_(Intercept)", "fp_(Intercept)")
  truth <- c(s$truth$pars[keep], p10_bar = cfg$p10)
  rows <- rbind(
    fit_rows("fp", {
      m <- fit_dyn(s$data, seed, fp = ~ 1)
      list(m = m, rows = summary_rows(m, truth, s$truth$z))
    }),
    fit_rows("no_fp", {
      m <- fit_dyn(s$data, seed, fp = NULL)
      list(m = m, rows = summary_rows(m, truth, s$truth$z))
    })
  )
  cbind(cfg, n_verified = length(unique(paste(s$data$verified$site,
                                              s$data$verified$season))), rows)
}

run_recovery <- function(cell, seed) {
  cfg <- recovery_cells[cell, ]
  s <- simulate_fp_dyn(nsite = cfg$nsite, nseason = cfg$nseason, nocc = 4,
                       psi1 = 0.5, b_psi1 = 0.5, col = 0.2, b_col = 1, ext = 0.2,
                       b_ext = -1, p11 = 0.6, p10 = 0.05, r11 = 0.8,
                       frac_verified = 0.2, nver = 2, seed = seed)
  rows <- fit_rows("fp", {
    m <- fit_dyn(s$data, seed, psi1 = ~ x, col = ~ rain, ext = ~ x, fp = ~ 1)
    list(m = m, rows = summary_rows(m, s$truth$pars, s$truth$z))
  })
  cbind(cfg, rows)
}

# parameters drawn from fp_priors(): intercepts Normal(0, 1.5), slopes
# Normal(0, 2.5), fp_gap Normal(2, 1.5) truncated at 0 (gamma0 = alpha -
# fp_gap), verified intercept Normal(0, 1.5)
run_sbc <- function(seed) {
  pr <- fp_priors()
  set.seed(seed)
  int <- function() stats::rnorm(1, 0, pr$intercept_sd)
  slope <- function() stats::rnorm(1, 0, pr$slope_sd)
  repeat {
    gap <- stats::rnorm(1, pr$gap_mean, pr$gap_sd)
    if (gap > 0) break
  }
  b <- list(psi1 = int(), b_psi1 = slope(), col = int(), b_col = slope(),
            ext = int(), b_ext = slope(), alpha = int(), delta0 = int())
  s <- simulate_fp_dyn(nsite = 100, nseason = 5, nocc = 4,
                       psi1 = stats::plogis(b$psi1), b_psi1 = b$b_psi1,
                       col = stats::plogis(b$col), b_col = b$b_col,
                       ext = stats::plogis(b$ext), b_ext = b$b_ext,
                       p11 = stats::plogis(b$alpha),
                       p10 = stats::plogis(b$alpha - gap),
                       r11 = stats::plogis(b$delta0), frac_verified = 0.2, nver = 2,
                       seed = seed)
  truth <- c(s$truth$pars, fp_gap = gap)
  t0 <- Sys.time()
  res <- tryCatch({
    m <- fit_dyn(s$data, seed, psi1 = ~ x, col = ~ rain, ext = ~ x, fp = ~ 1)
    dr <- posterior::as_draws_matrix(as_draws_df(m))
    thin <- dr[round(seq(1, posterior::ndraws(dr), length.out = 99)), names(truth)]
    list(ranks = colSums(sweep(unclass(thin), 2, truth, "<")),
         diag = m$diagnostics)
  }, error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !inherits(res, "error")
  data.frame(variable = names(truth), truth = unname(truth),
             rank = if (ok) unname(res$ranks[names(truth)]) else NA_real_,
             ok = ok, msg = if (ok) NA_character_ else conditionMessage(res),
             max_rhat = if (ok) res$diag$max_rhat else NA_real_,
             divergences = if (ok) res$diag$divergences else NA_integer_,
             seconds = secs)
}

run_one <- function(part, rep, cell, seed) {
  f <- file.path(job_dir, sprintf("%s-cell%02d-rep%04d.rds", part, cell, rep))
  if (file.exists(f)) return(invisible(NULL))
  out <- switch(part,
    ident = run_ident(cell, seed),
    recovery = run_recovery(cell, seed),
    sbc = run_sbc(seed))
  saveRDS(cbind(part = part, rep = rep, cell = cell, seed = seed, out), f)
  invisible(NULL)
}

crywolf::fp_compile()                          # once, before the workers
invisible(crywolf:::.compile_model("fp_occu_dyn"))
future::plan(future::multisession, workers = workers)
t0 <- Sys.time()
furrr::future_pwalk(tasks, run_one,
                    .options = furrr::furrr_options(seed = TRUE, chunk_size = 1))
message("done: ", nrow(tasks), " tasks in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 2), " h")
