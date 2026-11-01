# =============================================================================
# dev/dyn-kfold/run_kfold.R
# Step 8 of the dynamic model (dev/roadmap-v2.md): blocked cross-validation
# on the curiango data of the dynamic model, with the nights of each site as
# seasons (nightly use), against the model of independent site-nights.
# LOO is not an option: the log_lik of fp_occu_dyn() is one value per site,
# and curiango has seven sites.
#
# Data: custom classifier at threshold 0.40, one-minute occasions, verified
# minutes, recorder = site. Covariates standardized over the data: bio and
# aci per minute; temp and rh per site-night (on the minutes of the night);
# temp_d and rh_d, the means of temp and rh over the sites recorded on each
# date, as season covariates (every date from 2023-11-25 to 2025-02-26 has at
# least one site, so they exist for the transitions without data too).
#
# Models, both with random_recorder, det ~ bio + aci and fp ~ temp + rh:
#   dyn  fp_occu_dyn(psi1 = ~ 1, col = ~ temp_d + rh_d, ext = ~ temp_d + rh_d)
#   ind  fp_occu(occ = ~ temp_d + rh_d), sites = site-nights
#
# Folds: the dates are cut into blocks of 7 nights (weeks), and the weeks are
# dealt at random into K = 10 folds, the same for every site. The first and
# the last week stay in every training set, so that the training data span
# the same dates as the full data and fp_data() fills the same nightly grid.
# A fold leaves out all the minutes and verified minutes of its weeks.
#
# Held-out log density, for each site and fold, on every posterior draw of
# the training fit: for dyn, log p(y_full | theta) - log p(y_train | theta),
# from the forward algorithm (.dyn_loglik(), equal to the Stan log_lik), with
# the left-out site-nights as seasons without data in the training fit, so
# the states of those nights are predicted from the neighbouring nights; for
# ind, the sum of the log-likelihoods of the site's left-out site-nights.
# The elpd of a site and fold is the log of the mean over draws of
# exp(log density); summarise.R adds and compares them.
#
# Every fit: 4 chains of 1000 + 1000 iterations, adapt_delta = 0.95, seed
# 91000 + task. Each task is saved in dev/dyn-kfold/jobs/ and skipped on a
# rerun. Run from the package root, with crywolf installed from this
# checkout:
#   WORKERS=20 THREADS=1 Rscript dev/dyn-kfold/run_kfold.R
# TASKS (e.g. "dyn-1") runs only some of the 20 tasks, for a pilot. Do not
# edit this file while it runs: Rscript reads it as it goes.
# =============================================================================

library(crywolf)

workers <- as.integer(Sys.getenv("WORKERS", "20"))
threads <- as.integer(Sys.getenv("THREADS", "1"))
only <- Sys.getenv("TASKS", "")
job_dir <- file.path("dev", "dyn-kfold", "jobs")
dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)
K <- 10
block <- 7

zs <- function(v) as.numeric(scale(v))
m <- curiango$minutes
nt <- curiango$nights
daily <- stats::aggregate(cbind(temp, rh) ~ night, nt, mean)
daily <- data.frame(night = daily$night, temp_d = zs(daily$temp), rh_d = zs(daily$rh))
i <- match(m$site_night, nt$site_night)
occasions <- data.frame(file = paste(m$site_night, format(m$time, "%Y%m%d_%H%M")),
                        site = m$site, site_night = m$site_night, night = nt$night[i],
                        time = m$time, bio = zs(m$bio), aci = zs(m$aci),
                        temp = zs(nt$temp[i]), rh = zs(nt$rh[i]), recorder = m$site,
                        y = as.integer(m$conf_custom >= 0.40))
verified <- data.frame(site = curiango$verified$site,
                       site_night = curiango$verified$site_night,
                       night = nt$night[match(curiango$verified$site_night, nt$site_night)],
                       y = curiango$verified$y)
stopifnot(!anyNA(occasions$night), !anyNA(verified$night))
site_nights <- merge(nt[c("site_night", "site", "night")], daily, by = "night")

nights <- sort(unique(nt$night))
week <- as.integer(nights - nights[1]) %/% block + 1L
inner <- seq_len(max(week))[-c(1, max(week))]
set.seed(20261101)
fold_of_week <- integer(max(week))
fold_of_week[inner] <- sample(rep_len(seq_len(K), length(inner)))
fold_of <- function(night) fold_of_week[week[match(night, nights)]]

forms <- list(det = ~ bio + aci, fp = ~ temp + rh, occ = ~ temp_d + rh_d,
              psi1 = ~ 1, col = ~ temp_d + rh_d, ext = ~ temp_d + rh_d)

# the data without fold k, or all the data (k = NULL)
make_data <- function(model, k = NULL) {
  keep <- function(night) if (is.null(k)) rep(TRUE, length(night)) else fold_of(night) != k
  keep_occ <- keep(occasions$night)
  keep_ver <- keep(verified$night)
  suppressWarnings(
    if (model == "dyn") {
      fp_data(occasions[keep_occ, ], data.frame(site = sort(unique(occasions$site))),
              verified[keep_ver, ], site = "site", recorder = "recorder",
              season = "night", season_covs = daily)
    } else {
      fp_data(occasions[keep_occ, ], site_nights[keep(site_nights$night), ],
              verified[keep_ver, ],
              site = "site_night", recorder = "recorder")
    }
  )
}

fit_model <- function(model, d, seed) {
  args <- list(d, det = forms$det, fp = forms$fp, random_recorder = TRUE, seed = seed,
               chains = 4, parallel_chains = 4, threads_per_chain = threads,
               iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
               refresh = 0, show_messages = FALSE)
  if (model == "dyn") {
    do.call(fp_occu_dyn, c(args, list(psi1 = forms$psi1, col = forms$col, ext = forms$ext)))
  } else {
    do.call(fp_occu, c(args, list(occ = forms$occ)))
  }
}

# apply f to the draws of pars in chunks of n, binding the S x nsite results
by_draws <- function(pars, f, n = 200) {
  S <- nrow(pars$beta)
  chunks <- split(seq_len(S), ceiling(seq_len(S) / n))
  out <- lapply(chunks, function(k) {
    f(lapply(pars, function(p) {
      if (is.null(p)) NULL else if (is.null(dim(p))) p[k] else p[k, , drop = FALSE]
    }))
  })
  do.call(rbind, out)
}

# held-out log density, S draws x 7 sites (in the order of the dyn sites)
heldout_ll <- function(model, fit, d_full, k) {
  sites <- sort(unique(occasions$site))
  if (model == "dyn") {
    sd_tr <- fit$stan_data
    sd_full <- crywolf:::.build_dyn_stan_data(d_full, forms$psi1, forms$col, forms$ext,
                                              forms$det, forms$fp, TRUE, fit$priors)$stan_data
    stopifnot(identical(fit$data$seasons, d_full$seasons),
              identical(fit$data$sites$site, sites),
              identical(sd_tr$X, sd_full$X), identical(sd_tr$Xcol, sd_full$Xcol),
              identical(sd_tr$Xext, sd_full$Xext), identical(sd_tr$nrec, sd_full$nrec))
    pars <- crywolf:::.dyn_pars(fit$fit, sd_tr)
    by_draws(pars, function(p) {
      crywolf:::.dyn_loglik(sd_full, p) - crywolf:::.dyn_loglik(sd_tr, p)
    })
  } else {
    sd_full <- crywolf:::.build_stan_data(d_full, forms$occ, forms$det, forms$fp, TRUE,
                                          fit$priors)$stan_data
    stopifnot(identical(fit$stan_data$nrec, sd_full$nrec))
    pars <- crywolf:::.ppc_pars(fit$fit, sd_full, NULL)
    ll <- by_draws(pars, function(p) {
      crywolf:::.ppc_core(sd_full, p, "automated", "freeman-tukey", ll_only = TRUE)$ll
    })
    sn <- d_full$sites
    held <- fold_of(sn$night) == k
    out <- sapply(sites, function(s) rowSums(ll[, held & sn$site == s, drop = FALSE]))
    unname(out)
  }
}

run_task <- function(model, k, seed) {
  file <- file.path(job_dir, sprintf("%s-%02d.rds", model, k))
  if (file.exists(file)) return(invisible(file))
  t0 <- Sys.time()
  d_full <- make_data(model)
  fit <- fit_model(model, make_data(model, k), seed)
  minutes_fit <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ll <- heldout_ll(model, fit, d_full, k)
  held <- fold_of(occasions$night) == k
  out <- list(model = model, fold = k, seed = seed, ll = ll,
              sites = sort(unique(occasions$site)),
              n_heldout = c(table(factor(occasions$site[held],
                                         levels = sort(unique(occasions$site))))),
              summary = summary(fit)$table, diagnostics = fit$diagnostics,
              minutes_fit = minutes_fit,
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, file)
  invisible(file)
}

tasks <- expand.grid(fold = seq_len(K), model = c("dyn", "ind"), stringsAsFactors = FALSE)
tasks$seed <- 91000 + seq_len(nrow(tasks))
tasks$id <- paste(tasks$model, tasks$fold, sep = "-")
if (nzchar(only)) tasks <- tasks[tasks$id %in% strsplit(only, ",")[[1]], ]
message(nrow(tasks), " task(s), ", workers, " worker(s) x 4 chains x ", threads, " thread(s)")
future::plan(future::multisession, workers = min(workers, nrow(tasks)))
res <- furrr::future_pmap(tasks[c("model", "fold", "seed")],
                          function(model, fold, seed) {
                            tryCatch(run_task(model, fold, seed),
                                     error = function(e) conditionMessage(e))
                          },
                          .options = furrr::furrr_options(seed = NULL))
bad <- !vapply(res, function(r) is.character(r) && file.exists(r), logical(1))
if (any(bad)) {
  message("Failed: ", paste(tasks$id[bad], unlist(res[bad]), sep = ": ", collapse = "\n"))
  quit(status = 1)
}
message("done")
