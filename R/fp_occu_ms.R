#' Multispecies false-positive occupancy model
#'
#' Fits the false-positive occupancy model of [fp_occu()] to several species
#' at once, with the coefficients of each species drawn from distributions
#' of the community, so that rare species borrow strength from the others.
#' **This function is in development** (see [crywolf-api]): its arguments,
#' results and output names may change in a minor version, and it lacks the
#' methods of [fp_occu()] other than `print()`, `summary()` and
#' `as_draws_df()`.
#'
#' Each species \eqn{s} has its own occupancy coefficients
#' \eqn{\beta_s = \mu_\beta + \sigma_\beta z_s} (all the terms of `occ`, with
#' \eqn{z_s} standard normal), its own true-detection intercept
#' \eqn{\alpha_s = \alpha + \sigma_\alpha z_s}, and its own gap between the
#' true-detection and false-positive intercepts,
#' \eqn{g_s = \mu_g \exp(\sigma_g z_s)}, which keeps the false-positive rate
#' of every species below its true-detection rate at the reference. The
#' slopes of `det`, `fp` and `ver`, the intercept of the verified method and
#' the recorder effect are shared by the species. The false positives of a
#' species are independent of those of the others: the classifier errs in
#' labelling a species, not in confusing two species.
#'
#' The species are those of the data (`fp_data(species = )`), and a species
#' is never assumed to be absent from the list: there is no data
#' augmentation for species never detected.
#'
#' @param data An `fp_data` object with a species column; see the `species`
#'   argument of [fp_data()].
#' @param occ,det,fp,random_recorder,ver,priors As in [fp_occu()]. The
#'   standard deviations among species have the half-Normal priors of
#'   `fp_priors(species_sd = )`; the community means have the priors of the
#'   intercepts, slopes and gap of [fp_priors()].
#' @param chains,parallel_chains,threads_per_chain,iter_warmup,iter_sampling,adapt_delta,seed,refresh,...
#'   Passed to the sampler, as in [fp_occu()]. Runs of detections
#'   (`markov`) and site effects (`random_fp`, `random_det`) are not
#'   available.
#'
#' @return An object of class `fp_occu_ms`, a list with `fit`, `draws`
#'   (readable names: see Details), `stan_data`, `formulas`, `design`,
#'   `data`, `random_recorder`, `priors`, `sampler`, `diagnostics` and
#'   `call`.
#'
#' @details
#' **Output names.** The community means have the names of [fp_occu()]
#' (`occ_*`, `det_*`, `fp_*`, `fp_gap`, the median gap, and `p11_ref`,
#' `p10_ref` at the community means); the standard deviations among species
#' are `sd_occ_*`, `sd_det_(Intercept)` and `sd_log_fp_gap`; and the values
#' of each species carry its name in brackets: `occ_*[species]`,
#' `det_(Intercept)[species]`, `fp_gap[species]`, `p11_ref[species]` and
#' `p10_ref[species]`.
#' @seealso [fp_data()], [fp_occu()]
#' @export
#' @examples
#' \dontrun{
#' # name every species to keep (`schema.json` lists them)
#' x <- fp_read_sheepdog("occupancy_export",
#'                       species = c("ela_bic", "lep_fus", "phy_cuv"))
#' d <- fp_draws_data(x, n_draws = 1)[[1]]
#' m <- fp_occu_ms(d, fp = ~ 1)
#' m
#' }
fp_occu_ms <- function(data, occ = ~ 1, det = ~ 1, fp = ~ 1,
                       random_recorder = !is.null(data$recorder),
                       ver = ~ 1, priors = fp_priors(),
                       chains = 4, parallel_chains = chains, threads_per_chain = "auto",
                       iter_warmup = 1000, iter_sampling = 1000,
                       adapt_delta = 0.95, seed = NULL, refresh = NULL, ...) {
  if (!inherits(data, "fp_data") || is.null(data$species)) {
    stop("`data` must be an `fp_data` object with a species column; see ",
         "`fp_data(species = )`.", call. = FALSE)
  }
  na <- intersect(names(list(...)), c("markov", "random_fp", "random_det"))
  if (length(na)) {
    stop("`", paste(na, collapse = "`, `"), "` is not available in fp_occu_ms().",
         call. = FALSE)
  }
  no_ver <- is.null(data$verified) || nrow(data$verified) == 0
  no_conf <- is.null(data$confirmed) || all(is.na(data$auto[[data$confirmed[1]]]))
  if (no_ver && no_conf) {
    warning("No verified data: the false-positive rates are identified only by ",
            "the priors. Treat the estimates with caution.", call. = FALSE)
  }
  cores <- .sampler_cores(threads_per_chain, parallel_chains)
  threads_per_chain <- cores$threads_per_chain
  parallel_chains <- cores$parallel_chains
  sd_list <- .build_ms_stan_data(data, occ, det, fp, random_recorder, priors, ver = ver)
  if (is.null(refresh)) {
    refresh <- if (interactive()) max(1, round((iter_warmup + iter_sampling) / 10)) else 0
  }
  mod <- .compile_model("fp_occu_ms")
  run_sampler <- function(..., output_basename = .output_basename()) {
    mod$sample(..., output_basename = output_basename)
  }
  fit <- run_sampler(data = sd_list$stan_data, chains = chains,
                     parallel_chains = parallel_chains,
                     threads_per_chain = threads_per_chain,
                     iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                     adapt_delta = adapt_delta, seed = seed, refresh = refresh, ...)
  fit$draws()
  fit$sampler_diagnostics()
  draws <- .pretty_draws_ms(fit, sd_list)
  diag <- .diagnostics(fit, draws)
  out <- structure(
    list(fit = fit, draws = draws, stan_data = sd_list$stan_data,
         formulas = list(occ = occ, det = det, fp = fp, ver = ver),
         design = sd_list$design, data = data,
         random_recorder = isTRUE(sd_list$stan_data$estimate_re == 1),
         priors = priors,
         sampler = list(chains = chains, parallel_chains = parallel_chains,
                        threads_per_chain = threads_per_chain,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = adapt_delta),
         diagnostics = diag, call = match.call()),
    class = "fp_occu_ms"
  )
  .warn_diagnostics(diag)
  out
}

#' @rdname fp_occu_ms
#' @param x,object An `fp_occu_ms` object.
#' @param digits Number of digits to print.
#' @export
print.fp_occu_ms <- function(x, digits = 3, ...) {
  sd <- x$stan_data
  f <- x$formulas
  cat("Multispecies Bayesian false-positive occupancy model (",
      if (isTRUE(sd$obs_confirm == 1)) "observation-confirmation" else "multiple-methods",
      " design)\n\n", sep = "")
  cat("  occupancy:        ", .fmt_formula(f$occ), "by species\n")
  cat("  true detection:   ", .fmt_formula(f$det), "(intercept by species)",
      if (x$random_recorder) "+ (1 | recorder)" else "", "\n")
  cat("  false positive:   ",
      if (is.null(f$fp)) "none (standard occupancy model)" else
        paste(.fmt_formula(f$fp), "(gap by species)"), "\n")
  if (!is.null(f$ver) && length(all.vars(f$ver))) {
    cat("  verified:         ", .fmt_formula(f$ver), "\n")
  }
  cat(sprintf("  data:              %d species at %d sites, %d automated and %d verified occasions\n\n",
              sd$nspecies, sd$nsite, sum(sd$n2) + sd$N3, sd$N1))
  tab <- summary(x)$table
  print(as.data.frame(tab[!grepl("\\[", tab$variable), ]), digits = digits, row.names = FALSE)
  cat("\n  Values of each species in summary().\n")
  .print_diag(x$diagnostics)
  invisible(x)
}

#' @rdname fp_occu_ms
#' @param prob Width of the equal-tailed credible interval.
#' @return `summary()` returns an object of class `summary.fp_occu`, as for
#'   [fp_occu()].
#' @export
summary.fp_occu_ms <- function(object, prob = 0.95, ...) {
  summary.fp_occu(object, prob = prob)
}

#' @rdname fp_occu_ms
#' @export
as_draws_df.fp_occu_ms <- function(x, ...) {
  x$draws
}

# The data of fp_occu_ms.stan. Each site x species is a unit, numbered
# (site - 1) * nspecies + species; the units form an fp_data whose "sites"
# are the units, so that .build_stan_data() builds the observation model of
# each unit as fp_occu() builds that of each site, and the occupancy design
# of each unit from the covariates of its site.
.build_ms_stan_data <- function(data, occ, det, fp, random_recorder, priors, ver = ~ 1) {
  units <- .ms_units(data)
  sl <- .build_stan_data(units, occ, det, fp, random_recorder, priors, ver = ver)
  nsp <- length(data$species_names)
  sd <- sl$stan_data
  sd[c("estimate_fp_re", "fp_site_sd", "estimate_det_re", "det_site_sd")] <- NULL
  sd$nunit <- sd$nsite
  sd$nsite <- nrow(data$sites)
  sd$nspecies <- nsp
  sd$sp_unit <- rep(seq_len(nsp), sd$nsite)
  # priors saved before this one existed
  sd$species_sd <- if (is.null(priors$species_sd)) 1 else priors$species_sd
  sl$stan_data <- sd
  sl$design$species <- data$species_names
  sl
}

# the site x species units of a multispecies data set, as an fp_data with
# one "site" per unit (column .unit); the table of units carries the site
# covariates and the species
.ms_units <- function(data) {
  s <- data$site
  sites <- data$sites
  sp <- data$species_names
  nsp <- length(sp)
  unit_of <- function(df) {
    (match(df[[s]], sites[[s]]) - 1L) * nsp + match(as.character(df[[data$species]]), sp)
  }
  data$auto$.unit <- unit_of(data$auto)
  if (!is.null(data$verified)) data$verified$.unit <- unit_of(data$verified)
  usites <- sites[rep(seq_len(nrow(sites)), each = nsp), , drop = FALSE]
  rownames(usites) <- NULL
  usites$.unit <- seq_len(nrow(usites))
  usites[[data$species]] <- rep(sp, nrow(sites))
  data$sites <- usites
  data$site <- ".unit"
  # sites x species, unit by unit
  data$known_present <- as.vector(t(data$known_present))
  data
}

# readable draws: the names of fp_occu() for the community means, sd_* for
# the standard deviations among species and [species] for the values of
# each species
.pretty_draws_ms <- function(fit, sd_list) {
  sd <- sd_list$stan_data
  des <- sd_list$design
  map <- .name_map(des, sd$estimate_fp == 1, sd$estimate_re == 1,
                   has_verified = sd$N1 > 0)
  map <- map[names(map) != "p10_bar"]
  names(map) <- sub("^beta\\[", "mu_beta[", names(map))
  names(map)[names(map) == "fp_gap"] <- "mu_gap"
  sp <- des$species
  nsp <- length(sp)
  occ <- des$occ_names
  by_sp <- function(stan, name) {
    stats::setNames(sprintf("%s[%s]", name, sp), sprintf("%s[%d]", stan, seq_len(nsp)))
  }
  map <- c(map,
           stats::setNames(paste0("sd_occ_", occ), sprintf("sigma_beta[%d]", seq_along(occ))),
           stats::setNames(sprintf("occ_%s[%s]", rep(occ, nsp), rep(sp, each = length(occ))),
                           sprintf("beta_sp[%d,%d]", rep(seq_along(occ), nsp),
                                   rep(seq_len(nsp), each = length(occ)))),
           c(sigma_alpha = "sd_det_(Intercept)"),
           by_sp("alpha_sp", "det_(Intercept)"),
           by_sp("p11_sp", "p11_ref"),
           if (sd$estimate_fp == 1) {
             c(c(sigma_gap = "sd_log_fp_gap"), by_sp("gap_sp", "fp_gap"),
               by_sp("p10_sp", "p10_ref"))
           })
  dr <- posterior::as_draws_df(fit$draws(variables = unique(sub("\\[.*", "", names(map)))))
  dr <- posterior::subset_draws(dr, variable = names(map))
  posterior::variables(dr) <- unname(map[posterior::variables(dr)])
  dr
}
