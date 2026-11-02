#' Simulation-based calibration
#'
#' Checks that the sampler recovers the posterior correctly for a given
#' design and set of priors (Talts et al. 2018; Modrák et al. 2025). Each
#' replicate draws parameter values from the priors, simulates data with the
#' design of `data` (as [fp_prior_predictive()] does), fits the model, and
#' records the rank of the true value among the posterior draws. If the
#' computation is correct, the ranks are uniform. Besides the parameters, it
#' ranks the log-likelihood of the simulated data at the true values among
#' its values at the posterior draws (`log_lik`), a test quantity that
#' depends on the data and is sensitive to problems the parameter ranks
#' can miss (Modrák et al. 2025).
#'
#' @param data An `fp_data` object that gives the design.
#' @param occ,det,fp,ver,markov,random_recorder,random_fp,random_det,priors As in
#'   [fp_occu()]. The
#'   priors are both the ones the data are drawn from and the ones the model
#'   is fitted with, as the method requires.
#' @param reps Number of replicates.
#' @param ndraws Number of posterior draws, evenly thinned, among which the
#'   rank is computed; thinning reduces autocorrelation. It must not exceed
#'   the number of draws of each fit, and `ndraws + 1` should be a multiple
#'   of the number of bins used by [base::summary()] and [plot()] (10 by
#'   default).
#' @param base_seed,parallel,fit_args As in [fp_calibrate()].
#'
#' @details
#' This differs from [fp_calibrate()]: there, data come from fixed parameter
#' values and the question is whether a model's intervals cover them, which
#' tests the model against the data-generating process (including
#' misspecification). Here, data come from the priors of the model itself,
#' so the question is only whether the computation is right: a sampler
#' problem, an error in the model code, or a prior that leads to extreme
#' data shows up as ranks that are not uniform. A U-shaped histogram means
#' posteriors that are too narrow; a hump, too wide; a slope, a bias.
#'
#' Prior draws can produce extreme data (for instance, no detection at
#' all); such replicates are fitted as they are, and their divergent
#' transitions are recorded. The `log_lik` quantity is not computed with
#' `random_fp = TRUE` or `random_det = TRUE`, whose site effects are not
#' kept by the simulation.
#' [base::summary()] uses a chi-squared test of the binned ranks, which is
#' coarse; look at the histograms from [plot()] as well.
#'
#' @return A data frame of class `fp_sbc` with one row per replicate and
#'   parameter: `rep`, `parameter`, `truth`, `rank` (0 to `ndraws`),
#'   `divergences`, and `ok`. [base::summary()] tests the uniformity of the ranks
#'   of each parameter, and [plot()] draws their histograms with the band
#'   expected under uniformity.
#' @references
#' Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
#' Validating Bayesian inference algorithms with simulation-based
#' calibration. arXiv:1804.06788.
#'
#' Modrák, M., Moon, A. H., Kim, S., Bürkner, P., Huurre, N., Faltejsková,
#' K., Gelman, A., & Vehtari, A. (2025). Simulation-based calibration
#' checking for Bayesian computation: the choice of test quantities shapes
#' sensitivity. *Bayesian Analysis*, 20, 461–488. \doi{10.1214/23-BA1404}
#' @seealso [fp_calibrate()], [fp_prior_predictive()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 100, nocc = 10, seed = 1)
#' sbc <- fp_sbc(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
#'               reps = 100, parallel = TRUE)
#' summary(sbc)
#' plot(sbc)
#' }
fp_sbc <- function(data, occ = ~ 1, det = ~ 1, fp = ~ 1, ver = ~ 1, markov = NULL,
                   random_recorder = !is.null(data$recorder),
                   random_fp = FALSE, random_det = FALSE, priors = fp_priors(),
                   reps = 100,
                   ndraws = 99, base_seed = 1L, parallel = FALSE,
                   fit_args = list(chains = 4, parallel_chains = 1,
                                   iter_warmup = 1000, iter_sampling = 1000,
                                   adapt_delta = 0.95)) {
  .check_fp_data(data)
  sl <- .build_stan_data(data, occ, det, fp, random_recorder, priors,
                         random_fp = random_fp, ver = ver, markov = markov,
                         random_det = random_det)
  gen <- list(data = data, formulas = list(occ = occ, det = det, fp = fp, ver = ver),
              stan_data = sl$stan_data, design = sl$design,
              occasions = sl$occasions)
  fit_args <- utils::modifyList(
    list(show_messages = FALSE, show_exceptions = FALSE, refresh = 0,
         threads_per_chain = 1), fit_args)
  jobs <- data.frame(rep = seq_len(reps), seed = as.integer(base_seed) + seq_len(reps))
  n_draws <- utils::modifyList(as.list(formals(fp_occu)[c("chains", "iter_sampling")]),
                               fit_args)
  if (ndraws > n_draws$chains * n_draws$iter_sampling) {
    stop("`ndraws` exceeds the number of posterior draws of each fit ",
         "(chains x iter_sampling).", call. = FALSE)
  }
  fp_compile()

  run_one <- function(rep, seed) {
    .set_seed(seed)
    p <- .prior_draw(sl$stan_data, priors)
    sim <- .simulate_from(gen, p)
    truth <- .sbc_truth(p, sl$design, sl$stan_data)
    with_ll <- !isTRUE(sl$stan_data$estimate_fp_re == 1) &&
      !isTRUE(sl$stan_data$estimate_det_re == 1)
    if (with_ll) truth <- c(truth, log_lik = NA_real_)
    fit <- tryCatch(suppressWarnings(do.call(fp_occu, c(list(
      data = sim$data, occ = occ, det = det, fp = fp, ver = ver, markov = markov,
      random_recorder = random_recorder, random_fp = random_fp,
      random_det = random_det, priors = priors, seed = seed), fit_args))),
      error = function(e) NULL)
    if (is.null(fit)) {
      return(data.frame(rep = rep, parameter = names(truth), truth = unname(truth),
                        rank = NA_integer_, divergences = NA_integer_, ok = FALSE))
    }
    d <- posterior::as_draws_matrix(fit$draws)
    keep <- round(seq(1, nrow(d), length.out = ndraws))
    d <- d[keep, setdiff(names(truth), "log_lik"), drop = FALSE]
    if (with_ll) {
      # log-likelihood of the simulated data: at the truth and at each draw
      truth[["log_lik"]] <- sum(.ppc_core(fit$stan_data, .truth_pars(p, fit$stan_data),
                                          "automated", "freeman-tukey")$ll)
      ll <- rowSums(unclass(posterior::as_draws_matrix(fit$fit$draws("log_lik"))))[keep]
      d <- cbind(d, log_lik = ll)
    }
    data.frame(rep = rep, parameter = names(truth), truth = unname(truth),
               rank = purrr::map_int(names(truth), function(v) {
                 sum(d[, v] < truth[[v]])
               }),
               divergences = fit$diagnostics$divergences, ok = TRUE)
  }
  out <- purrr::list_rbind(.run_jobs(jobs, run_one, parallel))
  structure(tibble::as_tibble(out), class = c("fp_sbc", "tbl_df", "tbl", "data.frame"),
            ndraws = ndraws)
}

# true values of the parameters, with the output names of fp_occu()
.sbc_truth <- function(p, design, sd) {
  c(stats::setNames(p$beta, paste0("occ_", design$occ_names)),
    stats::setNames(p$alpha, paste0("det_", design$det_names)),
    if (sd$N1 > 0) c(`verified_(Intercept)` = p$delta0),
    if (isTRUE(sd$Kver_s > 0))
      stats::setNames(p$delta_s, paste0("verified_", design$ver_names[-1])),
    if (sd$estimate_fp == 1) c(`fp_(Intercept)` = p$gamma0,
                               fp_gap = p$alpha[1] - p$gamma0),
    if (sd$estimate_fp == 1 && sd$Kfp > 0)
      stats::setNames(p$gamma_fp, paste0("fp_", design$fp_names[-1])),
    if (sd$estimate_re == 1) c(sd_recorder = p$sigma_u),
    if (isTRUE(sd$estimate_fp_re == 1)) c(sd_fp_site = p$sigma_fp),
    if (isTRUE(sd$estimate_det_re == 1)) c(sd_det_site = p$sigma_det))
}

# true parameter values in the layout of .ppc_pars(), as one draw
.truth_pars <- function(p, sd) {
  row <- function(v) matrix(v, nrow = 1)
  list(beta = row(p$beta), alpha = row(p$alpha), delta0 = p$delta0,
       delta_s = row(if (isTRUE(sd$Kver_s > 0)) p$delta_s else numeric(0)),
       u = if (sd$estimate_re == 1) row(p$u),
       gamma0 = if (sd$estimate_fp == 1) p$gamma0,
       gamma_fp = row(if (sd$estimate_fp == 1 && sd$Kfp > 0) p$gamma_fp else numeric(0)))
}

.check_bins <- function(L, bins) {
  if ((L + 1) %% bins != 0) {
    stop("`bins` must divide the number of possible ranks, ndraws + 1 = ",
         L + 1, ", so that every bin is equally likely.", call. = FALSE)
  }
}

#' @export
summary.fp_sbc <- function(object, bins = 10, ...) {
  L <- attr(object, "ndraws")
  .check_bins(L, bins)
  ok <- object[object$ok, , drop = FALSE]
  out <- purrr::map(split(ok, ok$parameter), function(g) {
    counts <- tabulate(cut(g$rank, breaks = seq(-0.5, L + 0.5, length.out = bins + 1),
                           labels = FALSE), nbins = bins)
    test <- suppressWarnings(stats::chisq.test(counts))
    data.frame(parameter = g$parameter[1], n = nrow(g),
               mean_rank = mean(g$rank) / L, p_uniform = test$p.value,
               divergent = mean(g$divergences > 0))
  }) |> purrr::list_rbind()
  structure(out, class = c("summary.fp_sbc", "data.frame"),
            failed = sum(!object$ok[!duplicated(object$rep)]))
}

#' @export
print.summary.fp_sbc <- function(x, digits = 3, ...) {
  cat("Simulation-based calibration\n")
  cat("  mean_rank: near 0.5 when the posterior is unbiased; p_uniform: chi-squared\n",
      " test of uniform ranks (small values flag a problem).\n\n")
  print(as.data.frame(x), digits = digits, row.names = FALSE)
  if (attr(x, "failed") > 0) cat("\n", attr(x, "failed"), "replicate(s) failed to fit.\n")
  invisible(x)
}

#' @rdname fp_sbc
#' @param x An `fp_sbc` object.
#' @param bins Number of histogram bins.
#' @param ... Unused.
#' @export
plot.fp_sbc <- function(x, bins = 10, ...) {
  .require_pkg("ggplot2")
  L <- attr(x, "ndraws")
  .check_bins(L, bins)
  ok <- x[x$ok, , drop = FALSE]
  n <- length(unique(ok$rep))
  expected <- n / bins
  band <- stats::qbinom(c(0.005, 0.995), n, 1 / bins)
  ggplot2::ggplot(ok, ggplot2::aes(rank)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2],
                      fill = "grey85") +
    ggplot2::geom_hline(yintercept = expected, colour = "grey40") +
    ggplot2::geom_histogram(breaks = seq(-0.5, L + 0.5, length.out = bins + 1),
                            fill = "#2DE2E6", colour = "white", alpha = 0.9) +
    ggplot2::facet_wrap(~ parameter) +
    ggplot2::labs(x = "Rank of the true value among the posterior draws",
                  y = "Replicates",
                  caption = "Grey band: 99% range of each bin's count under uniform ranks") +
    ggplot2::theme_minimal()
}

utils::globalVariables("rank")
