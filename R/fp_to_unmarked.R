#' Convert data for unmarked::occuFP()
#'
#' Turns an `fp_data` object into the `unmarkedFrameOccuFP` that
#' `unmarked::occuFP()` fits by maximum likelihood, the frequentist
#' counterpart of [fp_occu()] (see `vignette("unmarked")`). With
#' `occasion = "occasion"` it is the reverse of [fp_data_wide()].
#'
#' @param data An `fp_data` object, for instance built with
#'   [fp_aggregate()] and [fp_data()] from [fp_read_sheepdog()].
#' @param obs_covs Optional names of numeric columns of `data$auto` to pass
#'   as observation covariates.
#' @param occasion Optional name of a column of `data$auto` (and of
#'   `data$verified`, if it has one) that numbers the occasions. Each
#'   occasion then goes to the column of its rank among all the values, so
#'   that occasions missing at a site leave `NA` in their place. Without it,
#'   a site's occasions are placed in the order of their rows.
#'
#' @details
#' The detection matrix has one row per site, in the order of `data$sites`,
#' and the verified occasions first: `type = c(n_verified, n_automated, 0)`,
#' where `n_verified` and `n_automated` are the numbers of verified and
#' automated occasion columns (see `occasion`); sites with fewer occasions
#' get `NA`.
#'
#' The observation covariate `METH` is `"verified"` on the verified columns
#' and `"auto"` on the automated ones (the reference level), because
#' `occuFP()` describes true detection by both methods with one formula:
#' use `detformula = ~ METH` for a separate detection probability of the
#' verified method, as [fp_occu()] has. The covariates in `obs_covs` exist
#' for the automated occasions only; they are 0 on the verified columns,
#' because `occuFP()` drops every observation whose covariate is missing.
#' Centre them, so that 0 is a typical value, and read a covariate in
#' `detformula` as acting on the automated method.
#'
#' Recorders are not passed on (`occuFP()` has no random effects), and data
#' of the observation-confirmation design (`confirmed` in [fp_data()]) are
#' not supported.
#'
#' @return An `unmarkedFrameOccuFP` object.
#' @seealso [fp_data_wide()], `vignette("unmarked")`
#' @section Experimental:
#' This function is experimental (new in 2.0.0, part of the integration with
#' 'sheepdog'): its arguments and results may change in a minor version.
#' See [crywolf-api].
#'
#' @export
#' @examplesIf requireNamespace("unmarked", quietly = TRUE)
#' sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
#' umf <- fp_to_unmarked(sim$data, occasion = "occasion")
#' unmarked::occuFP(detformula = ~ METH, FPformula = ~ 1, stateformula = ~ x,
#'                  data = umf)
fp_to_unmarked <- function(data, obs_covs = NULL, occasion = NULL) {
  .require_pkg("unmarked")
  .check_fp_data(data)
  if (!is.null(data$confirmed)) {
    stop("Data of the observation-confirmation design cannot be passed to ",
         "occuFP(); fit them with fp_occu().", call. = FALSE)
  }
  .check_df(data$auto, "data$auto", c(obs_covs, occasion))
  if ("METH" %in% obs_covs) {
    stop("`METH` is the method covariate made by fp_to_unmarked(); rename ",
         "the column.", call. = FALSE)
  }
  bad <- obs_covs[!purrr::map_lgl(obs_covs, function(v) is.numeric(data$auto[[v]]))]
  if (length(bad)) {
    stop("Observation covariates must be numeric: ", paste(bad, collapse = ", "),
         ".", call. = FALSE)
  }
  ids <- data$sites[[data$site]]
  # rows of one table as a sites x occasions matrix: occasions in the column
  # of their rank (with `occasion`) or in row order
  wide <- function(tb, values) {
    if (is.null(tb) || nrow(tb) == 0) return(matrix(NA, length(ids), 0))
    i <- match(tb[[data$site]], ids)
    pos <- if (!is.null(occasion) && occasion %in% names(tb)) {
      match(tb[[occasion]], sort(unique(tb[[occasion]])))
    } else {
      stats::ave(seq_along(i), i, FUN = seq_along)
    }
    if (anyDuplicated(cbind(i, pos))) {
      stop("Two occasions of a site have the same `", occasion, "`.", call. = FALSE)
    }
    m <- matrix(NA, length(ids), max(pos))
    m[cbind(i, pos)] <- values
    m
  }
  y_ver <- wide(data$verified, data$verified[[data$y]])
  y_auto <- wide(data$auto, data$auto[[data$y]])
  n_ver <- ncol(y_ver)
  if (n_ver == 0) {
    warning("No verified occasions: occuFP() cannot tell true from false ",
            "detections without them.", call. = FALSE)
  }
  # observation covariates in the long layout of unmarked: site by site,
  # occasions within sites
  long <- function(m) as.vector(t(m))
  meth <- rep(c("verified", "auto"), times = c(n_ver, ncol(y_auto)))
  covs <- data.frame(METH = factor(rep(meth, length(ids)), levels = c("auto", "verified")))
  covs[obs_covs] <- purrr::map(obs_covs, function(v) {
    long(cbind(matrix(0, length(ids), n_ver), wide(data$auto, data$auto[[v]])))
  })
  site_covs <- data$sites[setdiff(names(data$sites), data$site)]
  unmarked::unmarkedFrameOccuFP(
    y = cbind(y_ver, y_auto),
    siteCovs = if (ncol(site_covs)) site_covs else NULL,
    obsCovs = covs, type = c(n_ver, ncol(y_auto), 0))
}
