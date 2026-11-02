#' Detection histories over the draws of a classifier calibration
#'
#' Builds one version of the data per posterior draw of the calibration
#' exported by sheepdog (`calibration_draws.parquet`, the probability that
#' a prediction with a given score is a true detection, at each site and
#' species), for [fp_occu_mc()]. In each version an occasion is a detection
#' when at least one of its predictions is judged true under that draw.
#'
#' @param x An export read by [fp_read_sheepdog()], with calibration draws.
#'   Read it with `species` set to the target species, or to several species
#'   for the multispecies model (see Value).
#' @param rule How a prediction is judged true under a draw with
#'   probability `p` of being true: `"sample"` (default) draws it true with
#'   probability `p`, so that no threshold is involved; `"threshold"` takes
#'   it as true when `p >= p_min`, as `sheepdog::history_draws()` does.
#' @param p_min Threshold of `rule = "threshold"`.
#' @param n_draws Optional number of calibration draws to use (the first
#'   ones, in the order of their identifiers); all by default.
#' @param seed Optional seed of `rule = "sample"`. The random numbers of
#'   the session continue afterwards as if the call had drawn none.
#' @param occasions The occasions, one row per `site` and `occasion_id` of
#'   the export, with the column `site` names (the sites of the model, such
#'   as site-nights) and any detection covariates. Defaults to
#'   `x$occasions`. Every version carries a copy of it, so keep only the
#'   columns the model needs when there are many occasions and draws.
#' @param sites The sites table of [fp_data()], with the column `site`
#'   names. Defaults to `x$sites`.
#' @param site Column of `occasions` and `sites` that identifies the sites
#'   of the model.
#'
#' @details
#' **The rules.** `"sample"` judges each prediction on its own, so an
#' occasion with many low-scoring predictions is likely to hold at least one
#' judged true; it assumes that the predictions of an occasion are true or
#' false independently, given the calibration. `"threshold"` keeps only the
#' predictions the calibration is nearly sure of, so it gives fewer
#' detections: on the `curiango` data (custom classifier, minute
#' occasions), about 8% of the minutes with `"sample"` against about 1%
#' with `"threshold"`.
#'
#' When the calibration has keys beyond site and species
#' (`schema$calibration$keys`, such as `"month"`), a prediction takes the
#' value of its occasion: each such key is a column of the occasions of the
#' export and of `calibration_draws.parquet`. In an export of schema
#' version 1.0, without a `month` column, the month is taken from the
#' recording date of the occasion.
#'
#' A prediction takes the calibration value of the grid point nearest its
#' score, as in sheepdog: grid index
#' `round((score - grid_min) / grid_step) + 1`, within the grid. A prediction without calibration draws (no draws
#' for its site, species and keys, or a key missing for its occasion) cannot be scored: its occasion is
#' set to `NA` (not surveyed) in every version, with a warning, rather than
#' counted as a non-detection.
#'
#' **Histories drawn by sheepdog.** An export of schema 1.2 may hold the
#' histories sheepdog drew with `sheepdog::occasion_draws()`
#' (`occasion_draws.parquet`), by the same rules. They are read instead of
#' drawn again when they were made for the target species with the same
#' `rule`, the same `seed` (`"sample"`, which then needs `seed`) or
#' `p_min` (`"threshold"`), at least the draws asked for, and every
#' prediction has its occasion in `occasions`; the values are the same.
#'
#' **Calibrations with continuous covariates.** When the calibration has
#' continuous covariates (`schema$calibration$covariates`, schema 1.3, such
#' as a standardised temperature of the night), the probability that a
#' prediction is true depends on its occasion, not only on its score, and
#' sheepdog draws the histories when it exports; the export holds no
#' calibration draws. The histories are then read as they are: `rule`,
#' `p_min` and `seed` left at their defaults take the values of the export,
#' and other values, or more draws than were drawn, are an error (export
#' again). Occasions sheepdog could not score (without a value of a
#' covariate) are treated as not surveyed, with a warning, and occasions
#' not in `occasions` are left out.
#'
#' The recordings annotated by ear (`x$verified`) are the verified
#' occasions of every version: they do not depend on the calibration. The
#' reviewed clips (`x$validation`) are not used here; they are the data of
#' the calibration.
#'
#' Histories built this way hold few false positives, so a false-positive
#' process fitted to them (the default of [fp_occu()]) is informed mostly by
#' its priors; compare with `fp = NULL`.
#'
#' @return A list of [fp_data()] objects, one per calibration draw, with
#'   the attributes `draws` (the draw identifiers) and `rule`. When the
#'   detections hold several species, each version has every occasion once
#'   per species, with the column `species` (the short codes of sheepdog),
#'   as `fp_data(species = "species")`; the species are those with
#'   predictions.
#' @seealso [fp_occu_mc()], [fp_read_sheepdog()]
#' @section Experimental:
#' This function is experimental (new in 2.0.0, part of the integration with
#' 'sheepdog'): its arguments and results may change in a minor version.
#' See [crywolf-api].
#'
#' @export
#' @examples
#' \dontrun{
#' x <- fp_read_sheepdog("occupancy_export", species = "nyc_alb")
#' versions <- fp_draws_data(x, n_draws = 20, seed = 1)
#' mc <- fp_occu_mc(versions, occ = ~ 1, fp = NULL, threads_per_chain = 1)
#' summary(mc)
#' }
fp_draws_data <- function(x, rule = c("sample", "threshold"), p_min = 0.95,
                          n_draws = NULL, seed = NULL, occasions = x$occasions,
                          sites = x$sites, site = "site") {
  rule_given <- !missing(rule)
  p_min_given <- !missing(p_min)
  rule <- match.arg(rule)
  if (!inherits(x, "fp_sheepdog")) {
    stop("`x` must be an export read by fp_read_sheepdog().", call. = FALSE)
  }
  # continuous covariates of the calibration (schema 1.3): sheepdog drew the
  # histories when it exported, and there are no calibration draws to draw
  # them from here
  covariates <- unlist(x$schema$calibration$covariates)
  if (is.null(x$schema$calibration) ||
      (!length(covariates) && !"calibration_draws" %in% names(x$files))) {
    stop("The export has no calibration draws; write it with ",
         "sheepdog::export_occupancy_data(fit = ).", call. = FALSE)
  }
  .require_pkg("arrow")
  .check_df(occasions, "occasions", c("site", "occasion_id", site))
  det <- x$detections
  occ_key <- paste(occasions$site, occasions$occasion_id, sep = "\r")
  if (anyDuplicated(occ_key)) {
    stop("`occasions` must have one row per site and occasion_id.", call. = FALSE)
  }
  # several species: the occasions once per species, species after species;
  # occasion i of species s is row i + (s - 1) * n_occ
  sp <- sort(unique(det$species))
  n_occ <- nrow(occasions)
  auto0 <- occasions
  if (length(sp) > 1) {
    if ("species" %in% names(occasions)) {
      stop("`occasions` has a column `species`; with several species it is ",
           "added by fp_draws_data().", call. = FALSE)
    }
    auto0 <- occasions[rep(seq_len(n_occ), length(sp)), , drop = FALSE]
    auto0$species <- rep(sp, each = n_occ)
    rownames(auto0) <- NULL
  }
  row_of <- function(o, species) o + (match(species, sp) - 1L) * n_occ
  if (length(covariates)) {
    d <- .covariate_draws(x, if (rule_given) rule, if (p_min_given) p_min, seed,
                          n_draws, sp)
    ids <- d$ids
    rule <- d$rule
    saved <- d$saved
    o <- match(paste(saved$site, saved$occasion_id, sep = "\r"), occ_key)
    if (anyNA(o)) {
      warning(length(unique(paste(saved$site, saved$occasion_id)[is.na(o)])),
              " occasion(s) with predictions are not in `occasions` and were ",
              "left out.", call. = FALSE)
    }
    saved <- saved[!is.na(o), , drop = FALSE]
    n_na <- sum(is.na(saved$y[saved$draw == ids[1]]))
    if (n_na) {
      warning(n_na, " occasion(s) with predictions were not scored by sheepdog ",
              "(no value of a covariate or key of the calibration): they are ",
              "treated as not surveyed.", call. = FALSE)
    }
  } else {
    cal <- as.data.frame(arrow::read_parquet(x$files[["calibration_draws"]]))
    # the origin of the grid with the exact bits sheepdog used: the value
    # read back from schema.json can differ in the last digit and send a
    # score half-way between two points to the other one
    grid_min <- min(cal$confidence[cal$grid_index == 1])
    cal <- cal[cal$species %in% sp, , drop = FALSE]
    ids <- sort(unique(cal$draw))
    if (!is.null(n_draws)) ids <- utils::head(ids, n_draws)

    # keys of the calibration beyond site, species and the grid (such as
    # "month"): columns of the occasions of the export and of the parquet
    keys <- setdiff(unlist(x$schema$calibration$keys), c("site", "species", "grid_index"))
    # (a key column the schema does not announce shows up below as several
    # values per key and draw)
    if (!all(keys %in% names(cal))) {
      stop("schema.json and calibration_draws.parquet disagree on the keys of ",
           "the calibration: the parquet lacks ",
           paste(setdiff(keys, names(cal)), collapse = ", "), ".", call. = FALSE)
    }
    oe <- x$occasions
    # exports of schema 1.0 have no month column: the month of recording
    if ("month" %in% keys && !"month" %in% names(oe)) {
      oe$month <- substr(oe$recording_date, 1, 6)
    }
    if (!all(keys %in% names(oe))) {
      stop("The occasions of the export lack the calibration key(s) ",
           paste(setdiff(keys, names(oe)), collapse = ", "), ".", call. = FALSE)
    }
    e <- match(paste(det$site, det$occasion_id, sep = "\r"),
               paste(oe$site, oe$occasion_id, sep = "\r"))
    det_keys <- lapply(keys, function(k) as.character(oe[[k]][e]))
    # P(TP) as a matrix: one row per site x species (x keys) x grid point,
    # one column per draw
    key <- function(...) paste(..., sep = "\r")
    cal <- cal[cal$draw %in% ids, , drop = FALSE]
    cal_key <- do.call(key, c(list(cal$site, cal$species),
                              lapply(keys, function(k) as.character(cal[[k]])),
                              list(cal$grid_index)))
    rows <- unique(cal_key)
    if (anyDuplicated(cbind(cal_key, cal$draw))) {
      stop("calibration_draws.parquet has more than one value for some key and ",
           "draw; is a key column missing from schema.json?", call. = FALSE)
    }
    P <- matrix(NA_real_, length(rows), length(ids))
    P[cbind(match(cal_key, rows), match(cal$draw, ids))] <- cal$p_tp
    step <- x$schema$calibration$grid_step
    g <- round((det$confidence - grid_min) / step) + 1
    g <- pmin(pmax(g, 1), max(cal$grid_index))
    o <- row_of(match(paste(det$site, det$occasion_id, sep = "\r"), occ_key), det$species)
    k <- match(do.call(key, c(list(det$site, det$species), det_keys, list(g))), rows)
    # an occasion without a value of a key cannot be matched to a calibration
    k[Reduce(`|`, lapply(det_keys, is.na), FALSE)] <- NA
    if (any(is.na(o))) {
      warning(sum(is.na(o)), " prediction(s) are not in `occasions` and were ",
              "left out.", call. = FALSE)
    }
    # an occasion with a prediction that has no calibration cannot be scored:
    # it is treated as not surveyed (y = NA), not as a non-detection
    unscored <- unique(o[is.na(k) & !is.na(o)])
    if (length(unscored)) {
      warning(sum(is.na(k) & !is.na(o)), " prediction(s) of site x species",
              paste0(" x ", keys, collapse = ""), " without calibration draws",
              if (length(keys)) " or without a value of a key", ": their ",
              length(unscored), " occasion(s) are treated as not surveyed.",
              call. = FALSE)
    }
    # the histories sheepdog drew with the same settings, when it wrote them
    # (occasion_draws.parquet, schema 1.2): the same values, computed once
    saved <- .saved_occasion_draws(x, rule, p_min, seed, ids, sp,
                                   complete = !anyNA(o))
    keep <- !is.na(k) & !is.na(o) & !(o %in% unscored)
    if (!any(keep)) {
      stop("No prediction could be matched to the calibration draws.", call. = FALSE)
    }
    k <- k[keep]
    o <- o[keep]
  }

  verified <- x$verified
  if (!is.null(verified) && site != "site") {
    v <- match(paste(verified$site, verified$occasion_id, sep = "\r"), occ_key)
    verified[[site]] <- occasions[[site]][v]
    if (anyNA(v)) {
      warning(sum(is.na(v)), " annotated recording(s) are not in `occasions` and ",
              "were left out.", call. = FALSE)
    }
    verified <- verified[!is.na(v), , drop = FALSE]
  }
  if (!is.null(verified)) verified <- verified[c(site, setdiff(names(verified), site))]
  if (!is.null(saved)) {
    # occasions without predictions are not in the table: no detection
    Y <- matrix(0L, nrow(auto0), length(ids))
    Y[cbind(row_of(match(paste(saved$site, saved$occasion_id, sep = "\r"), occ_key),
                   saved$species),
            match(saved$draw, ids))] <- saved$y
  } else if (!is.null(seed) && rule == "sample") {
    # the user's random numbers continue as if this call had drawn none
    old <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    }
    .set_seed(seed)
    on.exit(if (is.null(old)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old, envir = globalenv())
    }, add = TRUE)
  }

  out <- purrr::imap(ids, function(id, j) {
    auto <- auto0
    if (is.null(saved)) {
      p <- P[k, j]
      true <- if (rule == "sample") stats::runif(length(p)) < p else p >= p_min
      auto$y <- as.integer(tabulate(o[true], nbins = nrow(auto)) > 0)
      auto$y[unscored] <- NA
    } else {
      auto$y <- Y[, j]
    }
    # the design checks of fp_data() speak once, for the first version
    species <- if (length(sp) > 1) "species"
    if (j == 1) return(fp_data(auto, sites, verified, site = site, species = species))
    withCallingHandlers(fp_data(auto, sites, verified, site = site, species = species),
                        crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
  })
  structure(out, draws = ids, rule = rule)
}

# The detection histories of occasion_draws.parquet, when sheepdog drew them
# as fp_draws_data() would here: same species, rule, p_min or seed, at least
# these draws, and every prediction in `occasions` (`complete`). NULL
# otherwise.
.saved_occasion_draws <- function(x, rule, p_min, seed, ids, species, complete) {
  od <- x$schema$occasion_draws
  if (is.null(od) || !"occasion_draws" %in% names(x$files) || !complete ||
      !identical(od$rule, rule) || !all(species %in% unlist(od$species))) {
    return(NULL)
  }
  same <- if (rule == "sample") {
    !is.null(seed) && !is.null(od$seed) && isTRUE(all.equal(od$seed, seed))
  } else {
    !is.null(od$p_min) && isTRUE(all.equal(od$p_min, p_min))
  }
  if (!same) return(NULL)
  saved <- as.data.frame(arrow::read_parquet(x$files[["occasion_draws"]]))
  saved <- saved[saved$species %in% species & saved$draw %in% ids, , drop = FALSE]
  if (!all(ids %in% saved$draw)) return(NULL)
  saved
}

# The detection histories sheepdog drew at export for a calibration with
# continuous covariates (schema 1.3), the only ones there are: `rule`,
# `p_min` and `seed` left NULL take the ones of the export, and other values
# are an error. A list with `saved` (occasion_draws.parquet, target species,
# the draws used), `ids` and `rule`.
.covariate_draws <- function(x, rule, p_min, seed, n_draws, species) {
  od <- x$schema$occasion_draws
  covariates <- paste(unlist(x$schema$calibration$covariates), collapse = ", ")
  again <- paste0("; export it again with sheepdog::export_occupancy_data() ",
                  "set to the histories needed.")
  if (is.null(od) || !"occasion_draws" %in% names(x$files)) {
    stop("The calibration of the export has continuous covariates (", covariates,
         ") and the export holds no histories drawn by sheepdog ",
         "(occasion_draws.parquet)", again, call. = FALSE)
  }
  differ <- function(what, asked, drawn) {
    stop("The calibration of the export has continuous covariates (", covariates,
         "), so the histories can only be those sheepdog drew, with ", what, " = ",
         drawn, ", not ", asked, again, call. = FALSE)
  }
  if (!all(species %in% unlist(od$species))) {
    differ("species", paste(species, collapse = ", "),
           paste(unlist(od$species), collapse = ", "))
  }
  if (is.null(rule)) rule <- od$rule
  if (!identical(rule, od$rule)) differ("rule", dQuote(rule, FALSE), dQuote(od$rule, FALSE))
  if (rule == "sample" && !is.null(seed) && !isTRUE(all.equal(od$seed, seed))) {
    differ("seed", seed, od$seed)
  }
  if (rule == "threshold" && !is.null(p_min) && !isTRUE(all.equal(od$p_min, p_min))) {
    differ("p_min", p_min, od$p_min)
  }
  saved <- as.data.frame(arrow::read_parquet(x$files[["occasion_draws"]]))
  saved <- saved[saved$species %in% species, , drop = FALSE]
  ids <- sort(unique(saved$draw))
  if (!is.null(n_draws)) {
    if (n_draws > length(ids)) differ("n_draws", n_draws, length(ids))
    ids <- utils::head(ids, n_draws)
  }
  list(saved = saved[saved$draw %in% ids, , drop = FALSE], ids = ids, rule = rule)
}
