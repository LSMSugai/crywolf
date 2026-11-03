#' Assemble data for a false-positive occupancy model
#'
#' Collects the automated detections, the verified detections, and the site
#' table into one object that [fp_occu()] can fit. All tables are in long
#' format, one row per occasion, which accommodates unequal numbers of
#' occasions per site without padding. Use [fp_data_wide()] if your data are
#' detection matrices, as in 'unmarked'.
#'
#' @param auto Data frame of the automated method (the one that can return
#'   false positives), one row per site and occasion. It needs a site column
#'   and a 0/1 detection column, and may hold observation covariates for the
#'   true-detection process and a recorder (or observer) identifier.
#' @param sites Data frame with one row per site: the site identifier plus any
#'   site covariates used in the occupancy (`occ`) or false-positive (`fp`)
#'   formulas of [fp_occu()]. Sites with no automated occasions are allowed and
#'   contribute only through their verified data, if any.
#' @param verified Optional data frame of the verified method (detections
#'   checked by an expert and therefore free of false positives), one row per
#'   site and occasion, with the same site and detection columns. Without it
#'   the false-positive rate is identified only by the priors, and
#'   [fp_occu()] warns about it.
#' @param site,y Names of the site and detection columns, shared by all
#'   tables.
#' @param recorder Optional name of a column of `auto` that identifies the
#'   recorder or observer, used for a random intercept on true detection.
#' @param confirmed Optional names of two columns of `auto` holding, for the
#'   occasions whose observations were checked by an expert, the number of
#'   observations confirmed as the target species and the number confirmed as
#'   something else, in that order (the observation-confirmation design of
#'   Chambert et al. 2015). Use `NA` in both for occasions that were not
#'   checked. All observations of a checked occasion must have been checked,
#'   so an occasion with a detection has at least one confirmed observation,
#'   and one without has none. See Details of [fp_occu()].
#' @param season Optional name of a column of `auto` (and of `verified`, if
#'   given) holding the season, or primary period, of each occasion, for the
#'   dynamic model. Seasons are ordered by the levels of a factor, and
#'   otherwise by sorting their values, so numbers, dates and ISO dates as
#'   text (such as `night` or `recording_date` of a 'sheepdog' export) all
#'   work. When the seasons are dates (`Date`, or text as `YYYY-MM-DD` or
#'   `YYYYMMDD`), the dates missing between the first and the last are added
#'   as seasons without data, in steps of the shortest interval between
#'   seasons, so that every transition spans the same time (one night, for
#'   nightly seasons, even when no site recorded on some nights). Other
#'   seasons (numbers, factors, other text) are taken as they come:
#'   consecutive seasons are one transition apart, whatever the time between
#'   them. A site without data in a season contributes to that season only
#'   through the transitions of its occupancy state. [fp_occu()] ignores
#'   this column.
#' @param season_covs Optional data frame of covariates of the seasons, for
#'   the colonization, extinction and initial occupancy formulas of
#'   [fp_occu_dyn()] (and its detection and false-positive formulas). It has
#'   the season column named by `season` and, for covariates that vary among
#'   sites too, the site column; one row per season, or per site and season.
#'   Seasons added between dates need rows only if a formula uses them.
#' @param species Optional name of a column of `auto` (and of `verified`, if
#'   given) holding the species, for the multispecies model
#'   [fp_occu_ms()], which is in development. The tables are then in long format by species too: one
#'   row per occasion and species, with every surveyed occasion present for
#'   every species (`y = 0` where the species was not detected), since the
#'   non-detections are data. The species are the levels of a factor, and
#'   otherwise the sorted values of `auto` and `verified`. It cannot be
#'   combined with `season` yet, and the single-species models refuse such
#'   data.
#'
#' @details
#' Besides errors for data that cannot be fitted, `fp_data()` warns about
#' data that can be fitted but are probably not what was meant: sites with
#' no automated occasion (often a failed join of `auto` and `sites`), very
#' unequal automated effort among sites (ten times or more), automated
#' occasions that are all detections or all non-detections, verified data
#' without any detection, and confirmed occasions without any false
#' observation. With `season`, it also warns about seasons without data at
#' any site (unused levels of a factor, or dates added between the first
#' and the last), about seasons without verified or confirmed data, in which
#' colonization and extinction are hard to tell from false positives, and
#' about seasons at unequal intervals that could not be filled in. With
#' `species`, it warns about species never detected by the automated method,
#' species without any verified detection, whose false-positive rates are
#' then informed only by the other species and the priors, and sites where
#' the species have unequal numbers of occasions (often occasions without a
#' detection left out for some species).
#' [fp_aggregate()] builds `auto` from classifier output.
#'
#' @return An object of class `fp_data`: a list with the tables `auto`,
#'   `verified` and `sites` (rows with a missing detection dropped), the
#'   column names `site`, `y`, `recorder`, `confirmed` and `season`, the
#'   vector `known_present`, which flags sites with at least one verified
#'   detection or confirmed true observation, and, with `season`, the
#'   seasons in time order in `seasons`, the step between dates in days in
#'   `season_step` (`NULL` for seasons that are not dates) and the table
#'   `season_covs`. With `species`, the species in `species_names`, and
#'   `known_present` is a matrix of sites by species.
#' @references
#' Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
#' positive detections in species occurrence data under different study
#' designs. *Ecology*, 96, 332–339. \doi{10.1890/14-1507.1}
#' @seealso [fp_data_wide()], [fp_occu()], [simulate_fp()]
#' @export
#' @examples
#' sim <- simulate_fp(nsite = 50, nocc = 5, seed = 1)
#' sim$data
#'
#' # the same structure built by hand
#' d <- fp_data(auto = sim$data$auto, sites = sim$data$sites,
#'              verified = sim$data$verified, recorder = "recorder")
#' d
fp_data <- function(auto, sites, verified = NULL, site = "site", y = "y",
                    recorder = NULL, confirmed = NULL, season = NULL,
                    season_covs = NULL, species = NULL) {
  .check_df(auto, "auto", c(site, y, recorder, confirmed, season, species))
  .check_df(sites, "sites", site)
  if (anyDuplicated(sites[[site]])) {
    stop("`sites` must have one row per site; duplicated identifiers found in `",
         site, "`.", call. = FALSE)
  }
  auto <- auto[!is.na(auto[[y]]), , drop = FALSE]
  .check_binary(auto[[y]], "auto")
  .check_sites_known(auto[[site]], sites[[site]], "auto", site)

  if (!is.null(verified)) {
    .check_df(verified, "verified", c(site, y, season, species))
    verified <- verified[!is.na(verified[[y]]), , drop = FALSE]
    .check_binary(verified[[y]], "verified")
    .check_sites_known(verified[[site]], sites[[site]], "verified", site)
  }
  if (!is.null(recorder) && anyNA(auto[[recorder]])) {
    stop("Column `", recorder, "` of `auto` has missing values.", call. = FALSE)
  }

  if (!is.null(confirmed)) .check_confirmed(auto, confirmed, y)
  species_names <- NULL
  if (!is.null(species)) {
    if (!is.null(season)) {
      stop("`species` and `season` cannot be combined yet.", call. = FALSE)
    }
    species_names <- .species_levels(auto, verified, species)
  }
  seasons <- season_step <- NULL
  if (!is.null(season)) {
    seasons <- .season_levels(auto, verified, season)
    # the levels of a factor are taken as they are
    if (!is.factor(auto[[season]])) seasons <- .fill_seasons(seasons)
    season_step <- attr(seasons, "step")
    attr(seasons, "step") <- NULL
  }
  if (!is.null(season_covs)) .check_season_covs(season_covs, season, seasons, sites, site)
  .warn_design(auto, sites, verified, site, y, confirmed)
  if (!is.null(season)) .warn_seasons(auto, verified, confirmed, season, seasons)
  if (!is.null(species)) .warn_species(auto, verified, site, y, species, species_names)

  # sites (by species) with a verified detection or a confirmed true
  # observation
  known_present <- matrix(0L, nrow(sites), max(1L, length(species_names)),
                          dimnames = list(NULL, species_names))
  mark <- function(d, rows) {
    sp <- if (is.null(species)) 1L else match(as.character(d[[species]][rows]), species_names)
    known_present[cbind(match(d[[site]][rows], sites[[site]]), sp)] <<- 1L
  }
  if (!is.null(verified) && nrow(verified) > 0) mark(verified, verified[[y]] == 1)
  if (!is.null(confirmed)) mark(auto, which(auto[[confirmed[1]]] > 0))
  if (is.null(species)) known_present <- known_present[, 1]

  structure(
    list(auto = auto, verified = verified, sites = sites,
         site = site, y = y, recorder = recorder, confirmed = confirmed,
         known_present = known_present, season = season, seasons = seasons,
         season_step = season_step, season_covs = season_covs,
         species = species, species_names = species_names),
    class = "fp_data"
  )
}

#' Assemble data from detection matrices
#'
#' Convenience wrapper around [fp_data()] for data stored as site-by-occasion
#' matrices, the layout used by 'unmarked'. Missing values (`NA`) mark
#' occasions that were not sampled, so sites may differ in effort.
#'
#' @param y_auto Matrix of automated detections (sites in rows, occasions in
#'   columns), coded 0/1/`NA`.
#' @param y_verified Optional matrix of verified detections with the same
#'   number of rows; unverified sites are rows of `NA`.
#' @param site_covs Optional data frame of site covariates, one row per site,
#'   in the same order as the rows of `y_auto`.
#' @param obs_covs Optional named list of matrices, each the same size as
#'   `y_auto`, holding observation covariates of the automated method.
#' @param recorder Optional vector, one value per site, identifying the
#'   recorder or observer.
#'
#' @return An object of class `fp_data`, as from [fp_data()].
#' @seealso [fp_data()]
#' @export
#' @examples
#' set.seed(1)
#' ya <- matrix(rbinom(40, 1, 0.3), 10, 4)
#' yv <- matrix(NA, 10, 2); yv[1:3, ] <- rbinom(6, 1, 0.5)
#' fp_data_wide(ya, yv, site_covs = data.frame(x = rnorm(10)))
fp_data_wide <- function(y_auto, y_verified = NULL, site_covs = NULL,
                         obs_covs = NULL, recorder = NULL) {
  y_auto <- as.matrix(y_auto)
  nsite <- nrow(y_auto)
  sites <- data.frame(site = seq_len(nsite))
  if (!is.null(site_covs)) {
    if (nrow(site_covs) != nsite) {
      stop("`site_covs` must have one row per row of `y_auto`.", call. = FALSE)
    }
    sites <- cbind(sites, site_covs)
  }
  to_long <- function(m) {
    data.frame(site = rep(seq_len(nrow(m)), times = ncol(m)),
               occasion = rep(seq_len(ncol(m)), each = nrow(m)),
               y = as.vector(m))
  }
  auto <- to_long(y_auto)
  if (!is.null(obs_covs)) {
    if (is.null(names(obs_covs)) || any(!nzchar(names(obs_covs)))) {
      stop("`obs_covs` must be a named list of matrices.", call. = FALSE)
    }
    purrr::iwalk(obs_covs, function(m, nm) {
      if (!identical(dim(as.matrix(m)), dim(y_auto))) {
        stop("Observation covariate `", nm, "` must have the same dimensions ",
             "as `y_auto`.", call. = FALSE)
      }
    })
    auto <- cbind(auto, purrr::map(obs_covs, function(m) as.vector(as.matrix(m))) |>
                        as.data.frame())
  }
  rec_col <- NULL
  if (!is.null(recorder)) {
    if (length(recorder) != nsite) {
      stop("`recorder` must have one value per site.", call. = FALSE)
    }
    auto$recorder <- recorder[auto$site]
    rec_col <- "recorder"
  }
  verified <- if (is.null(y_verified)) NULL else to_long(as.matrix(y_verified))
  fp_data(auto = auto, sites = sites, verified = verified, recorder = rec_col)
}

#' @export
print.fp_data <- function(x, ...) {
  ns <- nrow(x$sites)
  nver_sites <- if (is.null(x$verified)) 0L else length(unique(x$verified[[x$site]]))
  by <- list(x$auto[[x$site]])
  if (!is.null(x$species)) by <- c(by, list(x$auto[[x$species]]))
  # by site (and species): sites x species without occasions are NA
  naive_auto <- mean(tapply(x$auto[[x$y]], by, max), na.rm = TRUE)
  cat("<fp_data>\n")
  cat(sprintf("  sites:                  %d\n", ns))
  cat(sprintf("  automated occasions:    %d (naive occupancy %.2f)\n",
              nrow(x$auto), naive_auto))
  if (!is.null(x$species)) cat(sprintf("  species:                %d\n", length(x$species_names)))
  ver_det <- if (is.null(x$verified)) character(0) else
    unique(x$verified[[x$site]][x$verified[[x$y]] == 1])
  cat(sprintf("  verified sites:         %d (%d occasions; verified detections at %d sites)\n",
              nver_sites, if (is.null(x$verified)) 0L else nrow(x$verified),
              length(ver_det)))
  if (!is.null(x$confirmed)) {
    n_true <- x$auto[[x$confirmed[1]]]
    true_sites <- unique(x$auto[[x$site]][!is.na(n_true) & n_true > 0])
    cat(sprintf("  confirmed occasions:    %d (%d with a true observation, at %d sites)\n",
                sum(!is.na(n_true)), sum(n_true > 0, na.rm = TRUE), length(true_sites)))
  }
  cat(if (is.null(x$species)) "  sites known occupied:   " else "  site x species known occupied: ",
      sum(x$known_present), "\n", sep = "")
  if (!is.null(x$recorder)) {
    cat(sprintf("  recorders:              %d\n", length(unique(x$auto[[x$recorder]]))))
  }
  if (!is.null(x$season)) {
    n_sites <- .sites_per_season(x)
    cat(sprintf("  seasons:                %d%s (sites with data per season: %d to %d, median %g)\n",
                length(n_sites),
                if (is.null(x$season_step)) "" else sprintf(", %g day(s) apart", x$season_step),
                min(n_sites), max(n_sites), stats::median(n_sites)))
  }
  covs <- setdiff(names(x$sites), x$site)
  if (length(covs)) cat("  site covariates:       ", paste(covs, collapse = ", "), "\n")
  ocovs <- setdiff(names(x$auto), c(x$site, x$y, x$recorder, x$confirmed,
                                    x$season, x$species, "occasion"))
  if (length(ocovs)) cat("  observation covariates:", paste(ocovs, collapse = ", "), "\n")
  scovs <- setdiff(names(x$season_covs), c(x$site, x$season))
  if (length(scovs)) cat("  season covariates:     ", paste(scovs, collapse = ", "), "\n")
  if (!is.null(x$species)) {
    cat("  by species:\n")
    print(.species_table(x), row.names = FALSE)
  }
  invisible(x)
}

# validation helpers ----------------------------------------------------------

# the data of the single-species models
.check_fp_data <- function(data) {
  if (!inherits(data, "fp_data")) {
    stop("`data` must be an `fp_data` object; see `fp_data()`.", call. = FALSE)
  }
  if (!is.null(data$species)) {
    stop("`data` holds several species (`fp_data(species = )`); fit them with ",
         "fp_occu_ms(), or keep the rows of one species to fit it alone.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.check_df <- function(df, what, cols) {
  if (!is.data.frame(df)) stop("`", what, "` must be a data frame.", call. = FALSE)
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop("`", what, "` lacks column(s): ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

.check_binary <- function(v, what) {
  if (!all(v %in% c(0, 1))) {
    stop("Detections in `", what, "` must be coded 0/1.", call. = FALSE)
  }
  invisible(TRUE)
}

# warnings about data that fit but are unlikely to be what the user meant;
# their class lets the simulation functions silence them
.warn_design <- function(auto, sites, verified, site, y, confirmed) {
  warn <- function(...) {
    msg <- paste0(...)
    warning(structure(class = c("crywolf_design_warning", "warning", "condition"),
                      list(message = msg, call = NULL)))
  }
  n_occ <- table(factor(auto[[site]], levels = sites[[site]]))
  if (any(n_occ == 0)) {
    warn(sum(n_occ == 0), " site(s) have no automated occasion; they inform ",
         "only the occupancy model, through their covariates. Check the join ",
         "of `auto` and `sites` if this is unexpected.")
  }
  if (any(n_occ > 0) && max(n_occ) >= 10 * min(n_occ[n_occ > 0])) {
    warn("Automated effort is very unequal among sites (", min(n_occ[n_occ > 0]),
         " to ", max(n_occ), " occasions); check that occasions were not lost ",
         "or duplicated.")
  }
  if (nrow(auto) > 0 && length(unique(auto[[y]])) == 1) {
    warn("Every automated occasion is ",
         if (auto[[y]][1] == 1) "a detection" else "a non-detection",
         "; the detection and false-positive rates cannot be estimated.")
  }
  if (!is.null(verified) && nrow(verified) > 0 && all(verified[[y]] == 0)) {
    warn("No verified occasion is a detection, so the verified data never ",
         "prove occupancy and anchor the false-positive rate only weakly.")
  }
  if (!is.null(confirmed) && any(!is.na(auto[[confirmed[2]]])) &&
      sum(auto[[confirmed[2]]], na.rm = TRUE) == 0) {
    warn("No confirmed occasion holds a false observation; the false-positive ",
         "rate will be estimated close to zero. Check that the checked ",
         "occasions were chosen independently of the classifier's score.")
  }
  invisible(TRUE)
}

# the seasons in time order: the levels of a factor, otherwise the sorted
# values of `auto` and `verified` together
.season_levels <- function(auto, verified, season) {
  a <- auto[[season]]
  v <- if (is.null(verified)) NULL else verified[[season]]
  if (anyNA(a) || anyNA(v)) {
    stop("Column `", season, "` has missing values.", call. = FALSE)
  }
  if (is.factor(a)) {
    lv <- levels(a)
    bad <- setdiff(as.character(v), lv)
    if (length(bad)) {
      stop("Seasons in `verified` are not levels of `auto$", season, "` (e.g. ",
           paste(utils::head(bad, 3), collapse = ", "), ").", call. = FALSE)
    }
    return(lv)
  }
  if (!is.null(v) && !(is.numeric(a) && is.numeric(v)) &&
      !identical(class(a), class(v))) {
    stop("Column `", season, "` must have the same class in `auto` and ",
         "`verified`.", call. = FALSE)
  }
  sort(unique(c(a, v)))
}

# Seasons that are dates: the missing dates between the first and the last
# are added, in steps of the shortest interval, so that each transition
# spans the same time. Attribute "step": that interval in days (NULL when
# the seasons are not dates).
.fill_seasons <- function(lv) {
  d <- .season_dates(lv)
  gaps <- if (is.numeric(lv)) diff(lv) else numeric(0)
  numeric_gaps <- length(gaps) > 1 && diff(range(gaps)) > 1e-8 * max(abs(gaps))
  if (is.null(d) || length(d) < 2) {
    if (numeric_gaps) {
      .design_warning("Consecutive seasons are unequally spaced (", paste(utils::head(
        unique(gaps), 3), collapse = ", "), "); each counts as one transition.")
    }
    return(lv)
  }
  dn <- as.numeric(d)
  step <- min(diff(dn))
  if (any((dn - dn[1]) %% step != 0)) {
    .design_warning("The season dates are not evenly spaced by their shortest ",
                    "interval (", step, " day(s)); each pair of consecutive ",
                    "seasons counts as one transition, whatever its length.")
    return(lv)
  }
  full <- seq(d[1], d[length(d)], by = step)
  out <- if (inherits(lv, "Date")) full else format(full, attr(d, "format"))
  structure(out, step = step)
}

# the seasons as dates (Date, or text as YYYY-MM-DD or YYYYMMDD), with the
# text format as attribute; NULL if they are not dates
.season_dates <- function(lv) {
  if (inherits(lv, "Date")) return(structure(lv, format = NULL))
  if (!is.character(lv)) return(NULL)
  fmt <- if (all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", lv))) "%Y-%m-%d" else
    if (all(grepl("^[0-9]{8}$", lv))) "%Y%m%d" else return(NULL)
  d <- as.Date(lv, format = fmt)
  if (anyNA(d)) return(NULL)
  structure(d, format = fmt)
}

.check_season_covs <- function(season_covs, season, seasons, sites, site) {
  if (is.null(season)) {
    stop("`season_covs` needs `season`.", call. = FALSE)
  }
  .check_df(season_covs, "season_covs", season)
  keys <- intersect(c(site, season), names(season_covs))
  if (anyDuplicated(season_covs[keys])) {
    stop("`season_covs` must have one row per ", paste(keys, collapse = " and "),
         ".", call. = FALSE)
  }
  bad <- season_covs[[season]][is.na(.season_index(season_covs, season, seasons))]
  if (length(bad)) {
    stop("Seasons in `season_covs` are not seasons of the data (e.g. ",
         paste(utils::head(bad, 3), collapse = ", "), ").", call. = FALSE)
  }
  if (site %in% keys) .check_sites_known(season_covs[[site]], sites[[site]], "season_covs", site)
  invisible(TRUE)
}

.design_warning <- function(...) {
  warning(structure(class = c("crywolf_design_warning", "warning", "condition"),
                    list(message = paste0(...), call = NULL)))
}

# index of the season of each row of a table; numbers and dates are matched
# by value, factor levels and text as text
.season_index <- function(df, season, seasons) {
  x <- df[[season]]
  if (is.character(seasons)) x <- as.character(x)
  match(x, seasons)
}

# number of sites with automated or verified data in each season
.sites_per_season <- function(x) {
  t_all <- .season_index(x$auto, x$season, x$seasons)
  s_all <- as.character(x$auto[[x$site]])
  if (!is.null(x$verified)) {
    t_all <- c(t_all, .season_index(x$verified, x$season, x$seasons))
    s_all <- c(s_all, as.character(x$verified[[x$site]]))
  }
  first <- !duplicated(data.frame(s_all, t_all))
  tabulate(t_all[first], nbins = length(x$seasons))
}

# the species: the levels of a factor, otherwise the sorted values of `auto`
# and `verified` together, as text
.species_levels <- function(auto, verified, species) {
  a <- auto[[species]]
  v <- if (is.null(verified)) NULL else verified[[species]]
  if (anyNA(a) || anyNA(v)) {
    stop("Column `", species, "` has missing values.", call. = FALSE)
  }
  if (!is.factor(a)) return(sort(unique(as.character(c(a, v)))))
  bad <- setdiff(as.character(v), levels(a))
  if (length(bad)) {
    stop("Species in `verified` are not levels of `auto$", species, "` (e.g. ",
         paste(utils::head(bad, 3), collapse = ", "), ").", call. = FALSE)
  }
  levels(a)
}

# per species: automated occasions and detections, sites with an automated
# detection, naive occupancy, verified occasions and detections
.species_table <- function(x) {
  sp <- factor(as.character(x$auto[[x$species]]), levels = x$species_names)
  det <- x$auto[[x$y]] == 1
  n_sites <- function(rows) {
    vapply(split(x$auto[[x$site]][rows], sp[rows]), function(s) length(unique(s)), 1L)
  }
  vsp <- if (is.null(x$verified)) factor(character(0), levels = x$species_names) else
    factor(as.character(x$verified[[x$species]]), levels = x$species_names)
  vdet <- if (is.null(x$verified)) logical(0) else x$verified[[x$y]] == 1
  data.frame(species = x$species_names,
             occasions = as.vector(table(sp)),
             detections = as.vector(tapply(det, sp, sum, default = 0L)),
             sites_detected = n_sites(det),
             naive_occ = round(n_sites(det) / pmax(n_sites(rep(TRUE, length(sp))), 1), 2),
             verified = as.vector(table(vsp)),
             verified_det = as.vector(tapply(vdet, vsp, sum, default = 0L)))
}

.warn_species <- function(auto, verified, site, y, species, species_names) {
  warn <- .design_warning
  sp <- factor(as.character(auto[[species]]), levels = species_names)
  some <- function(nm) paste0(" (e.g. ", paste(utils::head(nm, 3), collapse = ", "), ")")
  never <- species_names[tapply(auto[[y]] == 1, sp, sum, default = 0L) == 0]
  if (length(never)) {
    warn(length(never), " species have no automated detection", some(never),
         "; their occupancy is informed mostly by the other species.")
  }
  vdet <- if (is.null(verified)) character(0) else
    as.character(verified[[species]][verified[[y]] == 1])
  unver <- setdiff(species_names, vdet)
  if (length(unver)) {
    warn(length(unver), " species have no verified detection", some(unver),
         "; their false-positive rates are informed only by the other species ",
         "and the priors.")
  }
  n <- table(auto[[site]], sp)
  uneven <- rownames(n)[apply(n, 1, function(r) length(unique(r)) > 1)]
  if (length(uneven)) {
    warn(length(uneven), " site(s) have unequal numbers of occasions among species",
         some(uneven), "; every surveyed occasion should be present for every ",
         "species, with y = 0 where it was not detected.")
  }
  invisible(TRUE)
}

.warn_seasons <- function(auto, verified, confirmed, season, seasons) {
  warn <- .design_warning
  nt <- length(seasons)
  t_auto <- .season_index(auto, season, seasons)
  t_ver <- if (is.null(verified)) integer(0) else .season_index(verified, season, seasons)
  sampled <- tabulate(c(t_auto, t_ver), nbins = nt) > 0
  if (any(!sampled)) {
    warn(sum(!sampled), " season(s) have no data at any site (e.g. ",
         paste(utils::head(seasons[!sampled], 3), collapse = ", "), "); they ",
         "enter only through the transitions of occupancy.")
  }
  t_check <- t_ver
  if (!is.null(confirmed)) t_check <- c(t_check, t_auto[!is.na(auto[[confirmed[1]]])])
  unchecked <- sampled & tabulate(t_check, nbins = nt) == 0
  if (any(unchecked)) {
    warn(sum(unchecked), " season(s) have no verified or confirmed data (e.g. ",
         paste(utils::head(seasons[unchecked], 3), collapse = ", "), "); there ",
         "colonization and extinction are hard to tell from false positives.")
  }
  invisible(TRUE)
}

.check_confirmed <- function(auto, confirmed, y) {
  if (length(confirmed) != 2) {
    stop("`confirmed` must name two columns: confirmed true and confirmed ",
         "false observations.", call. = FALSE)
  }
  n_true <- auto[[confirmed[1]]]
  n_false <- auto[[confirmed[2]]]
  if (!identical(is.na(n_true), is.na(n_false))) {
    stop("Columns `", confirmed[1], "` and `", confirmed[2], "` must be NA ",
         "on the same occasions.", call. = FALSE)
  }
  ok <- !is.na(n_true)
  if (any(n_true[ok] < 0 | n_false[ok] < 0 |
          n_true[ok] %% 1 != 0 | n_false[ok] %% 1 != 0)) {
    stop("Confirmed counts must be non-negative integers.", call. = FALSE)
  }
  if (any((n_true[ok] + n_false[ok] > 0) != (auto[[y]][ok] == 1))) {
    stop("On confirmed occasions, a detection must have at least one ",
         "confirmed observation and a non-detection none.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_sites_known <- function(v, known, what, site) {
  bad <- setdiff(unique(v), known)
  if (length(bad)) {
    stop(length(bad), " site(s) in `", what, "` are missing from `sites` (e.g. ",
         paste(utils::head(bad, 3), collapse = ", "), ").", call. = FALSE)
  }
  invisible(TRUE)
}
