# a small export with calibration draws, as sheepdog writes it: two sites,
# grid of scores from 0.1 in steps of 0.1, two draws; `edit` changes the
# calibration table before it is written (the parquet is written once: on
# Windows arrow keeps a file it has read mapped, so it cannot be rewritten)
write_cal_export <- function(dir, edit = identity) {
  dir.create(dir, showWarnings = FALSE)
  w <- function(d, name) utils::write.csv(d, file.path(dir, name), row.names = FALSE, na = "")
  w(data.frame(site = c("A", "A", "B", "B"), occasion_id = c("a1", "a2", "b1", "b2"),
               night = c("n1", "n2", "n1", "n1")), "occasions.csv")
  w(data.frame(site = c("A", "A", "B", "B"), occasion_id = c("a1", "a1", "b1", "b2"),
               species_code = "Sp one_sp-one", species = "sp_one",
               file_offset = c(0, 3, 0, 0), confidence = c(0.92, 0.18, 0.5, 0.31)),
    "detections.csv")
  w(data.frame(site = "B", occasion_id = "b2", species_code = "Sp one_sp-one",
               species = "sp_one", inclusion_prob = 0.1, verdict = "positive", present = 1),
    "recordings_verified.csv")
  w(data.frame(site = c("A", "B"), x = c(-1, 1)), "sites.csv")
  grid <- expand.grid(grid_index = 1:9, site = c("A", "B"), draw = c(3L, 8L),
                      stringsAsFactors = FALSE)
  # draw 3: every prediction true; draw 8: true only from a score of 0.9
  grid$p_tp <- ifelse(grid$draw == 3L, 1, as.numeric(grid$grid_index >= 9))
  grid$species <- "sp_one"
  grid$confidence <- 0.1 * grid$grid_index
  arrow::write_parquet(edit(grid), file.path(dir, "calibration_draws.parquet"))
  jsonlite::write_json(list(schema = "sheepdog-occupancy", schema_version = 1L,
                            sheepdog_version = "3.2.1.9000", occasion_unit = "file",
                            min_confidence = 0.1,
                            calibration = list(engine = "brms", n_draws = 2L,
                                               grid_step = 0.1, grid_min = 0.1)),
                       file.path(dir, "schema.json"), auto_unbox = TRUE)
  dir
}

test_that("each calibration draw gives one version of the detections", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export")))
  v <- fp_draws_data(x, rule = "threshold", p_min = 0.5)
  expect_length(v, 2)
  expect_equal(attr(v, "draws"), c(3L, 8L))
  expect_s3_class(v[[1]], "fp_data")
  # draw 3: every occasion with a prediction; draw 8: only a1 (score 0.92)
  expect_equal(v[[1]]$auto$y, c(1L, 0L, 1L, 1L))
  expect_equal(v[[2]]$auto$y, c(1L, 0L, 0L, 0L))
  expect_equal(v[[1]]$verified$y, 1L)
  # probabilities of 0 and 1 make the sampled versions equal to these
  s <- fp_draws_data(x, seed = 1)
  expect_equal(s[[1]]$auto$y, v[[1]]$auto$y)
  expect_equal(s[[2]]$auto$y, v[[2]]$auto$y)
  expect_length(fp_draws_data(x, n_draws = 1), 1)
})

test_that("a calibration by month is matched by the month of recording", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  # January: every prediction true; February: none
  by_month <- function(cal) rbind(transform(cal, month = "202401", p_tp = 1),
                                  transform(cal, month = "202402", p_tp = 0))
  dir <- write_cal_export(file.path(tempdir(), "cal-export-month"), by_month)
  occ <- utils::read.csv(file.path(dir, "occasions.csv"))
  occ$recording_date <- c("20240101", "20240201", "20240101", "20240201")
  utils::write.csv(occ, file.path(dir, "occasions.csv"), row.names = FALSE)
  sc <- jsonlite::read_json(file.path(dir, "schema.json"))
  sc$calibration$keys <- list("site", "species", "month", "grid_index")
  jsonlite::write_json(sc, file.path(dir, "schema.json"), auto_unbox = TRUE)
  x <- fp_read_sheepdog(dir)
  v <- fp_draws_data(x, rule = "threshold", p_min = 0.5)
  # a1 and b1 recorded in January, b2 in February (a2 has no prediction)
  expect_equal(v[[1]]$auto$y, c(1L, 0L, 1L, 0L))
})

test_that("predictions without calibration leave their occasions unsurveyed", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export-nocal"),
                                         function(cal) cal[cal$site == "A", ]))
  expect_warning(
    v <- withCallingHandlers(fp_draws_data(x, rule = "threshold", p_min = 0.5),
                             crywolf_design_warning = function(w) invokeRestart("muffleWarning")),
    "2 occasion\\(s\\) are treated as not surveyed")
  # fp_data() keeps surveyed occasions only: site B drops out of auto
  expect_equal(v[[1]]$auto$y, c(1L, 0L))
  expect_equal(unique(v[[1]]$auto$site), "A")
  dir <- write_cal_export(file.path(tempdir(), "cal-export-nocal2"),
                          function(cal) cal[cal$site == "C", ])
  expect_error(suppressWarnings(fp_draws_data(fp_read_sheepdog(dir))),
               "No prediction could be matched")
})

test_that("keys of schema.json and the parquet must agree", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  # a month column the schema does not announce
  dir <- write_cal_export(file.path(tempdir(), "cal-export-keys"),
                          function(cal) rbind(transform(cal, month = "202401"),
                                              transform(cal, month = "202402")))
  expect_error(fp_draws_data(fp_read_sheepdog(dir)), "more than one value")
  # a key the schema announces and the parquet lacks
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export-keys3")))
  x$schema$calibration$keys <- c("site", "species", "month", "grid_index")
  expect_error(fp_draws_data(x), "disagree.*lacks month")
  # a key column missing from both: two values per key and draw
  dir <- write_cal_export(file.path(tempdir(), "cal-export-keys2"),
                          function(cal) rbind(cal, cal))
  expect_error(fp_draws_data(fp_read_sheepdog(dir)), "more than one value")
})

test_that("the random numbers of the session are left untouched", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export-rng")))
  set.seed(42)
  ref <- stats::runif(1)
  set.seed(42)
  invisible(fp_draws_data(x, seed = 1))
  expect_identical(stats::runif(1), ref)
})

test_that("duplicated occasions are refused", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export-dup")))
  occ <- rbind(x$occasions, x$occasions[1, ])
  expect_error(fp_draws_data(x, occasions = occ), "one row per site")
})

test_that("the sites of the model can be other than the recorders", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export-nights")))
  occ <- x$occasions
  occ$site_night <- paste(occ$site, occ$night)
  nights <- data.frame(site_night = unique(occ$site_night))
  v <- fp_draws_data(x, rule = "threshold", occasions = occ, sites = nights,
                     site = "site_night")
  expect_equal(nrow(v[[1]]$sites), 3)
  expect_equal(v[[1]]$verified$site_night, "B n1")
})

test_that("an export without calibration draws is refused", {
  skip_if_not_installed("jsonlite")
  dir <- write_cal_export(file.path(tempdir(), "cal-export-none"))
  unlink(file.path(dir, "calibration_draws.parquet"))
  x <- fp_read_sheepdog(dir)
  expect_error(fp_draws_data(x), "no calibration draws")
})

test_that("the keys of the calibration are columns of the occasions", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  # a calibration by recorder: r1 every prediction true, r2 none
  by_rec <- function(cal) rbind(transform(cal, recorder = "r1", p_tp = 1),
                                transform(cal, recorder = "r2", p_tp = 0))
  dir <- write_cal_export(file.path(tempdir(), "cal-export-recorder"), by_rec)
  occ <- utils::read.csv(file.path(dir, "occasions.csv"))
  occ$recorder <- c("r1", "r1", "r2", NA)
  utils::write.csv(occ, file.path(dir, "occasions.csv"), row.names = FALSE, na = "")
  sc <- jsonlite::read_json(file.path(dir, "schema.json"))
  sc$calibration$keys <- list("site", "species", "recorder", "grid_index")
  jsonlite::write_json(sc, file.path(dir, "schema.json"), auto_unbox = TRUE)
  x <- fp_read_sheepdog(dir)
  # b2 has no recorder: it cannot be scored
  expect_warning(
    v <- withCallingHandlers(fp_draws_data(x, rule = "threshold", p_min = 0.5),
                             crywolf_design_warning = function(w) invokeRestart("muffleWarning")),
    "site x species x recorder without calibration draws or without a value of a key: their 1 occasion")
  expect_equal(v[[1]]$auto$y, c(1L, 0L, 0L))
  # a key the occasions do not have
  x$occasions$recorder <- NULL
  expect_error(fp_draws_data(x), "lack the calibration key\\(s\\) recorder")
})

test_that("an export of sheepdog with a calibration by month is read", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  # written by sheepdog (data-raw/make_export_fixture.R of sheepdog), schema
  # 1.2; the calibration lacks the scores above 0.5 at NHE in 202112
  dir <- test_path("fixtures", "sheepdog-export")
  x <- suppressWarnings(fp_read_sheepdog(dir, species = "phy_cuv"))
  x$files <- x$files[names(x$files) != "occasion_draws"]   # computed here
  expect_warning(
    v <- withCallingHandlers(fp_draws_data(x, rule = "threshold", p_min = 0.5),
                             crywolf_design_warning = function(w) invokeRestart("muffleWarning")),
    "14 occasion\\(s\\) are treated as not surveyed")
  expect_length(v, 20)
  # the same histories, computed by hand from the parquet
  cal <- as.data.frame(arrow::read_parquet(file.path(dir, "calibration_draws.parquet")))
  det <- x$detections
  occ <- paste(x$occasions$site, x$occasions$occasion_id)
  det$month <- x$occasions$month[match(paste(det$site, det$occasion_id), occ)]
  det$grid_index <- pmin(pmax(round((det$confidence - 0.1) / 0.05) + 1, 1),
                         max(cal$grid_index))
  d1 <- merge(det, cal[cal$draw == attr(v, "draws")[1], ],
              by = c("site", "species", "month", "grid_index"), all.x = TRUE)
  unscored <- unique(paste(d1$site, d1$occasion_id)[is.na(d1$p_tp)])
  hit <- unique(paste(d1$site, d1$occasion_id)[!is.na(d1$p_tp) & d1$p_tp >= 0.5])
  expect_length(unscored, 14)
  # fp_data() keeps the surveyed occasions only
  expect_identical(v[[1]]$auto$y, as.integer(setdiff(occ, unscored) %in% hit))
  expect_gt(sum(v[[1]]$auto$y), 0)
})

test_that("the histories drawn by sheepdog are those of fp_draws_data()", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  quiet <- function(e) suppressWarnings(withCallingHandlers(
    e, crywolf_design_warning = function(w) invokeRestart("muffleWarning")))
  dir <- test_path("fixtures", "sheepdog-export")
  x <- suppressWarnings(fp_read_sheepdog(dir, species = "phy_cuv"))
  computed <- x
  computed$files <- computed$files[names(computed$files) != "occasion_draws"]
  y <- function(v) lapply(v, function(d) d$auto$y)
  # rule "sample", seed 2026, 20 draws in occasion_draws.parquet
  expect_identical(y(quiet(fp_draws_data(x, seed = 2026))),
                   y(quiet(fp_draws_data(computed, seed = 2026))))
  expect_identical(y(quiet(fp_draws_data(x, seed = 2026, n_draws = 5))),
                   y(quiet(fp_draws_data(computed, seed = 2026)))[1:5])

  # the table is read only for the settings it was drawn with: a copy with
  # every scored occasion detected tells which way the values came
  copy <- file.path(tempdir(), "sheepdog-export-draws")
  dir.create(copy, showWarnings = FALSE)
  file.copy(list.files(dir, full.names = TRUE), copy, overwrite = TRUE)
  od <- as.data.frame(arrow::read_parquet(file.path(dir, "occasion_draws.parquet")))
  od$y[!is.na(od$y)] <- 1L
  unlink(file.path(copy, "occasion_draws.parquet"))
  arrow::write_parquet(od, file.path(copy, "occasion_draws.parquet"))
  x <- suppressWarnings(fp_read_sheepdog(copy, species = "phy_cuv"))
  # detections in the first draw: every scored occasion of the table
  all_one <- function(v) sum(v[[1]]$auto$y) == sum(!is.na(od$y[od$draw == min(od$draw)]))
  expect_true(all_one(quiet(fp_draws_data(x, seed = 2026))))
  expect_false(all_one(quiet(fp_draws_data(x, seed = 1))))
  expect_false(all_one(quiet(fp_draws_data(x))))
  expect_false(all_one(quiet(fp_draws_data(x, rule = "threshold"))))
  # an occasion with predictions left out of `occasions`: computed here
  drop <- x$occasions$occasion_id != x$detections$occasion_id[1]
  expect_false(all_one(quiet(fp_draws_data(x, seed = 2026, occasions = x$occasions[drop, ]))))
})

test_that("a calibration with continuous covariates gives sheepdog's histories", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  quiet <- function(e) suppressWarnings(withCallingHandlers(
    e, crywolf_design_warning = function(w) invokeRestart("muffleWarning")))
  # schema 1.3: tp ~ logit_conf + temp_z, 20 histories drawn by sheepdog with
  # rule "sample" and seed 2026, and no calibration_draws.parquet
  dir <- test_path("fixtures", "sheepdog-export-1.3")
  x <- suppressWarnings(fp_read_sheepdog(dir, species = "phy_cuv"))
  expect_false("calibration_draws" %in% names(x$files))
  expect_output(print(x), "temp_z")
  od <- as.data.frame(arrow::read_parquet(file.path(dir, "occasion_draws.parquet")))
  v <- quiet(fp_draws_data(x))
  expect_length(v, 20)
  expect_identical(attr(v, "draws"), sort(unique(od$draw)))
  expect_identical(attr(v, "rule"), "sample")
  # by hand: the value sheepdog drew, 0 for the occasions without predictions
  first <- od[od$draw == min(od$draw), ]
  m <- match(paste(x$occasions$site, x$occasions$occasion_id),
             paste(first$site, first$occasion_id))
  expect_identical(v[[1]]$auto$y, ifelse(is.na(m), 0L, first$y[m]))
  # the settings it was drawn with are the only ones there are
  expect_identical(lapply(quiet(fp_draws_data(x, seed = 2026, n_draws = 5)),
                          function(d) d$auto$y),
                   lapply(v[1:5], function(d) d$auto$y))
  expect_error(fp_draws_data(x, seed = 1), "seed = 2026, not 1")
  expect_error(fp_draws_data(x, rule = "threshold"), "rule = \"sample\"")
  expect_error(fp_draws_data(x, n_draws = 30), "n_draws = 20, not 30")
  without <- x
  without$files <- without$files[names(without$files) != "occasion_draws"]
  expect_error(fp_draws_data(without), "continuous covariates \\(temp_z\\)")

  # occasions sheepdog could not score are not surveyed; occasions left out
  # of `occasions` are dropped
  copy <- file.path(tempdir(), "sheepdog-export-1.3-na")
  dir.create(copy, showWarnings = FALSE)
  file.copy(list.files(dir, full.names = TRUE), copy, overwrite = TRUE)
  na_occ <- first$occasion_id[1]
  od$y[od$occasion_id == na_occ] <- NA
  unlink(file.path(copy, "occasion_draws.parquet"))
  arrow::write_parquet(od, file.path(copy, "occasion_draws.parquet"))
  x <- suppressWarnings(fp_read_sheepdog(copy, species = "phy_cuv"))
  design <- function(e) withCallingHandlers(
    e, crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
  expect_warning(v <- design(fp_draws_data(x)), "1 occasion\\(s\\) with predictions were not scored")
  # fp_data() keeps only the surveyed occasions
  expect_false(na_occ %in% v[[1]]$auto$occasion_id)
  expect_equal(nrow(v[[1]]$auto), nrow(x$occasions) - 1)
  keep <- !x$occasions$occasion_id %in% first$occasion_id[1:2]
  expect_warning(v <- design(fp_draws_data(x, occasions = x$occasions[keep, ])),
                 "2 occasion\\(s\\) with predictions are not in `occasions`")
  expect_equal(nrow(v[[1]]$auto), sum(keep))
})

test_that("an export of several species gives the data of all of them", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  # sp_two is never true under any draw
  dir <- write_cal_export(file.path(tempdir(), "cal-export-ms"), edit = function(g) {
    rbind(g, transform(g, species = "sp_two", p_tp = 0))
  })
  det <- utils::read.csv(file.path(dir, "detections.csv"))
  two <- transform(det[1:2, ], species_code = "Sp two_sp-two", species = "sp_two")
  utils::write.csv(rbind(det, two), file.path(dir, "detections.csv"), row.names = FALSE)
  x <- fp_read_sheepdog(dir, species = c("sp_one", "sp_two"))
  v <- suppressWarnings(fp_draws_data(x, rule = "threshold", p_min = 0.5))
  d <- v[[1]]
  expect_equal(d$species, "species")
  expect_equal(d$species_names, c("sp_one", "sp_two"))
  expect_equal(d$auto$species, rep(c("sp_one", "sp_two"), each = 4))
  # draw 3: every prediction of sp_one is true
  expect_equal(d$auto$y, c(1, 0, 1, 1, 0, 0, 0, 0))
  expect_equal(unname(d$known_present[, "sp_one"]), c(0L, 1L))
  expect_equal(unname(d$known_present[, "sp_two"]), c(0L, 0L))
})

test_that("the multispecies export of sheepdog gives its histories", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  # written by sheepdog (sheepdog#74), schema 1.3: three species, one of them
  # rare, calibration with a continuous covariate, 20 histories drawn
  dir <- test_path("fixtures", "sheepdog-export-1.3-ms")
  schema <- jsonlite::read_json(file.path(dir, "schema.json"), simplifyVector = TRUE)
  x <- suppressWarnings(fp_read_sheepdog(dir, species = schema$species))
  expect_equal(sort(unique(x$detections$species)), c("ela_bic", "lep_fus", "phy_cuv"))
  v <- suppressWarnings(fp_draws_data(x, n_draws = 2))
  expect_length(v, 2)
  d <- v[[2]]
  expect_equal(d$species_names, c("ela_bic", "lep_fus", "phy_cuv"))
  expect_equal(nrow(d$auto), 3 * nrow(x$occasions))
  # the histories are those of occasion_draws.parquet, 0 where absent
  od <- as.data.frame(arrow::read_parquet(file.path(dir, "occasion_draws.parquet")))
  od <- od[od$draw == attr(v, "draws")[2], ]
  key <- function(t) paste(t$site, t$occasion_id, t$species)
  expect_equal(d$auto$y[match(key(od), key(d$auto))], od$y)
  expect_equal(sum(d$auto$y), sum(od$y, na.rm = TRUE))
  # the rare species is found only at ARA and BAI
  tab <- crywolf:::.species_table(d)
  expect_lte(tab$sites_detected[tab$species == "lep_fus"], 2)
})
