# an export of the sheepdog-occupancy schema, version 1, written by hand as
# sheepdog::export_occupancy_data() writes it (empty cells for NA)
pauraque <- "Nyctidromus albicollis_Common Pauraque"
write_export <- function(dir, version = 1L, validation = NULL) {
  dir.create(dir, showWarnings = FALSE)
  w <- function(d, name) utils::write.csv(d, file.path(dir, name), row.names = FALSE, na = "")
  w(data.frame(
    site = c("01", "01", "02", "02"), occasion_id = c("a1", "a2", "b1", "b2"),
    source_audio = c("/r/a1.wav", "/r/a2.wav", "/r/b1.wav", NA),
    recording_date = c("20240101", "20240101", "20240101", NA),
    recording_hour = c("20", "21", "20", NA), night = c("20240101", "20240101", "20240101", NA),
    start = c(0, 0, 0, NA), duration = NA), "occasions.csv")
  w(data.frame(
    site = c("01", "01", "02", "02", "02"), occasion_id = c("a1", "a1", "b1", "b1", "b2"),
    species_code = c(pauraque, pauraque, pauraque, "Other bird_Other", pauraque),
    species = c("compau", "compau", "compau", "other", "compau"),
    file_offset = c(0, 3, 6, 0, 0), confidence = c(0.9, 0.4, 0.7, 0.8, 0.3)),
    "detections.csv")
  if (is.null(validation)) {
    validation <- data.frame(
      site = c("01", "02", "02", "02"), occasion_id = c("a1", "b1", "b1", "b2"),
      species_code = pauraque, species = "compau", clip_start = c(0, 6, 6, 0),
      clip_end = c(3, 9, 9, 3), confidence = c(0.9, 0.7, 0.7, 0.3),
      sampling_mode = c("random_detections", "random_detections", "top", "random_detections"),
      inclusion_prob = c(0.5, 0.5, NA, 0.5), n_detections_occasion = c(2, 1, 1, 1),
      verdict = c("positive", "negative", "negative", "disagreement"), tp = c(1, 0, 0, NA),
      reviewer = c("ana", "ana", "ana", "ana; bia"), n_reviewers = c(1, 1, 1, 2))
  }
  w(validation, "validation.csv")
  w(data.frame(
    site = c("01", "02", "02"), occasion_id = c("a2", "b1", "b2"),
    species_code = c(pauraque, pauraque, "Other bird_Other"),
    species = c("compau", "compau", "other"), inclusion_prob = 0.1,
    verdict = c("positive", "not_reviewed", "negative"), present = c(1, NA, 0)),
    "recordings_verified.csv")
  w(data.frame(site = c("01", "02"), temp = c(20, 25)), "sites.csv")
  jsonlite::write_json(list(schema = "sheepdog-occupancy", schema_version = version,
                            sheepdog_version = "3.2.1.9000", occasion_unit = "file",
                            min_confidence = 0.1),
                       file.path(dir, "schema.json"), auto_unbox = TRUE)
  dir
}

test_that("an export is read into the tables of fp_aggregate() and fp_data()", {
  skip_if_not_installed("jsonlite")
  dir <- write_export(file.path(tempdir(), "sheepdog-export"))
  expect_warning(x <- suppressMessages(fp_read_sheepdog(dir, species = "compau")),
                 "sampling_mode = \"top\"")
  expect_s3_class(x, "fp_sheepdog")
  expect_identical(x$occasions$site, c("01", "01", "02", "02"))   # identifiers as text
  expect_true(is.na(x$occasions$source_audio[4]))   # an occasion with site and id only
  expect_equal(nrow(x$detections), 4)               # the other species dropped
  expect_true("start" %in% names(x$detections))
  expect_equal(nrow(x$validation), 3)               # the "top" clip dropped
  expect_identical(x$validation$verdict, c(TRUE, FALSE, NA))   # "disagreement" is NA
  expect_identical(x$validation$verdict_label[3], "disagreement")
  expect_equal(x$verified$y, 1L)                    # state unknown: dropped
  expect_output(print(x), "reviewed clips")
  # matched by `species` as well as by the full label
  expect_equal(nrow(suppressWarnings(suppressMessages(
    fp_read_sheepdog(dir, species = pauraque)))$detections), 4)
  auto <- fp_aggregate(x$occasions, x$detections, threshold = 0.5,
                       by = c("site", "occasion_id"), validation = x$validation)
  expect_equal(auto$y, c(1L, 0L, 1L, 0L))
  expect_equal(auto$n_true, c(1L, NA, 0L, NA))
  d <- fp_data(auto, x$sites, x$verified, confirmed = c("n_true", "n_false"))
  expect_s3_class(d, "fp_data")
})

test_that("verdicts in both tables of an export are warned about", {
  skip_if_not_installed("jsonlite")
  dir <- write_export(file.path(tempdir(), "sheepdog-export-overlap"))
  x <- suppressWarnings(suppressMessages(fp_read_sheepdog(dir, species = "compau")))
  auto <- fp_aggregate(x$occasions, x$detections, threshold = 0.5,
                       by = c("site", "occasion_id"), validation = x$validation)
  ver <- x$verified
  ver$occasion_id <- "a1"      # a reviewed clip, annotated as a recording too
  expect_warning(fp_data(auto, x$sites, ver, confirmed = c("n_true", "n_false")),
                 "1 verified occasion", class = "crywolf_design_warning")
})

test_that("clips of the top design are dropped, or kept once", {
  skip_if_not_installed("jsonlite")
  dir <- write_export(file.path(tempdir(), "sheepdog-export-top"))
  w <- character()
  k <- withCallingHandlers(suppressMessages(fp_read_sheepdog(dir, species = "compau",
                                                             keep_top = TRUE)),
                           warning = function(e) {
                             w <<- c(w, conditionMessage(e))
                             invokeRestart("muffleWarning")
                           })
  expect_match(w, "Keeping 1", all = FALSE)
  # the top clip is the same clip as a random one: kept once
  expect_equal(nrow(k$validation), 3)
  top_only <- data.frame(site = "01", occasion_id = "a1", species_code = pauraque,
                         species = "compau", clip_start = 0, clip_end = 3,
                         confidence = 0.9, sampling_mode = "top", inclusion_prob = NA,
                         n_detections_occasion = 2, verdict = "positive", tp = 1)
  dir2 <- write_export(file.path(tempdir(), "sheepdog-export-toponly"),
                       validation = top_only)
  t <- suppressWarnings(suppressMessages(fp_read_sheepdog(dir2, species = "compau")))
  expect_null(t$validation)
})

test_that("an export of several species read without `species` warns", {
  skip_if_not_installed("jsonlite")
  dir <- write_export(file.path(tempdir(), "sheepdog-export-species"))
  w <- character()
  withCallingHandlers(suppressMessages(fp_read_sheepdog(dir)),
                      warning = function(e) {
                        w <<- c(w, conditionMessage(e))
                        invokeRestart("muffleWarning")
                      })
  expect_match(w, "several species", all = FALSE)
})

test_that("other schemas and versions are refused", {
  skip_if_not_installed("jsonlite")
  dir <- write_export(file.path(tempdir(), "sheepdog-export-v2"), version = 2L)
  expect_error(fp_read_sheepdog(dir), "Version 2")
  # minor versions only add to the schema
  dir <- write_export(file.path(tempdir(), "sheepdog-export-v1.1"), version = 1.1)
  expect_s3_class(suppressWarnings(fp_read_sheepdog(dir)), "fp_sheepdog")
  expect_error(fp_read_sheepdog(tempdir()), "No schema.json")
})

test_that("the export of sheepdog itself is read", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("sheepdog")
  skip_if_not(exists("export_occupancy_data", envir = asNamespace("sheepdog")),
              "this sheepdog has no export_occupancy_data()")
  ex <- system.file("extdata", "example", package = "sheepdog")
  skip_if(!nzchar(ex), "sheepdog has no example data")
  out <- file.path(tempdir(), "sheepdog-real-export")
  unlink(out, recursive = TRUE)
  det <- suppressMessages(sheepdog::read_selection_tables(file.path(ex, "birdnet")))
  suppressMessages(sheepdog::export_occupancy_data(
    out, det, recordings = file.path(ex, "recordings.csv"),
    validation = list.files(file.path(ex, "validation"), full.names = TRUE)))
  sp <- utils::read.csv(file.path(out, "detections.csv"))$species[1]
  x <- suppressWarnings(suppressMessages(fp_read_sheepdog(out, species = sp)))
  expect_s3_class(x, "fp_sheepdog")
  expect_true(all(x$detections$species == sp))
  expect_true(is.character(x$occasions$occasion_id))
  auto <- fp_aggregate(x$occasions, x$detections, threshold = 0.5,
                       by = c("site", "occasion_id"), validation = x$validation)
  expect_equal(nrow(auto), nrow(x$occasions))
})

test_that("an export of schema 1.2 brings occasion covariates and text keys", {
  skip_if_not_installed("jsonlite")
  # written by sheepdog (data-raw/make_export_fixture.R of sheepdog);
  # validation.csv trimmed to the target species, to keep the package small
  dir <- test_path("fixtures", "sheepdog-export")
  expect_warning(x <- fp_read_sheepdog(dir, species = "phy_cuv"), "sampling_mode = \"top\"")
  expect_equal(x$schema$schema_version, 1.2)
  expect_true(all(c("calibration_draws", "occasion_draws") %in% names(x$files)))
  expect_identical(x$schema$occasion_covariates, c("hour", "month"))
  expect_true(is.numeric(x$occasions$hour))
  expect_type(x$occasions$month, "character")   # a key of the calibration
  expect_identical(sort(unique(x$occasions$month)),
                   c("202110", "202111", "202112", "202201", "202202", "202203"))
  expect_identical(unique(x$detections$species), "phy_cuv")
})
