test_that("fp_data validates its inputs", {
  sites <- data.frame(site = 1:3, x = c(-1, 0, 1))
  auto <- data.frame(site = rep(1:3, each = 2), y = c(0, 1, 0, 0, 1, 1))
  expect_s3_class(fp_data(auto, sites), "fp_data")
  expect_error(fp_data(transform(auto, y = y * 2), sites), "0/1")
  expect_error(fp_data(auto, rbind(sites, sites)), "one row per site")
  expect_error(fp_data(transform(auto, site = site + 1), sites), "missing from")
  expect_error(fp_data(auto[, "site", drop = FALSE], sites), "lacks column")
})

test_that("missing detections are dropped and known_present is set", {
  sites <- data.frame(site = 1:3)
  auto <- data.frame(site = 1:3, y = c(1, NA, 0))
  ver <- data.frame(site = c(1, 1, 3), y = c(0, 1, 0))
  # site 2 keeps no automated occasion once its missing detection is dropped
  expect_warning(d <- fp_data(auto, sites, verified = ver), "no automated occasion")
  expect_equal(nrow(d$auto), 2)
  expect_equal(d$known_present, c(1L, 0L, 0L))
})

test_that("fp_data_wide matches the long layout", {
  ya <- matrix(c(1, 0, NA, 1, 0, 0), 3, 2)
  yv <- matrix(c(1, NA, NA, 0, NA, NA), 3, 2)
  d <- fp_data_wide(ya, yv, site_covs = data.frame(x = 1:3),
                    obs_covs = list(bio = matrix(1:6, 3, 2)),
                    recorder = c("a", "b", "a"))
  expect_equal(nrow(d$auto), 5)
  expect_equal(nrow(d$verified), 2)
  expect_true(all(c("bio", "recorder") %in% names(d$auto)))
  expect_equal(d$known_present, c(1L, 0L, 0L))
  expect_error(fp_data_wide(ya, obs_covs = list(bio = matrix(1, 2, 2))),
               "same dimensions")
})

test_that("confirmed occasions are validated and mark sites as occupied", {
  auto <- data.frame(site = c(1, 1, 2, 2, 3), y = c(1, 0, 1, 0, 1),
                     n_true = c(1, 0, 0, NA, NA), n_false = c(2, 0, 1, NA, NA))
  sites <- data.frame(site = 1:3)
  d <- fp_data(auto, sites, confirmed = c("n_true", "n_false"))
  expect_equal(d$known_present, c(1L, 0L, 0L))
  expect_output(print(d), "confirmed occasions:    3 \\(1 with a true")

  bad <- auto; bad$n_false[4] <- 0
  expect_error(fp_data(bad, sites, confirmed = c("n_true", "n_false")), "same occasions")
  bad <- auto; bad$n_true[2] <- 1
  expect_error(fp_data(bad, sites, confirmed = c("n_true", "n_false")), "detection must")
  bad <- auto; bad$n_true[1] <- -1
  expect_error(fp_data(bad, sites, confirmed = c("n_true", "n_false")), "non-negative")
  expect_error(fp_data(auto, sites, confirmed = "n_true"), "two columns")
})

test_that("an occasion both verified and confirmed is warned about", {
  sites <- data.frame(site = 1:2)
  conf <- c("n_true", "n_false")
  auto <- data.frame(site = c(1, 1, 2, 2), occasion_id = c("a", "b", "a", "b"),
                     y = c(1, 0, 1, 0), n_true = c(1, NA, 0, NA),
                     n_false = c(0, NA, 1, NA))
  expect_warning(fp_data(auto, sites, confirmed = conf,
                         verified = data.frame(site = 1, occasion_id = "a", y = 1)),
                 "1 verified occasion", class = "crywolf_design_warning")

  # the same identifier at another site, or an unchecked occasion, is no overlap
  expect_no_warning(fp_data(auto, sites, confirmed = conf,
    verified = data.frame(site = c(1, 2), occasion_id = c("b", "b"), y = 1)))
  # without identifiers in both tables, or without `confirmed`, nothing to compare
  expect_no_warning(fp_data(auto, sites, confirmed = conf,
                            verified = data.frame(site = 1, y = 1)))
  expect_no_warning(fp_data(auto[names(auto) != "occasion_id"], sites, confirmed = conf,
    verified = data.frame(site = 1, occasion_id = "a", y = 1)))
  expect_no_warning(fp_data(auto, sites,
    verified = data.frame(site = 1, occasion_id = "a", y = 1)))
  # missing identifiers do not match each other
  na_id <- transform(auto, occasion_id = NA_character_)
  expect_no_warning(fp_data(na_id, sites, confirmed = conf,
    verified = data.frame(site = 1, occasion_id = NA_character_, y = 1)))
})

test_that("the overlap of verified and confirmed occasions is by species", {
  auto <- data.frame(site = rep(1:2, each = 2), occasion_id = "a",
                     species = rep(c("s1", "s2"), 2), y = c(1, 0, 0, 1),
                     n_true = c(1, NA, NA, 1), n_false = c(1, NA, NA, 0))
  warned <- function(sp) {
    w <- character()
    withCallingHandlers(
      fp_data(auto, data.frame(site = 1:2), species = "species",
              confirmed = c("n_true", "n_false"),
              verified = data.frame(site = 1, occasion_id = "a", species = sp, y = 1)),
      crywolf_design_warning = function(e) {
        w <<- c(w, conditionMessage(e))
        invokeRestart("muffleWarning")
      })
    any(grepl("1 verified occasion", w))
  }
  expect_true(warned("s1"))
  expect_false(warned("s2"))   # same site and occasion, unchecked for this species
})

test_that("the overlap of verified and confirmed occasions is by season", {
  # occasion identifiers that restart in each season
  auto <- data.frame(site = rep(1:2, each = 2), season = rep(1:2, 2), occasion_id = 1,
                     y = c(1, 0, 0, 1), n_true = c(1, NA, NA, 1),
                     n_false = c(1, NA, NA, 0))
  warned <- function(s) {
    w <- character()
    withCallingHandlers(
      fp_data(auto, data.frame(site = 1:2), season = "season",
              confirmed = c("n_true", "n_false"),
              verified = data.frame(site = 1, season = s, occasion_id = 1, y = 1)),
      crywolf_design_warning = function(e) {
        w <<- c(w, conditionMessage(e))
        invokeRestart("muffleWarning")
      })
    any(grepl("1 verified occasion", w))
  }
  expect_true(warned(1))
  expect_false(warned(2))
})

test_that("print separates verified sites, confirmed occasions, and known occupancy", {
  s <- simulate_fp(nsite = 60, nocc = 8, frac_confirmed = 0.4, frac_verified = 0, seed = 1)
  out <- capture.output(print(s$data))
  expect_true(any(grepl("verified sites: +0 \\(0 occasions; verified detections at 0 sites\\)", out)))
  n_true <- s$data$auto$n_true
  true_sites <- length(unique(s$data$auto$site[!is.na(n_true) & n_true > 0]))
  expect_true(any(grepl(paste0("at ", true_sites, " sites\\)"), out)))
  expect_true(any(grepl(paste0("sites known occupied: +", sum(s$data$known_present)), out)))
})

test_that("fp_data holds several species in long format", {
  sites <- data.frame(site = 1:3)
  auto <- expand.grid(site = 1:3, occasion = 1:2, sp = c("b", "a"),
                      stringsAsFactors = FALSE)
  auto$y <- c(1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  ver <- data.frame(site = c(1, 2), sp = c("b", "a"), y = c(1, 0))
  w <- character()
  d <- withCallingHandlers(fp_data(auto, sites, ver, species = "sp"),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  expect_match(w, "1 species have no automated detection \\(e.g. a\\)", all = FALSE)
  expect_match(w, "1 species have no verified detection \\(e.g. a\\)", all = FALSE)
  expect_equal(d$species_names, c("a", "b"))
  expect_equal(d$known_present,
               matrix(c(0L, 0L, 0L, 1L, 0L, 0L), 3, dimnames = list(NULL, c("a", "b"))))
  tab <- crywolf:::.species_table(d)
  expect_equal(tab$detections, c(0, 3))
  expect_equal(tab$sites_detected, c(0, 2))
  expect_equal(tab$verified_det, c(0, 1))
  expect_output(print(d), "by species")
  # a factor sets the species and their order
  auto$sp <- factor(auto$sp, levels = c("b", "a", "c"))
  d <- suppressWarnings(fp_data(auto, sites, ver, species = "sp"))
  expect_equal(colnames(d$known_present), c("b", "a", "c"))
})

test_that("fp_data with species checks the layout and refuses other uses", {
  sites <- data.frame(site = 1:2)
  auto <- data.frame(site = c(1, 1, 2, 2, 1), sp = c("a", "b", "a", "b", "a"),
                     y = c(1, 0, 1, 1, 0), night = 1)
  ver <- data.frame(site = 1:2, sp = c("a", "b"), y = 1)
  expect_warning(d <- fp_data(auto, sites, ver, species = "sp"),
                 "unequal numbers of occasions among species \\(e.g. 1\\)")
  expect_error(fp_data(auto, sites, species = "sp", season = "night"), "cannot be combined")
  expect_error(fp_data(auto, sites, ver[, c("site", "y")], species = "sp"), "lacks column")
  expect_error(fp_data(transform(auto, sp = NA), sites, species = "sp"), "missing values")
  expect_error(fp_data(transform(auto, sp = factor(sp)), sites,
                       transform(ver, sp = "z"), species = "sp"), "not levels")
  expect_error(fp_occu_dyn(d), "season column")
  expect_error(fp_occu(d), "several species")
  expect_error(fp_prior_predictive(d), "several species")
  expect_error(fp_occu_mc(list(d, d)), "several species")
})
