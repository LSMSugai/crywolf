# Refreshes every precomputed result of the vignettes, in vignettes/*.rds.
# Run after a change to the Stan program or to fp_occu() that alters the
# posterior; for a change that does not (a refactor or speed-up that passes
# tests/testthat/test-regression.R), dev/restamp-precomputed.R is enough.
#
# Takes several hours (the curiango fits are the longest). Most scripts use
# the installed package, so install the current version first; then run from
# the package root: Rscript dev/precompute-all.R
# Uses at most 10 cores; set OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and
# VECLIB_MAXIMUM_THREADS=1 when other heavy jobs share the machine.
# vignettes/site-effects-scenarios.rds is not refreshed here: it comes from
# dev/random-det/simulate-scenarios.R (1600 fits, about 11 hours), whose
# prototype Stan program is independent of the package.
scripts <- c("precompute-vignette.R", "precompute-calibration.R",
             "precompute-calibration-fit.R",
             "precompute-sbc.R", "precompute-sbc-random-det.R",
             "precompute-fp-variation.R", "precompute-markov.R",
             "precompute-curiango-random-det.R", "precompute-threshold.R",
             "precompute-threshold-fp.R", "precompute-dynamic.R")
run <- function(script, args = character()) {
  message("== ", script, " ", paste(args, collapse = " "), " (", format(Sys.time()), ")")
  status <- system2("Rscript", c(file.path("dev", script), args))
  if (status != 0) stop(script, " failed.", call. = FALSE)
}
purrr::walk(scripts, run)
# the curiango fits are run per classifier, then combined; partial outputs
# of earlier runs in dev/curiango-fits/ would be combined too, so clear them
unlink(list.files(file.path("dev", "curiango-fits"), pattern = "[.]rds$", full.names = TRUE))
purrr::walk(c("custom", "birdnet", "combine"),
            function(a) run("precompute-curiango.R", a))
message("All precomputed results refreshed (", format(Sys.time()), ").")
