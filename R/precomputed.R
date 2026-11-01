# Stamps for the precomputed results of the vignettes ---------------------------
#
# Results that take hours to compute are saved in vignettes/*.rds by the
# scripts in dev/. Each carries the MD5 of the Stan program that made it (or
# that the posterior regression test showed to give the same posterior), so
# that a vignette can say when its results may be out of date.

# of the program and of the functions it includes
.stan_md5 <- function() {
  files <- c(fp_stan_file(),
             file.path(dirname(fp_stan_file()), "fp_emissions.stanfunctions"))
  paste(tools::md5sum(files), collapse = "-")
}

.stamp <- function(x) {
  attr(x, "crywolf_stan_md5") <- .stan_md5()
  attr(x, "crywolf_version") <- as.character(utils::packageVersion("crywolf"))
  x
}

# a note for the vignette when the stamp does not match the Stan program in
# use; an empty string otherwise
.stale_note <- function(x, file) {
  if (identical(attr(x, "crywolf_stan_md5"), .stan_md5())) return("")
  paste0("*The precomputed results on this page (`", file, "`) were made with ",
         "another version of the Stan program than the one in this version of ",
         "crywolf and may be out of date. They are refreshed by ",
         "`dev/precompute-all.R` in the source repository.*\n\n")
}
