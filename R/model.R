# Compilation and caching of the Stan model -----------------------------------

.fp_env <- new.env(parent = emptyenv())

#' Compile (or load) the Stan model
#'
#' Compiles the false-positive occupancy model the first time it is needed and
#' keeps the executable in the user cache directory
#' (`tools::R_user_dir("crywolf", "cache")`), so later sessions load it
#' without recompiling. [fp_occu()] calls this function itself; call it
#' directly to compile ahead of time, for example before launching parallel
#' workers in [fp_calibrate()].
#'
#' The model is compiled with threading enabled (`stan_threads`), so that
#' [fp_occu()] can split the likelihood of each chain among several threads
#' (argument `threads_per_chain`).
#'
#' @param force Logical. Recompile even if a current executable exists.
#' @param quiet Logical. Suppress compiler output.
#'
#' @return A `CmdStanModel` object (invisibly).
#' @export
#' @examples
#' \dontrun{
#' fp_compile()
#' }
fp_compile <- function(force = FALSE, quiet = TRUE) {
  .compile_model("fp_occu", force = force, quiet = quiet)
}

# Compiles (or loads) one of the Stan programs of the package, `fp_occu` or
# `fp_occu_dyn`; both #include fp_emissions.stanfunctions
.compile_model <- function(name, force = FALSE, quiet = TRUE) {
  if (!force && !is.null(.fp_env[[name]])) return(invisible(.fp_env[[name]]))
  .check_cmdstan()
  # the executable compiled when the package was installed from source (if
  # it is missing, stan_package_model() would compile one without threads)
  src <- .stan_file(name)
  exe <- file.path(dirname(src),
                   if (.Platform$OS.type == "windows") paste0(name, ".exe") else name)
  mod <- if (!force && file.exists(exe)) {
    tryCatch(instantiate::stan_package_model(name = name, package = "crywolf"),
             error = function(e) NULL)
  }
  if (is.null(mod) || !.has_threads(mod)) {
    # binary packages (such as those on CRAN), or an installation without a
    # recent enough CmdStan, come without it: compile once into the user cache
    dir <- tools::R_user_dir("crywolf", which = "cache")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    # the program and the functions it includes; a change to either
    # recompiles
    files <- c(src, file.path(dirname(src), "fp_emissions.stanfunctions"))
    changed <- FALSE
    for (f in files) {
      dst_f <- file.path(dir, basename(f))
      if (!file.exists(dst_f) || !identical(readLines(f), readLines(dst_f))) {
        file.copy(f, dst_f, overwrite = TRUE)
        changed <- TRUE
      }
    }
    dst <- file.path(dir, basename(src))
    threads <- list(stan_threads = TRUE)
    mod <- suppressWarnings(
      cmdstanr::cmdstan_model(dst, quiet = quiet, force_recompile = force || changed,
                              cpp_options = threads)
    )
    if (!.has_threads(mod)) {
      mod <- cmdstanr::cmdstan_model(dst, quiet = quiet, force_recompile = TRUE,
                                     cpp_options = threads)
    }
  }
  .fp_env[[name]] <- mod
  invisible(mod)
}

#' Path to the Stan program
#'
#' @return The path to the Stan file shipped with the package, for users who
#'   want to read the model or modify a copy of it. The likelihood of the
#'   observations of each site is in `fp_emissions.stanfunctions`, in the
#'   same folder, which the program includes.
#' @export
#' @examples
#' fp_stan_file()
fp_stan_file <- function() {
  .stan_file("fp_occu")
}

.stan_file <- function(name) {
  file <- paste0(name, ".stan")
  path <- system.file("bin", "stan", file, package = "crywolf")
  # under devtools::load_all() the program is still in the sources
  if (!nzchar(path)) {
    path <- system.file("src", "stan", file, package = "crywolf", mustWork = TRUE)
  }
  path
}

.has_threads <- function(mod) {
  opts <- mod$cpp_options()
  v <- unlist(opts[tolower(names(opts)) == "stan_threads"])
  isTRUE(as.logical(v[1]))
}

.min_cmdstan <- "2.40.0"

# the version of CmdStan that cmdstanr finds, or NULL
.cmdstan_found <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(NULL)
  path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  if (!nzchar(path)) return(NULL)
  tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NULL)
}

.cmdstan_recent <- function(v) utils::compareVersion(v, .min_cmdstan) >= 0

# a CmdStan that crywolf can use: found, and recent enough (tests and
# vignettes skip the fits otherwise)
.has_cmdstan <- function() {
  v <- .cmdstan_found()
  !is.null(v) && .cmdstan_recent(v)
}

.check_cmdstan <- function() {
  v <- .cmdstan_found()
  if (is.null(v)) {
    stop("CmdStan was not found. Install it once with ",
         "`cmdstanr::install_cmdstan()` and try again.", call. = FALSE)
  }
  if (!.cmdstan_recent(v)) {
    stop("crywolf needs CmdStan 2.40 or later (found ", v, "). Update it ",
         "with `cmdstanr::install_cmdstan()`.", call. = FALSE)
  }
  invisible(TRUE)
}
