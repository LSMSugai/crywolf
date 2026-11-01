libs <- file.path(R_PACKAGE_DIR, "libs", R_ARCH)
dir.create(libs, recursive = TRUE, showWarnings = FALSE)
for (file in c("symbols.rds", Sys.glob(paste0("*", SHLIB_EXT)))) {
  if (file.exists(file)) {
    file.copy(file, file.path(libs, file))
  }
}
inst_stan <- file.path("..", "inst", "stan")
if (dir.exists(inst_stan)) {
  warning(
    "Stan models in inst/stan/ are deprecated in {instantiate} ",
    ">= 0.0.4.9001 (2024-01-03). Please put them in src/stan/ instead."
  )
  if (file.exists("stan")) {
    warning("src/stan/ already exists. Not copying models from inst/stan/.")
  } else {
    message("Copying inst/stan/ to src/stan/.")
    fs::dir_copy(path = inst_stan, new_path = "stan")
  }
}
bin <- file.path(R_PACKAGE_DIR, "bin")
if (!file.exists(bin)) {
  dir.create(bin, recursive = TRUE, showWarnings = FALSE)
}
bin_stan <- file.path(bin, "stan")
fs::dir_copy(path = "stan", new_path = bin_stan)
callr::r(
  func = function(bin_stan) {
    # crywolf needs CmdStan 2.40 or later (tuples as reduce_sum() arguments).
    # Without it, install anyway: fp_compile() compiles the model on first
    # use after CmdStan is installed or updated, and says how otherwise.
    v <- if (requireNamespace("cmdstanr", quietly = TRUE)) {
      tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NULL)
    }
    # devtools::load_all() installs the package into a temporary library
    # each time; it uses the model in the user cache instead
    if (grepl("devtools_install_", bin_stan, fixed = TRUE)) {
      return(invisible(NULL))
    }
    if (is.null(v) || utils::compareVersion(v, "2.40.0") < 0) {
      message("crywolf: CmdStan 2.40 or later not found (found ",
              if (is.null(v)) "none" else v, "); the Stan model was not ",
              "compiled now and will be compiled on first use after CmdStan ",
              "is installed or updated with cmdstanr::install_cmdstan().")
      return(invisible(NULL))
    }
    # a failed compilation does not stop the installation: fp_compile()
    # compiles the model on first use and reports any error then
    tryCatch(
      instantiate::stan_package_compile(
        models = instantiate::stan_package_model_files(path = bin_stan),
        cpp_options = list(stan_threads = TRUE)
      ),
      error = function(e) {
        message("crywolf: the Stan model could not be compiled now (",
                conditionMessage(e), "); it will be compiled on first use.")
      }
    )
  },
  args = list(bin_stan = bin_stan),
  show = TRUE,
  stderr = "2>&1"
)
