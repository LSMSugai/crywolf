# The arguments of the stable functions are frozen from 2.0.0
# (help("crywolf-api")): names, order and defaults change only in a major
# version. A deliberate change updates the snapshot with
# testthat::snapshot_accept("api").
test_that("the arguments of the stable functions do not change", {
  experimental <- c("fp_occu_mc", "fp_draws_data", "fp_to_unmarked", "fp_occu_ms")
  ns <- asNamespace("crywolf")
  # re-exports (as_draws_df(), loo()) belong to their own packages
  reexports <- c("as_draws_df", "loo")
  exports <- sort(setdiff(getNamespaceExports("crywolf"), c(experimental, reexports)),
                  method = "radix")
  # S3 methods by the name of the function that implements them
  s3 <- getNamespaceInfo(ns, "S3methods")
  methods <- sort(ifelse(is.na(s3[, 3]), paste(s3[, 1], s3[, 2], sep = "."), s3[, 3]),
                  method = "radix")
  methods <- methods[!grepl("fp_occu_(mc|ms)$", methods)]
  sig <- function(name, f) {
    args <- vapply(formals(f), function(a) {
      if (identical(a, quote(expr = ))) return("")
      paste0(" = ", gsub("\\s+", " ", paste(deparse(a), collapse = " ")))
    }, "")
    paste0(name, "(", paste0(names(args), args, collapse = ", "), ")")
  }
  out <- c(vapply(exports, function(n) sig(n, getExportedValue("crywolf", n)), ""),
           vapply(methods, function(n) sig(n, get(n, envir = ns)), ""))
  expect_snapshot(cat(out, sep = "\n"))
})
