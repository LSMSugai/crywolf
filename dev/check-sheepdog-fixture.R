# Checks that the sheepdog fixtures of the tests are still identical to those
# of sheepdog (tests/testthat/fixtures/ on its main branch), which sheepdog's
# contract test (test-contract.R) checks against the schema:
#   tests/testthat/fixtures/sheepdog-export      <- export-1.2
#   tests/testthat/fixtures/sheepdog-export-1.3  <- export-1.3
#   tests/testthat/fixtures/sheepdog-export-1.3-ms <- export-1.3-ms
# Files are compared by their git blob SHA, as listed by the GitHub API, so
# nothing is downloaded unless the copies differ.
# Run from the package root: Rscript dev/check-sheepdog-fixture.R [--update]
# With --update, the files that differ are downloaded and the ones sheepdog no
# longer has are removed; then rerun the tests of fp_read_sheepdog() and
# fp_draws_data().
update <- "--update" %in% commandArgs(trailingOnly = TRUE)
fixtures <- c("sheepdog-export" = "export-1.2",
              "sheepdog-export-1.3" = "export-1.3",
              "sheepdog-export-1.3-ms" = "export-1.3-ms")
differ <- FALSE
for (i in seq_along(fixtures)) {
  local <- file.path("tests/testthat/fixtures", names(fixtures)[i])
  api <- paste0("https://api.github.com/repos/diogoprov/sheepdog/contents/",
                "tests/testthat/fixtures/", fixtures[i], "?ref=main")
  remote <- jsonlite::fromJSON(api)
  remote <- remote[remote$type == "file", c("name", "sha", "download_url")]
  files <- list.files(local)
  sha <- if (length(files)) system2("git", c("hash-object", file.path(local, files)),
                                    stdout = TRUE) else character()
  changed <- remote$name[!remote$name %in% files |
                           sha[match(remote$name, files)] != remote$sha]
  gone <- setdiff(files, remote$name)
  if (!length(changed) && !length(gone)) {
    message(local, " is identical to sheepdog's ", fixtures[i], " (",
            length(files), " files).")
    next
  }
  differ <- TRUE
  if (length(changed)) message(local, " differs from sheepdog: ", toString(changed))
  if (length(gone)) message(local, " has files not in sheepdog: ", toString(gone))
  if (update) {
    dir.create(local, showWarnings = FALSE)
    for (f in changed) {
      utils::download.file(remote$download_url[remote$name == f], file.path(local, f),
                           mode = "wb", quiet = TRUE)
    }
    unlink(file.path(local, gone))
  }
}
if (differ && !update) {
  message("Rerun with --update to copy sheepdog's fixtures.")
  quit(status = 1)
}
if (differ) message("Updated. Rerun devtools::test(filter = \"read_sheepdog|draws_data\").")
