# Prints crywolf-cheatsheet.html to PDF with Chrome, for the download on the
# site. Output: pkgdown/assets/cheatsheet/crywolf-cheatsheet.pdf
# Run from the package root. Needs Google Chrome; set CHROME to its path if it
# is elsewhere.
chrome <- Sys.getenv("CHROME",
                     "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
source_file <- normalizePath("data-raw/cheatsheet/crywolf-cheatsheet.html")
dir.create("pkgdown/assets/cheatsheet", recursive = TRUE, showWarnings = FALSE)
output <- file.path(getwd(), "pkgdown/assets/cheatsheet/crywolf-cheatsheet.pdf")
status <- system2(chrome, c("--headless", "--disable-gpu",
                            "--no-pdf-header-footer",
                            shQuote(paste0("--print-to-pdf=", output)),
                            shQuote(paste0("file://", source_file))))
stopifnot(status == 0, file.exists(output))
