# Reruns the cell of random_det with sd = 1, on the same data and seeds, with
# 2000 + 2000 iterations: 7 of its 100 fits had R-hat of 1.01 to 1.02 with
# 1000 + 1000 and missed criterion 1 (dev/dyn-val/README.md).
# Run from the package root: WORKERS=30 Rscript dev/dyn-val/run_randet_long.R
Sys.setenv(PARTS = "randet", CELLS = "2", ITER = "2000",
           JOB_DIR = file.path("dev", "dyn-val", "jobs-long"))
source(file.path("dev", "dyn-val", "run_val.R"))
