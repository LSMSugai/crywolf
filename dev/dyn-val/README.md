# Validation of the dynamic model before it stops being experimental

`run_val.R` runs the designs that `dev/dyn-sim/` left out; its header gives
them. `dev/dyn-kfold/run_kfold.R` adds, on the curiango data, the models with
`markov` and the model of independent site-nights with `random_det`. The
criteria below were set on 04/11/2026, before these results existed.

## Criteria for `fp_occu_dyn()` to stop being experimental

Coverage is that of 95% intervals. Bias is on the probability scale.

| # | Criterion | Threshold |
|---|---|---|
| 1 | Convergence in simulation | in each cell, at least 95% of fits with R-hat < 1.01 and no divergences; no fit ends in error |
| 2 | Calibration, `psi1`, `col`, `ext`, `p10_bar`, `psi_fs` | pooled over the cells of each p11, verification and season level: coverage at least 0.90 and absolute bias at most 0.05, with p11 >= 0.30 and verification >= 10% |
| 3 | Designs outside criterion 2 | p11 = 0.15 or no verification: results reported in the help page as limits; they do not block |
| 4 | Observation-confirmation | `calib`, cells with confirmed occasions and `b_fp = 0`: as criterion 2 |
| 5 | False positives that vary | `calib`, `b_fp = 1`: `fp ~ x` meets criterion 2, and the help page says what `fp ~ 1` costs |
| 6 | SBC | new seeds, 1000 replicates: chi-square p > 0.005 for each of the 10 parameters (0.05 with Bonferroni) |
| 7 | `random_det` and `markov` | with the option on: coverage at least 0.90 for `sd_det_site` and `det_y_prev`, and criterion 2 for the other quantities. An option that fails stays experimental, without blocking the function |
| 8 | Real data | a second data set besides curiango, with PPC, sensitivity to the priors of `fp_gap` and a comparison with the model of independent seasons |
| 9 | Interface | one release (2.2.0) with no change to arguments, classes or output names, then the names added to the API snapshot test |
| 10 | Documentation | help page with the limits of criteria 3 and 7, and `fp_calibrate()` either accepting dynamic designs or saying that it does not |

Criteria 1 to 7 come from these runs. Criteria 8 to 10 do not need AWS. The
comparison of `dyn`, `dyn_mk`, `ind_mk` and `ind_rd` by cross-validation
informs the vignette and is not a criterion.

## Outcome of criteria 1 to 7

The run of 04-05/11/2026 (c7a.8xlarge spot instances in us-east-2, from commit
3a338a9) finished all 6,320 tasks, and no fit ended in error.
Tables are in `results.md`.

| # | Outcome | Detail |
|---|---|---|
| 1 | Not met as written | 9 of the 144 `lowp` cells fall below 95% (lowest 0.60), all with p11 = 0.15 and no verification. `random_det` with sd = 1 reaches 0.93: 7 fits have R-hat of 1.01 to 1.02, mostly in the detection intercept, with no divergences. With 2000 + 2000 iterations it reaches 0.99. All other cells reach at least 0.967 |
| 2 | Met | 45 of 45 groups; lowest coverage 0.906, largest absolute bias 0.024 |
| 3 | Limits to report | p11 = 0.15 with no verification: bias of `col` from 0.12 to 0.17 and of `ext` from 0.09 to 0.13. With 10% verified: bias of `col` and `ext` above 0.05 with 3 or 6 seasons |
| 4 | Met | Both models, 10% and 30% of occasions confirmed |
| 5 | Met | `fp ~ x` in all three verification forms. `fp ~ 1` costs coverage of `p10_bar` (0.72 to 0.89) and a bias of up to 0.018 in `psi_fs` |
| 6 | Met | Lowest p = 0.039 (`verified_(Intercept)`); `ext_(Intercept)` has p = 0.99 |
| 7 | `markov` met; `random_det` met, with longer chains for large site effects | Coverage of `sd_det_site` 0.93 and 0.97, of `det_y_prev` 0.94 and 0.90. `random_det` with sd = 1 meets criterion 1 only with 2000 + 2000 iterations |

`run_randet_long.R` reran the cell of `random_det` with sd = 1 on the same data
and seeds, after these results were known. One of its 100 fits has R-hat of 1.01,
and none has divergences. The estimates did not change: the bias and coverage of
every quantity stayed within 0.01 of the first run.

Slopes such as `fp_x` are not among the quantities of criterion 2. Its coverage
was 0.89 to 0.97. Coverage from 100 replicates has a standard error near 0.02,
so values at 0.90 do not separate from the threshold.

## Run on AWS

With the recipe of `dev/dyn-sim/aws/` (see `dev/dyn-sim/README.md`). Commit
first, then:

```bash
export AWS_PROFILE=pcps JOBS=dev/dyn-val/jobs
B=pcps-sim-560246495553
dev/dyn-sim/aws/pack_code.sh $B crywolf-dyn-val
ID=$(dev/dyn-sim/aws/launch_spot.sh $B crywolf-dyn-val dev/dyn-val/run_val.R c7a.16xlarge 62 pcps-sim-ec2 8)
dev/dyn-sim/aws/supervise.sh $B crywolf-dyn-val dev/dyn-val/run_val.R 62 pcps-sim-ec2 8 $ID c7a.16xlarge c7i.16xlarge m7a.16xlarge
aws s3 sync s3://$B/fpcalib/runs/crywolf-dyn-val/jobs dev/dyn-val/jobs
```

The rerun of `random_det` goes the same way, with its own run name, with
`JOBS=dev/dyn-val/jobs-long`, the script `dev/dyn-val/run_randet_long.R` and
30 workers on a c7a.8xlarge. Then:

```bash
Rscript dev/dyn-val/summarise.R
```
