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
