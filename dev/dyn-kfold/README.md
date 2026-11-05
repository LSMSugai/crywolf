# Blocked cross-validation of the dynamic model (step 8)

`run_kfold.R` compares, on the curiango data, `fp_occu_dyn()` with the nights
of each site as seasons against `fp_occu()` with independent site-nights, with
and without `markov`, and the latter with `random_det`, by 10-fold
cross-validation over blocks of 7 nights; its header gives the design. LOO is not an option, because the log_lik of the dynamic model is
one value per site and curiango has seven sites. Each of the 50 tasks (5
models x 10 folds) is saved in `jobs/` (ignored by git) and skipped on a
rerun. `summarise.R` reads the job files and writes `results.rds` and the
tables of `results.md`.

A fit takes about 46 minutes with 4 chains on one core each (local pilot,
fold 1, Mac). On AWS, with the recipe of `dev/dyn-sim/aws/` (see
`dev/dyn-sim/README.md`) and `JOBS` pointing at this job folder. Commit
first, then:

```bash
export AWS_PROFILE=pcps JOBS=dev/dyn-kfold/jobs
B=pcps-sim-560246495553
dev/dyn-sim/aws/pack_code.sh $B crywolf-kfold
ID=$(dev/dyn-sim/aws/launch_spot.sh $B crywolf-kfold dev/dyn-kfold/run_kfold.R c7a.24xlarge 20 pcps-sim-ec2 4)
dev/dyn-sim/aws/supervise.sh $B crywolf-kfold dev/dyn-kfold/run_kfold.R 20 pcps-sim-ec2 4 $ID c7a.24xlarge c7i.24xlarge m7a.24xlarge
aws s3 sync s3://$B/fpcalib/runs/crywolf-kfold/jobs dev/dyn-kfold/jobs
Rscript dev/dyn-kfold/summarise.R
```

## Run of 01/11/2026

On a c7a.24xlarge spot instance, in 2 h (fits of 47 to 66 minutes), from
commit 4cb8ec7. The training fit of `ind`, fold 8 (seed 91018), did not
converge: one chain stayed in another mode (p11 about 0.14 and p10 about
0.07, against 0.06 and 0.004 in the other chains and folds), with R-hat 1.54
and bulk ESS 7. It was rerun on the Mac with `SEED_OFFSET=100` (seed 91118),
and converged (R-hat 1.008, bulk ESS 809). Its elpd changed little (-805.9
against -806.1 with the stuck chain). The other 19 fits have R-hat <= 1.011
and bulk ESS >= 610; one has a divergence. Results in `results.md`.

## Run of 04-05/11/2026

The 30 new fits (`dyn_mk`, `ind_mk`, `ind_rd`) ran on c7a.8xlarge spot instances
in us-east-2, from commit 3a338a9, as run `crywolf-kfold2`, which reused the 20
job files of the run above. All 50 fits have R-hat <= 1.011 and bulk ESS >= 547;
three have one divergence each.

`dyn_mk` has the highest elpd: 206.9 (SE 33.4) above `ind_mk`, 361.5 (SE 91.5)
above `ind_rd` and 481.3 (SE 120.1) above `dyn`. It is ahead of the other models
at all sites and in all folds. `random_det` gains 311.3 (SE 100.1) over `ind`,
and `markov` 465.9 over `ind`. These standard errors treat the units as
independent. The elpd of `ind_rd` comes from 400 draws, against 4,000 for the
other models. This comparison informs the vignette and is not a criterion of
`dev/dyn-val/README.md`.
