# Simulations of the dynamic model (steps 6 and 7)

`run_sims.R` runs the identification study (step 6), the recovery study and
the simulation-based calibration (step 7) of `fp_occu_dyn()`; its header
gives the design. Each task is saved in `jobs/` (ignored by git) and skipped
on a rerun. `summarise.R` reads the job files and writes `results.rds` and
the tables of `results.md`.

Locally, with crywolf installed from this checkout:

```bash
REPS=1 WORKERS=8 Rscript dev/dyn-sim/run_sims.R    # test, one replicate per cell
```

On AWS (spot), with the recipe of the compendium fpcalib (`aws/`), adapted in
`aws/`: the account, bucket, instance profile and quota are those of pcps-sim
and fpcalib (profile `pcps`, region us-east-1, prefix `fpcalib/runs/`, profile
`pcps-sim-ec2`, 192 spot vCPUs). Commit first, then:

```bash
export AWS_PROFILE=pcps
B=pcps-sim-560246495553
dev/dyn-sim/aws/pack_code.sh $B crywolf-dyn-test
dev/dyn-sim/aws/launch_spot.sh $B crywolf-dyn-test dev/dyn-sim/run_sims.R c7a.2xlarge 8 pcps-sim-ec2 2 1
dev/dyn-sim/aws/pack_code.sh $B crywolf-dyn
ID=$(dev/dyn-sim/aws/launch_spot.sh $B crywolf-dyn dev/dyn-sim/run_sims.R c7a.48xlarge 190 pcps-sim-ec2 6)
dev/dyn-sim/aws/supervise.sh $B crywolf-dyn dev/dyn-sim/run_sims.R 190 pcps-sim-ec2 6 $ID c7a.48xlarge c7i.48xlarge m7a.48xlarge
aws s3 sync s3://$B/fpcalib/runs/crywolf-dyn/jobs dev/dyn-sim/jobs
Rscript dev/dyn-sim/summarise.R
```
