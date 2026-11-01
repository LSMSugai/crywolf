# crywolf <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/diogoprov/crywolf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/diogoprov/crywolf/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23385045.svg)](https://doi.org/10.5281/zenodo.23385045)
[![R-universe](https://diogoprov.r-universe.dev/crywolf/badges/version)](https://diogoprov.r-universe.dev/crywolf)
[![Codecov test coverage](https://codecov.io/gh/diogoprov/crywolf/graph/badge.svg)](https://app.codecov.io/gh/diogoprov/crywolf)
<!-- badges: end -->

Bayesian false-positive occupancy models for data from automated detectors,
with simulation-based checks of whether the correction is calibrated for
your design.

Acoustic classifiers and other automated detectors return false positives.
**crywolf** fits the multiple-methods (site-confirmation) false-positive
occupancy model of Miller et al. (2011) in Stan: the classifier output may be
wrong, and a verified subset of detections is certain. It also fits the
observation-confirmation design of Chambert et al. (2015), in which an expert
checks the observations of some automated occasions. The false-positive
rate can be constant or depend on site covariates, which matters when false
positives track the same environmental gradient as the species.

## Installation

crywolf is distributed through
[R-universe](https://diogoprov.r-universe.dev), which provides binaries for
Windows, macOS and Linux. The command below also installs 'cmdstanr' from
the Stan R-universe; CmdStan itself is installed once with
`install_cmdstan()`:

```r
install.packages("crywolf", repos = c("https://diogoprov.r-universe.dev",
                                      "https://stan-dev.r-universe.dev",
                                      getOption("repos")))
cmdstanr::install_cmdstan()
```

crywolf needs CmdStan 2.40 or later (since crywolf 1.7.0). If you
installed CmdStan before September 2026, update it with
`cmdstanr::install_cmdstan()`; `cmdstanr::cmdstan_version()` shows the
version in use.

The development version can be installed from GitHub with
`remotes::install_github("diogoprov/crywolf")`.

## Example

```r
library(crywolf)

# simulated data in which false positives track a gradient correlated
# with the occupancy covariate
sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
                   conf_cor = 0.7, seed = 2026)

# the simulated recorders have no effect, so leave the recorder term out
m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)  # constant rate
m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)  # rate depends on x

coef(m0, "occ"); coef(m1, "occ")                # true slope: -0.8
loo::loo_compare(loo(m0), loo(m1))

# is a constant rate calibrated for this design, and does verifying more help?
fp_verification_curve(nsite = 300, nocc = 20, conf_effect = 1.5,
                      conf_cor = 0.7,
                      models = list(constant = ~ 1, fp_on_x = ~ x),
                      reps = 50)
```

## Main functions

| Function | Purpose |
|---|---|
| `fp_read_birdnet()`, `fp_read_validation()`, `fp_aggregate()` | read BirdNET selection tables and reviewed clips, and build the automated occasions |
| `fp_threshold_sweep()` | fit the model at several classifier thresholds |
| `fp_data()`, `fp_data_wide()` | assemble automated, verified, and confirmed detections and site covariates |
| `fp_occu()` | fit the model; methods `print`, `summary`, `coef`, `predict`, `loo`, `as_draws_df` |
| `fp_priors()`, `fp_gap_prior()`, `fp_prior_predictive()` | prior settings, a prior from validated detections, and what the priors imply for data like yours |
| `fp_ppc()` | posterior predictive check of the detections per site |
| `plot()`, `fp_effects()`, `fp_sites()` | coefficients and diagnostics, marginal effects, and the probability that each site is occupied |
| `simulate_fp()` | simulate data under homogeneous, unstructured, or structured false positives |
| `fp_calibrate()` | coverage and bias of the credible intervals for a design, by simulation |
| `fp_verification_curve()` | coverage against the fraction of sites verified |
| `fp_calibrate_fit()`, `fp_sbc()` | calibration on data simulated from a fit, and simulation-based calibration of the computation |
| `loo()`, `fp_loo_weights()` | leave-one-site-out cross-validation (exact refits for unreliable sites) and stacking weights |
| `plot()` | coverage or bias of a calibration run |
| `curiango` | acoustic monitoring of the Common Pauraque in the Pantanal (Piatti et al. 2026) |

## Related tools

- [fpcalib](https://github.com/diogoprov/fpcalib) is the research compendium
  of the simulation study behind the calibration checks
  (https://doi.org/10.5281/zenodo.23379100). `simulate_fp(seed = s)`
  reproduces its datasets exactly.
- [sheepdog](https://github.com/diogoprov/sheepdog)
  (https://provetelab.org/sheepdog/,
  https://doi.org/10.5281/zenodo.22696210) is an R package that cuts audio
  clips from BirdNET predictions for manual review and records the verdict
  on each clip, which can feed the verified or confirmed data of
  `fp_data()`. It is distributed through the same R-universe as crywolf:
  `install.packages("sheepdog", repos = c("https://diogoprov.r-universe.dev",
  getOption("repos")))`. From version 2.0.0 it writes `clip_manifest.csv` and
  `validation_by_clip.csv` with English column names;
  `fp_read_validation()` reads these and the files of earlier versions
  (`validacao_por_clipe.csv`). The models in
  crywolf treat the checked recordings as chosen independently of the
  classifier's confidence score. Clips picked because they have the highest
  scores contain fewer false positives than a random sample and would
  understate the false-positive rate, so draw the clips at random for this
  use.

## Citation

If you use crywolf, please cite the package (`citation("crywolf")`;
https://doi.org/10.5281/zenodo.23385045, which always resolves to the latest
version) and
Miller et al. (2011), or Chambert et al. (2015) for the observation-confirmation
design.

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C., Bailey,
L. L., & Weir, L. A. (2011). Improving occupancy estimation when two types of
observational error occur: non-detection and species misidentification.
*Ecology*, 92, 1422–1428. https://doi.org/10.1890/10-1396.1

Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332–339. https://doi.org/10.1890/14-1507.1

## Why "crywolf"?

The name comes from Aesop's fable *The Boy Who Cried Wolf* (Perry Index
210), in which a shepherd boy raises the alarm about a wolf so many times
when there is none that nobody believes him when the wolf finally comes
([Wikipedia](https://en.wikipedia.org/wiki/The_Boy_Who_Cried_Wolf)). An
automated detector is that boy: it sometimes cries "species!" when the
species is absent. Ignoring the false alarms biases occupancy estimates, and
discarding every alarm wastes the true ones. crywolf takes the middle path:
it uses a verified subset of detections to estimate how often the detector
cries wolf, and corrects for it.

The logo is a nod to the 1980s and to
["Cry Wolf"](https://en.wikipedia.org/wiki/Cry_Wolf_(A-ha_song)), the 1986
single by the Norwegian band a-ha from the album *Scoundrel Days*. It borrows
the synthwave look of the decade: a neon magenta frame on deep purple, a
striped setting sun, and a rounded display typeface. The wolf howls in sound
waves over a cyan audio equalizer, standing for the acoustic recordings
behind the package. The package has no connection with
the band, and the logo uses none of its artwork or lyrics.

## License

GPL (>= 3). The wolf silhouette in the logo is *Canis lupus monstrabilis* by
Tracy A. Heath, from PhyloPic, in the public domain (CC0 1.0).
