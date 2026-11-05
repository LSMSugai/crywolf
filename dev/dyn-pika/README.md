# Second data set for the dynamic model (criterion 8)

`run.R` fits `fp_occu_dyn()` to the American pika surveys of Goldman et al.
(2023, *Ecology and Evolution* 13: e10019, doi:10.1002/ece3.10019). Its header
gives the design. `post_process.csv` comes unchanged from Dryad
(doi:10.5061/dryad.ttdz08m2n, CC0). The data were collected by the OSU Cascades
Human and Ecosystem Resilience and Sustainability Lab with support from the
National Park Service, and uses of them should attribute that source.

The data have 771 plots in 5 park units and 5,126 surveys from 2010 to 2021. Each
fit takes about 4 minutes on 8 cores, so it runs locally:

```bash
Rscript dev/dyn-pika/run.R
```

## Run of 05/11/2026

All five fits have 0 divergences, R-hat of at most 1.004 and bulk ESS of at least
1,109.

| Requirement | Result |
|---|---|
| PPC, automated | p = 0.78 (`dyn`); p = 0.21 (`ind`) |
| PPC, verified | p = 0.16 (`dyn`) |
| Priors of `fp_gap` | means of 1, 2 and 3 and an SD of 3 change no reported estimate by more than 0.023 on its own scale |
| Independent seasons | `dyn` is ahead of `ind` by 178.7 elpd (SE 42.6), with no Pareto k above 0.7 |

Colonization is 0.048 and extinction 0.071 per year, and `fp_gap` is 3.14.
The article reports about 0.06 for both with a different model, with a false
detection rate for each sign. The comparison with it is only of order of
magnitude.

Two simplifications matter. The verified occasion shares its survey with the
automated ones, and unambiguous scat is a subset of fresh scat, so the two methods
are not fully independent given the state. And the sign enters `det` and `fp` as
a covariate, whereas the article gives each sign its own model.

Goldman et al. treat calls as ambiguous, although the false calls they expect
come from pikas just outside the plot, not from other species.
